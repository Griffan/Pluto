//
// Created by Fan Zhang on 8/6/15.
//

#include "PBWTHaplotyper.h"
#include "MemoryAllocators.h"

//debug related
#define DEBUG false

static const float UNDERFLOW_MIN = std::numeric_limits<float>::min() * 1e2;

static void printLeftMatrix(float *probability, int numStates) {
    for (int i = 0; i < numStates; ++i) {
        for (int j = 0; j <= i; ++j, probability++) {
            fprintf(stderr, "(%d,%d):%9.9f\t", i, j, *probability);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}

//initiation
PBWTHaplotyper::PBWTHaplotyper(int nhaps, int nsnps) {
    markers = nsnps;
    individuals = nhaps / 2.0;
    nSampleCopy = 0;//additional, the original haps not included
}

PBWTHaplotyper::PBWTHaplotyper() {//TODO:change the number
    nSampleCopy = 0;//additional, the original haps not included
}

#ifdef HETERSITE
void PBWTHaplotyper::InitAuxillary() {

    tmpHaps = AllocateCharMatrix(individuals * 2 + (individuals - phased) * nSampleCopy * 2, markers);

    tmpGeno = AllocateCharMatrix(individuals, markers * 3);

    tmpPenetrance = new float[markers * 9];
}
#endif

PBWTHaplotyper::~PBWTHaplotyper() {
#ifdef HETERSITE
    if (tmpHaps != nullptr) {
        for (int i = 0; i < 2 * individuals + (individuals - phased) * nSampleCopy * 2; ++i) {
            delete[] tmpHaps[i];
        }
        delete[] tmpHaps;
    }
    if (tmpGeno != nullptr) {
        for (int i = 0; i < individuals; ++i) {
            delete[] tmpGeno[i];
        }
        delete[] tmpGeno;
    }
    if (tmpPenetrance != nullptr)
        delete[] tmpPenetrance;
#endif

    if (sampledHaps != nullptr) {
        for (int l = 0; l < nSampleCopy * GetUnphasedNum() * 2; ++l) {
            delete[] sampledHaps[l];
        }
        delete[] sampledHaps;
    }

    if (Wrapper != nullptr)
        delete Wrapper;

    if (genotypes) {
        for (int i = 0; i < GetUnphasedNum(); i++) {
            if (genotypes[i]) delete[] genotypes[i];
        }
        delete[] genotypes;
        genotypes = nullptr;
    }

    if (genoProbs) {
        for (int i = 0; i < GetUnphasedNum(); i++) {
            if (genoProbs[i]) delete[] genoProbs[i];
        }
        delete[] genoProbs;
        genoProbs = nullptr;
    }

    if (haplotypes) {
        int nTarget = runningModel & PHASE ? GetUnphasedNum() : individuals;
        for (int i = 0; i < nTarget; i++) {
            if (haplotypes[i * 2]) delete[] haplotypes[i * 2];
            if (haplotypes[i * 2] + 1) delete[] haplotypes[i * 2 + 1];
        }
        delete[] haplotypes;
        haplotypes = nullptr;
    }

    DestroyPvalueMatrix();
}


//random setup
void PBWTHaplotyper::InitialSampleCopy(Random *rand) {

    if (rand == nullptr)
        rand = &globalRandom;
    CalculatePhred2Prob();

    if (nSampleCopy == 0) return;
    sampledHaps = new char *[nSampleCopy * GetUnphasedNum() * 2];
    for (int l = 0; l < nSampleCopy * GetUnphasedNum() * 2; ++l) {
        sampledHaps[l] = new char[markers];
    }
    int nTarget = GetUnphasedNum();

    for (int j = 0; j < markers; j++) {
        double mac = 0;
        int markerIndex = 3 * j;

//        double hyperprior11 = freq1s[j] * freq1s[j];
//        double hyperprior12 = 2.0 * freq1s[j] * (1.0 - freq1s[j]);
//        double hyperprior22 = (1.0 - freq1s[j]) * (1.0 - freq1s[j]);
//
//
//        for (int i = 0; i < nTarget; i++) {
//            double post11 = hyperprior11 * phred2prob[(size_t) genotypes[i][markerIndex]];
//            double post12 = hyperprior12 * phred2prob[(size_t) genotypes[i][markerIndex + 1]];
//            double post22 = hyperprior22 * phred2prob[(size_t) genotypes[i][markerIndex + 2]];
//            double sumpost = post11 + post12 + post22;
//            post11 /= sumpost;
//            post12 /= sumpost;
//            post22 /= sumpost;
//
//            // estimated counts of AL2
//            mac += post12 + 2 * post22;
//        }
//
//        //here, each person contributes two alleles
//        double freq = 0.5 * mac / (double) nTarget;
//
//        double prior_11 = (1.0 - freq) * (1.0 - freq);
//        double prior_12 = 2.0 * freq * (1.0 - freq);
//        double prior_22 = freq * freq;

        double prior_11 = (1.0 - freq1s[j]) * (1.0 - freq1s[j]);
        double prior_12 = 2.0 * freq1s[j] * (1.0 - freq1s[j]);
        double prior_22 = freq1s[j] * freq1s[j];

        for (int i = 0; i < nTarget; i++) {

            double posterior_11 = prior_11 * phred2prob[(size_t) genotypes[i][markerIndex]];
            double posterior_12 = prior_12 * phred2prob[(size_t) genotypes[i][markerIndex + 1]];
            double posterior_22 = prior_22 * phred2prob[(size_t) genotypes[i][markerIndex + 2]];
            double sum = posterior_11 + posterior_12 + posterior_22;

            if (sum == 0)
                printf("Problem!\n");

            posterior_11 /= sum;
            posterior_12 /= sum;

            for (int k = 0; k < nSampleCopy; ++k) {
                double r = rand->Next();
                int index = i * nSampleCopy + k;
                if (r < posterior_11) {
                    sampledHaps[index * 2][j] = 0;
                    sampledHaps[index * 2 + 1][j] = 0;
                } else if (r < posterior_11 + posterior_12) {
                    bool bit = rand->Binary();

                    sampledHaps[index * 2][j] = bit;
                    sampledHaps[index * 2 + 1][j] = bit ^ 1;
                } else {
                    sampledHaps[index * 2][j] = 1;
                    sampledHaps[index * 2 + 1][j] = 1;
                }
            }
        }
    }
}

//phasing

void PBWTHaplotyper::SwapIndividuals(int a, int b) {
    // if (b < 0 || b >= individuals)
    //   printf("Bad Swap!");

    Swap(genotypes[a], genotypes[b]);
    Swap(haplotypes[a * 2], haplotypes[b * 2]);
    Swap(haplotypes[a * 2 + 1], haplotypes[b * 2 + 1]);

    if (nSampleCopy > 0) {
        Swap(sampledHaps[a * 2], sampledHaps[(GetUnphasedNum() - 1) * 2]);
        Swap(sampledHaps[a * 2 + 1], sampledHaps[(GetUnphasedNum() - 1) * 2 + 1]);
    }
    if (diseaseCount) {
        Swap(diseaseStatus[a], diseaseStatus[b]);
        Swap(diseaseScores[a], diseaseScores[b]);
    }

    if (weights != nullptr) {
        float temp = weights[a];
        weights[a] = weights[b];
        weights[b] = temp;
    }
}

void PBWTHaplotyper::PrepareRefSetPBWTWrapper() {
    if (Wrapper != nullptr) {
        delete Wrapper;
        Wrapper = nullptr;
    }
    Wrapper = new PBWTWrapper(2 * phased, markers, PvalueMatrix, 0);
    Wrapper->SetHaps(haplotypes, 2 * (individuals - phased), 2 * individuals, nullptr, 0, 0, thetas, 0);
    Wrapper->CursorBackwards();//calculate backwards order of suffix
    Wrapper->CursorForwards();
}

#ifdef NAIVE
int PBWTHaplotyper::LoopThroughChromosomesHighPrecision() {

    ResetCrossovers();

    if (isRev) ReverseInput();
    clock_t t = clock();
    if (Wrapper != nullptr) {
        delete Wrapper;
        Wrapper = nullptr;
    }
    printf("[HighPrecision]build model start...\n");
    Wrapper = new PBWTWrapper(2 * individuals + (individuals - phased) * nSampleCopy * 2, markers, PvalueMatrix, prefixLength);
    Wrapper->SetHaps(haplotypes, 0, 2 * individuals, sampledHaps, 0, (individuals - phased) * nSampleCopy * 2, thetas);
    Wrapper->CursorBackwards();//calculate backwards order of suffix
    Wrapper->CursorForwards();
    clock_t t1 = clock();
    printf("[HighPrecision]build model end time:%.2f sec\n", (float) (t1 - t) / CLOCKS_PER_SEC);
    for (int i = individuals - 1; i >= 0; i--) {

        if (i < individuals - phased) {
            SwapIndividuals(i, individuals - 1);

            fprintf(stderr, "[HighPrecision]phasing individual %d...\n", i);
            FwdBwdLocalParameter localParameter;
            InitialFwdValues(individuals - 1, localParameter);
            ForwardAlgorithm(individuals - 1, localParameter);
            BackwardSampling(&globalRandom, individuals - 1, haplotypes, localParameter);
            for (int j = 0; j < nSampleCopy; ++j) {//n copy per individual
                BackwardSampling(&globalRandom, j + i * nSampleCopy, sampledHaps, localParameter);
            }
//            t1=clock();
//            printf("sampling time:%.2f sec\n", (float) (t1-t) / CLOCKS_PER_SEC);

#ifdef _DEBUG
            if (!SanityCheck())
               {
               printf("\nProblems above occurred haplotyping individual %d\n\n", i);
               Print();
               }
#endif
            SwapIndividuals(i, individuals - 1);
        }
    }
    t = clock();
    printf("[HighPrecision]forward algorithm and sampling time:%.2f sec\n", (float) (t - t1) / CLOCKS_PER_SEC);
    if (isRev) ReverseInput();
    return 0;
}
#endif

#include <exception>
#include <Error.h>

class SamplingException : public std::exception {
    virtual const char *what() const throw() {
        return "Current sample encountered unexpected sampling space!\n";
    }
} samplingException;

#ifdef _OPENMP
#include <omp.h>
#endif

int PBWTHaplotyper::LoopThroughChromosomesPhasing(Pedigree &ped) {

    ResetCrossovers();
    if (isRev and not(runningModel & PHASE)) ReverseInput();
    ConstructGraph();
    clock_t t1 = clock();

    int nTarget = GetUnphasedNum();
#ifdef _OPENMP
    omp_set_num_threads(nThread);
#pragma omp parallel for
#endif
    for (int i = nTarget - 1; i >= 0; i--) {
//            SwapIndividuals(i, individuals - 1);
        fprintf(stderr, "[%s]phasing individual %d:%s...\n\n", __FUNCTION__, i, ped[i].pid.c_str());
        LocalForwardBackwardSampling(i);

#ifdef _DEBUG
        if (!SanityCheck())
           {
           printf("\nProblems above occurred haplotyping individual %d\n\n", i);
           Print();
           }
#endif
//            SwapIndividuals(i, individuals - 1);
    }
    clock_t t = clock();
    fprintf(stderr, "[%s]forward algorithm and sampling time:%.2f sec\n\n", __FUNCTION__,
            (float) (t - t1) / CLOCKS_PER_SEC);
    if (isRev and not(runningModel & PHASE)) ReverseInput();
    return 0;
}

int PBWTHaplotyper::LoopThroughChromosomesGenotyping(Pedigree &ped) {

    ResetCrossovers();
    if (isRev and not(runningModel & PHASE)) ReverseInput();
    ConstructGraph();
    clock_t t1 = clock();

    int nTarget = GetUnphasedNum();
#ifdef _OPENMP
    omp_set_num_threads(nThread);
#pragma omp parallel for
#endif
    for (int i = nTarget - 1; i >= 0; i--) {
//            SwapIndividuals(i, individuals - 1);
        fprintf(stderr, "[%s]phasing individual %d:%s...\n\n", __FUNCTION__, i, ped[i].pid.c_str());
        LocalForwardBackward(i);

#ifdef _DEBUG
        if (!SanityCheck())
           {
           printf("\nProblems above occurred haplotyping individual %d\n\n", i);
           Print();
           }
#endif
//            SwapIndividuals(i, individuals - 1);
    }
    clock_t t = clock();
    fprintf(stderr, "[%s]forward-backward algorithm time:%.2f sec\n\n", __FUNCTION__,
            (float) (t - t1) / CLOCKS_PER_SEC);
    if (isRev and not(runningModel & PHASE)) ReverseInput();
    return 0;
}

void PBWTHaplotyper::ConstructGraph() {
    if (loadGraph != "Empty") {//read graph
        clock_t t = clock();
        if (Wrapper != nullptr) {
            delete Wrapper;
            Wrapper = nullptr;
        }
        Wrapper = new PBWTWrapper(2 * phased, markers);
        Wrapper->SetHaps(haplotypes, 0, 0, nullptr, 0, 0, thetas, phased);
        if (loadGraph.find(".DAG") != std::string::npos)
            Wrapper->Graph.ReadDAG(loadGraph);
        else if (loadGraph.find(".json") != std::string::npos)
            Wrapper->Graph.FromJson(loadGraph);
        else
            warning("unknown input graph format!");
        Wrapper->PrintSummary();
        //for debug
//        loadGraph="reference.panel.DAG";
//        Wrapper->Graph.WriteDAG(loadGraph);
//        Wrapper->ResetWrapper();
//        Wrapper->Graph.ReadDAG(loadGraph);
//
//        loadGraph="reference.panel.DAG2";
//        Wrapper->Graph.WriteDAG(loadGraph);
        //for debug
        clock_t t1 = clock();
        fprintf(stderr, "Done loading graph in time:%.2f sec\n\n", (float) (t1 - t) / CLOCKS_PER_SEC);
    } else {//construct graph, all are phased samples
        clock_t t = clock();
        if (Wrapper != nullptr) {
            delete Wrapper;
            Wrapper = nullptr;
        }
//        Wrapper = new PBWTWrapper(2 * phased, markers, PvalueMatrix, prefixLength);
//        Wrapper->SetHaps(haplotypes, 2 * (individuals - phased), 2 * individuals, nullptr, 0, 0, thetas);
        Wrapper = new PBWTWrapper(2 * individuals /*+ (individuals - phased) * nSampleCopy*/, markers, PvalueMatrix,
                                  prefixLength);
        Wrapper->SetHaps(haplotypes, 0, 2 * individuals, 0, 0, 0, thetas, phased);
        Wrapper->CursorBackwards();//calculate backwards order of suffix
        Wrapper->CursorForwards();
        Wrapper->Graph.WriteDAG(graphFilePrefix + ".DAG");
//        Wrapper->Graph.ToJson(graphFilePrefix + ".json");
        clock_t t1 = clock();
        printf("Done building graph in time:%.2f sec\n", (float) (t1 - t) / CLOCKS_PER_SEC);
    }
}

bool PBWTHaplotyper::ReverseInput() {//ensure this function only called by ITERATIVE
    int begin = 0;
    int end = markers - 1;
    for (; begin < end; ++begin, --end) {
        int nTarget = runningModel & PHASE ? GetUnphasedNum() : individuals;
        for (int i = 0; i < nTarget; ++i) {
            //haplotypes
            std::swap(haplotypes[i * 2][begin], haplotypes[i * 2][end]);
            std::swap(haplotypes[i * 2 + 1][begin], haplotypes[i * 2 + 1][end]);
        }
        for (int i = 0; i < GetUnphasedNum(); ++i) {
            //genotypes
            std::swap(genotypes[i][begin * 3], genotypes[i][end * 3]);
            std::swap(genotypes[i][begin * 3 + 1], genotypes[i][end * 3 + 1]);
            std::swap(genotypes[i][begin * 3 + 2], genotypes[i][end * 3 + 2]);
        }
        for (int i = 0; i < GetUnphasedNum() * nSampleCopy; ++i) {
            //haplotypes
            std::swap(sampledHaps[i * 2][begin], sampledHaps[i * 2][end]);
            std::swap(sampledHaps[i * 2 + 1][begin], sampledHaps[i * 2 + 1][end]);
        }
        //penetrance
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                std::swap(Penetrance(begin, j, k), Penetrance(end, j, k));
            }
        }
        std::swap(error_models[begin], error_models[end]);
        std::swap(freq1s[begin], freq1s[end]);
    }

    begin = 0;
    end = markers - 2;
    for (; begin < end; ++begin, --end) {
        std::swap(thetas[begin], thetas[end]);
    }
    return true;
}

void
PBWTHaplotyper::ImputeAlleles(int marker, int state1, int state2, Random *rand, int currentIndividual, char **haps) {


    int currentHap1 = 2 * currentIndividual;
    int currentHap2 = currentHap1 + 1;

    int copied1 = GetAllele(marker, state1);//Wrapper->clusterAllele[marker][state1];//haplotypes[state1][marker];
    int copied2 = GetAllele(marker, state2); //Wrapper->clusterAllele[marker][state2];//haplotypes[state2][marker];
//    fprintf(stdout,"marker %d copied genotype: %d|%d\n",marker,copied1,copied2);

    int markerIndex = marker * 3;

    int ph11 = (unsigned char) genotypes[currentIndividual][markerIndex];
    int ph12 = (unsigned char) genotypes[currentIndividual][markerIndex + 1];
    int ph22 = (unsigned char) genotypes[currentIndividual][markerIndex + 2];

    CalculatePhred2Prob();

    double posterior_11 = Penetrance(marker, copied1 + copied2, 0) * phred2prob[ph11];
    double posterior_12 = Penetrance(marker, copied1 + copied2, 1) * phred2prob[ph12];
    double posterior_22 = Penetrance(marker, copied1 + copied2, 2) * phred2prob[ph22];
    double sum = posterior_11 + posterior_12 + posterior_22;

    posterior_11 /= sum;
    posterior_22 /= sum;
    if (sum == 0)
        printf("Problem!\n");

    double r = rand->Next();
//    fprintf(stdout, "[Choice]individual: %d,Marker:%d copied (%d,%d) and orginal gl:(%d,%d,%d) and posterior:(%g,%g,%g)\tChoice:%g\n",
//            currentIndividual, marker, copied1, copied2, ph11, ph12, ph22, posterior_11, posterior_12, posterior_22, r);
    if (r < posterior_11)//homo ref alleles
    {
        if (copied1 != copied2 || copied1 != 0)
            fprintf(stdout,
                    "[posterior_11]individual: %d,Homo ref Marker:%d from else (%d,%d) and orginal gl:(%d,%d,%d), posterior:(%g,%g,%g),choice:%g\n",
                    currentIndividual, marker, copied1, copied2, ph11, ph12, ph22, posterior_11, posterior_12,
                    posterior_22, r);

        haps[currentHap1][marker] = 0;
        haps[currentHap2][marker] = 0;
    } else if (r < posterior_11 + posterior_22)//home alt alleles
    {
        if (copied1 != copied2 || copied1 != 1)
            fprintf(stdout,
                    "[posterior_22]individual: %d,Homo alt Marker:%d from else (%d,%d) and orginal gl:(%d,%d,%d), posterior:(%g,%g,%g),choice:%g\n",
                    currentIndividual, marker, copied1, copied2, ph11, ph12, ph22, posterior_11, posterior_12,
                    posterior_22, r);
        haps[currentHap1][marker] = 1;
        haps[currentHap2][marker] = 1;
    } else if (copied1 != copied2)//heter states and heter alleles
    {
        double rate = GetErrorRate(marker);

        if (rand->Next() < rate * rate / ((rate * rate) + (1.0 - rate) * (1.0 - rate)))//if both alleles mutated
        {
            copied1 = !copied1;
            copied2 = !copied2;
        }

        haps[currentHap1][marker] = copied1;
        haps[currentHap2][marker] = copied2;
    } else//heter alleles but copied homo states
    {
        fprintf(stdout,
                "[posterior_12]individual: %d, Heter Marker:%d from else (%d,%d) and orginal gl:(%d,%d,%d), posterior:(%g,%g,%g),choice:%g\n",
                currentIndividual, marker, copied1, copied2, ph11, ph12, ph22, posterior_11, posterior_12, posterior_22,
                r);
        bool bit = rand->Binary();
        haps[currentHap1][marker] = bit;
        haps[currentHap2][marker] = bit ^ 1;
    }

    int imputed1 = haps[currentHap1][marker];
    int imputed2 = haps[currentHap2][marker];
    //fprintf(stdout,"imputed genotype: %d|%d\n",imputed1,imputed2);
    //int differences = abs(copied1 - imputed1) + abs(copied2 - imputed2);
    int differences = abs(copied1 + copied2 - imputed1 - imputed2);

    error_models[marker].matches += 2 - differences;
    error_models[marker].mismatches += differences;
}

void
PBWTHaplotyper::ImputeAllelesRaw(int marker, int state1, int state2, Random *rand, int currentIndividual, char **haps) {


    int currentHap1 = 2 * currentIndividual;
    int currentHap2 = currentHap1 + 1;

    int copied1 = GetAllele(marker, state1);
    int copied2 = GetAllele(marker, state2);

    haps[currentHap1][marker] = copied1;
    haps[currentHap2][marker] = copied2;

    int imputed1 = haps[currentHap1][marker];
    int imputed2 = haps[currentHap2][marker];

    int differences = abs(copied1 + copied2 - imputed1 - imputed2);

    error_models[marker].matches += 2 - differences;
    error_models[marker].mismatches += differences;
}

void PBWTHaplotyper::ImputeAllele(int haplotype, int marker, int state, char **haps) {
    // if (updateDiseaseScores) UpdateDiseaseScores(marker, state);

    haps[haplotype][marker] = GetAllele(marker, state);
}

void PBWTHaplotyper::FillPath(int haplotype, int fromMarker, int toMarker, int state, char **haps) {
    fromMarker++;

    while (fromMarker < toMarker)
        ImputeAllele(haplotype, fromMarker++, state, haps);
}

//HMM version two

float PBWTHaplotyper::GetRecombRate(int from, int to) {
    if (from < 0 or from >= markers - 1) {
        fprintf(stderr, "site out of range in GetRecombRate\n");
        exit(EXIT_FAILURE);
    }
    if (isRev) {
        from = Wrapper->nMarkers - 1 - from;
    }
    return Wrapper->recomRate[from];
}

void PBWTHaplotyper::RandomSetup(Random *rand) {
    if (rand == NULL)
        rand = &globalRandom;

    CalculatePhred2Prob();

    UpdateFreq();

    for (int j = 0; j < markers; j++) {
        int markerIndex = 3 * j;

        double prior_11 = freq1s[j] * freq1s[j];
        double prior_12 = 2.0 * freq1s[j] * (1.0 - freq1s[j]);
        double prior_22 = (1.0 - freq1s[j]) * (1.0 - freq1s[j]);

        for (int i = 0; i < GetUnphasedNum(); i++) {

            double posterior_11 = prior_11 * phred2prob[genotypes[i][markerIndex]];
            double posterior_12 = prior_12 * phred2prob[genotypes[i][markerIndex + 1]];
            double posterior_22 = prior_22 * phred2prob[genotypes[i][markerIndex + 2]];
            double sum = posterior_11 + posterior_12 + posterior_22;

            if (sum == 0)
                printf("Problem!\n");

            posterior_11 /= sum;
            posterior_12 /= sum;
            //experiment
//            posterior_22 /= sum;
//            SetGP(i,j,0,0,posterior_11);
//            SetGP(i,j,0,1,posterior_12);
//            SetGP(i,j,1,1,posterior_22);
            //experiment end

            double r = rand->Next();
            if (r < posterior_11) {
                haplotypes[i * 2][j] = 0;
                haplotypes[i * 2 + 1][j] = 0;
            } else if (r < posterior_11 + posterior_12) {
                bool bit = rand->Binary();

                haplotypes[i * 2][j] = bit;
                haplotypes[i * 2 + 1][j] = bit ^ 1;
            } else {
                haplotypes[i * 2][j] = 1;
                haplotypes[i * 2 + 1][j] = 1;
            }
//            int min = std::min(std::min(posterior_11, posterior_12), posterior_22);//phred score
//            if (min == posterior_11) {
//                haplotypes[i * 2][j] = 0;
//                haplotypes[i * 2 + 1][j] = 0;
//            } else if (min == posterior_12) {
//                bool bit = rand->Binary();
//
//                haplotypes[i * 2][j] = bit;
//                haplotypes[i * 2 + 1][j] = bit ^ 1;
//            } else {
//                haplotypes[i * 2][j] = 1;
//                haplotypes[i * 2 + 1][j] = 1;
//            }

        }
    }
}

//without recombination
struct EdgePair {
    int childNode1;
    int childNode2;
    int parentNode1;
    int parentNode2;
    float fwd;

    EdgePair(int a, int b, int c, int d, float e) {
        childNode1 = a;
        childNode2 = b;
        parentNode1 = c;
        parentNode2 = d;
        fwd = e;
    }
};

inline bool EdgePaircomparator(const EdgePair &lhs, const EdgePair &rhs) {
    return lhs.fwd > rhs.fwd;
}

int PBWTHaplotyper::InitialFwdValues(int sampleIndex, FwdBwdLocalParameter &localParameter) {
    localParameter.isRec.assign(markers, false);
    localParameter.fwdValueSum.assign(markers, 0.f);
    std::unordered_map<int, float> dummy;
    dummy.reserve(HASH_RESERVE);
    localParameter.fwdValueNode1Sum.assign(markers, dummy);
    localParameter.fwdValueNode2Sum.assign(markers, dummy);

    localParameter.states = GetStateNumFrom(0);//actually only 1 state

    float prior = 1.f / (localParameter.states * localParameter.states);
    float tmpFwdValue(0.f), gl(0.f);

    double posterior[3];
    CalculatePosteriorGL(0, sampleIndex, posterior);

    for (StateIndex i = 0; i < localParameter.states; ++i) {
        for (StateIndex j = 0; j < localParameter.states; ++j) {
            char allele1 = GetAllele(0, i);
            char allele2 = GetAllele(0, j);
//            gl = GetEmissionProb(0, sampleIndex, i, j);
            gl = GetGL(sampleIndex, 0, allele1, allele2);
//            gl = posterior[allele1+allele2];
            if (gl > genoThresh) {
                tmpFwdValue = prior * gl;
                if (0)
                    fprintf(stderr,
                            "normal debug in (%d,%d)prior:%g\t%g(%d,%d)\tfreq:%f\n",
                            i, j, prior, gl, allele1, allele2, freq1s[0]);
//                if(i == j) tmpFwdValue *= 0.5;
                if (tmpFwdValue < UNDERFLOW_MIN) {
//                    tmpFwdValue = UNDERFLOW_MIN;
                    continue;
                }
                localParameter.fwdValueSum[0] += tmpFwdValue;
                localParameter.FillParentsNodeVec(0, i, j, 0, 0, tmpFwdValue);
            }
        }
    }

    tmpFwdValue = 0.f;
    for (auto kv:localParameter.parentsNodeVec[0]) {
        StateIndex tmpNode1 = localParameter.GetFirst(kv.first);
        StateIndex tmpNode2 = localParameter.GetSecond(kv.first);

        //below is debug code
        if (DEBUG)
            for (int j = 0; j < localParameter.GetFwdVec(0, kv.second).size(); ++j) {
                float &tmpFwd = localParameter.GetFwdVec(0, kv.second).at(j);
                tmpFwd /= localParameter.fwdValueSum[0];
                tmpFwdValue += tmpFwd;
                int first = localParameter.GetFirst(localParameter.GetSourceVec(0, kv.second).at(j));
                int second = localParameter.GetSecond(localParameter.GetSourceVec(0, kv.second).at(j));
                fprintf(stderr,
                        "Init forward marker:%d\tfwdValueSum:%g\tpair:%lld(%d,%d)=(%d,%d)\tcurrentFwd:%g\tnormalized currentFwd:%g\tfrom:(%d,%d)\n",
                        0, localParameter.fwdValueSum[0], kv.first, tmpNode1, tmpNode2,
                        GetAllele(0, tmpNode1), GetAllele(0, tmpNode2), localParameter.GetFwdVec(0, kv.second).at(j),
                        tmpFwd, first, second);
            }
        else
            for (auto &tmpFwd: localParameter.GetFwdVec(0, kv.second)) {
                tmpFwd /= localParameter.fwdValueSum[0];
                tmpFwdValue += tmpFwd;
            }
        localParameter.fwdValueNode1Sum[0][tmpNode1] += tmpFwdValue;
        localParameter.fwdValueNode2Sum[0][tmpNode2] += tmpFwdValue;
    }
    return 0;
}

int PBWTHaplotyper::ResetFwdValues(FwdBwdLocalParameter &localParameter) {
    localParameter.fwdValueSum.clear();
    localParameter.fwdValueNode1Sum.clear();
    localParameter.fwdValueNode2Sum.clear();

    for (int k = 0; k < markers; ++k) {
//        localParameter.genuienParents[k].clear();
        localParameter.parentsNodeVec[k].clear();
    }
    return 0;
}

int PBWTHaplotyper::FindRecSite(vector<bool> &siteVec, int sampleIndex) {
    for (int i = 0; i < markers; ++i) {
        if (haplotypes[sampleIndex * 2][i] + haplotypes[sampleIndex * 2 + 1][i] == 1)
            siteVec[i] = true;
    }
    return 0;
}

int PBWTHaplotyper::ForwardAlgorithmRec(int sampleIndex, FwdBwdLocalParameter &localParameter, bool greedyMode) {
    Random *rand = &globalRandom;
    fprintf(stderr, "[%s] Process begin ...\n", __FUNCTION__);
    float prevFwdValue(0.f);
    float tmpFwdValue(0.f), lowestFwd(0.f);
    float gl(0.f), glParents(0.f);

    int fitPair = 0;
    int numCrediablePair = 0;

    std::priority_queue<EdgePair, std::vector<EdgePair>,
            std::function<bool(const EdgePair &, const EdgePair &)> > EdgePairList(EdgePaircomparator);

    char allele1(0), allele2(0);
    StateIndex childNode1(0), childNode2(0), parentNode1(0), parentNode2(0);

    float recRate(0.f);
    float baseProb(0.f);

    bool reenter(false);
    int fiberCnt(0), thinFiberCnt(0);
    bool thinFiberPair(false);
//    availablePair.ResetMarkerIndexAt(0);//TODO: use improved ibs methods maybe
//    availablePair.NextMarker();//marker 0 doesn't have parents
    for (int i = 1; i < markers; i++) {
        fiberCnt = 0;
        thinFiberCnt = 0;
        thinFiberPair = false;
        localParameter.fwdValueSum[i] = 0.f;
        numCrediablePair = 0;
        fitPair = 0;
        reenter = false;


        REENTRY:
        for (auto kv: localParameter.parentsNodeVec[i - 1]) {

            parentNode1 = localParameter.GetFirst(kv.first);//current round parents, last round children
            parentNode2 = localParameter.GetSecond(kv.first);
//            if(parentNode1 > parentNode2) continue;
            //sum over all the parents of current parents that lead to a child pair
            prevFwdValue = localParameter.GetSumValueFromContainer(localParameter.GetFwdVec(i - 1, kv.second));
            //i child of i-1 child state
            for (allele1 = 0; allele1 < 2; ++allele1) {
                childNode1 = GetChildNode(i - 1, parentNode1, allele1);
                if (childNode1 == -1) continue;
                //i child of i-1 child state
                for (allele2 = 0; allele2 < 2; ++allele2) {
                    childNode2 = GetChildNode(i - 1, parentNode2, allele2);
                    if (childNode2 == -1) continue;
//                    if (parentNode2 == parentNode1 and childNode2 > childNode1) continue;
//                    gl = GetEmissionProb(i, sampleIndex, childNode1, childNode2);
                    gl = GetGL(sampleIndex, i, allele1, allele2);
                    if (greedyMode and gl > genoThresh and i > 20)//and i % step != 0)
                    {
                        tmpFwdValue = prevFwdValue * GetTransitionProb(i - 1, parentNode1, childNode1) *
                                      GetTransitionProb(i - 1, parentNode2, childNode2) * gl;
//                        if(childNode1 == childNode2) tmpFwdValue *= 0.5;

                        if (0)//&& (i >= 9173 && i <= 9175 ))
                            fprintf(stderr,
                                    "normal debug marker %d from (%d,%d) to (%d,%d)=[%d,%d] prevFwdValue:%g\ttp1:%g(%g)\ttp2:%g(%g)\t%g(%d,%d)[%d,%d,%d]\tfreq:%f\n",
                                    i, parentNode1, parentNode2,
                                    childNode1, childNode2, GetAllele(i, childNode1), GetAllele(i, childNode2),
                                    prevFwdValue,
                                    GetTransitionProb(i - 1, parentNode1, childNode1),
                                    GetTransitionFreq(i - 1, parentNode1, childNode1),
                                    GetTransitionProb(i - 1, parentNode2, childNode2),
                                    GetTransitionFreq(i - 1, parentNode2, childNode2),
                                    gl, allele1, allele2,
                                    i, childNode1, childNode2,
                                    freq1s[i]);
//                        if(allele1 != GetAllele(i,childNode1)) continue;//TODO:reconsider
//                        if(allele2 != GetAllele(i,childNode2)) continue;
                        assert(allele1 == GetAllele(i, childNode1));
                        assert(allele2 == GetAllele(i, childNode2));

                        fiberCnt += 2;
                        int tmpThinFiberCnt = GetTransitionFreq(i - 1, parentNode1, childNode1) < 3 ? 1 :
                                              0 + GetTransitionFreq(i - 1, parentNode2, childNode2) < 3 ? 1 : 0;
                        if (tmpThinFiberCnt == 2) thinFiberPair = true;
                        thinFiberCnt += tmpThinFiberCnt;

                        if (tmpFwdValue < UNDERFLOW_MIN) {
//                            tmpFwdValue = UNDERFLOW_MIN;
                            continue;
                        }
                        //from (parentNode1, parentNode2) to (childNode1, childNode2), childNodes are present in conditional graph, but parentNode are not necessarily present
//                        if(localParameter.fwdValueSum[i] > 0 and tmpFwdValue/localParameter.fwdValueSum[i] < genoThresh) continue;
                        localParameter.fwdValueSum[i] += tmpFwdValue;
                        localParameter.FillParentsNodeVec(i, childNode1, childNode2, parentNode1, parentNode2,
                                                          tmpFwdValue);
                        fitPair++;
                    } else if (reenter)// and gl > 0)
                    {
                        localParameter.isRec[i - 1] = false;
                        tmpFwdValue = prevFwdValue *
                                      GetTransitionProb(i - 1, parentNode1, childNode1) *
                                      GetTransitionProb(i - 1, parentNode2, childNode2) *
                                      gl;//i fwdValueSum
//                        if(childNode1 == childNode2) tmpFwdValue *= 0.5;

                        if (0 && (i >= 9173 && i <= 9175))
                            fprintf(stderr,
                                    "mutate debug from (%d,%d) to (%d,%d) prevFwdValue:%g\ttp1:%g(%g)\ttp2:%g(%g)\t%g(%d,%d)[%d,%d,%d]\tfreq:%f\n",
                                    parentNode1,
                                    parentNode2, childNode1, childNode2,
                                    prevFwdValue, GetTransitionProb(i - 1, parentNode1, childNode1),
                                    GetTransitionFreq(i - 1, parentNode1, childNode1),
                                    GetTransitionProb(i - 1, parentNode2, childNode2),
                                    GetTransitionFreq(i - 1, parentNode2, childNode2),
                                    gl, allele1, allele2, genotypes[sampleIndex][i * 3],
                                    genotypes[sampleIndex][i * 3 + 1], genotypes[sampleIndex][i * 3 + 2], freq1s[i]);

//                        if(localParameter.fwdValueSum[i] > 0 and tmpFwdValue/localParameter.fwdValueSum[i] < genoThresh) continue;

                        if (tmpFwdValue < UNDERFLOW_MIN) {
//                            tmpFwdValue = UNDERFLOW_MIN;
                            continue;
                        }
                        fitPair++;

                        if (numCrediablePair < 50) {
                            EdgePairList.push(
                                    EdgePair(childNode1, childNode2, parentNode1, parentNode2, tmpFwdValue));
                            numCrediablePair++;
                            if (tmpFwdValue < lowestFwd) {
                                lowestFwd = tmpFwdValue;
                            }
                        } else if (tmpFwdValue >
                                   lowestFwd)// EdgePairList full and should be added into List, pop out lowest
                        {
                            EdgePairList.pop();
                            EdgePairList.push(
                                    EdgePair(childNode1, childNode2, parentNode1, parentNode2, tmpFwdValue));
                            lowestFwd = EdgePairList.top().fwd;
                        }
                    }
                }
            }
        }
        //fix singleton
        {
            if ((thinFiberCnt * 2 >= fiberCnt or thinFiberPair) && !reenter)//if more than 25% singletons
            {
                fitPair = 0;
                localParameter.ClearParentsNodeVec(i);
            }
        }

        if (fitPair == 0 && !reenter)//process orphan nodes
        {
            fprintf(stderr, "[%s] sample %d recombined at marker %d\n", __FUNCTION__, sampleIndex, i);
            recRate = GetRecombRate(i - 1, i);//start from marker 1 but store at index 0
            localParameter.isRec[i - 1] = true;
            localParameter.states = GetStateNumFrom(i);
            for (childNode1 = 0; childNode1 < localParameter.states; ++childNode1) {//dest nodeA
                for (childNode2 = 0; childNode2 < localParameter.states; ++childNode2) {//dest nodeB
                    allele1 = GetAllele(i, childNode1);
                    allele2 = GetAllele(i, childNode2);
                    gl = GetGL(sampleIndex, i, allele1, allele2);
                    if (gl < genoThresh) {
                        continue;
                    }
                    //TODO:expand to all possible parent nodes
                    for (auto tmpParentNode1:GetParentNodes(i, childNode1)) {//each parent of nodeA
                        for (auto tmpParentNode2:GetParentNodes(i, childNode2)) {//each parent of nodeB
                            allele1 = GetAllele(i - 1, tmpParentNode1);
                            allele2 = GetAllele(i - 1, tmpParentNode2);
                            glParents = GetGL(sampleIndex, i - 1, allele1, allele2);
                            if (glParents < genoThresh) continue;//genotype of all possible parentNode pair of childNode
//                            if (parentNode2 > parentNode1 and childNode2 == childNode1) continue;
                            if (localParameter.parentsNodeVec[i - 1].find(
                                    localParameter.MakePair(tmpParentNode1, tmpParentNode2)) !=
                                localParameter.parentsNodeVec[i - 1].end()) {
                                prevFwdValue = localParameter.GetSumValueFromContainer(
                                        localParameter.GetFwdVec(i - 1, tmpParentNode1, tmpParentNode2));
                            } else
                                prevFwdValue = 0.f;

                            if (DEBUG || 0 && (i >= 9750 && i <= 9750))
                                fprintf(stderr,
                                        "normal debug from (%d,%d) to (%d,%d) prevFwdValue:%g\ttp1:%g(%g)\ttp2:%g(%g)\t%g(%d,%d)\tfreq:%f\n",
                                        tmpParentNode1,
                                        tmpParentNode2, childNode1, childNode2,
                                        prevFwdValue, GetTransitionProb(i - 1, tmpParentNode1, childNode1),
                                        GetTransitionFreq(i - 1, tmpParentNode1, childNode1),
                                        GetTransitionProb(i - 1, tmpParentNode2, childNode2),
                                        GetTransitionFreq(i - 1, tmpParentNode2, childNode2),
                                        gl, GetAllele(i, childNode1), GetAllele(i, childNode2), freq1s[i]);
                            fitPair++;

                            baseProb = GetTransitionProb(i - 1, tmpParentNode1, childNode1) *
                                       GetTransitionProb(i - 1, tmpParentNode2, childNode2) * gl;

                            tmpFwdValue = prevFwdValue * baseProb * (1 - recRate) * (1 - recRate);

                            if (localParameter.fwdValueNode1Sum[i - 1].find(tmpParentNode1) !=
                                localParameter.fwdValueNode1Sum[i - 1].end())
                                tmpFwdValue += localParameter.fwdValueNode1Sum[i - 1][tmpParentNode1] *
                                               GetHapProbAt(i - 1, tmpParentNode2) *
                                               baseProb * (1 - recRate) * recRate;

                            if (localParameter.fwdValueNode2Sum[i - 1].find(tmpParentNode2) !=
                                localParameter.fwdValueNode2Sum[i - 1].end())
                                tmpFwdValue += localParameter.fwdValueNode2Sum[i - 1][tmpParentNode2] *
                                               GetHapProbAt(i - 1, tmpParentNode1) *
                                               baseProb * (1 - recRate) * recRate;

                            tmpFwdValue += GetHapProbAt(i - 1, tmpParentNode2) *
                                           GetHapProbAt(i - 1, tmpParentNode1) * baseProb * recRate * recRate;
                            if (DEBUG)
                                fprintf(stderr,
                                        "2 normal debug from (%d,%d) to (%d,%d) tmpFwdValue:%g\ttp1:%g(%g)\ttp2:%g(%g)\t%g(%d,%d)\tfreq:%f\n",
                                        tmpParentNode1,
                                        tmpParentNode2, childNode1, childNode2,
                                        tmpFwdValue, GetTransitionProb(i - 1, tmpParentNode1, childNode1),
                                        GetTransitionFreq(i - 1, tmpParentNode1, childNode1),
                                        GetTransitionProb(i - 1, tmpParentNode2, childNode2),
                                        GetTransitionFreq(i - 1, tmpParentNode2, childNode2),
                                        gl, GetAllele(i, childNode1), GetAllele(i, childNode2), freq1s[i]);
//                            if(childNode1 == childNode2) tmpFwdValue *= 0.5;
                            if (tmpFwdValue < UNDERFLOW_MIN) {
//                                tmpFwdValue = UNDERFLOW_MIN;
                                continue;

                            }
                            if (numCrediablePair < 1024) {
                                EdgePairList.push(
                                        EdgePair(childNode1, childNode2, tmpParentNode1, tmpParentNode2, tmpFwdValue));
                                numCrediablePair++;
                                if (tmpFwdValue < lowestFwd) {
                                    lowestFwd = tmpFwdValue;
                                }
                            } else if (tmpFwdValue >
                                       lowestFwd)// EdgePairList full and should be added into List, pop out lowest
                            {
                                EdgePairList.pop();
                                EdgePairList.push(
                                        EdgePair(childNode1, childNode2, tmpParentNode1, tmpParentNode2, tmpFwdValue));
                                lowestFwd = EdgePairList.top().fwd;
                            }
                        }
                    }
                }
            }
        }

        if (fitPair == 0 && !reenter) {
            //fprintf(stderr,"fatal error at marker %d\n",i);
            //exit(EXIT_FAILURE);
            fprintf(stderr, "[%s] sample %d genotype mutated at marker %d\n", __FUNCTION__, sampleIndex, i);
            reenter = true;
            goto REENTRY;
        }
        if (fitPair == 0) {
            fprintf(stderr, "[%s] sample %d mutate and recombine failed at marker %d\n", __FUNCTION__, sampleIndex, i);
            exit(EXIT_FAILURE);
        }


        while (not EdgePairList.empty()) {
            EdgePair tmpEdgePair = EdgePairList.top();
            EdgePairList.pop();
            if (tmpEdgePair.fwd > 0.f) {

//                if(localParameter.fwdValueSum[i] > 0 and tmpFwdValue/localParameter.fwdValueSum[i] < genoThresh) continue;
                localParameter.fwdValueSum[i] += tmpEdgePair.fwd;
                localParameter.FillParentsNodeVec(i, tmpEdgePair.childNode1, tmpEdgePair.childNode2,
                                                  tmpEdgePair.parentNode1, tmpEdgePair.parentNode2, tmpEdgePair.fwd);
            }
        }

        for (auto &kv: localParameter.parentsNodeVec[i]) {
            StateIndex tmpNode1 = localParameter.GetFirst(kv.first);
            StateIndex tmpNode2 = localParameter.GetSecond(kv.first);
            float tmp(0.f), tmp2(0.f);

//below is debug code
            if (DEBUG)
                for (int j = 0; j < localParameter.GetFwdVec(i, kv.second).size(); ++j) {
                    float &tmpFwd = localParameter.GetFwdVec(i, kv.second).at(j);
                    float notNormalized = tmpFwd;
                    tmpFwd /= localParameter.fwdValueSum[i];
                    tmp += tmpFwd;
                    int first = localParameter.GetFirst(localParameter.GetSourceVec(i, kv.second).at(j));
                    int second = localParameter.GetSecond(localParameter.GetSourceVec(i, kv.second).at(j));
                    fprintf(stderr,
                            "forward marker:%d\tfwdValueSum:%g\tpair:(%d,%d)=(%d,%d)\tcurrentFwd:%g\tnormalized currentFwd:%g\tfrom:(%d,%d)\n",
                            i, localParameter.fwdValueSum[i], tmpNode1, tmpNode2,
                            GetAllele(i, tmpNode1), GetAllele(i, tmpNode2),
                            notNormalized, tmpFwd, first, second);
                }
            else
                for (auto &tmpFwd: localParameter.GetFwdVec(i, kv.second)) {
                    tmpFwd /= localParameter.fwdValueSum[i];
                    tmp += tmpFwd;
                }

            localParameter.fwdValueNode1Sum[i][tmpNode1] += tmp;
            localParameter.fwdValueNode2Sum[i][tmpNode2] += tmp;
        }
//        availablePair.NextMarker();//last one is empty
    }
//    exit(EXIT_FAILURE);
    fprintf(stderr, "[%s] Process end ...\n", __FUNCTION__);
    return 0;
}

int PBWTHaplotyper::ForwardAlgorithmRestrict(int sampleIndex, FwdBwdLocalParameter &localParameter) {
    fprintf(stderr, "[%s] Process begin ...\n", __FUNCTION__);
    float prevFwdValue(0.f);
    float tmpFwdValue(0.f), lowestFwd(0.f);
    float gl(0.f), glParents(0.f);

    int fitPair = 0;
    int numCrediablePair = 0;

    std::priority_queue<EdgePair, std::vector<EdgePair>,
            std::function<bool(const EdgePair &, const EdgePair &)> > EdgePairList(EdgePaircomparator);

    char allele1(0), allele2(0);
    StateIndex childNode1(0), childNode2(0), parentNode1(0), parentNode2(0);

    float recRate(0.f);
    float baseProb(0.f);


    for (int i = 1; i < markers; i++) {
        localParameter.fwdValueSum[i] = 0.f;
        numCrediablePair = 0;
//        fprintf(stderr, "[%s] sample %d recombined at marker %d\n", __FUNCTION__, sampleIndex, i);
        recRate = GetRecombRate(i - 1, i);//start from marker 1 but store at index 0
        localParameter.isRec[i - 1] = true;
        for (auto nodePair:localParameter.GetStatePair(i)) {
            childNode1 = localParameter.GetFirst(nodePair);
            childNode2 = localParameter.GetSecond(nodePair);
            allele1 = GetAllele(i, childNode1);
            allele2 = GetAllele(i, childNode2);
            gl = GetGL(sampleIndex, i, allele1, allele2);
            if (gl < genoThresh) {
                continue;
            }
            //TODO:expand to all possible parent nodes
            for (auto tmpParentNode1:GetParentNodes(i, childNode1)) {//each parent of nodeA
                for (auto tmpParentNode2:GetParentNodes(i, childNode2)) {//each parent of nodeB
                    allele1 = GetAllele(i - 1, tmpParentNode1);
                    allele2 = GetAllele(i - 1, tmpParentNode2);
                    glParents = GetGL(sampleIndex, i - 1, allele1, allele2);
                    if (glParents < genoThresh) continue;//genotype of all possible parentNode pair of childNode
//                            if (parentNode2 > parentNode1 and childNode2 == childNode1) continue;
                    if (localParameter.parentsNodeVec[i - 1].find(
                            localParameter.MakePair(tmpParentNode1, tmpParentNode2)) !=
                        localParameter.parentsNodeVec[i - 1].end()) {
                        prevFwdValue = localParameter.GetSumValueFromContainer(
                                localParameter.GetFwdVec(i - 1, tmpParentNode1, tmpParentNode2));
                    } else
                        prevFwdValue = 0.f;

                    if (DEBUG || 0 && (i >= 9750 && i <= 9750))
                        fprintf(stderr,
                                "normal debug from (%d,%d) to (%d,%d) prevFwdValue:%g\ttp1:%g(%g)\ttp2:%g(%g)\t%g(%d,%d)\tfreq:%f\n",
                                tmpParentNode1,
                                tmpParentNode2, childNode1, childNode2,
                                prevFwdValue, GetTransitionProb(i - 1, tmpParentNode1, childNode1),
                                GetTransitionFreq(i - 1, tmpParentNode1, childNode1),
                                GetTransitionProb(i - 1, tmpParentNode2, childNode2),
                                GetTransitionFreq(i - 1, tmpParentNode2, childNode2),
                                gl, GetAllele(i, childNode1), GetAllele(i, childNode2), freq1s[i]);
                    fitPair++;

                    baseProb = GetTransitionProb(i - 1, tmpParentNode1, childNode1) *
                               GetTransitionProb(i - 1, tmpParentNode2, childNode2) * gl;

                    tmpFwdValue = prevFwdValue * baseProb * (1 - recRate) * (1 - recRate);

                    if (localParameter.fwdValueNode1Sum[i - 1].find(tmpParentNode1) !=
                        localParameter.fwdValueNode1Sum[i - 1].end())
                        tmpFwdValue += localParameter.fwdValueNode1Sum[i - 1][tmpParentNode1] *
                                       GetHapProbAt(i - 1, tmpParentNode2) *
                                       baseProb * (1 - recRate) * recRate;

                    if (localParameter.fwdValueNode2Sum[i - 1].find(tmpParentNode2) !=
                        localParameter.fwdValueNode2Sum[i - 1].end())
                        tmpFwdValue += localParameter.fwdValueNode2Sum[i - 1][tmpParentNode2] *
                                       GetHapProbAt(i - 1, tmpParentNode1) *
                                       baseProb * (1 - recRate) * recRate;

                    tmpFwdValue += GetHapProbAt(i - 1, tmpParentNode2) *
                                   GetHapProbAt(i - 1, tmpParentNode1) * baseProb * recRate * recRate;
                    if (DEBUG)
                        fprintf(stderr,
                                "2 normal debug from (%d,%d) to (%d,%d) tmpFwdValue:%g\ttp1:%g(%g)\ttp2:%g(%g)\t%g(%d,%d)\tfreq:%f\n",
                                tmpParentNode1,
                                tmpParentNode2, childNode1, childNode2,
                                tmpFwdValue, GetTransitionProb(i - 1, tmpParentNode1, childNode1),
                                GetTransitionFreq(i - 1, tmpParentNode1, childNode1),
                                GetTransitionProb(i - 1, tmpParentNode2, childNode2),
                                GetTransitionFreq(i - 1, tmpParentNode2, childNode2),
                                gl, GetAllele(i, childNode1), GetAllele(i, childNode2), freq1s[i]);
//                            if(childNode1 == childNode2) tmpFwdValue *= 0.5;
                    if (tmpFwdValue < UNDERFLOW_MIN) {
//                                tmpFwdValue = UNDERFLOW_MIN;
                        continue;

                    }
                    if (numCrediablePair < 50 * 1024) {
                        EdgePairList.push(
                                EdgePair(childNode1, childNode2, tmpParentNode1, tmpParentNode2, tmpFwdValue));
                        numCrediablePair++;
                        if (tmpFwdValue < lowestFwd) {
                            lowestFwd = tmpFwdValue;
                        }
                    } else if (tmpFwdValue >
                               lowestFwd)// EdgePairList full and should be added into List, pop out lowest
                    {
                        EdgePairList.pop();
                        EdgePairList.push(
                                EdgePair(childNode1, childNode2, tmpParentNode1, tmpParentNode2, tmpFwdValue));
                        lowestFwd = EdgePairList.top().fwd;
                    }
                }
            }
        }

        if (fitPair == 0) {
            fprintf(stderr, "[%s] sample %d mutate and recombine failed at marker %d\n", __FUNCTION__, sampleIndex, i);
            exit(EXIT_FAILURE);
        }


        while (not EdgePairList.empty()) {
            EdgePair tmpEdgePair = EdgePairList.top();
            EdgePairList.pop();
            if (tmpEdgePair.fwd > 0.f) {

//                if(localParameter.fwdValueSum[i] > 0 and tmpFwdValue/localParameter.fwdValueSum[i] < genoThresh) continue;
                localParameter.fwdValueSum[i] += tmpEdgePair.fwd;
                localParameter.FillParentsNodeVec(i, tmpEdgePair.childNode1, tmpEdgePair.childNode2,
                                                  tmpEdgePair.parentNode1, tmpEdgePair.parentNode2, tmpEdgePair.fwd);
            }
        }

        for (auto &kv: localParameter.parentsNodeVec[i]) {
            StateIndex tmpNode1 = localParameter.GetFirst(kv.first);
            StateIndex tmpNode2 = localParameter.GetSecond(kv.first);
            float tmp(0.f), tmp2(0.f);

//below is debug code
            if (DEBUG)
                for (int j = 0; j < localParameter.GetFwdVec(i, kv.second).size(); ++j) {
                    float &tmpFwd = localParameter.GetFwdVec(i, kv.second).at(j);
                    float notNormalized = tmpFwd;
                    tmpFwd /= localParameter.fwdValueSum[i];
                    tmp += tmpFwd;
                    int first = localParameter.GetFirst(localParameter.GetSourceVec(i, kv.second).at(j));
                    int second = localParameter.GetSecond(localParameter.GetSourceVec(i, kv.second).at(j));
                    fprintf(stderr,
                            "forward marker:%d\tfwdValueSum:%g\tpair:(%d,%d)=(%d,%d)\tcurrentFwd:%g\tnormalized currentFwd:%g\tfrom:(%d,%d)\n",
                            i, localParameter.fwdValueSum[i], tmpNode1, tmpNode2,
                            GetAllele(i, tmpNode1), GetAllele(i, tmpNode2),
                            notNormalized, tmpFwd, first, second);
                }
            else
                for (auto &tmpFwd: localParameter.GetFwdVec(i, kv.second)) {
                    tmpFwd /= localParameter.fwdValueSum[i];
                    tmp += tmpFwd;
                }

            localParameter.fwdValueNode1Sum[i][tmpNode1] += tmp;
            localParameter.fwdValueNode2Sum[i][tmpNode2] += tmp;
        }
    }
    fprintf(stderr, "[%s] Process end ...\n", __FUNCTION__);
    return 0;
}

int PBWTHaplotyper::BackwardSamplingRec(Random *rand, int sampleIndex, char **sampledHaps,
                                        FwdBwdLocalParameter &localParameter) {

    double choice(0.);
    double sum(0.), subSum(0.);
    float gl(0.f);
    float baseProb(0.f), transProb0(0.f), transProb1(0.f);
    float np1(0.f), np2(0.f);
    float edgeProb0(0.f), edgeProb1(0.f);
    StateIndex sampledParent0(0), sampledChild0(0);
    StateIndex sampledParent1(0), sampledChild1(0);
    StateIndex tmpSampledParent0(0), tmpSampledParent1(0);
    double sampledFwd(0.f);//not normalized fwd value
    float recRate(0.f);
    int sampled = 0;
    unsigned long long sampledChildPair(0), sampledParentPair(0);//state pair holder

    choice = rand->Uniform(0, 1);
    //print last marker sates for debug
    int count = 0;
    for (auto iter = localParameter.parentsNodeVec[markers - 1].begin();
         iter != localParameter.parentsNodeVec[markers - 1].end(); ++iter) {
        count++;
        sampledFwd = localParameter.GetSumValueFromContainer(localParameter.GetFwdVec(markers - 1, iter->second));
        if (DEBUG)
            fprintf(stderr, "(%d,%d)\tvalue:%g\tcount:%d\n", localParameter.GetFirst(iter->first),
                    localParameter.GetSecond(iter->first), sampledFwd, count);

    }

    //sample last marker
    for (auto kv: localParameter.parentsNodeVec[markers - 1]) {
        int destIndex = kv.second;//child nodes pair index
        sampledChildPair = kv.first;//child nodes pair
        int sourceVecSize = localParameter.GetSourceVec(markers - 1, destIndex).size();
        for (int sourceIndex = 0; sourceIndex < sourceVecSize; ++sourceIndex) {//actually sampling edges
            sum += localParameter.GetFwd(markers - 1, destIndex, sourceIndex);
            if (DEBUG)
                fprintf(stderr, "sampling last marker:(%d,%d) with sum:%g, recRate:%g\n",
                        localParameter.GetFirst(kv.first), localParameter.GetSecond(kv.first), sum, recRate);
            if (sum > choice) {
                sampledChild0 = localParameter.GetFirst(sampledChildPair);
                sampledChild1 = localParameter.GetSecond(sampledChildPair);
                sampledParentPair = localParameter.GetSource(markers - 1, destIndex, sourceIndex);
                sampledParent0 = localParameter.GetFirst(sampledParentPair);
                sampledParent1 = localParameter.GetSecond(sampledParentPair);
                edgeProb0 = GetEdgeProbAt(markers - 2, sampledParent0, GetAllele(markers - 1, sampledChild0));
                edgeProb1 = GetEdgeProbAt(markers - 2, sampledParent1, GetAllele(markers - 1, sampledChild1));
                sampledFwd = localParameter.GetFwd(markers - 1, destIndex, sourceIndex) *
                             localParameter.fwdValueSum[markers - 1];
                goto INIT_SAMPLE_BREAK;
            }
        }
    }
    INIT_SAMPLE_BREAK:

    ImputeAlleles(markers - 1, sampledChild0, sampledChild1, rand, sampleIndex, sampledHaps);//fill in last marker

    float testSum(0.f);

    for (int i = markers - 2; i > 0; --i) {

//        if (i == 4102 || i == 4103 || i == 4104) {
//        fprintf(stderr,
//                "marker:%d\ttestSum:%g\tsum:%g\tchoice:%g\tprev(%d,%d) to (%d,%d)\tsampledFwd:%g\tavailablePair.Size:%d\trecRate:%g\tgl:%g\tprobEdge0:%g\tprobEdge1:%g\tfwdValueSum:%g\n",
//                i, testSum, sum,
//                choice, sampledParent0, sampledParent1, sampledChild0, sampledChild1,  sampledFwd,
//                availablePair.Size(),
//                recRate, gl0, edgeProb0, edgeProb1, fwdValueSum[i]);

        gl = GetGL(sampleIndex, i + 1, GetAllele(i + 1, sampledChild0), GetAllele(i + 1, sampledChild1));
//        gl = GetEmissionProb(i + 1, sampleIndex, sampledChild0, sampledChild1);
        if (DEBUG)
            fprintf(stderr, "marker:%d\tsum:%g\tchoice:%g\t(%d,%d)\tvalue:%g\tgl:%g[%d,%d](%g,%g,%g)\tstateNum:%d\n", i,
                    sum, choice,
                    sampledParent0,
                    sampledParent1, sampledFwd, gl, GetAllele(i + 1, sampledChild0), GetAllele(i + 1, sampledChild1),
                    GetGL(sampleIndex, i + 1, 0, 0),
                    GetGL(sampleIndex, i + 1, 0, 1),
                    GetGL(sampleIndex, i + 1, 1, 1),
                    GetStateNumFrom(i));

        sum = 0.f;
        testSum = 0.f;

        ImputeAlleles(i, sampledParent0, sampledParent1, rand, sampleIndex, sampledHaps);
        choice = rand->Uniform(0, sampledFwd);


        if (!localParameter.isRec[i]) {//normal hap match without rec
            baseProb = GetTransitionProb(i, sampledParent0, sampledChild0) *
                       GetTransitionProb(i, sampledParent1, sampledChild1) * gl;
            //we are actually sampling i-1's states
            int destIndex = localParameter.GetDestIndex(i, sampledParent0, sampledParent1);
            int sourceVecSize = localParameter.GetSourceVec(i, sampledParent0, sampledParent1).size();
            for (int grandParentsIndex = 0; grandParentsIndex < sourceVecSize; ++grandParentsIndex) {
                sum += baseProb * localParameter.GetFwd(i, destIndex, grandParentsIndex);//kv.second is normalized
                if (0 && i == 9174) {
                    fprintf(stderr,
                            "sampling marker:%d\tfrom:(%d,%d)\tto pair:(%d,%d)=(%d,%d)\tFwd:%g\ttp1:%g\ttp2:%g\tgl:%gsourceVecSize:%d\n",
                            i, sampledParent0, sampledParent1, sampledChild0, sampledChild1,
                            GetAllele(i + 1, sampledChild0), GetAllele(i + 1, sampledChild1),
                            localParameter.GetFwd(i, destIndex, grandParentsIndex),
                            GetTransitionProb(i, sampledParent0, sampledChild0),
                            GetTransitionProb(i, sampledParent1, sampledChild1), gl, sourceVecSize);
                }
                if (sum > choice) {
                    sampledChild0 = sampledParent0;
                    sampledChild1 = sampledParent1;
                    sampledParentPair = localParameter.GetSource(i, destIndex, grandParentsIndex);
                    sampledParent0 = localParameter.GetFirst(sampledParentPair);
                    sampledParent1 = localParameter.GetSecond(sampledParentPair);
                    edgeProb0 = GetEdgeProbAt(i - 1, sampledParent0, GetAllele(i, sampledChild0));
                    edgeProb1 = GetEdgeProbAt(i - 1, sampledParent1, GetAllele(i, sampledChild1));
                    sampledFwd = localParameter.GetFwd(i, destIndex, grandParentsIndex) * localParameter.fwdValueSum[i];
                    break;
                }
            }
        } else {
            baseProb = edgeProb0 * edgeProb1 * gl;//site i
            recRate = GetRecombRate(i, i + 1);//start from marker 1 but store at index 0
            for (auto kv:localParameter.parentsNodeVec[i]) {// node pair at site i, sampled (sampledParent0,sampledParent1) for i
                uint64_t parentNodePair = kv.first;
                np1 = GetHapProbAt(i, localParameter.GetFirst(parentNodePair));
                np2 = GetHapProbAt(i, localParameter.GetSecond(parentNodePair));
                subSum = 0.;

                tmpSampledParent0 = localParameter.GetFirst(parentNodePair);
                tmpSampledParent1 = localParameter.GetSecond(parentNodePair);

                bool noJump0 = tmpSampledParent0 == sampledParent0;
                bool noJump1 = tmpSampledParent1 == sampledParent1;
                if (noJump0 && noJump1)//no recomb
                    subSum += (1 - recRate) * (1 - recRate) * baseProb / (np1 * np2);
                if (noJump0)//noJump for first
                    subSum += (1 - recRate) * recRate * baseProb / np1;
                if (noJump1)//noJump for second
                    subSum += (1 - recRate) * recRate * baseProb / np2;
                subSum += recRate * recRate * baseProb;

                int destIndex = kv.second;
                int sourceVecSize = localParameter.GetSourceVec(i, destIndex).size();
                //all grandParents ending at parentNodePair.first and parentNodePair.second
                for (int grandParentsIndex = 0; grandParentsIndex < sourceVecSize; ++grandParentsIndex) {
                    sum += subSum * localParameter.GetFwd(i, destIndex, grandParentsIndex);
                    //GetFwd is the fwd value at site i-1, contributed by grandParents at site i-2
                    if (sum > choice) {
                        sampledChild0 = tmpSampledParent0;//previous correct parents rejuvenate
                        sampledChild1 = tmpSampledParent1;
                        //override no-rec Imputation
                        ImputeAlleles(i, sampledChild0, sampledChild1, rand, sampleIndex, sampledHaps);
                        parentNodePair = localParameter.GetSource(i, destIndex, grandParentsIndex);
                        sampledParent0 = localParameter.GetFirst(parentNodePair);
                        sampledParent1 = localParameter.GetSecond(parentNodePair);
                        edgeProb0 = GetEdgeProbAt(i - 1, sampledParent0, GetAllele(i, sampledChild0));
                        edgeProb1 = GetEdgeProbAt(i - 1, sampledParent1, GetAllele(i, sampledChild1));
                        sampledFwd =
                                localParameter.GetFwd(i, destIndex, grandParentsIndex) * localParameter.fwdValueSum[i];
                        goto SAMPLE_BREAK;
                    }
                }
            }
            SAMPLE_BREAK:;
        }
        if (sum < choice) {
            fprintf(stderr,
                    "inside marker:%d\ttestSum:%g\tsum:%g\tchoice:%g\t(%d,%d)\ttotalValue:%g\tsampled:%d\tsubSum:%g\trecRate:%g\tbaseProb:%g\tfwdValueSum:%g\n",
                    i, testSum, sum, choice, sampledParent0, sampledParent1, sampledFwd, sampled, subSum,
                    recRate,
                    baseProb, localParameter.fwdValueSum[i]);
            count = 0;
            for (auto iter = localParameter.parentsNodeVec[i].begin();
                 iter != localParameter.parentsNodeVec[i].end(); ++iter) {
                count++;
                sampledFwd = localParameter.GetSumValueFromContainer(localParameter.GetFwdVec(i - 1, iter->second));
                if (DEBUG)
                    fprintf(stderr, "(%d,%d)\tvalue:%g\tcount:%d\n", localParameter.GetFirst(iter->first),
                            localParameter.GetSecond(iter->first), sampledFwd, count);

            }
            throw samplingException;
        }
    }

    ImputeAlleles(0, sampledParent0, sampledParent1, rand, sampleIndex, sampledHaps);
    return 0;
}


//auxilary functions for forward-backward method
float PBWTHaplotyper::MatureBwdValue(FwdBwdLocalParameter &localParameter, int site) {
    float recRate = GetRecombRate(site, site + 1);
    float sumBaseBwd = 0.f;
    for (auto kv:localParameter.parentsNodeVec[site]) {
        uint64_t sampledChildPair = kv.first;
        StateIndex sampledChild1 = localParameter.GetFirst(sampledChildPair);//parent of previous marker
        StateIndex sampledChild2 = localParameter.GetSecond(sampledChildPair);
        float np1 = GetHapProbAt(site, sampledChild1);
        float np2 = GetHapProbAt(site, sampledChild2);

        float bwdValue = (1 - recRate) * (1 - recRate) * localParameter.bwdValue[sampledChildPair] / (np1 * np2);
        bwdValue += (1 - recRate) * recRate * localParameter.bwdValueNode1Sum[sampledChild1] / np1;
        bwdValue += (1 - recRate) * recRate * localParameter.bwdValueNode2Sum[sampledChild2] / np2;
        bwdValue += recRate * recRate * localParameter.bwdValueSum;

        localParameter.finalBwdValue[sampledChildPair] = bwdValue;//marker i, same as sampleChild
        sumBaseBwd += bwdValue;
    }
    return sumBaseBwd;
}

void PBWTHaplotyper::SetNextBwdValue(FwdBwdLocalParameter &localParameter, float gl, uchar allele1, uchar allele2,
                                     int destIndex, float bwdValue, int site) {
    for (int sourceIndex = 0; sourceIndex < localParameter.GetSourceVec(site, destIndex).size(); ++sourceIndex) {
        uint64_t sampledParentPair = localParameter.GetSource(site, destIndex, sourceIndex);
        StateIndex sampledParent1 = localParameter.GetFirst(sampledParentPair);
        StateIndex sampledParent2 = localParameter.GetSecond(sampledParentPair);

        float edgeProb1 = GetEdgeProbAt(site - 1, sampledParent1, allele1);
        float edgeProb2 = GetEdgeProbAt(site - 1, sampledParent2, allele2);

        float nextBaseBwd = bwdValue * edgeProb1 * edgeProb2 * gl;

        localParameter.bwdValue[sampledParentPair] += nextBaseBwd;//markers - 2 tmp bwdValue
        localParameter.bwdValueNode1Sum[sampledParent1] += nextBaseBwd;
        localParameter.bwdValueNode2Sum[sampledParent2] += nextBaseBwd;
        localParameter.bwdValueSum += nextBaseBwd;
    }
}

void
PBWTHaplotyper::ProcessFwdBwd(int sampleIndex, FwdBwdLocalParameter &localParameter, float sumBaseBwd,
                              float *tmpGL, int site) {
    for (auto kv:localParameter.parentsNodeVec[site]) {
        uint64_t sampledChildPair = kv.first;
        int destIndex = kv.second;
        StateIndex sampledChild1 = localParameter.GetFirst(sampledChildPair);
        StateIndex sampledChild2 = localParameter.GetSecond(sampledChildPair);
        char allele1 = GetAllele(site, sampledChild1);
        char allele2 = GetAllele(site, sampledChild2);
        float gl = GetGL(sampleIndex, site, allele1, allele2);
        float fwdValue = localParameter.GetSumFwdValueFrom(site, sampledChild1, sampledChild2);//TODO:double check
        localParameter.finalBwdValue[sampledChildPair] /= sumBaseBwd;
        float bwdValue = localParameter.finalBwdValue[sampledChildPair];
//        if(site == 9750)
//        {
//            fprintf(stderr,
//                    "site:%d\tsampled:(%d,%d) gl:%g(%d,%d)\tfwd:%g\tbwd:%g\n",
//                    site, sampledChild1, sampledChild2, gl, GetAllele(site, sampledChild1), GetAllele(site, sampledChild2),
//                    fwdValue, localParameter.finalBwdValue[sampledChildPair]);
//        }
        tmpGL[allele1 + allele2] += fwdValue * bwdValue;//in decimal format
        if (site != 0) SetNextBwdValue(localParameter, gl, allele1, allele2, destIndex, bwdValue, site);
    }
}

int PBWTHaplotyper::BackwardAlgorithmRec(Random *rand, int sampleIndex, char **sampledHaps,
                                         FwdBwdLocalParameter &localParameter) {

    float sumBaseBwd(0.f);//for normalization
    uint64_t sampledChildPair(0);//state pair holder
    float tmpGL[3] = {0.f, 0.f, 0.f};

    //calculate last marker
    for (auto kv:localParameter.parentsNodeVec[markers - 1]) {
        sampledChildPair = kv.first;
        localParameter.finalBwdValue[sampledChildPair] = 1;
        sumBaseBwd += 1;
    }
    ProcessFwdBwd(sampleIndex, localParameter, sumBaseBwd, tmpGL, markers - 1);
    FillGL(sampleIndex, markers - 1, tmpGL);
    //end of last marker
    for (int i = markers - 2; i >= 0; --i) {
        tmpGL[0] = tmpGL[1] = tmpGL[2] = 0;
        sumBaseBwd = MatureBwdValue(localParameter, i);
        localParameter.ResetBwdValue();//clear tmp bwdValues
        ProcessFwdBwd(sampleIndex, localParameter, sumBaseBwd, tmpGL, i);
        FillGL(sampleIndex, i, tmpGL);
    }
    return 0;
}

int PBWTHaplotyper::LocalForwardBackwardSampling(int sampleIndex) {
    FwdBwdLocalParameter localParameter(markers);
    try {
        InitialFwdValues(sampleIndex, localParameter);
        ForwardAlgorithmRec(sampleIndex, localParameter, true);
        BackwardSamplingRec(&globalRandom, sampleIndex, haplotypes, localParameter);
//        for (int j = 0; j < nSampleCopy; ++j) {//n copy per individual
//            BackwardSamplingRec(&globalRandom, j + sampleIndex * nSampleCopy, sampledHaps, localParameter);
//        }
    }
    catch (std::exception &e) {
        fprintf(stderr, "%dth individual phasing failed!\n", sampleIndex);
        fprintf(stderr, e.what());
    }
    return 0;
}

int PBWTHaplotyper::LocalForwardBackward(int sampleIndex) {
    FwdBwdLocalParameter localParameter(markers);
    try {
        InitialFwdValues(sampleIndex, localParameter);
        ForwardAlgorithmRec(sampleIndex, localParameter, false);
        BackwardSamplingRec(&globalRandom, sampleIndex, haplotypes, localParameter);
        BackwardAlgorithmRec(&globalRandom, sampleIndex, haplotypes, localParameter);
        CorrectGenotype();
//        localParameter.ResetParentsNodeVec();
//        FindAccessibleStates(sampleIndex, localParameter);
//        localParameter.ResetParentsNodeVec();
//        InitialFwdValues(sampleIndex, localParameter);
//        ForwardAlgorithmRestrict(sampleIndex, localParameter);
//        BackwardAlgorithmRec(&globalRandom, sampleIndex, haplotypes, localParameter);
//        CorrectGenotype();
//        for (int j = 0; j < nSampleCopy; ++j) {//n copy per individual
//            BackwardAlgorithmRec(&globalRandom, j + sampleIndex * nSampleCopy, sampledHaps, localParameter);
//        }
    }
    catch (std::exception &e) {
        fprintf(stderr, "%dth individual phasing failed!\n", sampleIndex);
        fprintf(stderr, e.what());
    }
//    exit(EXIT_FAILURE);
    return 0;
}

int PBWTHaplotyper::CorrectGenotype() {
    Random *rand = &globalRandom;
    for (int i = 0; i < individuals - phased; ++i) {
        for (int j = 0; j < markers; ++j) {
            float minGL = std::numeric_limits<float>::min();
            int geno = 4;
            for (int k = 0; k < 3; ++k) {
                if (genoProbs[i][j * 3 + k] > minGL) {
                    minGL = genoProbs[i][j * 3 + k];
                    geno = k;
                }
            }
            if (geno != haplotypes[2 * i][j] + haplotypes[2 * i + 1][j]) {
//                fprintf(stdout, "[CorrectGenotype]individual: %d, Heter Marker:%d from else (%d,%d) and gp:(%g,%g,%g)\n",
//                        i, j, haplotypes[2*i][j], haplotypes[2*i+1][j], genoProbs[i][j*3 + 0], genoProbs[i][j*3 + 1], genoProbs[i][j*3 + 2]);
                if (geno == 0) {
                    haplotypes[2 * i][j] = 0;
                    haplotypes[2 * i + 1][j] = 0;
                } else if (geno == 2) {
                    haplotypes[2 * i][j] = 1;
                    haplotypes[2 * i + 1][j] = 1;
                } else//heter
                {
                    bool bit = rand->Binary();
                    haplotypes[2 * i][j] = bit;
                    haplotypes[2 * i + 1][j] = bit ^ 1;
                }

            }
        }
    }
    return 0;
}


int PBWTHaplotyper::FindAccessibleStates(int sampleIndex, FwdBwdLocalParameter &localParameter) {
    StateIndex parentNode1(-1), parentNode2(-1), childNode1(-1), childNode2(-1);
    char allele1(-1), allele2(-1);

    int geno = haplotypes[sampleIndex * 2][0] + haplotypes[sampleIndex * 2 + 1][0];

    StateIndex numStates = GetStateNumFrom(0);//actually only 1 state

    for (StateIndex i = 0; i < numStates; ++i) {
        for (StateIndex j = 0; j < numStates; ++j) {
            char allele1 = GetAllele(0, i);
            char allele2 = GetAllele(0, j);
            if (geno == 0 or geno == 2) {
                if (allele1 + allele2 == geno) {
                    localParameter.AddStatePair(0, i, j);
                    localParameter.FillParentsNodeVec(0, i, j, 0, 0, 0);
                }
            } else {
                localParameter.AddStatePair(0, i, j);
                localParameter.FillParentsNodeVec(0, i, j, 0, 0, 0);
            }
        }
    }

    for (int i = 1; i < markers; i++) {
        int geno = haplotypes[sampleIndex * 2][i] + haplotypes[sampleIndex * 2 + 1][i];
        for (auto kv: localParameter.parentsNodeVec[i - 1]) {
            parentNode1 = localParameter.GetFirst(kv.first);//current round parents, last round children
            parentNode2 = localParameter.GetSecond(kv.first);
            //i child of i-1 child state
            for (allele1 = 0; allele1 < 2; ++allele1) {
                childNode1 = GetChildNode(i - 1, parentNode1, allele1);
                if (childNode1 == -1) continue;
                //i child of i-1 child state
                for (allele2 = 0; allele2 < 2; ++allele2) {
                    childNode2 = GetChildNode(i - 1, parentNode2, allele2);
                    if (childNode2 == -1) continue;
                    if (geno == 0 or geno == 2) {
                        if (allele1 + allele2 == geno) {
                            localParameter.AddStatePair(i, childNode1, childNode2);
                            localParameter.FillParentsNodeVec(i, childNode1, childNode2, parentNode1, parentNode2, 0);
                        }
                    } else {
                        localParameter.AddStatePair(i, childNode1, childNode2);
                        localParameter.FillParentsNodeVec(i, childNode1, childNode2, parentNode1, parentNode2, 0);
                    }
                }
            }
        }
    }

    return 0;
}

//memory management
bool PBWTHaplotyper::SetErrorAndTheta(std::vector<float> &holderError, std::vector<float> &holderTheta) {
    for (int i = 0; i < markers - 1; i++) {
        thetas[i] = holderTheta[i + 1];
        SetErrorRate(i, holderError[i]);
//	fprintf(stderr,"out of %d, marker %d:%f\n",markers,i,holderError[i]);
    }
    SetErrorRate(markers - 1, holderError[markers - 1]);
//    fprintf(stderr,"got here!\n");
    holderError.clear();
    holderTheta.clear();
    return true;
}

bool PBWTHaplotyper::AllocateMemory(int persons, int m)//for GraphSeperated, only alloc mem for genotypes when phasing
{
    individuals = persons;
    markers = m;
    if ((runningModel & PHASE) || (runningModel & ITERATIVE)) {//enter from phase by ref
        thetas = new float[markers - 1];
        for (int i = 0; i < markers - 1; i++)
            thetas[i] = 0.01;
        genotypes = AllocateCharMatrix(GetUnphasedNum(), markers * 3);
        genoProbs = AllocateFloatMatrix(GetUnphasedNum(), markers * 3);
        penetrances = new float[markers * 9];
        error_models = new Errors[markers];
        SetErrorRate(0.01);
        crossovers = new int[markers - 1];
    } else//enter from graph construct
    {
        thetas = nullptr;
        genotypes = nullptr;
        genoProbs = nullptr;
        penetrances = nullptr;
        error_models = nullptr;
        crossovers = nullptr;
    }
    freq1s = new double[markers];

    if (runningModel & PHASE)
        haplotypes = AllocateCharMatrix(GetUnphasedNum() * 2, markers);
    else
        haplotypes = AllocateCharMatrix(individuals * 2, markers);

    Wrapper = nullptr;
    marginals = nullptr;

    leftMatrices = nullptr;
    leftProbabilities = nullptr;

    memoryBlock = nullptr;
    smallMemoryBlock = nullptr;
    smallFree = 0;

    stack = nullptr;
    stackPtr = -1;


    orderedGenotypeFlags = nullptr;

    return readyForUse = true;
}
