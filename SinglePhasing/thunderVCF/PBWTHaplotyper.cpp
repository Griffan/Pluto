//
// Created by Fan Zhang on 8/6/15.
//

#include "PBWTHaplotyper.h"
#include "MemoryAllocators.h"

//debug related
#define DEBUG false
static const float UNDERFLOW_MIN = std::numeric_limits<float>::min() * 100;

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
        for (int l = 0; l < nSampleCopy * (individuals - phased) * 2; ++l) {
            delete[] sampledHaps[l];
        }
        delete[] sampledHaps;
    }

    if (Wrapper != nullptr)
        delete Wrapper;

    ReleaseMemoryBlock();
    DestroyPvalueMatrix();
}

//memory management
void PBWTHaplotyper::ReleaseMemoryBlock() {
    for (std::unordered_map<int, std::vector<float *> >::iterator iter = memoryBlockList.begin();
         iter != memoryBlockList.end(); ++iter) {
        for (size_t i = 0; i < iter->second.size(); i++) {
            if (iter->second[i] != nullptr) delete[] iter->second[i];
        }
    }
}

void PBWTHaplotyper::GetMemoryBlock(int marker) {
    if (!economyMode || marker == 0 || marker > stack[stackPtr] + gridSize) {
        stack[++stackPtr] = marker;
        leftMatrices[marker] = GetLargeBlock();

        ResetReuseablePool();
    } else
        leftMatrices[marker] = GetReuseableBlock();
}

float *PBWTHaplotyper::GetLargeBlock() {
    int blockSize = orderedGenotypes ? states * states : states * (states + 1) / 2;
    if (numInUse.find(blockSize) == numInUse.end()) {
        numInUse[blockSize] = 0;
        memoryBlockList[blockSize] = std::vector<float *>(0, nullptr);
    }
    if (numInUse[blockSize] < (int) memoryBlockList[blockSize].size()) {
        numInUse[blockSize]++;
        return memoryBlockList[blockSize][numInUse[blockSize] - 1];
    } else {
        memoryBlockList[blockSize].push_back(AllocateMemoryBlock());
        numInUse[blockSize]++;
        return memoryBlockList[blockSize][numInUse[blockSize] - 1];
    }
}

float *PBWTHaplotyper::GetReuseableBlock() {
    int blockSize = orderedGenotypes ? states * states : states * (states + 1) / 2;
    if (numInUse.find(blockSize) == numInUse.end()) {
        numInUse[blockSize] = 0;
        memoryBlockList[blockSize] = std::vector<float *>(0, nullptr);
    }
    if (numInUse[blockSize] < (int) memoryBlockList[blockSize].size()) {
        numInUse[blockSize]++;//TODO::reset needed
        return memoryBlockList[blockSize][numInUse[blockSize] - 1];
    } else {
        memoryBlockList[blockSize].push_back(AllocateMemoryBlock());
        numInUse[blockSize]++;
        return memoryBlockList[blockSize][numInUse[blockSize] - 1];
    }

}

void PBWTHaplotyper::ResetMemoryPool() {
    nextAvailable = nextSmallAvailable = 0;
    nextReuseable = markers - 1;
    stackPtr = -1;
    for (std::unordered_map<int, int>::iterator iter = numInUse.begin(); iter != numInUse.end(); ++iter) {
        iter->second = 0;
    }
}

void PBWTHaplotyper::ResetReuseablePool() {
    nextReuseable = markers - 1;
}

/*
void PBWTHaplotyper::RetrieveMemoryBlock(int marker) {
    if (stack[stackPtr] <= marker) {
//        fprintf(stderr, "%d out from RetrieveMemory\n",marker);
        return;
    } else {
        ResetReuseablePool();

        float *from = leftMatrices[stack[--stackPtr]];

        for (int i = stack[stackPtr] + 1; i <= marker; i++) {
            int markerindex = i * 3;
            {
                leftMatrices[i] = GetReuseableBlock();

                Transpose(i, from, leftMatrices[i]);
                ConditionOnData(leftMatrices[i], i, genotypes[states / 2][markerindex],
                                genotypes[states / 2][markerindex + 1], genotypes[states / 2][markerindex + 2]);

                from = leftMatrices[i];
            }
        }
    }
}
*/
bool PBWTHaplotyper::ForceMemoryAllocation() {
    // Cycle through individuals, with the exact same steps as the actual
    // haplotyper and request memory ... by requesting all memory upfront,
    // we force crashes to happen early.
    for (int i = 0; i < individuals - phased; i++) {
        ResetMemoryPool();
        GetMemoryBlock(0);

        if (leftMatrices[0] == nullptr)
            return false;

        int skipped = 0;
        for (int j = 1; j < markers; j++)
            //if (genotypes[i][j] != GENOTYPE_MISSING || j == markers - 1)
        {
            GetMemoryBlock(j);

            if (leftMatrices[j] == nullptr)
                return false;
        }
        //else
        //skipped++;

        if (skipped == 0) break;
    }

    if (!phased)
        return true;

    ResetMemoryPool();
    for (int j = 0; j < markers; j++) {
        GetSmallMemoryBlock(j);

        if (leftMatrices[j] == nullptr)
            return false;
    }

    return true;
}

//random setup
void PBWTHaplotyper::InitialSampleCopy(Random *rand) {

    if (rand == nullptr)
        rand = &globalRandom;
    CalculatePhred2Prob();

    if (nSampleCopy == 0) return;
    sampledHaps = new char *[nSampleCopy * (individuals - phased) * 2];
    for (int l = 0; l < nSampleCopy * (individuals - phased) * 2; ++l) {
        sampledHaps[l] = new char[markers];
    }

    for (int j = 0; j < markers; j++) {
        double mac = 0;
        int markerindex = 3 * j;

        double hyperprior11 = freq1s[j] * freq1s[j];
        double hyperprior12 = 2.0 * freq1s[j] * (1.0 - freq1s[j]);
        double hyperprior22 = (1.0 - freq1s[j]) * (1.0 - freq1s[j]);

        for (int i = 0; i < individuals; i++) {
            double post11 = hyperprior11 * phred2prob[(size_t) genotypes[i][markerindex]];
            double post12 = hyperprior12 * phred2prob[(size_t) genotypes[i][markerindex + 1]];
            double post22 = hyperprior22 * phred2prob[(size_t) genotypes[i][markerindex + 2]];
            double sumpost = post11 + post12 + post22;
            post11 /= sumpost;
            post12 /= sumpost;
            post22 /= sumpost;

            // estimated counts of AL2
            mac += post12 + 2 * post22;
        }

        //here, each person contributes two alleles
        double freq = 0.5 * mac / (double) individuals;

        double prior_11 = (1.0 - freq) * (1.0 - freq);
        double prior_12 = 2.0 * freq * (1.0 - freq);
        double prior_22 = freq * freq;


        for (int i = 0; i < individuals - phased; i++) {

            double posterior_11 = prior_11 * phred2prob[(size_t) genotypes[i][markerindex]];
            double posterior_12 = prior_12 * phred2prob[(size_t) genotypes[i][markerindex + 1]];
            double posterior_22 = prior_22 * phred2prob[(size_t) genotypes[i][markerindex + 2]];
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
        Swap(sampledHaps[a * 2], sampledHaps[(individuals - phased - 1) * 2]);
        Swap(sampledHaps[a * 2 + 1], sampledHaps[(individuals - phased - 1) * 2 + 1]);
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
    Wrapper->SetHaps(haplotypes, 2 * (individuals - phased), 2 * individuals, nullptr, 0, 0, thetas);
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

class SamplingException : public std::exception {
    virtual const char *what() const throw() {
        return "Current sample encountered unexpected sampling space!\n";
    }
} samplingException;

#ifdef _OPENMP

#include <omp.h>

#endif

int PBWTHaplotyper::LoopThroughChromosomesRecomb(Pedigree &ped) {

    ResetCrossovers();

    ConstructGraph();

    clock_t t1 = clock();
#ifdef _OPENMP
    omp_set_num_threads(nThread);
#pragma omp parallel for
#endif
    for (int i = individuals - 1; i >= 0; i--) {
//            SwapIndividuals(i, individuals - 1);
        fprintf(stderr, "[%s]phasing individual %d:%s...\n\n", __FUNCTION__, i, ped[i].pid.c_str());
        LocalForwadBackWard(i);
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
    if (isRev) ReverseInput();
    return 0;
}

void PBWTHaplotyper::ConstructGraph() {
    if (isRev) ReverseInput();

    if (loadGraph != "Empty") {
        clock_t t = clock();
        if (Wrapper != nullptr) {
            delete Wrapper;
            Wrapper = nullptr;
        }
        //Wrapper = new PBWTWrapper(2 * phasedForByRef, markers, PvalueMatrix, prefixLength);
        //Wrapper->ReleaseWrapperMemory();
        Wrapper = new PBWTWrapper(2 * phased, markers);
        Wrapper->SetHaps(haplotypes, 0, 0, nullptr, 0, 0, thetas);
        Wrapper->Graph.ReadContainer(loadGraph);
        //for debug
//        loadGraph="reference.panel.DAG";
//        Wrapper->Graph.WriteContainer(loadGraph);
//        Wrapper->ResetWrapper();
//        Wrapper->Graph.ReadContainer(loadGraph);
//
//        loadGraph="reference.panel.DAG2";
//        Wrapper->Graph.WriteContainer(loadGraph);
        //for debug
        clock_t t1 = clock();
        fprintf(stderr, "Done loading graph in time:%.2f sec\n\n", (float) (t1 - t) / CLOCKS_PER_SEC);
    } else {//construct graph, all are phased samples
        clock_t t = clock();
        if (Wrapper != nullptr) {
            delete Wrapper;
            Wrapper = nullptr;
        }
        Wrapper = new PBWTWrapper(2 * phased, markers, PvalueMatrix, prefixLength);
        Wrapper->SetHaps(haplotypes, 0, 2 * individuals, nullptr, 0, 0, thetas);
        Wrapper->CursorBackwards();//calculate backwards order of suffix
        Wrapper->CursorForwards();
        Wrapper->Graph.WriteContainer(outputPrefix + ".DAG");
        clock_t t1 = clock();
        printf("Done building graph in time:%.2f sec\n", (float) (t1 - t) / CLOCKS_PER_SEC);
    }
}

#ifdef HETERSITE
int PBWTHaplotyper::LoopThroughChromosomesViaPBWTWithHeterOnly() {

    ResetCrossovers();

    for (int i = individuals - 1; i >= 0; i--) {

        if (i < individuals - phased) {
            indexBeingSampled = i;
            SwapIndividuals(i, individuals - 1);
            if (isRev) ReverseInput();
            clock_t t = clock();
            ExtractHeterSites(individuals - 1);
            Wrapper->SetHaps(haplotypes, 0, 2 * individuals, sampledHaps, 0, (individuals - phased) * nSampleCopy * 2,
                             thetas);
            Wrapper->CursorBackwards();//calculate backwards order of suffix
            Wrapper->CursorForwards();

//#ifdef DEBUG
//            {
//                Wrapper->PrintHap(tmpHaps, Wrapper->a[0]);
//
//                // Wrapper->PrintHap(tmpHaps,Wrapper->a[6]);
//                Wrapper->PrintHap(tmpHaps, Wrapper->a[Wrapper->N - 1]);
//                // Wrapper->PrintMatrix(Wrapper->a,"a array matrix");
//                Wrapper->PrintMatrix(Wrapper->d, "d array");
//                //Wrapper->PrintVector(Wrapper->a[Wrapper->N-7],"last a array");
//            }
//#endif
            printf("%d markers used for individual %d\n", markers, i);
            clock_t t1 = clock();
            printf("build model time:%.2f sec\n", (float) (t1 - t) / CLOCKS_PER_SEC);
            //exit(EXIT_SUCCESS);
//            if (weights != NULL)
//                ScaleWeights();
//
//            if (updateDiseaseScores)
//                ScoreNPL();

            //if (i < individuals - phased) {
            fprintf(stderr, "phasing individual %d...\n", i);
//            ScoreLeftConditional();
            ForwardAlgorithm(0, <#initializer#>);
            t = clock();
            fprintf(stderr,"forward algorithm time:%.2f sec\n", (float) (t - t1) / CLOCKS_PER_SEC);

//            SampleChromosomes(&globalRandom);
            BackwardSampling(&globalRandom, individuals - 1, haplotypes, <#initializer#>);

            for (int j = 0; j < nSampleCopy; ++j) {//n copy per individual
                BackwardSampling(&globalRandom, j + i * nSampleCopy, sampledHaps, <#initializer#>);
            }
            t1 = clock();
            fprintf(stderr,"sampling time:%.2f sec\n", (float) (t1 - t) / CLOCKS_PER_SEC);
            //exit(EXIT_SUCCESS);

//            if (updateDiseaseScores && diseaseCount)
//                IntegrateNPL();

#ifdef _DEBUG
            if (!SanityCheck())
               {
               printf("\nProblems above occurred haplotyping individual %d\n\n", i);
               Print();
               }
#endif
            FillHeterSitesBack(individuals - 1);
            SwapIndividuals(i, individuals - 1);
            //Wrapper->ResetWrapper();
            if (isRev) ReverseInput();
        }


    }

    return 0;
}
#endif

bool PBWTHaplotyper::ReverseInput() {
    int begin = 0;
    int end = markers - 1;
    for (; begin < end; ++begin, --end) {
        for (int i = 0; i < individuals; ++i) {
            //haplotypes
            std::swap(haplotypes[i * 2][begin], haplotypes[i * 2][end]);
            std::swap(haplotypes[i * 2 + 1][begin], haplotypes[i * 2 + 1][end]);
            //genotypes
            std::swap(genotypes[i][begin * 3], genotypes[i][end * 3]);
            std::swap(genotypes[i][begin * 3 + 1], genotypes[i][end * 3 + 1]);
            std::swap(genotypes[i][begin * 3 + 2], genotypes[i][end * 3 + 2]);
        }

        for (int i = 0; i < (individuals - phased) * nSampleCopy; ++i) {
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

    int markerindex = marker * 3;
//    int ph11 = (unsigned char) genotypes[states / 2][markerindex];
//    int ph12 = (unsigned char) genotypes[states / 2][markerindex + 1];
//    int ph22 = (unsigned char) genotypes[states / 2][markerindex + 2];
    int ph11 = (unsigned char) genotypes[currentIndividual][markerindex];
    int ph12 = (unsigned char) genotypes[currentIndividual][markerindex + 1];
    int ph22 = (unsigned char) genotypes[currentIndividual][markerindex + 2];

    CalculatePhred2Prob();

    double posterior_11 = Penetrance(marker, copied1 + copied2, 0) * phred2prob[ph11];
    double posterior_12 = Penetrance(marker, copied1 + copied2, 1) * phred2prob[ph12];
    double posterior_22 = Penetrance(marker, copied1 + copied2, 2) * phred2prob[ph22];
    double sum = posterior_11 + posterior_12 + posterior_22;

    posterior_11 /= sum;
    posterior_22 /= sum;

    double r = rand->Next();

    if (r < posterior_11)//homo ref alleles
    {
        if (copied1 != copied2 || copied1 != 0)
            fprintf(stdout, "individidual: %d,Homo ref Marker:%d from else (%d,%d) and orginal gl:(%d,%d,%d)\n",
                    currentIndividual, marker, copied1, copied2, ph11, ph12, ph22);

        haps[currentHap1][marker] = 0;
        haps[currentHap2][marker] = 0;
    } else if (r < posterior_11 + posterior_22)//home alt alleles
    {
        if (copied1 != copied2 || copied1 != 1)
            fprintf(stdout, "individual: %d,Homo alt Marker:%d from else (%d,%d) and orginal gl:(%d,%d,%d)\n",
                    currentIndividual, marker, copied1, copied2, ph11, ph12, ph22);
        haps[currentHap1][marker] = 1;
        haps[currentHap2][marker] = 1;
    } else if (copied1 != copied2)//heter states and heter alleles
    {
//        double rate = GetErrorRate(marker);

//        if (rand->Next() < rate * rate / ((rate * rate) + (1.0 - rate) * (1.0 - rate)))//if both alleles mutated
//        {
//            copied1 = !copied1;
//            copied2 = !copied2;
//        }

        haps[currentHap1][marker] = copied1;
        haps[currentHap2][marker] = copied2;
    } else//homo states but heter alleles
    {
        fprintf(stdout, "individual: %d, Heter Marker:%d from else (%d,%d) and orginal gl:(%d,%d,%d)\n",
                currentIndividual, marker, copied1, copied2, ph11, ph12, ph22);
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

double PBWTHaplotyper::GetGL(int individual, int marker, uchar allele1, uchar allele2) {
    return phred2prob[(size_t) genotypes[individual][3 * marker + allele1 + allele2]];
}

float PBWTHaplotyper::GetRecombRate(int marker) {
    if (isRev) {
        marker = Wrapper->nMarkers - 1 - marker;
    }
    return Wrapper->recomRate[marker];
}

void PBWTHaplotyper::RandomSetup(Random *rand) {
    if (rand == NULL)
        rand = &globalRandom;

    CalculatePhred2Prob();

    for (int j = 0; j < markers; j++) {
        int markerindex = 3 * j;
        for (int i = 0; i < individuals; i++) {

            int posterior_11 = genotypes[i][markerindex];
            int posterior_12 = genotypes[i][markerindex + 1];
            int posterior_22 = genotypes[i][markerindex + 2];
            int min = std::min(std::min(posterior_11, posterior_12), posterior_22);//phred score
            if (min == posterior_11) {
                haplotypes[i * 2][j] = 0;
                haplotypes[i * 2 + 1][j] = 0;
            } else if (min == posterior_12) {
                bool bit = rand->Binary();

                haplotypes[i * 2][j] = bit;
                haplotypes[i * 2 + 1][j] = bit ^ 1;
            } else {
                haplotypes[i * 2][j] = 1;
                haplotypes[i * 2 + 1][j] = 1;
            }

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

#ifdef NAIVE
int PBWTHaplotyper::ForwardAlgorithm(int sampleIndex, FwdBwdLocalParameter &localParameter) {
    float prevFwdValue(0.f);
    float tmpFwdValue(0.f);
    float gl(0.f);
    float lowestFwd(0.f);//100th smallest fwd value
    int numCrediablePair(0);
//    int numInitialHetSite(0);

    int fitPair = 0;
    int totalPair = 0;
    int noChildPair = 0;

    std::priority_queue<EdgePair, std::vector<EdgePair>, std::function<bool(const EdgePair &,
                                                                            const EdgePair &)> > EdgePairList(
            EdgePaircomparator);

    uchar allele1, allele2;
    StateIndex childNode1, childNode2, parentNode1, parentNode2;
    for (int i = 1; i < markers; i++) {
        fitPair = 0;
        totalPair = 0;
        noChildPair = 0;
        localParameter.fwdValueSum[i] = 0.f;
        lowestFwd = 0.f;
        numCrediablePair = 0;

        for (auto iter = localParameter.genuienParents[i - 1].begin();
             iter != localParameter.genuienParents[i - 1].end(); ++iter)//all states at site i-1, parentNode1: hap1
        {
            parentNode1 = iter->first;
            for (auto iter2 = iter->second.begin(); iter2 !=
                                                    iter->second.end(); ++iter2)// parentNode2:hap2; iter2->second:all the source states to current state
            {
                parentNode2 = iter2->first;
                prevFwdValue = localParameter.SumFwdValueFromOriginVec(
                        iter2->second);//sum over all the parents of current parents that lead to a child pair
                totalPair++;
                for (allele1 = 0; allele1 < 2; ++allele1) {
                    childNode1 = GetChildNode(i - 1, parentNode1, allele1);//i child of i-1 child state
                    if (childNode1 == -1) continue;
                    for (allele2 = 0; allele2 < 2; ++allele2) {
                        childNode2 = GetChildNode(i - 1, parentNode2, allele2);//i child of i-1 child state
                        if (childNode2 == -1) continue;
                        gl = GetGL(sampleIndex, i, allele1, allele2);//i gl
//                        fprintf(stderr,"site:%d,prev(%d,%d) to current(%d,%d) : allell1:%d\tallele2:%d\tedgeNumHap1:%g\tedgeNumHap2:%g\tgl:%g\n",
//                                    i,parentNode1,parentNode2,childNode1,childNode2,allele1,allele2,GetTransitionProb(i-1,parentNode1,childNode1),GetTransitionProb(i-1,parentNode2,childNode2),gl);
                        if (gl > 1e-1 || i < 20) {//if fwd algorithm broke, relex genotype constraint
                            fitPair++;
                            tmpFwdValue = prevFwdValue *
                                          GetTransitionProb(i - 1, parentNode1, childNode1) *
                                          GetTransitionProb(i - 1, parentNode2, childNode2) *
                                          gl;//i fwdValueSum
                            if (tmpFwdValue < UNDERFLOW_MIN && prevFwdValue > 0) {
                                tmpFwdValue = UNDERFLOW_MIN;
                            }
//                            fprintf(stderr,"site:%d,prev(%d,%d) to current(%d,%d) : %g and prevFwd:%g\ttp1:%g\ttp2:%g\tgl:%g fwdValueSum:%g\n",
//                                    i,parentNode1,parentNode2,childNode1,childNode2,tmpFwdValue,prevFwdValue,
//                                    GetTransitionProb(i-1,parentNode1,childNode1),
//                                    GetTransitionProb(i-1,parentNode2,childNode2),gl,fwdValueSum[i]);
//                            (*nextFwdValuePtr)[childNode1][childNode2].push_back(origin(parentNode1,parentNode2,tmpFwdValue));
                            localParameter.fwdValueSum[i] += tmpFwdValue;
                            localParameter.genuienParents[i][childNode1][childNode2][std::make_pair(parentNode1,
                                                                                     parentNode2)] = tmpFwdValue;
                        } else {
                            tmpFwdValue = prevFwdValue *
                                          GetTransitionProb(i - 1, parentNode1, childNode1) *
                                          GetTransitionProb(i - 1, parentNode2, childNode2) *
                                          gl;//i fwdValueSum
                            if (tmpFwdValue < UNDERFLOW_MIN && prevFwdValue > 0) {
                                tmpFwdValue = UNDERFLOW_MIN;
                            }

                            if (numCrediablePair < 200) {
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
        }

        while (not EdgePairList.empty()) {
            EdgePair tmpEdgePair = EdgePairList.top();
            EdgePairList.pop();
            localParameter.genuienParents[i][tmpEdgePair.childNode1][tmpEdgePair.childNode2][std::make_pair(tmpEdgePair.parentNode1,
                                                                                             tmpEdgePair.parentNode2)] = tmpEdgePair.fwd;
            localParameter.fwdValueSum[i] += tmpEdgePair.fwd;
        }
//        if (fitPair == 0) {
//            fprintf(stderr, "[Warning]marker %d broken! %d totalPair, %d noChildPair\n", i, totalPair, noChildPair);
////            exit(EXIT_FAILURE);
//            brokenList[i] = true;
//            goto REENTRY;
//            UpdateStateNum(GetStateNumFrom(i));
//            prevFwdValue = 1.f / (states * states);
//            for (int j = 0; j < states; ++j)
//                for (int k = 0; k < states; ++k) {
//                    allele1 = GetAllele(i, j);
//                    allele2 = GetAllele(i, k);
//                    gl = GetGL(sampleIndex, i, allele1, allele2);//i gl
//                    if (gl > 1e-1) {
//                        tmpFwdValue = prevFwdValue * gl;//i fwdValueSum
//                        if (tmpFwdValue < UNDERFLOW_MIN && prevFwdValue > 0) {
//                            tmpFwdValue = UNDERFLOW_MIN;
//                        }
//                        fwdValueSum[i] += tmpFwdValue;
//                        genuienParents[i][j][k][std::make_pair(0, 0)] = tmpFwdValue;
//                    }
//                }
//        }
//        fprintf(stderr, "Effective rate %f out of totalPair %d at site %d...\n", fitPair / double(totalPair), totalPair,i);
//        fprintf(stderr,"site:%d report overall fwd:%g\n",i,fwdValueSum[i]);
        for (auto iter = localParameter.genuienParents[i].begin(); iter != localParameter.genuienParents[i].end(); ++iter)
            for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); ++iter2)
                for (auto iter3 = iter2->second.begin(); iter3 != iter2->second.end(); ++iter3) {
                    iter3->second /= localParameter.fwdValueSum[i];
                }
    }
    return 0;
}

int PBWTHaplotyper::BackwardSampling(Random *rand, int sampleIndex, char **sampledHaps,
                                     FwdBwdLocalParameter localParameter) {
    double choice(0.);
    double sum(0.);
    double gl0(0.f);
    float transProb0(0.f), transProb1(0.f);
//    float prevFwdValue(0.f);
    StateIndex first0(0), second0(0);

    StateIndex sampledFirst(0);
    StateIndex sampledSecond(0);
    double sampledFwd(0.f);

    choice = rand->Uniform(0, 1);
    int count=0;
    for(auto iter=localParameter.genuienParents[markers-1].begin();iter!=localParameter.genuienParents[markers-1].end();++iter)
        for(auto iter2=iter->second.begin();iter2!=iter->second.end();++iter2)
        {
            count++;
            sampledFwd = localParameter.SumFwdValueFromOriginVec(iter2->second);
            fprintf(stderr,"(%d,%d)\tvalue:%g\tcount:%d\n",iter->first,iter2->first,sampledFwd,count);


        }
    for (auto iter = localParameter.genuienParents[markers - 1].begin(); iter != localParameter.genuienParents[markers - 1].end(); ++iter)
        for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); ++iter2) {
            sampledFwd = localParameter.SumFwdValueFromOriginVec(iter2->second);
            sum += sampledFwd;
            if (sum > choice) {
                sampledFirst = iter->first;
                sampledSecond = iter2->first;
                goto SAMPLE_BREAK;
            }
        }
    SAMPLE_BREAK:
//    first0 = sampledFirst;
//    second0 = sampledSecond;
//    ImputeAlleles(markers - 1, first0, second0, rand, sampleIndex, sampledHaps);
////    fprintf(stderr,"marker:%d\tsum:%g\tchoice:%g\t(%d,%d)\tvalue:%g\toverallFwd:%g\n",markers-1,sum,choice,sampledFirst,sampledSecond,prevFwdValue,fwdValueSum[markers-1]);
//    sum = 0.;
//    choice = rand->Uniform(0, prevFwdValue);
//    for (auto kv: genuienParents[markers - 1][first0][second0]) {
////            transProb0 = GetTransitionProb(markers - 2, kv.first.first, sampledFirst);
////            transProb1 = GetTransitionProb(markers - 2, kv.first.second, sampledSecond);
////            prevFwdValue = kv.second * fwdValueSum[markers - 1] / (transProb0 * transProb1 * gl0);
//        sum += kv.second;
//        if (sum > choice) {
//            sampledFirst = kv.first.first;
//            sampledSecond = kv.first.second;
//            gl0 = GetGL(individuals - 1, markers - 1, GetAllele(markers - 1, first0), GetAllele(markers - 1, second0));
//            transProb0 = GetTransitionProb(markers - 2, sampledFirst, first0);
//            transProb1 = GetTransitionProb(markers - 2, sampledSecond, second0);
//            sampledFwd = kv.second * fwdValueSum[markers - 1] / (transProb0 * transProb1 * gl0);
//            break;
//        }
//    }


    for (int i = markers - 1; i > 0; --i) {
        fprintf(stderr,"marker:%d\tsum:%g\tchoice:%g\t(%d,%d)\tvalue:%g\n",i,sum,choice,sampledFirst,sampledSecond,sampledFwd);
//        fprintf(stderr,"marker:%d\tsum:%g\n",i,sum);

        sum = 0.;
        first0 = sampledFirst;
        second0 = sampledSecond;

        ImputeAlleles(i, first0, second0, rand, sampleIndex, sampledHaps);

        {
            if (fabs(sampledFwd - UNDERFLOW_MIN) < std::numeric_limits<double>::epsilon())
                sampledFwd = localParameter.SumFwdValueFromOriginVec(localParameter.genuienParents[i][first0][second0]);

            choice = rand->Uniform(0, sampledFwd);
            if (localParameter.genuienParents[i].find(first0) == localParameter.genuienParents[i].end() ||
                    localParameter.genuienParents[i][first0].find(second0) == localParameter.genuienParents[i][first0].end()) {
                fprintf(stderr, "marker %d does not have state %d or %d\n", i, first0, second0);
                exit(EXIT_FAILURE);
            }

            for (auto kv:localParameter.genuienParents[i][first0][second0])//we are actually sampling i-1's states
            {
                sum += kv.second;//kv.second is normalized
                if (sum > choice) {
                    sampledFirst = kv.first.first;
                    sampledSecond = kv.first.second;
                    gl0 = GetGL(sampleIndex, i, GetAllele(i, first0), GetAllele(i, second0));
                    transProb0 = GetTransitionProb(i - 1, sampledFirst, first0);
                    transProb1 = GetTransitionProb(i - 1, sampledSecond, second0);
                    sampledFwd = kv.second * localParameter.fwdValueSum[i] / (transProb0 * transProb1 * gl0);
                    break;
                }
            }
        }
//        FillPath(sampleIndex * 2, i, i0 + 1, sampledFirst,haplotypes);
//        FillPath(sampleIndex * 2 + 1, i, i0 + 1, sampledSecond,haplotypes);
    }

    ImputeAlleles(0, sampledFirst, sampledSecond, rand, sampleIndex, sampledHaps);
    return 0;
}

#endif
// with recombination

//beagle version
#ifdef RECBEAGLE
int PBWTHaplotyper::ForwardAlgorithmRecBeagle() {
    int SampleIndex = individuals - 1;
    UpdateStateNum(GetStateNumFrom(0));
    InitialFwdValues(<#initializer#>);
    float prevFwdValue(0.f);
    float tmpFwdValue(0.f),nodePairFwd(0.f),lowestFwd(0.f);
    float gl(0.f);


    int totalPair = 0;
    int noChildPair = 0;
    int numCrediablePair =0;

    char allele1(0), allele2(0);
//    StateIndex childNode1(0), childNode2(0), parentNode1(0), parentNode2(0);

    float recRate(0.f);
    float baseProb(0.f);
    availablePair.ResetMarkerIndexAt(0);//TODO: use improved ibs methods maybe
    availablePair.NextMarker();//marker 0 doesn't have parents
    for (int i = 1; i < markers; i++) {
        totalPair = 0;
        noChildPair = 0;
        fwdValueSum[i] = 0;
        recRate = GetRecombRate(i - 1);//start from marker 1 but store at index 0
        numCrediablePair=0;
        //process orphan nodes
        UpdateStateNum(GetStateNumFrom(i));
        for (StateIndex childNode1 = 0; childNode1 <states; ++childNode1) {//dest nodeA
            for (StateIndex childNode2 = 0; childNode2 <states; ++childNode2) {//dest nodeB
                gl = GetGL(SampleIndex, i, GetAllele(i,childNode1), GetAllele(i,childNode2));
                if ( gl < 1e-1) continue;
                nodePairFwd = 0.f;
                for (auto parentNode1:GetParentNodes(i,childNode1)) {
                    for (auto parentNode2:GetParentNodes(i,childNode2))
                    {
                        if(genuienParents[i-1].find(parentNode1)!=genuienParents[i-1].end() && genuienParents[i-1][parentNode1].find(parentNode2)!= genuienParents[i-1][parentNode1].end())
                            prevFwdValue = SumFwdValueFromOriginVec(genuienParents[i-1][parentNode1][parentNode2]);//sum over all the parents of current parents that lead to a child pair
                        else
                            prevFwdValue =0.f;

                        baseProb = GetTransitionProb(i - 1, parentNode1, childNode1) *
                                   GetTransitionProb(i - 1, parentNode2, childNode2) * gl;

                        tmpFwdValue = prevFwdValue * baseProb * (1 - recRate) * (1 - recRate);
//                        fprintf(stderr,"(%d,%d) to (%d,%d) tmp output1 tmpFwdValue:%g\tprevFwdValue:%g\tbaseProb:%g\n",parentNode1,parentNode2,childNode1,childNode2,tmpFwdValue,prevFwdValue,baseProb);
                        if(fwdValueNode1Sum[i - 1].find(parentNode1)!= fwdValueNode1Sum[i - 1].end())
                        tmpFwdValue += fwdValueNode1Sum[i - 1][parentNode1] * GetHapProbAt(i - 1, parentNode2) *
                                       baseProb * (1 - recRate) * recRate;
//                        fprintf(stderr,"tmp output2 tmpFwdValue:%g\tfwdValueNode1Sum[i - 1][parentNode1]:%g\n",tmpFwdValue,fwdValueNode1Sum[i - 1][parentNode1]);
                        if(fwdValueNode2Sum[i - 1].find(parentNode2)!= fwdValueNode2Sum[i - 1].end())
                            tmpFwdValue += fwdValueNode2Sum[i - 1][parentNode2] * GetHapProbAt(i - 1, parentNode1) *
                                       baseProb * (1 - recRate) * recRate;
//                        fprintf(stderr,"tmp output3 tmpFwdValue:%g\tfwdValueNode2Sum[i - 1][parentNode2]:%g\n",tmpFwdValue,fwdValueNode2Sum[i - 1][parentNode2]);
                        tmpFwdValue += GetHapProbAt(i - 1, parentNode2) *
                                       GetHapProbAt(i - 1, parentNode1) * baseProb * recRate * recRate;
//                        fprintf(stderr,"tmp output4 tmpFwdValue:%g\n",tmpFwdValue);
//                        if (tmpFwdValue < UNDERFLOW_MIN && prevFwdValue > 0) {
//                            tmpFwdValue = UNDERFLOW_MIN;
//                        }


                        if (numCrediablePair < 100) {
//                                fprintf(stderr,"site:%d push overall fwd:%g\n",i,tmpFwdValue);
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

        while (not EdgePairList.empty()) {
            EdgePair tmpEdgePair = EdgePairList.top();
            EdgePairList.pop();
            if (tmpEdgePair.fwd > 0.f) {
//                fprintf(stderr,"site:%d report overall fwd:%g\n",i,tmpEdgePair.fwd);
                genuienParents[i][tmpEdgePair.childNode1][tmpEdgePair.childNode2][std::make_pair(
                        tmpEdgePair.parentNode1,
                        tmpEdgePair.parentNode2)] = tmpEdgePair.fwd;
                fwdValueSum[i] += tmpEdgePair.fwd;
                availablePair.FillNextAvailableStatePair(
                        std::make_pair(tmpEdgePair.childNode1, tmpEdgePair.childNode2));
            }
        }






        for (auto iter = genuienParents[i].begin(); iter != genuienParents[i].end(); ++iter)
            for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); ++iter2) {
                float tmp(0.f),tmp2(0.f);
                for (auto iter3 = iter2->second.begin(); iter3 != iter2->second.end(); ++iter3) {
                    tmp2 += iter3->second;
                    iter3->second /= fwdValueSum[i];
                    tmp += iter3->second;
                    fprintf(stderr, "foward marker:%d\tfwdValueSum:%g\tpair:(%d,%d)=(%d,%d)\tcurrentFwd:%g\tfrom:(%d,%d)\n",
                            i, fwdValueSum[i],iter->first,iter2->first,
                            GetAllele(i,iter->first),GetAllele(i,iter2->first),tmp2,iter3->first.first,iter3->first.second);
                }
                fwdValueNode1Sum[i][iter->first] += tmp;
                fwdValueNode2Sum[i][iter2->first] += tmp;

            }
//        if(i==9604) fprintf(stderr, "sampled:%d\tlalala marker:%d\tfwdValueSum:%g\n", sampled,i, fwdValueSum[i]);
        availablePair.NextMarker();//last one is empty
    }
    return 0;
}

int PBWTHaplotyper::BackwardSamplingRecBeagle(Random *rand, int SampleIndex, char **sampledHaps) {

    double choice(0.);
    double sum(0.), subSum(0.);
    float gl0(0.f);
    float baseProb(0.f);
    float np1(0.f), np2(0.f);

    int first0(0), second0(0);
    float edgeProbFirst0(0.f), edgeProbSecond0(0.f);

    int sampledGrandParent0(0),sampledGrandParent1(0);

    int sampledParent0(0);
    int sampledParent1(0);
    double sampledFwd(0.f);

    float recRate(0.f);

    int sampled = 0;

    choice = rand->Uniform(0, 1);
    availablePair.ResetMarkerIndexAt(markers - 1);//rewind to last actual point
    for (auto iter = genuienParents[markers - 1].begin(); iter != genuienParents[markers - 1].end(); ++iter)
        for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); ++iter2) {
            for (auto grandParents : iter2->second) {
                sum += grandParents.second;
                fprintf(stderr, "sampling last marker:(%d,%d) with sum:%g, recRate:%g\n",iter->first,iter2->first,sum,recRate);
                if (sum > choice) {
                    sampledParent0 = iter->first;
                    sampledParent1 = iter2->first;
                    sampledGrandParent0 = grandParents.first.first;
                    sampledGrandParent1 = grandParents.first.second;
                    edgeProbFirst0 = GetEdgeProbAt(markers - 2, sampledGrandParent0,
                                                   GetAllele(markers - 1, sampledParent0));
                    edgeProbSecond0 = GetEdgeProbAt(markers - 2, sampledGrandParent1,
                                                    GetAllele(markers - 1, sampledParent1));
                    sampledFwd = grandParents.second * fwdValueSum[markers - 1];
                    goto INIT_SAMPLE_BREAK;
                }
            }
        }
    INIT_SAMPLE_BREAK:
    availablePair.PrevMarker();
    float testSum(0.f);
    int tmpSampled=0;
    for (int i = markers - 1; i > 1; --i) {

//        if (i == 4102 || i == 4103 || i == 4104) {
        fprintf(stderr,
                "marker:%d\ttestSum:%g\tsum:%g\tchoice:%g\tprev(%d,%d) to (%d,%d)\tsampledFwd:%g\tavailablePair.Size:%d\trecRate:%g\tgl:%g\tprobEdge0:%g\tprobEdge1:%g\tfwdValueSum:%g\ttmpSampled:%d\n",
                i, testSum, sum,
                choice, sampledGrandParent0, sampledGrandParent1, sampledParent0, sampledParent1, sampledFwd, availablePair.Size(),
                recRate, gl0, edgeProbFirst0,edgeProbSecond0, fwdValueSum[i],tmpSampled);

        for (auto kv:genuienParents[i][sampledParent0][sampledParent1]) {
            fprintf(stderr, "marker:%d\tpair:(%d,%d)\ttestSum:%g\n", i, kv.first.first, kv.first.second, kv.second);
        }
//        }
        sum = 0.f;
        testSum=0.f;

        first0 = sampledParent0;
        second0 = sampledParent1;

        ImputeAlleles(i, first0, second0, rand, SampleIndex, sampledHaps);
        sampled = 0;
        choice = rand->Uniform(0, sampledFwd);

        if (genuienParents[i].find(first0) == genuienParents[i].end() ||
            genuienParents[i][first0].find(second0) == genuienParents[i][first0].end()) {
            fprintf(stderr, "marker %d does not have state %d or %d\n", i, first0, second0);
            exit(EXIT_FAILURE);
        }

        gl0 = GetGL(individuals - 1, i, GetAllele(i, first0), GetAllele(i, second0));
        recRate = GetRecombRate(i - 1);//start from marker 1 but store at index 0
        baseProb = edgeProbFirst0 * edgeProbSecond0 * gl0;//site i

//
//        if (genuienParents[i][first0][second0].size() > 0)//has perfect match parent
//        {
//            subSum =0.;
//            for (auto parents:genuienParents[i][first0][second0]) {
//                np1 = GetHapProbAt(i - 1, parents.first.first);
//                np2 = GetHapProbAt(i - 1, parents.first.second);
//                subSum += (1 - recRate) * (1 - recRate) * baseProb / (np1 * np2);
//                subSum += (1 - recRate) * recRate * baseProb / np1;
//                subSum += (1 - recRate) * recRate * baseProb / np2;
//                subSum += recRate * recRate * baseProb;
//                sum += subSum * parents.second;//grandParents.second is the fwd value at site i-1, contributed by grandParents at site i-2
//                if (sum > choice) {
//                    sampledParent0 = parentNodePair.first;
//                    sampledParent1 = parentNodePair.second;
//                    sampledGrandParent0 = grandParents.first.first;
//                    sampledGrandParent1 = grandParents.first.second;
//                    edgeProbFirst0 = GetEdgeProbAt(i - 2, sampledGrandParent0,
//                                                   GetAllele(i - 1, sampledParent0));
//                    edgeProbSecond0 = GetEdgeProbAt(i - 2, sampledGrandParent1,
//                                                    GetAllele(i - 1, sampledParent1));
//                    sampledFwd = grandParents.second * fwdValueSum[i - 1];
//                    testSum = grandParents.second;
//                    goto SAMPLE_BREAK;
//                }
//            }
//        }
//        else//no perfect match, then we consider other haplotypes
        while (!availablePair.IsEnd()) {

            std::pair<int, int> parentNodePair = availablePair.GetNextAvailableStatePair();
            np1 = GetHapProbAt(i - 1, parentNodePair.first);
            np2 = GetHapProbAt(i - 1, parentNodePair.second);
            subSum = 0.;
            tmpSampled=0;

//            if (genuienParents[i][first0][second0].find(
//                    std::make_pair(parentNodePair.first, parentNodePair.second)) !=
//                genuienParents[i][first0][second0].end())//no recomb
//
//            if (GetChildNode(i - 1, parentNodePair.first, GetAllele(i, first0)) == first0)//noJump for first
//                subSum += (1 - recRate) * recRate * baseProb / np1;
//
//            if (GetChildNode(i - 1, parentNodePair.second, GetAllele(i, second0)) == second0)//noJump for second
//                subSum += (1 - recRate) * recRate * baseProb / np2;

            if(parentNodePair.first == sampledGrandParent0 && parentNodePair.second == sampledGrandParent1)//no recomb
            {
                subSum += (1 - recRate) * (1 - recRate) * baseProb / (np1 * np2);
                tmpSampled+=4;
            }
            if (parentNodePair.first == sampledGrandParent0)//noJump for first
            {
                subSum += (1 - recRate) * recRate * baseProb / np1;
                tmpSampled+=2;
            }

            if (parentNodePair.second == sampledGrandParent1)//noJump for second
            {
                subSum += (1 - recRate) * recRate * baseProb / np2;
                tmpSampled+=1;
            }

            subSum += recRate * recRate * baseProb;
            if(tmpSampled > sampled) sampled=tmpSampled;

            fprintf(stderr, "sampling:(%d,%d)[%g,%g] to (%d,%d) with subSum:%g, recRate:%g, tmpSampled:%d\n",parentNodePair.first,parentNodePair.second,np1,np2,first0,second0,subSum,recRate,tmpSampled);
            if(genuienParents[i-1].find(parentNodePair.first)!= genuienParents[i-1].end() && genuienParents[i-1][parentNodePair.first].find(parentNodePair.second)!= genuienParents[i-1][parentNodePair.first].end())
            {
                for (auto grandParents:genuienParents[i - 1][parentNodePair.first][parentNodePair.second])
                {//all grandParents ending at parentNodePair.first and parentNodePair.second
                    sum += subSum * grandParents.second;//grandParents.second is the fwd value at site i-1, contributed by grandParents at site i-2
                    if (sum > choice) {
                        sampledParent0 = parentNodePair.first;
                        sampledParent1 = parentNodePair.second;
                        sampledGrandParent0 = grandParents.first.first;
                        sampledGrandParent1 = grandParents.first.second;
                        edgeProbFirst0 = GetEdgeProbAt(i - 2, sampledGrandParent0,
                                                       GetAllele(i - 1, sampledParent0));
                        edgeProbSecond0 = GetEdgeProbAt(i - 2, sampledGrandParent1,
                                                        GetAllele(i - 1, sampledParent1));
                        sampledFwd = grandParents.second * fwdValueSum[i - 1];
                        testSum = grandParents.second;
                        goto SAMPLE_BREAK;
                    }
                }
            } else
            {
                fprintf(stderr,"shouldn't happen here at marker %d\n",i-1);
            }
        }
        SAMPLE_BREAK:
        if (sum < choice)
            fprintf(stderr, "inside marker:%d\ttestSum:%g\tsum:%g\tchoice:%g\t(%d,%d)\ttotalValue:%g\tsampled:%d\tsubSum:%g\trecRate:%g\tbaseProb:%g\tfwdValueSum:%g\n",
                    i-1, testSum,sum, choice, sampledParent0, sampledParent1, sampledFwd, sampled, subSum, recRate, baseProb,fwdValueSum[i - 1]);
            availablePair.PrevMarker();
    }

    ImputeAlleles(1, sampledParent0, sampledParent1, rand, SampleIndex, sampledHaps);
    ImputeAlleles(0, sampledGrandParent0, sampledGrandParent1, rand, SampleIndex, sampledHaps);
    return 0;
}
#endif
//end of beagle version


//beagle variant version
#ifdef RECBEAGLEVARIANT
int PBWTHaplotyper::ForwardAlgorithmRecNew() {
    int SampleIndex = individuals - 1;
    UpdateStateNum(GetStateNumFrom(0));
    InitialFwdValues(<#initializer#>);
    float prevFwdValue(0.f);
    float tmpFwdValue(0.f);
    float gl(0.f);
    float lowestFwd(0.f);//100th smallest fwd value
    int numCrediablePair(0);

    int fitPair = 0;
    int totalPair = 0;
    int noChildPair = 0;

    char allele1(0), allele2(0);
    StateIndex childNode1(0), childNode2(0), parentNode1(0), parentNode2(0);

    float recRate(0.f);
    float baseProb(0.f);
    availablePair.ResetMarkerIndexAt(0);//TODO: use improved ibs methods maybe
    availablePair.NextMarker();//marker 0 doesn't have parents
    for (int i = 1; i < markers; i++) {
        fitPair = 0;
        totalPair = 0;
        noChildPair = 0;
        fwdValueSum[i] = 0;
        numCrediablePair = 0;
        recRate = GetRecombRate(i - 1);//start from marker 1 but store at index 0

        for (auto iter = genuienParents[i - 1].begin();
             iter != genuienParents[i - 1].end(); ++iter)//all states at site i-1, parentNode1: hap1
        {
            parentNode1 = iter->first;
            for (auto iter2 = iter->second.begin();
                 iter2 != iter->second.end(); ++iter2)// parentNode2:hap2; iter2->second:all the source states to current state
            {
                parentNode2 = iter2->first;

                prevFwdValue = SumFwdValueFromOriginVec(
                        iter2->second);//sum over all the parents of current parents that lead to a child pair
                totalPair++;
                for (allele1 = 0; allele1 < 2; ++allele1) {
                    childNode1 = GetChildNode(i - 1, parentNode1, allele1);//i child of i-1 child state
                    if (childNode1 == -1) continue;
                    for (allele2 = 0; allele2 < 2; ++allele2) {
                        childNode2 = GetChildNode(i - 1, parentNode2, allele2);//i child of i-1 child state
                        if (childNode2 == -1) continue;
                        gl = GetGL(SampleIndex, i, allele1, allele2);
                        if (gl > 1e-1) {
                            fitPair++;
                            baseProb = GetTransitionProb(i - 1, parentNode1, childNode1) *
                                       GetTransitionProb(i - 1, parentNode2, childNode2) * gl;

                            tmpFwdValue = prevFwdValue * baseProb * (1 - recRate) * (1 - recRate);

                            tmpFwdValue += fwdValueNode1Sum[i - 1][parentNode1] * GetHapProbAt(i - 1, parentNode2) *
                                           baseProb * (1 - recRate) * recRate;

                            tmpFwdValue += fwdValueNode2Sum[i - 1][parentNode2] * GetHapProbAt(i - 1, parentNode1) *
                                           baseProb * (1 - recRate) * recRate;

                            tmpFwdValue += GetHapProbAt(i - 1, parentNode2) *
                                           GetHapProbAt(i - 1, parentNode1) * baseProb * recRate * recRate;
//                            if(i==813)
//                            {
//                                fprintf(stderr,"special watchpoint--tp1:%g\ttp2:%g\tgl:%g\tpreFwdValue:%g\tfwdValueNode1Sum:%g\tfwdValueNode2Sum:%g\tnp1:%g\tnp2:%g\tfrom (%d,%d) to (%d,%d)\n",
//                                        GetTransitionProb(i - 1, parentNode1, childNode1),GetTransitionProb(i - 1, parentNode2, childNode2),gl,prevFwdValue,fwdValueNode1Sum[i - 1][parentNode1],fwdValueNode2Sum[i - 1][parentNode2],
//                                        GetHapProbAt(i - 1, parentNode1),GetHapProbAt(i - 1, parentNode2),parentNode1,parentNode2,childNode1,childNode2);
//                            }
                            if (tmpFwdValue < UNDERFLOW_MIN && prevFwdValue > 0) {
                                tmpFwdValue = UNDERFLOW_MIN;
                            }

                            fwdValueSum[i] += tmpFwdValue;
                            genuienParents[i][childNode1][childNode2][std::make_pair(parentNode1,
                                                                                     parentNode2)] = tmpFwdValue;
                            availablePair.FillNextAvailableStatePair(std::make_pair(parentNode1, parentNode2));

                        }
//                        else {
//                            baseProb = GetTransitionProb(i - 1, parentNode1, childNode1) *
//                                       GetTransitionProb(i - 1, parentNode2, childNode2) * gl;
//
//                            tmpFwdValue = prevFwdValue * baseProb * (1 - recRate) * (1 - recRate);
//
//                            tmpFwdValue += fwdValueNode1Sum[i - 1][parentNode1] * GetHapProbAt(i - 1, parentNode2) *
//                                           baseProb * (1 - recRate) * recRate;
//
//                            tmpFwdValue += fwdValueNode2Sum[i - 1][parentNode2] * GetHapProbAt(i - 1, parentNode1) *
//                                           baseProb * (1 - recRate) * recRate;
//
//                            tmpFwdValue += GetHapProbAt(i - 1, parentNode2) *
//                                           GetHapProbAt(i - 1, parentNode1) * baseProb * recRate * recRate;
//
//                            if (tmpFwdValue < UNDERFLOW_MIN && prevFwdValue > 0) {
//                                tmpFwdValue = UNDERFLOW_MIN;
//                            }
//
//                            if (numCrediablePair < 200) {
////                                fprintf(stderr,"site:%d push overall fwd:%g\n",i,tmpFwdValue);
//                                EdgePairList.push(
//                                        EdgePair(childNode1, childNode2, parentNode1, parentNode2, tmpFwdValue));
//                                numCrediablePair++;
//                                if (tmpFwdValue < lowestFwd) {
//                                    lowestFwd = tmpFwdValue;
//                                }
//                            } else if (tmpFwdValue >
//                                       lowestFwd)// EdgePairList full and should be added into List, pop out lowest
//                            {
//                                EdgePairList.pop();
//                                EdgePairList.push(
//                                        EdgePair(childNode1, childNode2, parentNode1, parentNode2, tmpFwdValue));
//                                lowestFwd = EdgePairList.top().fwd;
//                            }
//                        }
                    }
                }
//                availablePair.FillNextAvailableStatePair(std::make_pair(parentNode1, parentNode2));
            }
        }

//        while (not EdgePairList.empty()) {
//            EdgePair tmpEdgePair = EdgePairList.top();
//            EdgePairList.pop();
//            if (tmpEdgePair.fwd > 0.f) {
////                fprintf(stderr,"site:%d report overall fwd:%g\n",i,tmpEdgePair.fwd);
//                genuienParents[i][tmpEdgePair.childNode1][tmpEdgePair.childNode2][std::make_pair(
//                        tmpEdgePair.parentNode1,
//                        tmpEdgePair.parentNode2)] = tmpEdgePair.fwd;
//                fwdValueSum[i] += tmpEdgePair.fwd;
//                availablePair.FillNextAvailableStatePair(
//                        std::make_pair(tmpEdgePair.parentNode1, tmpEdgePair.parentNode2));
//            }
//        }

//       if(fitPair == 0 ) fprintf(stderr,"site:%d report overall fwd:%g\tfitPair:%d\n",i,fwdValueSum[i],fitPair);
//        int sampled(0);
        for (auto iter = genuienParents[i].begin(); iter != genuienParents[i].end(); ++iter)
            for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); ++iter2) {
                float tmp(0.f),tmp2(0.f);
                for (auto iter3 = iter2->second.begin(); iter3 != iter2->second.end(); ++iter3) {
                    tmp2 += iter3->second;
                    iter3->second /= fwdValueSum[i];
                    tmp += iter3->second;
                    fprintf(stderr, "foward marker:%d\tfwdValueSum:%g\tpair:(%d,%d)=(%d,%d)\tcurrentFwd:%g\tfrom:(%d,%d)\n",
                            i, fwdValueSum[i],iter->first,iter2->first,
                            GetAllele(i,iter->first),GetAllele(i,iter2->first),tmp2,iter3->first.first,iter3->first.second);
                }
                fwdValueNode1Sum[i][iter->first] += tmp;
                fwdValueNode2Sum[i][iter2->first] += tmp;

            }
//        if(i==9604) fprintf(stderr, "sampled:%d\tlalala marker:%d\tfwdValueSum:%g\n", sampled,i, fwdValueSum[i]);
        availablePair.NextMarker();//last one is empty
    }
    return 0;
}

int PBWTHaplotyper::BackwardSamplingRecNew(Random *rand, int SampleIndex, char **sampledHaps) {

    double choice(0.);
    double sum(0.), subSum(0.);
    float gl0(0.f);
    float baseProb(0.f);
    float np1(0.f), np2(0.f);

    int first0(0), second0(0);
    float edgeProbFirst0(0.f), edgeProbSecond0(0.f);

    int sampledGrandParent0(0),sampledGrandParent1(0);

    int sampledParent0(0);
    int sampledParent1(0);
    double sampledFwd(0.f);

    float recRate(0.f);

    int sampled = 0;

    choice = rand->Uniform(0, 1);
    availablePair.ResetMarkerIndexAt(markers - 1);//rewind to last actual point
    for (auto iter = genuienParents[markers - 1].begin(); iter != genuienParents[markers - 1].end(); ++iter)
        for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); ++iter2) {
            for (auto grandParents : iter2->second) {
                sum += grandParents.second;
                fprintf(stderr, "sampling last marker:(%d,%d) with sum:%g, recRate:%g\n",iter->first,iter2->first,sum,recRate);
                if (sum > choice) {
                    sampledParent0 = iter->first;
                    sampledParent1 = iter2->first;
                    sampledGrandParent0 = grandParents.first.first;
                    sampledGrandParent1 = grandParents.first.second;
                    edgeProbFirst0 = GetEdgeProbAt(markers - 2, sampledGrandParent0,
                                                   GetAllele(markers - 1, sampledParent0));
                    edgeProbSecond0 = GetEdgeProbAt(markers - 2, sampledGrandParent1,
                                                    GetAllele(markers - 1, sampledParent1));
                    sampledFwd = grandParents.second * fwdValueSum[markers - 1];
                    goto INIT_SAMPLE_BREAK;
                }
            }
        }
    INIT_SAMPLE_BREAK:
    float testSum(0.f);
    int tmpSampled=0;
    for (int i = markers - 1; i > 1; --i) {

//        if (i == 4102 || i == 4103 || i == 4104) {
            fprintf(stderr,
                    "marker:%d\ttestSum:%g\tsum:%g\tchoice:%g\tprev(%d,%d) to (%d,%d)\tsampledFwd:%g\tavailablePair.Size:%d\trecRate:%g\tgl:%g\tprobEdge0:%g\tprobEdge1:%g\tfwdValueSum:%g\ttmpSampled:%d\n",
                    i, testSum, sum,
                    choice, sampledGrandParent0, sampledGrandParent1, sampledParent0, sampledParent1, sampledFwd, availablePair.Size(),
                    recRate, gl0, edgeProbFirst0,edgeProbSecond0, fwdValueSum[i],tmpSampled);

            for (auto kv:genuienParents[i][sampledParent0][sampledParent1]) {
                fprintf(stderr, "marker:%d\tpair:(%d,%d)\ttestSum:%g\n", i, kv.first.first, kv.first.second, kv.second);
            }
//        }
        sum = 0.f;
        testSum=0.f;

        first0 = sampledParent0;
        second0 = sampledParent1;

        ImputeAlleles(i, first0, second0, rand, SampleIndex, sampledHaps);
        sampled = 0;
        choice = rand->Uniform(0, sampledFwd);

        if (genuienParents[i].find(first0) == genuienParents[i].end() ||
            genuienParents[i][first0].find(second0) == genuienParents[i][first0].end()) {
            fprintf(stderr, "marker %d does not have state %d or %d\n", i, first0, second0);
            exit(EXIT_FAILURE);
        }

        gl0 = GetGL(individuals - 1, i, GetAllele(i, first0), GetAllele(i, second0));
        recRate = GetRecombRate(i - 1);//start from marker 1 but store at index 0
        baseProb = edgeProbFirst0 * edgeProbSecond0 * gl0;//site i

//
//        if (genuienParents[i][first0][second0].size() > 0)//has perfect match parent
//        {
//            subSum =0.;
//            for (auto parents:genuienParents[i][first0][second0]) {
//                np1 = GetHapProbAt(i - 1, parents.first.first);
//                np2 = GetHapProbAt(i - 1, parents.first.second);
//                subSum += (1 - recRate) * (1 - recRate) * baseProb / (np1 * np2);
//                subSum += (1 - recRate) * recRate * baseProb / np1;
//                subSum += (1 - recRate) * recRate * baseProb / np2;
//                subSum += recRate * recRate * baseProb;
//                sum += subSum * parents.second;//grandParents.second is the fwd value at site i-1, contributed by grandParents at site i-2
//                if (sum > choice) {
//                    sampledParent0 = parentNodePair.first;
//                    sampledParent1 = parentNodePair.second;
//                    sampledGrandParent0 = grandParents.first.first;
//                    sampledGrandParent1 = grandParents.first.second;
//                    edgeProbFirst0 = GetEdgeProbAt(i - 2, sampledGrandParent0,
//                                                   GetAllele(i - 1, sampledParent0));
//                    edgeProbSecond0 = GetEdgeProbAt(i - 2, sampledGrandParent1,
//                                                    GetAllele(i - 1, sampledParent1));
//                    sampledFwd = grandParents.second * fwdValueSum[i - 1];
//                    testSum = grandParents.second;
//                    goto SAMPLE_BREAK;
//                }
//            }
//        }
//        else//no perfect match, then we consider other haplotypes
        while (!availablePair.IsEnd()) {

            std::pair<int, int> parentNodePair = availablePair.GetNextAvailableStatePair();
            np1 = GetHapProbAt(i - 1, parentNodePair.first);
            np2 = GetHapProbAt(i - 1, parentNodePair.second);
            subSum = 0.;
            tmpSampled=0;

//            if (genuienParents[i][first0][second0].find(
//                    std::make_pair(parentNodePair.first, parentNodePair.second)) !=
//                genuienParents[i][first0][second0].end())//no recomb
//
//            if (GetChildNode(i - 1, parentNodePair.first, GetAllele(i, first0)) == first0)//noJump for first
//                subSum += (1 - recRate) * recRate * baseProb / np1;
//
//            if (GetChildNode(i - 1, parentNodePair.second, GetAllele(i, second0)) == second0)//noJump for second
//                subSum += (1 - recRate) * recRate * baseProb / np2;

            if(parentNodePair.first == sampledGrandParent0 && parentNodePair.second == sampledGrandParent1)//no recomb
            {
                subSum += (1 - recRate) * (1 - recRate) * baseProb / (np1 * np2);
                tmpSampled+=4;
            }
            if (parentNodePair.first == sampledGrandParent0)//noJump for first
            {
                subSum += (1 - recRate) * recRate * baseProb / np1;
                tmpSampled+=2;
            }

            if (parentNodePair.second == sampledGrandParent1)//noJump for second
            {
                subSum += (1 - recRate) * recRate * baseProb / np2;
                tmpSampled+=1;
            }

            subSum += recRate * recRate * baseProb;
            if(tmpSampled > sampled) sampled=tmpSampled;

            fprintf(stderr, "sampling:(%d,%d)[%g,%g] to (%d,%d) with subSum:%g, recRate:%g, tmpSampled:%d\n",parentNodePair.first,parentNodePair.second,np1,np2,first0,second0,subSum,recRate,tmpSampled);
            if(genuienParents[i-1][parentNodePair.first][parentNodePair.second].size()>0)
            {
                for (auto grandParents:genuienParents[i - 1][parentNodePair.first][parentNodePair.second])
                {//all grandParents ending at parentNodePair.first and parentNodePair.second
                    sum += subSum * grandParents.second;//grandParents.second is the fwd value at site i-1, contributed by grandParents at site i-2
                    if (sum > choice) {
                        sampledParent0 = parentNodePair.first;
                        sampledParent1 = parentNodePair.second;
                        sampledGrandParent0 = grandParents.first.first;
                        sampledGrandParent1 = grandParents.first.second;
                        edgeProbFirst0 = GetEdgeProbAt(i - 2, sampledGrandParent0,
                                                       GetAllele(i - 1, sampledParent0));
                        edgeProbSecond0 = GetEdgeProbAt(i - 2, sampledGrandParent1,
                                                        GetAllele(i - 1, sampledParent1));
                        sampledFwd = grandParents.second * fwdValueSum[i - 1];
                        testSum = grandParents.second;
                        goto SAMPLE_BREAK;
                    }
                }
            } else
            {
                fprintf(stderr,"shouldn't happen here at marker %d\n",i-1);
            }
        }
        SAMPLE_BREAK:
//        sampledFwd = SumFwdValueFromOriginVec(genuienParents[i-1][sampledParent0][sampledParent1]) * fwdValueSum[i - 1];

        if (sum < choice)
            fprintf(stderr, "inside marker:%d\ttestSum:%g\tsum:%g\tchoice:%g\t(%d,%d)\ttotalValue:%g\tsampled:%d\tsubSum:%g\trecRate:%g\tbaseProb:%g\tfwdValueSum:%g\n",
                    i-1, testSum,sum, choice, sampledParent0, sampledParent1, sampledFwd, sampled, subSum, recRate, baseProb,fwdValueSum[i - 1]);
        availablePair.PrevMarker();
    }

    ImputeAlleles(1, sampledParent0, sampledParent1, rand, SampleIndex, sampledHaps);
    ImputeAlleles(0, sampledGrandParent0, sampledGrandParent1, rand, SampleIndex, sampledHaps);
    return 0;
}
#endif
//end of beagle variant version

//hybrid version

int PBWTHaplotyper::InitialFwdValues(int sampleIndex, FwdBwdLocalParameter &localParameter) {
    localParameter.isRec.assign(markers, false);
    localParameter.fwdValueSum.assign(markers, 0.f);
    std::unordered_map<int, float> dummy;
    dummy.reserve(HASH_RESERVE);
    localParameter.fwdValueNode1Sum.assign(markers, dummy);
    localParameter.fwdValueNode2Sum.assign(markers, dummy);

    localParameter.states = GetStateNumFrom(0);

    float prior = 1.f / (localParameter.states * localParameter.states);
    float tmpFwdValue(0.f), gl(0.f);

    for (StateIndex i = 0; i < localParameter.states; ++i) {
        for (StateIndex j = 0; j < localParameter.states; ++j) {
            int allele1 = GetAllele(0, i);
            int allele2 = GetAllele(0, j);
            gl = GetGL(sampleIndex, 0, allele1, allele2);
            if (gl > 1e-1) {
                tmpFwdValue = prior * gl;
                if (tmpFwdValue < UNDERFLOW_MIN) {
                    tmpFwdValue = UNDERFLOW_MIN;
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
                        "foward marker:%d\tfwdValueSum:%g\tpair:%d(%d,%d)=(%d,%d)\tcurrentFwd:%g\tnormalized currentFwd:%g\tfrom:(%d,%d)\n",
                        0, localParameter.fwdValueSum[0], kv.first,tmpNode1, tmpNode2,
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

int PBWTHaplotyper::ForwardAlgorithmRec(int sampleIndex, FwdBwdLocalParameter &localParameter) {
    float prevFwdValue(0.f);
    float tmpFwdValue(0.f), lowestFwd(0.f);
    float gl(0.f);

    int fitPair = 0;
    int numCrediablePair = 0;


    std::priority_queue<EdgePair, std::vector<EdgePair>, std::function<bool(const EdgePair &,
                                                                            const EdgePair &)> > EdgePairList(
            EdgePaircomparator);

    char allele1(0), allele2(0);
    StateIndex childNode1(0), childNode2(0), parentNode1(0), parentNode2(0);

    float recRate(0.f);
    float baseProb(0.f);

    bool reenter(false);
//    availablePair.ResetMarkerIndexAt(0);//TODO: use improved ibs methods maybe
//    availablePair.NextMarker();//marker 0 doesn't have parents
    for (int i = 1; i < markers; i++) {
        localParameter.fwdValueSum[i] = 0.f;
        numCrediablePair = 0;
        fitPair = 0;
        reenter = false;
        //all states at site i-1, parentNode1: hap1
        REENTRY:
        for (auto kv: localParameter.parentsNodeVec[i - 1]) {

            parentNode1 = localParameter.GetFirst(kv.first);//current round parents, last round children
            parentNode2 = localParameter.GetSecond(kv.first);
            //sum over all the parents of current parents that lead to a child pair
            prevFwdValue = localParameter.SumFwdValueFromOriginVec(localParameter.GetFwdVec(i - 1, kv.second));
            //i child of i-1 child state
            for (allele1 = 0; allele1 < 2; ++allele1) {
                childNode1 = GetChildNode(i - 1, parentNode1, allele1);
                if (childNode1 == -1) continue;
                //i child of i-1 child state
                for (allele2 = 0; allele2 < 2; ++allele2) {
                    childNode2 = GetChildNode(i - 1, parentNode2, allele2);
                    if (childNode2 == -1) continue;
                    gl = GetGL(sampleIndex, i, allele1, allele2);
                    if (gl > 1e-1 || i < 20) {
                        fitPair++;
                        tmpFwdValue = prevFwdValue * GetTransitionProb(i - 1, parentNode1, childNode1) *
                                      GetTransitionProb(i - 1, parentNode2, childNode2) * gl;
                        if (DEBUG)fprintf(stderr, "debug (%d,%d) prevFwdValue:%g\ttp1:%g\ttp2:%g\t%g\n", parentNode1, parentNode2,
                                prevFwdValue, GetTransitionProb(i - 1, parentNode1, childNode1),
                                GetTransitionProb(i - 1, parentNode2, childNode2), gl);
                        localParameter.fwdValueSum[i] += tmpFwdValue;
                        //from (parentNode1, parentNode2) to (childNode1, childNode2), childNodes are present in conditional graph, but parentNode are not necessarily present
//                        localParameter.parentsNodeVec[i][std::make_pair(childNode1, childNode2)][std::make_pair(parentNode1,
//                                                                                                 parentNode2)] = tmpFwdValue;
                        localParameter.FillParentsNodeVec(i, childNode1, childNode2, parentNode1, parentNode2,
                                                          tmpFwdValue);
                    } else if (reenter) {
                        tmpFwdValue = prevFwdValue *
                                      GetTransitionProb(i - 1, parentNode1, childNode1) *
                                      GetTransitionProb(i - 1, parentNode2, childNode2) *
                                      gl;//i fwdValueSum
                        if (tmpFwdValue < UNDERFLOW_MIN && prevFwdValue > 0) {
                            tmpFwdValue = UNDERFLOW_MIN;
                        }

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

        if (fitPair == 0 && !reenter)//process orphan nodes
        {
            fprintf(stderr, "[%s] sample %d recombined at marker %d\n", __FUNCTION__, sampleIndex, i);
            recRate = GetRecombRate(i - 1);//start from marker 1 but store at index 0
            localParameter.isRec[i - 1] = true;
            localParameter.states = GetStateNumFrom(i);
            for (childNode1 = 0; childNode1 < localParameter.states; ++childNode1) {//dest nodeA
                for (childNode2 = 0; childNode2 < localParameter.states; ++childNode2) {//dest nodeB
                    gl = GetGL(sampleIndex, i, GetAllele(i, childNode1),
                               GetAllele(i, childNode2));//genotype of each childNode pair
                    if (gl < 1e-1) continue;
                    for (auto tmpParentNode1:GetParentNodes(i, childNode1)) {
                        for (auto tmpParentNode2:GetParentNodes(i, childNode2)) {
                            if (GetGL(sampleIndex, i - 1, GetAllele(i - 1, tmpParentNode1),
                                      GetAllele(i - 1, tmpParentNode2)) < 1e-1)
                                continue;//genotype of all possible parentNode pair of childNode
                            fitPair++;
                            if (localParameter.parentsNodeVec[i - 1].find(
                                    localParameter.MakePair(tmpParentNode1, tmpParentNode2)) !=
                                localParameter.parentsNodeVec[i - 1].end()) {
                                prevFwdValue = localParameter.SumFwdValueFromOriginVec(
                                        localParameter.GetFwdVec(i - 1, tmpParentNode1, tmpParentNode2));
//                                        (localParameter.parentsNodeVec[i - 1][localParameter.MakePair(tmpParentNode1,
//                                                                             tmpParentNode2)]);//sum over all the parents of current parents that lead to a child pair
                            } else
                                prevFwdValue = 0.f;

                            baseProb = GetTransitionProb(i - 1, tmpParentNode1, childNode1) *
                                       GetTransitionProb(i - 1, tmpParentNode2, childNode2) * gl;

                            tmpFwdValue = prevFwdValue * baseProb * (1 - recRate) * (1 - recRate);
//                        fprintf(stderr,"(%d,%d) to (%d,%d) tmp output1 tmpFwdValue:%g\tprevFwdValue:%g\tbaseProb:%g\n",tmpParentNode1,tmpParentNode2,childNode1,childNode2,tmpFwdValue,prevFwdValue,baseProb);
                            if (localParameter.fwdValueNode1Sum[i - 1].find(tmpParentNode1) !=
                                localParameter.fwdValueNode1Sum[i - 1].end())
                                tmpFwdValue += localParameter.fwdValueNode1Sum[i - 1][tmpParentNode1] *
                                               GetHapProbAt(i - 1, tmpParentNode2) *
                                               baseProb * (1 - recRate) * recRate;
//                        fprintf(stderr,"tmp output2 tmpFwdValue:%g\tfwdValueNode1Sum[i - 1][tmpParentNode1]:%g\n",tmpFwdValue,localParameter.fwdValueNode1Sum[i - 1][tmpParentNode1]);
                            if (localParameter.fwdValueNode2Sum[i - 1].find(tmpParentNode2) !=
                                localParameter.fwdValueNode2Sum[i - 1].end())
                                tmpFwdValue += localParameter.fwdValueNode2Sum[i - 1][tmpParentNode2] *
                                               GetHapProbAt(i - 1, tmpParentNode1) *
                                               baseProb * (1 - recRate) * recRate;
//                        fprintf(stderr,"tmp output3 tmpFwdValue:%g\tfwdValueNode2Sum[i - 1][tmpParentNode2]:%g\n",tmpFwdValue,localParameter.fwdValueNode2Sum[i - 1][tmpParentNode2]);
                            tmpFwdValue += GetHapProbAt(i - 1, tmpParentNode2) *
                                           GetHapProbAt(i - 1, tmpParentNode1) * baseProb * recRate * recRate;
//                        fprintf(stderr,"tmp output4 tmpFwdValue:%g\n",tmpFwdValue);

                            if (numCrediablePair < 50) {
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


        while (not EdgePairList.empty()) {
            EdgePair tmpEdgePair = EdgePairList.top();
            EdgePairList.pop();
            if (tmpEdgePair.fwd > 0.f) {
//                localParameter.parentsNodeVec[i][std::make_pair(tmpEdgePair.childNode1, tmpEdgePair.childNode2)][std::make_pair(
//                        tmpEdgePair.parentNode1,
//                        tmpEdgePair.parentNode2)] = tmpEdgePair.fwd;
                localParameter.FillParentsNodeVec(i, tmpEdgePair.childNode1, tmpEdgePair.childNode2,
                                                  tmpEdgePair.parentNode1, tmpEdgePair.parentNode2, tmpEdgePair.fwd);
                localParameter.fwdValueSum[i] += tmpEdgePair.fwd;
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
                    tmpFwd /= localParameter.fwdValueSum[i];
                    tmp += tmpFwd;
                    int first = localParameter.GetFirst(localParameter.GetSourceVec(i, kv.second).at(j));
                    int second = localParameter.GetSecond(localParameter.GetSourceVec(i, kv.second).at(j));
                    fprintf(stderr,
                            "foward marker:%d\tfwdValueSum:%g\tpair:(%d,%d)=(%d,%d)\tcurrentFwd:%g\tnormalized currentFwd:%g\tfrom:(%d,%d)\n",
                            i, localParameter.fwdValueSum[i], tmpNode1, tmpNode2,
                            GetAllele(i, tmpNode1), GetAllele(i, tmpNode2),
                            localParameter.GetFwdVec(i, kv.second).at(j), tmpFwd, first, second);
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
    return 0;
}

int PBWTHaplotyper::BackwardAlgorithmRec(int sampleIndex, char **sampledHaps,
                                         FwdBwdLocalParameter &localParameter) {

    double choice(0.);
    double sum(0.), subSum(0.);
    float gl0(0.f);
    float baseProb(0.f), transProb0(0.f), transProb1(0.f);
    float np1(0.f), np2(0.f);

    int first0(0), second0(0);
    float edgeProb0(0.f), edgeProb1(0.f);

    StateIndex sampledGrandParent0(0), sampledGrandParent1(0);

    StateIndex sampledParent0(0), sampledChild0(0);
    StateIndex sampledParent1(0), sampledChild1(0);
    double sampledFwd(0.f);

    float recRate(0.f);
    int sampled = 0;

//    choice = rand->Uniform(0, 1);
    int count = 0;
    for (auto iter = localParameter.parentsNodeVec[markers - 1].begin();
         iter != localParameter.parentsNodeVec[markers - 1].end(); ++iter) {
        count++;
        sampledFwd = localParameter.SumFwdValueFromOriginVec(localParameter.GetFwdVec(markers - 1, iter->second));
//            fprintf(stderr,"(%d,%d)\tvalue:%g\tcount:%d\n",iter->first.first,iter->first.second,sampledFwd,count);

    }

    for (auto kv: localParameter.parentsNodeVec[markers - 1]) {
        int destIndex = kv.second;
        for (int sourceIndex = 0;
             sourceIndex < localParameter.GetSourceVec(markers - 1, destIndex).size(); ++sourceIndex) {
            sum += localParameter.GetFwd(markers - 1, destIndex, sourceIndex);
//            fprintf(stderr, "sampling last marker:(%d,%d) with sum:%g, recRate:%g\n",
//                    kv.first.first, kv.first.second, sum, recRate);
            if (sum > choice) {
                sampledChild0 = localParameter.GetFirst(kv.first);
                sampledChild1 = localParameter.GetSecond(kv.first);
                sampledParent0 = localParameter.GetFirst(localParameter.GetSource(markers - 1, destIndex, sourceIndex));
                sampledParent1 = localParameter.GetSecond(
                        localParameter.GetSource(markers - 1, destIndex, sourceIndex));
                edgeProb0 = GetEdgeProbAt(markers - 2, sampledParent0, GetAllele(markers - 1, sampledChild0));
                edgeProb1 = GetEdgeProbAt(markers - 2, sampledParent1, GetAllele(markers - 1, sampledChild1));
//                if(!isRec[markers - 1])
//                    sampledFwd = parent.second * fwdValueSum[markers - 1]/(GetTransitionProb(markers - 2, sampledParent0, first0) *
//                                                                           GetTransitionProb(markers - 2, sampledParent1, second0) *
//                                                                           GetGL(individuals - 1, markers - 1, GetAllele(markers - 1, first0), GetAllele(markers - 1, second0)));//tranProb not available if recombined
//                else
                sampledFwd = localParameter.GetFwd(markers - 1, destIndex, sourceIndex) *
                             localParameter.fwdValueSum[markers - 1];//tranProb not available if recombined

                goto INIT_SAMPLE_BREAK;
            }
        }
    }
    INIT_SAMPLE_BREAK:

//    ImputeAlleles(markers - 1, sampledChild0, sampledChild1, rand, sampleIndex, sampledHaps);

    float testSum(0.f);

    for (int i = markers - 2; i > 0; --i) {

//        if (i == 4102 || i == 4103 || i == 4104) {
//        fprintf(stderr,
//                "marker:%d\ttestSum:%g\tsum:%g\tchoice:%g\tprev(%d,%d) to (%d,%d)\tsampledFwd:%g\tavailablePair.Size:%d\trecRate:%g\tgl:%g\tprobEdge0:%g\tprobEdge1:%g\tfwdValueSum:%g\n",
//                i, testSum, sum,
//                choice, sampledParent0, sampledParent1, sampledChild0, sampledChild1,  sampledFwd,
//                availablePair.Size(),
//                recRate, gl0, edgeProb0, edgeProb1, fwdValueSum[i]);
////
//
//        for (auto kv:parentsNodeVec[i][std::make_pair(sampledParent0, sampledParent1)]) {
//            fprintf(stderr, "marker:%d\tpair:(%d,%d)\ttestSum:%g\n", i, kv.first.first, kv.first.second, kv.second);
//        }
//        }

//        fprintf(stderr,"marker:%d\tsum:%g\tchoice:%g\t(%d,%d)\tvalue:%g\n",i,sum,choice,sampledParent0,sampledParent1,sampledFwd);

        sum = 0.f;
        testSum = 0.f;

//        ImputeAlleles(i, sampledParent0, sampledParent1, rand, sampleIndex, sampledHaps);
//        choice = rand->Uniform(0, sampledFwd);

        gl0 = GetGL(sampleIndex, i + 1, GetAllele(i + 1, sampledChild0), GetAllele(i + 1, sampledChild1));

        if (!localParameter.isRec[i]) {//normal hap match without rec
            double choice0 = choice / (GetTransitionProb(i, sampledParent0, sampledChild0) *
                                       GetTransitionProb(i, sampledParent1, sampledChild1) * gl0);
//            for (auto grandParents:localParameter.parentsNodeVec[i][localParameter.MakePair(sampledParent0, sampledParent1)])//we are actually sampling i-1's states
            int destIndex = localParameter.GetDestIndex(i, sampledParent0, sampledParent1);
            for (int grandParentsIndex = 0; grandParentsIndex < localParameter.GetSourceVec(i, sampledParent0,
                                                                                            sampledParent1).size(); ++grandParentsIndex) {

                sum += localParameter.GetFwd(i, destIndex, grandParentsIndex);//kv.second is normalized
                if (sum > choice0) {
                    sampledChild0 = sampledParent0;
                    sampledChild1 = sampledParent1;
                    sampledParent0 = localParameter.GetFirst(localParameter.GetSource(i, destIndex, grandParentsIndex));
                    sampledParent1 = localParameter.GetSecond(
                            localParameter.GetSource(i, destIndex, grandParentsIndex));
                    edgeProb0 = GetEdgeProbAt(i - 1, sampledParent0, GetAllele(i, sampledChild0));
                    edgeProb1 = GetEdgeProbAt(i - 1, sampledParent1, GetAllele(i, sampledChild1));
                    sampledFwd = localParameter.GetFwd(i, destIndex, grandParentsIndex) * localParameter.fwdValueSum[i];
                    break;
                }
            }
        } else {

            baseProb = edgeProb0 * edgeProb1 * gl0;//site i
            recRate = GetRecombRate(i);//start from marker 1 but store at index 0
            for (auto kv:localParameter.parentsNodeVec[i]) {// node pair at site i - 1, sampled (first0,second0) for i

                uint64_t parentNodePair = kv.first;
                np1 = GetHapProbAt(i, localParameter.GetFirst(parentNodePair));
                np2 = GetHapProbAt(i, localParameter.GetSecond(parentNodePair));
                subSum = 0.;

                if (localParameter.GetFirst(parentNodePair) == sampledParent0 &&
                    localParameter.GetSecond(parentNodePair) == sampledParent1)//no recomb
                    subSum += (1 - recRate) * (1 - recRate) * baseProb / (np1 * np2);

                if (localParameter.GetFirst(parentNodePair) == sampledParent0)//noJump for first
                    subSum += (1 - recRate) * recRate * baseProb / np1;

                if (localParameter.GetSecond(parentNodePair) == sampledParent1)//noJump for second
                    subSum += (1 - recRate) * recRate * baseProb / np2;

                subSum += recRate * recRate * baseProb;

                int destIndex = kv.second;
                for (int grandParentsIndex = 0; grandParentsIndex < localParameter.GetSourceVec(i,
                                                                                                destIndex).size(); ++grandParentsIndex) {//all grandParents ending at parentNodePair.first and parentNodePair.second
//                    sum += subSum *
//                           grandParents.second;//grandParents.second is the fwd value at site i-1, contributed by grandParents at site i-2
                    sum += subSum *
                           localParameter.GetFwd(i, destIndex,
                                                 grandParentsIndex);//grandParents.second is the fwd value at site i-1, contributed by grandParents at site i-2
                    if (sum > choice) {
                        sampledChild0 = localParameter.GetFirst(parentNodePair);//previous correct parents rejuvenate
                        sampledChild1 = localParameter.GetSecond(parentNodePair);
//                        ImputeAlleles(i, sampledChild0, sampledChild1, rand, sampleIndex,
//                                      sampledHaps);//override no-rec Imputation
                        sampledParent0 = localParameter.GetFirst(
                                localParameter.GetSource(i, destIndex, grandParentsIndex));
                        sampledParent1 = localParameter.GetSecond(
                                localParameter.GetSource(i, destIndex, grandParentsIndex));
                        edgeProb0 = GetEdgeProbAt(i - 1, sampledParent0, GetAllele(i, sampledChild0));
                        edgeProb1 = GetEdgeProbAt(i - 1, sampledParent1, GetAllele(i, sampledChild1));
                        sampledFwd =
                                localParameter.GetFwd(i, destIndex, grandParentsIndex) * localParameter.fwdValueSum[i];
                        testSum = localParameter.GetFwd(i, destIndex, grandParentsIndex);
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
            throw samplingException;
        }
    }

//    ImputeAlleles(0, sampledParent0, sampledParent1, rand, sampleIndex, sampledHaps);
    return 0;
}


int PBWTHaplotyper::BackwardSamplingRec(Random *rand, int sampleIndex, char **sampledHaps,
                                        FwdBwdLocalParameter &localParameter) {

    double choice(0.);
    double sum(0.), subSum(0.);
    float gl0(0.f);
    float baseProb(0.f), transProb0(0.f), transProb1(0.f);
    float np1(0.f), np2(0.f);

    int first0(0), second0(0);
    float edgeProb0(0.f), edgeProb1(0.f);

    StateIndex sampledGrandParent0(0), sampledGrandParent1(0);

    StateIndex sampledParent0(0), sampledChild0(0);
    StateIndex sampledParent1(0), sampledChild1(0);
    double sampledFwd(0.f);

    float recRate(0.f);
    int sampled = 0;

    choice = rand->Uniform(0, 1);
    int count = 0;
    for (auto iter = localParameter.parentsNodeVec[markers - 1].begin();
         iter != localParameter.parentsNodeVec[markers - 1].end(); ++iter) {
        count++;
        sampledFwd = localParameter.SumFwdValueFromOriginVec(localParameter.GetFwdVec(markers - 1, iter->second));
        if (DEBUG)
            fprintf(stderr, "(%d,%d)\tvalue:%g\tcount:%d\n", localParameter.GetFirst(iter->first),
                    localParameter.GetSecond(iter->first), sampledFwd, count);

    }

    for (auto kv: localParameter.parentsNodeVec[markers - 1]) {
        int destIndex = kv.second;
        for (int sourceIndex = 0;
             sourceIndex < localParameter.GetSourceVec(markers - 1, destIndex).size(); ++sourceIndex) {
            sum += localParameter.GetFwd(markers - 1, destIndex, sourceIndex);
            if (DEBUG)
                fprintf(stderr, "sampling last marker:(%d,%d) with sum:%g, recRate:%g\n",
                        localParameter.GetFirst(kv.first), localParameter.GetSecond(kv.first), sum, recRate);
            if (sum > choice) {
                sampledChild0 = localParameter.GetFirst(kv.first);
                sampledChild1 = localParameter.GetSecond(kv.first);
                sampledParent0 = localParameter.GetFirst(localParameter.GetSource(markers - 1, destIndex, sourceIndex));
                sampledParent1 = localParameter.GetSecond(
                        localParameter.GetSource(markers - 1, destIndex, sourceIndex));
                edgeProb0 = GetEdgeProbAt(markers - 2, sampledParent0, GetAllele(markers - 1, sampledChild0));
                edgeProb1 = GetEdgeProbAt(markers - 2, sampledParent1, GetAllele(markers - 1, sampledChild1));
//                if(!isRec[markers - 1])
//                    sampledFwd = parent.second * fwdValueSum[markers - 1]/(GetTransitionProb(markers - 2, sampledParent0, first0) *
//                                                                           GetTransitionProb(markers - 2, sampledParent1, second0) *
//                                                                           GetGL(individuals - 1, markers - 1, GetAllele(markers - 1, first0), GetAllele(markers - 1, second0)));//tranProb not available if recombined
//                else
                sampledFwd = localParameter.GetFwd(markers - 1, destIndex, sourceIndex) *
                             localParameter.fwdValueSum[markers - 1];//tranProb not available if recombined

                goto INIT_SAMPLE_BREAK;
            }
        }
    }
    INIT_SAMPLE_BREAK:

    ImputeAlleles(markers - 1, sampledChild0, sampledChild1, rand, sampleIndex, sampledHaps);

    float testSum(0.f);

    for (int i = markers - 2; i > 0; --i) {

//        if (i == 4102 || i == 4103 || i == 4104) {
//        fprintf(stderr,
//                "marker:%d\ttestSum:%g\tsum:%g\tchoice:%g\tprev(%d,%d) to (%d,%d)\tsampledFwd:%g\tavailablePair.Size:%d\trecRate:%g\tgl:%g\tprobEdge0:%g\tprobEdge1:%g\tfwdValueSum:%g\n",
//                i, testSum, sum,
//                choice, sampledParent0, sampledParent1, sampledChild0, sampledChild1,  sampledFwd,
//                availablePair.Size(),
//                recRate, gl0, edgeProb0, edgeProb1, fwdValueSum[i]);
////
//
//        for (auto kv:parentsNodeVec[i][std::make_pair(sampledParent0, sampledParent1)]) {
//            fprintf(stderr, "marker:%d\tpair:(%d,%d)\ttestSum:%g\n", i, kv.first.first, kv.first.second, kv.second);
//        }
//        }
        if (DEBUG)
            fprintf(stderr, "marker:%d\tsum:%g\tchoice:%g\t(%d,%d)\tvalue:%g\tgl:%g\n", i, sum, choice, sampledParent0,
                    sampledParent1, sampledFwd, gl0);

        sum = 0.f;
        testSum = 0.f;

        ImputeAlleles(i, sampledParent0, sampledParent1, rand, sampleIndex, sampledHaps);
        choice = rand->Uniform(0, sampledFwd);

        gl0 = GetGL(sampleIndex, i + 1, GetAllele(i + 1, sampledChild0), GetAllele(i + 1, sampledChild1));

        if (!localParameter.isRec[i]) {//normal hap match without rec
            double choice0 = choice / (GetTransitionProb(i, sampledParent0, sampledChild0) *
                                       GetTransitionProb(i, sampledParent1, sampledChild1) * gl0);
//            for (auto grandParents:localParameter.parentsNodeVec[i][localParameter.MakePair(sampledParent0, sampledParent1)])//we are actually sampling i-1's states
            int destIndex = localParameter.GetDestIndex(i, sampledParent0, sampledParent1);
            for (int grandParentsIndex = 0; grandParentsIndex < localParameter.GetSourceVec(i, sampledParent0,
                                                                                            sampledParent1).size(); ++grandParentsIndex) {

                sum += localParameter.GetFwd(i, destIndex, grandParentsIndex);//kv.second is normalized
                if (sum > choice0) {
                    sampledChild0 = sampledParent0;
                    sampledChild1 = sampledParent1;
                    sampledParent0 = localParameter.GetFirst(localParameter.GetSource(i, destIndex, grandParentsIndex));
                    sampledParent1 = localParameter.GetSecond(
                            localParameter.GetSource(i, destIndex, grandParentsIndex));
                    edgeProb0 = GetEdgeProbAt(i - 1, sampledParent0, GetAllele(i, sampledChild0));
                    edgeProb1 = GetEdgeProbAt(i - 1, sampledParent1, GetAllele(i, sampledChild1));
                    sampledFwd = localParameter.GetFwd(i, destIndex, grandParentsIndex) * localParameter.fwdValueSum[i];
                    break;
                }
            }
        } else {

            baseProb = edgeProb0 * edgeProb1 * gl0;//site i
            recRate = GetRecombRate(i);//start from marker 1 but store at index 0
            for (auto kv:localParameter.parentsNodeVec[i]) {// node pair at site i - 1, sampled (first0,second0) for i

                uint64_t parentNodePair = kv.first;
                np1 = GetHapProbAt(i, localParameter.GetFirst(parentNodePair));
                np2 = GetHapProbAt(i, localParameter.GetSecond(parentNodePair));
                subSum = 0.;

                if (localParameter.GetFirst(parentNodePair) == sampledParent0 &&
                    localParameter.GetSecond(parentNodePair) == sampledParent1)//no recomb
                    subSum += (1 - recRate) * (1 - recRate) * baseProb / (np1 * np2);

                if (localParameter.GetFirst(parentNodePair) == sampledParent0)//noJump for first
                    subSum += (1 - recRate) * recRate * baseProb / np1;

                if (localParameter.GetSecond(parentNodePair) == sampledParent1)//noJump for second
                    subSum += (1 - recRate) * recRate * baseProb / np2;

                subSum += recRate * recRate * baseProb;

                int destIndex = kv.second;
                for (int grandParentsIndex = 0; grandParentsIndex < localParameter.GetSourceVec(i,
                                                                                                destIndex).size(); ++grandParentsIndex) {//all grandParents ending at parentNodePair.first and parentNodePair.second
//                    sum += subSum *
//                           grandParents.second;//grandParents.second is the fwd value at site i-1, contributed by grandParents at site i-2
                    sum += subSum *
                           localParameter.GetFwd(i, destIndex,
                                                 grandParentsIndex);//grandParents.second is the fwd value at site i-1, contributed by grandParents at site i-2
                    if (sum > choice) {
                        sampledChild0 = localParameter.GetFirst(parentNodePair);//previous correct parents rejuvenate
                        sampledChild1 = localParameter.GetSecond(parentNodePair);
                        ImputeAlleles(i, sampledChild0, sampledChild1, rand, sampleIndex,
                                      sampledHaps);//override no-rec Imputation
                        sampledParent0 = localParameter.GetFirst(
                                localParameter.GetSource(i, destIndex, grandParentsIndex));
                        sampledParent1 = localParameter.GetSecond(
                                localParameter.GetSource(i, destIndex, grandParentsIndex));
                        edgeProb0 = GetEdgeProbAt(i - 1, sampledParent0, GetAllele(i, sampledChild0));
                        edgeProb1 = GetEdgeProbAt(i - 1, sampledParent1, GetAllele(i, sampledChild1));
                        sampledFwd =
                                localParameter.GetFwd(i, destIndex, grandParentsIndex) * localParameter.fwdValueSum[i];
                        testSum = localParameter.GetFwd(i, destIndex, grandParentsIndex);
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
            throw samplingException;
        }
    }

    ImputeAlleles(0, sampledParent0, sampledParent1, rand, sampleIndex, sampledHaps);
    return 0;
}

int PBWTHaplotyper::LocalForwadBackWard(int sampleIndex) {
    FwdBwdLocalParameter localParameter(individuals, markers);
    InitialFwdValues(sampleIndex, localParameter);
    try {
        ForwardAlgorithmRec(sampleIndex, localParameter);
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

//end of hybrid version

// end with recombination

#ifdef HETERSITE
int PBWTHaplotyper::ExtractHeterSites(int individualToProcess) {//apply after swap individualToProcess to the back

    if (Wrapper != nullptr) {
        delete Wrapper;
        Wrapper = nullptr;
    }
    absoluteIndexToRelative.clear();
    relativeIndexToAbsolute.clear();
    tmpMarkers = 0;

    if (onlyHeterSite) {

        std::vector<bool> HeterIndex(markers, false);
        for (int i = 0; i < markers; ++i) {
//            fprintf(stderr,"through heter partat marker %d: %d\t%d\t%d\n",i,genotypes[individualToProcess][i*3],genotypes[individualToProcess][i*3+1],genotypes[individualToProcess][i*3+2]);
            if (genotypes[individualToProcess][i * 3 + 1] < genotypes[individualToProcess][i * 3] &&
                genotypes[individualToProcess][i * 3 + 1] < genotypes[individualToProcess][i * 3 + 2]) {
                HeterIndex[i] = true;
                absoluteIndexToRelative[i] = tmpMarkers++;
                relativeIndexToAbsolute.push_back(i);
            }
        }

        for (int i = 0; i < tmpMarkers; ++i) {
            int markerAbsoluteNow = relativeIndexToAbsolute[i];
            if (!HeterIndex[markerAbsoluteNow]) {
                continue;
            }//homo
            for (int j = 0; j < individuals; ++j) {

                tmpHaps[2 * j][i] = haplotypes[2 * j][markerAbsoluteNow];
                tmpHaps[2 * j + 1][i] = haplotypes[2 * j + 1][markerAbsoluteNow];
                tmpGeno[j][i * 3] = genotypes[j][markerAbsoluteNow * 3];
                tmpGeno[j][i * 3 + 1] = genotypes[j][markerAbsoluteNow * 3 + 1];
                tmpGeno[j][i * 3 + 2] = genotypes[j][markerAbsoluteNow * 3 + 2];
            }
            for (int j = 0; j < (individuals - phased) * nSampleCopy; ++j) {
                tmpHaps[individuals * 2 + 2 * j][i] = sampledHaps[2 * j][markerAbsoluteNow];
                tmpHaps[individuals * 2 + 2 * j + 1][i] = sampledHaps[2 * j + 1][markerAbsoluteNow];
            }
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    tmpPenetrance[i * 9 + j * 3 + k] = Penetrance(markerAbsoluteNow, j, k);
                }
            }
        }
        SwapTempHaps();
    } else {
        tmpMarkers = markers;
    }

    if (tmpMarkers == 0) {
        fprintf(stderr, "found 0 markers available...abort!\n");
        abort();
    }
    Wrapper = new PBWTWrapper(2 * individuals + (individuals - phased) * nSampleCopy * 2, tmpMarkers, PvalueMatrix, 0);

    return 0;
}

int PBWTHaplotyper::FillHeterSitesBack(int individualToProcess) {

    if (onlyHeterSite) {
        SwapTempHaps();
        int markerAbsoluteNow = 0;
        for (int i = 0; i < tmpMarkers; ++i) {
            markerAbsoluteNow = relativeIndexToAbsolute[i];
            haplotypes[2 * individualToProcess][markerAbsoluteNow] = tmpHaps[2 * individualToProcess][i];
            haplotypes[2 * individualToProcess + 1][markerAbsoluteNow] = tmpHaps[2 * individualToProcess + 1][i];
        }
    }
    for (auto &parents: genuienParents) {
        parents.clear();
    }
    return 0;
}
#endif

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
    if (runningModel & PHASE) {//enter from phase by ref
        thetas = new float[markers - 1];
        for (int i = 0; i < markers - 1; i++)
            thetas[i] = 0.01;
        genotypes = AllocateCharMatrix(individuals, markers * 3);
        penetrances = new float[markers * 9];
        error_models = new Errors[markers];
        freq1s = new double[markers];
        crossovers = new int[markers - 1];
    } else//enter from graph construct
    {
        thetas = nullptr;
        genotypes = nullptr;
        penetrances = nullptr;
        error_models = nullptr;
        freq1s = nullptr;
        crossovers = nullptr;
    }
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
