//
// Created by Fan Zhang on 8/6/15.
//

#ifndef PLUTO_PBWTHAPLOTYPER_H
#define PLUTO_PBWTHAPLOTYPER_H

#include <unordered_set>
#include "ShotgunHaplotyper.h"
#include "PBWTWrapper.h"
#include "GeneticDistanceMap.h"
#include "GzipFileType.h"
#include "FwdBwdLocalParameter.h"

#define INDEX 1
#define PHASE 2
#define ITERATIVE 4

class PBWTHaplotyper : public ShotgunHaplotyper {
public:
    std::string loadGraph = "Empty";//indicate if build graph from previously built graph
    std::string graphFilePrefix = "Empty";
    bool geneticMapAvailable = false;
    int runningModel = 0;
    int prefixLength = 0;
//    int phasedForByRef = 0;
    int nThread = 1;
    double genoThresh;


    GeneticDistanceMap GDMap;
/*!PBWTHaplotyper
 * Consturctor
 * @param nhaps number of haplotypes
 * @param nsnps number of snp markers
 */
    PBWTHaplotyper(int nhaps, int nsnps);

    PBWTHaplotyper();

    ~PBWTHaplotyper();

    /*!InitialHaplotypeByConsensus
     * copy consencus haplotypes into haplotypes char[][], and reset consensus
     * @param consensus
     */
    void InitialHaplotypeByConsensus(ConsensusBuilder &consensus) {
        for (int i = 0; i < (individuals - phased) * 2; ++i) {
            for (int j = 0; j < markers; ++j) {
                haplotypes[i][j] = consensus.consensus[i][j];
            }
        }
        consensus.stored = 0;
    }

    void ConstructGraph();

    void InitialSampleCopy(Random *rand);

    void SwapIndividuals(int a, int b);

    void PrepareRefSetPBWTWrapper();


    int LoopThroughChromosomesHighPrecision();

    int LoopThroughChromosomesSingleRound();

    int LoopThroughChromosomesViaPBWTWithHeterOnly();

    /*!LoopThroughChromosomesGenotyping
     * Phasing Algorithm
     * @param ped pedigree information
     * @return
     */
    int LoopThroughChromosomesGenotyping(Pedigree &ped);

    int LoopThroughChromosomesPhasing(Pedigree &ped);

    int CorrectGenotype();

    virtual void RandomSetup(Random *rand = NULL);

    virtual void ImputeAlleles(int marker, int state1, int state2, Random *rand, int currentIndividual, char **haps);

    virtual void ImputeAllelesRaw(int marker, int state1, int state2, Random *rand, int currentIndividual, char **haps);

    virtual void ImputeAllele(int haplotype, int marker, int state, char **haps);

    virtual void FillPath(int haplotype, int fromMarker, int toMarker, int state, char **haps);


    void SetUseRev(bool useOrNot) { isRev = useOrNot; }

    bool ReverseInput();

    void SetOnlyGT(bool onlyOrNot) { onlyGT = onlyOrNot; }

    bool GetOnlyGT() { return onlyGT; }

    void FillGenotypeLikelihood(long Phred11, long Phred12, long Phred22, int idv, int genoIndex)
    {
        long minValue = 127;
        long minPhred = std::min(Phred11, std::min(Phred12, Phred22));
        if( (Phred11 - minPhred) > 30) Phred11 = minValue;
        else if((Phred12 - minPhred) > 30) Phred12 = minValue;
        else if((Phred22 - minPhred) > 30) Phred22 = minValue;

        if(Phred11 > minValue) Phred11 =minValue;
        if(Phred12 > minValue) Phred12 =minValue;
        if(Phred22 > minValue) Phred22 =minValue;

//        long minPhred = std::min(Phred11, std::min(Phred12, Phred22));
//        if (minPhred == Phred12) {
//            if ((Phred11 - minPhred) > 30) Phred11 = 127;
//            if ((Phred22 - minPhred) > 30) Phred22 = 127;
//        } else {
//            if (Phred11 == minPhred) Phred22 = 127;
//            if (Phred22 == minPhred) Phred11 = 127;
//        }
        genotypes[idv][genoIndex] = static_cast<char>(Phred11);
        genotypes[idv][genoIndex + 1] = static_cast<char>(Phred12);
        genotypes[idv][genoIndex + 2] = static_cast<char>(Phred22);
    }

    int CalculatePosteriorGL(int markerIndex, int sampleIndex, double *posterior)//for RandomSetup
    {
        double prior[3];
        prior[0] = freq1s[markerIndex] * freq1s[markerIndex];
        prior[1] = 2.0 * freq1s[markerIndex] * (1.0 - freq1s[markerIndex]);
        prior[2] = (1.0 - freq1s[markerIndex]) * (1.0 - freq1s[markerIndex]);
        posterior[0] = prior[0] * phred2prob[genotypes[sampleIndex][markerIndex*3]];
        posterior[1] = prior[1] * phred2prob[genotypes[sampleIndex][markerIndex*3+1]];
        posterior[2] = prior[2] * phred2prob[genotypes[sampleIndex][markerIndex*3+2]];
        double sum = posterior[0] + posterior[1] + posterior[2];

        posterior[0] /= sum;
        posterior[1] /= sum;
        posterior[2] /= sum;

        return 0;
    }

//    double GetEmissionProb(int markerIndex, int sampleIndex, StateIndex A, StateIndex B)
//    {
//        int geno = GetAllele(markerIndex, A) + GetAllele(markerIndex, B);
//        double emmit[3]={0,0,0};
//        emmit[0] = Penetrance(markerIndex, geno, 0) * phred2prob[genotypes[sampleIndex][markerIndex*3]];
//        emmit[1] = Penetrance(markerIndex, geno, 1) * phred2prob[genotypes[sampleIndex][markerIndex*3+1]];
//        emmit[2] = Penetrance(markerIndex, geno, 2) * phred2prob[genotypes[sampleIndex][markerIndex*3+2]];
//        double sum = emmit[0] + emmit[1] + emmit[2];
//
//        return emmit[geno]/sum;
//    }

    int UpdateFreq() {
        if(runningModel & ITERATIVE) {
            for (int i = 0; i < markers; ++i) {
                float cnt = 0;
                for (int j = 0; j < individuals; ++j) {
                    cnt += haplotypes[j * 2][i] + haplotypes[j * 2 + 1][i];
                }
                freq1s[i] = cnt / (2 * individuals);
            }
        }
        return 0;
    }

    //genotype related


    float GetGL(int individual, int marker, char allele1, char allele2) {
        return (float)phred2prob[(size_t) genotypes[individual][3 * marker + allele1 + allele2]];
    }

    void SetGP(int individual, int marker, char allele1, char allele2, float p) {
        p = p > std::numeric_limits<float>::min() ? p : std::numeric_limits<float>::min();
//        int tmpPL = static_cast<int>(log10(p));
//        genoProbs[individual][3 * marker + allele1 + allele2] = tmpPL < 127 ? tmpPL:127;
        genoProbs[individual][3 * marker + allele1 + allele2] = p;
    }

    void ResetGL(int individual, int marker) {
        genotypes[individual][3 * marker] = 127;
        genotypes[individual][3 * marker + 1] = 127;
        genotypes[individual][3 * marker + 2] = 127;
    }

    void NormalizeGP(int individual, int marker, float *tmpGP) {
        double sum =  tmpGP[0] +  tmpGP[1] +  tmpGP[2];
        tmpGP[0] /= sum;
        tmpGP[1] /= sum;
        tmpGP[2] /= sum;
    }

    void FillGL(int individual, int marker, float * tmpGP)
    {

        NormalizeGP(individual, marker, tmpGP);
        SetGP(individual, marker, 0, 0, tmpGP[0]);
        SetGP(individual, marker, 0, 1, tmpGP[1]);
        SetGP(individual, marker, 1, 1, tmpGP[2]);
    }
    //end of genotype related

    int GetUnphasedNum()
    {
        int nTarget;
        if(runningModel & INDEX) nTarget = individuals;
        else if(runningModel & PHASE) nTarget = individuals - phased;
        else if(runningModel & ITERATIVE) nTarget = individuals - phased;
        else
        {
            fprintf(stderr,"Should not reach this code!\n");
            exit(EXIT_FAILURE);
        }

        return nTarget;
    }

    //Memory management functions
    bool SetErrorAndTheta(std::vector<float> &holderError, std::vector<float> &holderTheta);

    bool AllocateMemory(int nIndividuals, int nMarkers);

    //inline section
    inline float GetTransitionProb(int site, StateIndex from, StateIndex to) {
        return Wrapper->Graph.GetTransitionProb(site, from, to);//site is for parent
    }

    inline float GetTransitionFreq(int site, StateIndex from, StateIndex to) {
        return Wrapper->Graph.GetTransitionFreq(site, from, to);
    }

    inline float GetHapProbAt(int site, int index) {
        return Wrapper->GetHapProbAt(site, index);
    }

    inline float GetEdgeProbAt(int site, StateIndex from, char allele) {
        return Wrapper->Graph.GetEdgeProbFromParentNode(site, from, allele);//nhaps of allele / total nhaps
    }

    inline char GetAllele(int site, StateIndex state) const {
        return Wrapper->GetAllele(site, state);
    }

    inline StateIndex GetChildNode(int site, StateIndex state, char allele) {
        return Wrapper->Graph.GetChildNode(site, state, allele);
    }

    inline ParentSet GetParentNodes(int site, StateIndex state) {
//        return Wrapper->Graph.StateNodeMat[site][state]->GetParentIndexSet();
        return Wrapper->Graph.GetParentSet(site, state);
    }

    inline int GetStateNumFrom(int site) {
        return Wrapper->GetNumStates(site);
    }

    inline int GetCurrentIndividualState(int site, int chrom) {
        return Wrapper->haplotypeCluster[site][2 * (individuals - 1) + chrom];
    }

    inline void SetCurrentIndividualState(int site, int chrom, int state) {
        Wrapper->haplotypeCluster[site][2 * (individuals - 1) + chrom] = state;
    }

    inline void UpdateStateNum(int num) {
        states = num;
    }

    float GetRecombRate(int from, int to);
    //HMM version two

//    struct PairHash {
//    public:
//        template<typename T, typename U>
////        long operator()(const std::pair<T, U> &x) const//Cantor pairing function:
////        {
////            int a = std::hash<T>()(x.first);
////            int b = std::hash<U>()(x.second);
////            unsigned long A = (unsigned long) (a >= 0 ? 2 * (long) a : -2 * (long) a - 1);
////            unsigned long B = (unsigned long) (b >= 0 ? 2 * (long) b : -2 * (long) b - 1);
////            long C = (long) ((A >= B ? A * A + A + B : A + B * B) / 2);
////            return a < 0 && b < 0 || a >= 0 && b >= 0 ? C : -C - 1;
////        }
//        long operator()(const std::pair<T, U> &x) const//Cantor pairing function:
//        {
//            return x.first << ((sizeof(StateIndex) -1) * 32+31) | x.second;
//        }
//    };

    /*
    class AvailableParentStatePair {
//        std::vector<NodePair > nextAvailableStatePair;//marker/availablePair
        NodePair::iterator nextAvailableStatePairIter;
//        int MarkerIndex;
    public:
        int MarkerIndex;
        std::vector<NodePair> nextAvailableStatePair;//marker/availablePair
        AvailableParentStatePair(int nmarker) : nextAvailableStatePair(nmarker, NodePair()), MarkerIndex(0) {}

        ~AvailableParentStatePair() {
            nextAvailableStatePair.clear();
        }

        void FillNextAvailableStatePair(std::pair<int, int> p) {
            nextAvailableStatePair[MarkerIndex].insert(p);
        }

        void ResetMarkerIndexAt(int index) {
            MarkerIndex = index;
            nextAvailableStatePairIter = nextAvailableStatePair[MarkerIndex].begin();
        }

        std::pair<int, int> GetNextAvailableStatePair() {
            return *(nextAvailableStatePairIter++);
        }

        bool IsEnd() {
            return nextAvailableStatePairIter == nextAvailableStatePair[MarkerIndex].end();
        }

        void NextMarker() { MarkerIndex++; }

        void PrevMarker() {
//            nextAvailableStatePair[MarkerIndex].clear();
            MarkerIndex--;
            nextAvailableStatePairIter = nextAvailableStatePair[MarkerIndex].begin();
        }

        int Size() {
            return nextAvailableStatePair[MarkerIndex].size();
        }

        int Assign(NodePair &PairHash) {
            nextAvailableStatePair[MarkerIndex] = PairHash;
        }
    };
*/

    int InitialFwdValues(int sampleIndex, FwdBwdLocalParameter &localParameter);

    int ResetFwdValues(FwdBwdLocalParameter &localParameter);

    int FindRecSite(vector<bool> &siteVec, int sampleIndex);

    int LocalForwardBackward(int sampleIndex);

    int LocalForwardBackwardSampling(int sampleIndex);

        //without recombination
    int ForwardAlgorithm(int sampleIndex, FwdBwdLocalParameter &localParameter);

    int BackwardSampling(Random *rand, int sampleIndex, char **sampledHaps, FwdBwdLocalParameter localParameter);

    //with recombination
#ifdef RECBEAGLE
    /*!ForwadAlgorithmRec()
     * forward algorithm with recombination model added
     * @return
     */
    int ForwardAlgorithmRecBeagle();
    /*!BackwardSamplingRec()
     * sampling haplotype while backward algorithm with recombination model added
     * @param rand random number generator
     * @param SampleIndex index of sample about to sample for
     * @param sampledHaps char[][] array to store sampled haplotype
     * @return
     */
    int BackwardSamplingRecBeagle(Random *rand, int SampleIndex, char **sampledHaps);
#endif

    int FindAccessibleStates(int sampleIndex, FwdBwdLocalParameter &localParameter);

    /*!ForwardAlgorithmRec
     * Calculate Forward Marginal Probability
     * @param sampleIndex
     * @param localParameter
     * @return
     */
    int ForwardAlgorithmRec(int sampleIndex, FwdBwdLocalParameter &localParameter, bool greedyMode);
    int ForwardAlgorithmRestrict(int sampleIndex, FwdBwdLocalParameter &localParameter);

    /*!BackwardAlgorithmRec
     * Calculate Backward Marginal Probability
     * @param sampleIndex
     * @param sampledHaps
     * @param localParameter
     * @return
     */
    int BackwardAlgorithmRec(Random *rand, int sampleIndex, char **sampledHaps, FwdBwdLocalParameter &localParameter);

    /*!BackwardSamplingRec()
     * sampling haplotype while backward algorithm with recombination model added
     * @param rand random number generator
     * @param SampleIndex index of sample about to sample for
     * @param sampledHaps char[][] array to store sampled haplotype
     * @return
     */
    int
    BackwardSamplingRec(Random *rand, int sampleIndex, char **sampledHaps, FwdBwdLocalParameter &localParameter);

#ifdef RECBEAGLEVARIANT
    int ForwardAlgorithmRecNew();
    /*!BackwardSamplingRec()
     * sampling haplotype while backward algorithm with recombination model added
     * @param rand random number generator
     * @param SampleIndex index of sample about to sample for
     * @param sampledHaps char[][] array to store sampled haplotype
     * @return
     */
    int BackwardSamplingRecNew(Random *rand, int SampleIndex, char** sampledHaps);
#endif
    char **sampledHaps = nullptr;
    int nSampleCopy;

    float **genoProbs;


    //KS D value related

    float ***PvalueMatrix = nullptr;//10k X 10k
    inline int CalculatePvalueMatrix() {
        std::cerr << "Enter CalculatePvalueMatrix() " << std::endl;
        PvalueMatrix = new float **[101];
        for (int i = 1; i <= 100; ++i) {
            PvalueMatrix[i] = new float *[(int) floor(10000.0 / i) + 1];
            std::cerr << "Calculating " << i << " thousand" << std::endl;
            for (int j = i; j < (int) floor(10000.0 / i) + 1; ++j) {
                PvalueMatrix[i][j] = new float[1000];

                for (int D = 1000; D > 0; D--) {
                    PvalueMatrix[i][j][D - 1] = 1. - psmirnov2x(double(D) / 1000.0, i, j);
//                    std::cerr<<i<<"\t"<<j<<"\t"<<D<<"\t"<<PvalueMatrix[i][j][size_t(D*1000)-1]<<std::endl;
                }
            }
        }
        std::cerr << "Exit CalculatePvalueMatrix() " << std::endl;
        return 0;
    }

    inline int WritePvalueMatrix(std::string fileName) {
        //std::fstream fout("/Users/fanzhang/Downloads/PlutoTest/PvalueMatrix",std::ios_base::binary|std::ios_base::out);
        GzipFileType fout(fileName.c_str(), "w");
        if (!fout.isOpen()) {
            std::cerr << "open file " << fileName << " failed!" << std::endl;
            exit(EXIT_FAILURE);
        }
        for (int i = 1; i <= 100; ++i) {
            for (int j = i; j < (int) floor(10000.0 / i) + 1; ++j) {

                for (int D = 1000; D > 0; D--) {
                    fout.write((char *) &(PvalueMatrix[i][j][D - 1]), sizeof(float));
//                    std::cerr<<"write:"<<i<<"\t"<<j<<"\t"<<D<<"\t"<<PvalueMatrix[i][j][size_t(D*1000)-1]<<std::endl;

                }
            }
        }
        fout.close();
        return 0;
    }

    inline int ReadPvalueMatrix(std::string fileName) {
        //std::fstream fin("/Users/fanzhang/Downloads/PlutoTest/PvalueMatrix",std::ios_base::binary|std::ios_base::in);
        GzipFileType fin(fileName.c_str(), "r");
        if (!fin.isOpen()) {
            std::cerr << "open file " << fileName << " failed!" << std::endl;
            std::cerr << "please specify --calPvalueMatrix to obtain PvalueMatrix file!" << std::endl;
            exit(EXIT_FAILURE);
        }
        PvalueMatrix = new float **[101];
        for (int i = 1; i <= 100; ++i) {
            PvalueMatrix[i] = new float *[(int) floor(10000.0 / i) + 1];
            for (int j = i; j < (int) floor(10000.0 / i) + 1; ++j) {
                PvalueMatrix[i][j] = new float[1000];
                for (int D = 1000; D > 0; D--) {
                    fin.read((char *) &(PvalueMatrix[i][j][D - 1]), sizeof(float));
//                    std::cerr<<"read:"<<i<<"\t"<<j<<"\t"<<D<<"\t"<<PvalueMatrix[i][j][size_t(D*1000)-1]<<std::endl;

                }
            }
        }
        fin.close();
        return 0;
    }

    inline void DestroyPvalueMatrix() {
        //KS table
        if (PvalueMatrix) {
            for (int i = 1; i <= 100; ++i) {
                for (int j = i; j < (int) floor(10000.0 / i) + 1; ++j) {
                    delete[] PvalueMatrix[i][j];
                }
                delete[] PvalueMatrix[i];
            }
            delete[] PvalueMatrix;
        }
    }


protected:
    PBWTWrapper *Wrapper = nullptr;

    bool isRev = false;
    bool onlyGT;

    void SetNextBwdValue(FwdBwdLocalParameter &localParameter, float gl, uchar allele1, uchar allele2, int destIndex,
                         float bwdValue, int site);

    float MatureBwdValue(FwdBwdLocalParameter &localParameter, int site);

    void
    ProcessFwdBwd(int sampleIndex, FwdBwdLocalParameter &localParameter, float sumBaseBwd, float *tmpGL, int site);
};

#endif //PLUTO_PBWTHAPLOTYPER_H


