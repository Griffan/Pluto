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

#define INDEX 1
#define PHASE 2
#define ITERATIVE 4
#define HASH_RESERVE 4096
class PBWTHaplotyper : public ShotgunHaplotyper {
public:
    std::string loadGraph = "Empty";//indicate if build graph from previously built graph
    std::string outputPrefix = "Empty";
    bool geneticMapAvailable = false;
    int runningModel = 0;
    int prefixLength = 0;
//    int phasedForByRef = 0;
    int nThread = 0;

    GeneticDistanceMap GDMap;

    PBWTHaplotyper(int nhaps, int nsnps);

    PBWTHaplotyper();

//	void InitAuxillary();
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

//    void RandomSetup(Random * rand);
    void SwapIndividuals(int a, int b);

    void PrepareRefSetPBWTWrapper();

    int LoopThroughChromosomesHighPrecision();

    int LoopThroughChromosomesRecomb(Pedigree &ped);

    int LoopThroughChromosomesSingleRound();


    int LoopThroughChromosomesViaPBWTWithHeterOnly();

//	void Transpose(int site, float * source, float * dest);

    virtual void RandomSetup(Random *rand = NULL);
//    virtual void ScoreLeftConditional();

//    virtual void ConditionOnData(float * matrix, int marker, char phred11, char phred12, char phred22);

    virtual void ImputeAlleles(int marker, int state1, int state2, Random *rand, int currentIndividual, char **haps);

    virtual void ImputeAllelesRaw(int marker, int state1, int state2, Random *rand, int currentIndividual, char **haps);

    virtual void ImputeAllele(int haplotype, int marker, int state, char **haps);

    virtual void FillPath(int haplotype, int fromMarker, int toMarker, int state, char **haps);


    void SetUseRev(bool useOrNot) { isRev = useOrNot; }

    bool ReverseInput();

    void SetOnlyGT(bool onlyOrNot) { onlyGT = onlyOrNot; }

    bool GetOnlyGT() { return onlyGT; }

#ifdef HETERSITE
    bool onlyHeterSite=false;
    //subset markers related
    char** tmpHaps= nullptr;
    char ** tmpGeno= nullptr;
    float* tmpPenetrance = nullptr;
    int tmpMarkers;
    double max_num;
    std::vector<int> relativeIndexToAbsolute;
    std::unordered_map<int,int> absoluteIndexToRelative;
    inline int SwapTempHaps()
    {
        char** tmpH = haplotypes;
        haplotypes=tmpHaps;
        tmpHaps=tmpH;

        tmpH=genotypes;
        genotypes=tmpGeno;
        tmpGeno=tmpH;

        float * tmpP=penetrances;
        penetrances=tmpPenetrance;
        tmpPenetrance=tmpP;

        int tmpM=markers;
        markers=tmpMarkers;
        tmpMarkers=tmpM;
        return 0;
    }

    int ExtractHeterSites(int individualToProcess);
    int FillHeterSitesBack(int individualToProcess);
#endif

    //Memory management functions
    bool SetErrorAndTheta(std::vector<float> &holderError, std::vector<float> &holderTheta);

    bool AllocateMemory(int nIndividuals, int nMarkers);
//    virtual void EstimateMemoryInfo(int Individuals, int Markers, int States, bool Compact, bool Phased);
//    virtual void RetrieveMemoryBlock(int marker);



    virtual bool ForceMemoryAllocation();

    //inline section
    inline float GetTransitionProb(int site, StateIndex from, StateIndex to) {
        uchar allele = GetAllele(site + 1, to);
        if ((int) Wrapper->Graph.StateNodeMat.size() <= site) {
            fprintf(stderr, "site %d doesn't exist!\n", site);
            abort();
        }
        if ((int) Wrapper->Graph.StateNodeMat[site].size() <= from) {
            fprintf(stderr, "site:%d, from:%d states too large!\n", site, from);
            return 0.f;
        }
        if (Wrapper->Graph.StateNodeMat[site][from]->GetChildNodeIndex(allele) == -1 ||
            Wrapper->Graph.StateNodeMat[site][from]->GetChildNodeIndex(allele) != to) {
            fprintf(stderr, "site:%d from:%d to:%d states too large!\n", site, from, to);
            return 0.f;
        }
        return Wrapper->Graph.GetProbToCurrentNodeConditionalOnParentNode(site, from, allele);//site is for parent
    }

    inline float GetHapProbAt(int site, int index) {
        return Wrapper->GetHapProbAt(site, index);
    }

    inline float GetEdgeProbAt(int site, StateIndex from, uchar allele) {
        return Wrapper->Graph.GetEdgeProbFromParentNode(site, from, allele);//
    }

    inline uchar GetAllele(int site, StateIndex state)const {
        return Wrapper->GetAllele(site, state);
    }

    inline StateIndex GetChildNode(int site, StateIndex state, uchar allele) {
        return Wrapper->Graph.StateNodeMat[site][state]->GetChildNodeIndex(allele);
    }

    inline ParentSet &GetParentNodes(int site, StateIndex state) {
        return Wrapper->Graph.StateNodeMat[site][state]->GetParentIndexSet();
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

    double GetGL(int individual, int marker, uchar allele1, uchar allele2);

    float GetRecombRate(int marker);
    //HMM version two

    struct PairHash {
    public:
        template<typename T, typename U>
//        long operator()(const std::pair<T, U> &x) const//Cantor pairing function:
//        {
//            int a = std::hash<T>()(x.first);
//            int b = std::hash<U>()(x.second);
//            unsigned long A = (unsigned long) (a >= 0 ? 2 * (long) a : -2 * (long) a - 1);
//            unsigned long B = (unsigned long) (b >= 0 ? 2 * (long) b : -2 * (long) b - 1);
//            long C = (long) ((A >= B ? A * A + A + B : A + B * B) / 2);
//            return a < 0 && b < 0 || a >= 0 && b >= 0 ? C : -C - 1;
//        }
        long operator()(const std::pair<T, U> &x) const//Cantor pairing function:
        {
            return x.first<<(sizeof(StateIndex)*32)|x.second;
        }
    };

//    typedef std::unordered_map<std::pair<StateIndex, StateIndex>, float, PairHash> Source;//(nodeA,nodeB)->fwd
//    typedef std::unordered_map<std::pair<StateIndex, StateIndex>, Source, PairHash> DestToSource;

//    typedef std::unordered_map<u_int64_t, float> Source;//(nodeA,nodeB)->fwd
//    typedef std::unordered_map<u_int64_t, Source> DestToSource;

    typedef u_int64_t NodePair;
    typedef std::vector<u_int64_t> SourceVec;
    typedef std::vector<std::vector<u_int64_t> > Index2SourceVec;
    typedef std::vector<float> FwdVec;
    typedef std::vector<std::vector<float> > Index2FwdVec;
    typedef std::unordered_map<u_int64_t, int> Dest2SourceVecIndex;

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
    struct FwdBwdLocalParameter {
        int states;
        //from (parentNode1, parentNode2) to (childNode1, childNode2), childNodes are present in conditional graph, but parentNode are not necessarily present
//        std::vector<DestToSource> parentsNodeVec;

        std::vector<Dest2SourceVecIndex> parentsNodeVec;
        std::vector<Index2SourceVec > megaSourceVec;
        std::vector<Index2FwdVec> megaFwdVec;

        std::vector<float> fwdValueSum;
        std::vector<std::unordered_map<int, float> > fwdValueNode1Sum;
        std::vector<std::unordered_map<int, float> > fwdValueNode2Sum;
        std::vector<bool> isRec;

        FwdBwdLocalParameter(int individuals, int markers) {
//            DestToSource dummy;
//            dummy.reserve(HASH_RESERVE);
//            parentsNodeVec.assign(markers, dummy);
            Dest2SourceVecIndex dummy;
            dummy.reserve(HASH_RESERVE);
            parentsNodeVec.assign(markers,dummy);

            Index2SourceVec dummy2;
            dummy2.reserve(HASH_RESERVE);
            megaSourceVec.assign(markers, dummy2);

            Index2FwdVec dummy3;
            dummy3.reserve(HASH_RESERVE);
            megaFwdVec.assign(markers, dummy3);

        }
        inline u_int64_t MakePair(StateIndex first, StateIndex second)
        {
            return (u_int64_t)first<<32|second;
        }

        inline StateIndex GetFirst(u_int64_t pair)
        {
            return (StateIndex)(pair>>32);
        }

        inline StateIndex GetSecond(u_int64_t pair)
        {
            return (StateIndex)(pair&0xffff);
        }
        inline int FillParentsNodeVec(int i, StateIndex childNode1, StateIndex childNode2, StateIndex parentNode1, StateIndex parentNode2, float tmpFwdValue)
        {
//            parentsNodeVec[i][MakePair(childNode1, childNode2)][MakePair(parentNode1, parentNode2)] = tmpFwdValue;
            u_int64_t index = std::numeric_limits<u_int64_t >::max();
            if(parentsNodeVec[i].find(MakePair(childNode1, childNode2)) != parentsNodeVec[i].end())//Dest already exists
            {
                index = parentsNodeVec[i][MakePair(childNode1, childNode2)];
                megaFwdVec[i][index].push_back(tmpFwdValue);
                megaSourceVec[i][index].push_back(MakePair(parentNode1, parentNode2));
            }
            else
            {
                index = megaFwdVec[i].size();
                parentsNodeVec[i][MakePair(childNode1, childNode2)] = index;

                FwdVec dummyFwdVec(1, tmpFwdValue);
                dummyFwdVec.reserve(HASH_RESERVE/2.0);
                megaFwdVec[i].push_back(dummyFwdVec);

                SourceVec dummySourceVec(1, MakePair(parentNode1,parentNode2));
                dummySourceVec.reserve(HASH_RESERVE/2.0);
                megaSourceVec[i].push_back(dummySourceVec);
            }
            return 0;
        }

        inline FwdVec & GetFwdVec(int i, int destIndex)
        {
            return megaFwdVec[i][destIndex];
        }

        inline FwdVec & GetFwdVec(int i, StateIndex A, StateIndex B)
        {
            return GetFwdVec(i,parentsNodeVec[i][MakePair(A,B)]);
        }

        inline float GetFwd(int i, int destIndex, int sourceIndex)
        {
            return megaFwdVec[i][destIndex][sourceIndex];
        }

        inline SourceVec & GetSourceVec(int i,  int destIndex)
        {
            return megaSourceVec[i][destIndex];
        }

        inline SourceVec & GetSourceVec(int i, StateIndex A, StateIndex B)
        {
            return GetSourceVec(i,parentsNodeVec[i][MakePair(A,B)]);
        }

        inline u_int64_t GetSource(int i, int destIndex, int sourceIndex)
        {
            return megaSourceVec[i][destIndex][sourceIndex];
        }

        inline int GetDestIndex(int i, StateIndex A, StateIndex B)
        {
            return parentsNodeVec[i][MakePair(A,B)];
        }

        inline float SumFwdValueFromOriginVec(const FwdVec &a) const {
            if (a.size() == 0) return 0.f;
            float sum(0.f);
            for (auto kv:a) {
                sum += kv;
            }
            return sum;
        }
    };

    int InitialFwdValues(int sampleIndex, FwdBwdLocalParameter &localParameter);

    int ResetFwdValues(FwdBwdLocalParameter &localParameter);

    int LocalForwadBackWard(int sampleIndex);

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

    int ForwardAlgorithmRec(int sampleIndex, FwdBwdLocalParameter &localParameter);

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

    //memory management related
    std::unordered_map<int, std::vector<float *> > memoryBlockList;//size and list of address
    std::unordered_map<int, int> numInUse;//size and list of address
    int totalBlockNum;

    void GetMemoryBlock(int marker);//revise GetMemoryBlock function
    float *GetReuseableBlock();//revise GetReuseableBlock function
    float *GetLargeBlock();

    void ResetMemoryPool();

    void ResetReuseablePool();

    void ReleaseMemoryBlock();

};

#endif //PLUTO_PBWTHAPLOTYPER_H


