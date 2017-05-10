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

class PBWTHaplotyper : public ShotgunHaplotyper{
public:
	bool onlyHeterSite;
    bool geneticMapAvailable;

    GeneticDistanceMap GDMap;

    PBWTHaplotyper(int nhaps, int nsnps);
	PBWTHaplotyper();
	void InitAuxillary();
    ~PBWTHaplotyper();

    /*!InitialHaplotypeByConsensus
     * copy consencus haplotypes into haplotypes char[][], and reset consensus
     * @param consensus
     */
    void InitialHaplotypeByConsensus(ConsensusBuilder& consensus)
    {
        for (int i = 0; i < (individuals-phased)*2; ++i) {
            for (int j = 0; j < markers ; ++j) {
                haplotypes[i][j] = consensus.consensus[i][j];
            }
        }
        consensus.stored=0;
    }

	void InitialSampleCopy(Random * rand);
//    void RandomSetup(Random * rand);
    void SwapIndividuals(int a, int b);
    void PrepareRefSetPBWTWrapper();

    int LoopThroughChromosomesHighPrecision();
    int LoopThroughChromosomesRecomb();
    int LoopThroughChromosomesSingleRound();


	int LoopThroughChromosomesViaPBWTWithHeterOnly();

	void Transpose(int site, float * source, float * dest);

	virtual void RandomSetup(Random * rand = NULL);
    virtual void ScoreLeftConditional();

    virtual void ConditionOnData(float * matrix, int marker, char phred11, char phred12, char phred22);

    virtual void ImputeAlleles(int marker, int state1, int state2, Random *rand, int currentIndividual, char** haps);
	virtual void ImputeAllelesRaw(int marker, int state1, int state2, Random *rand, int currentIndividual, char** haps);
    virtual void ImputeAllele(int haplotype, int marker, int state, char** haps);
    virtual void FillPath(int haplotype, int fromMarker, int toMarker, int state, char** haps);
    virtual void SampleChromosomes(Random * rand);


    void SetUseRev(bool useOrNot){isRev=useOrNot;}
	bool ReverseInput();

    void SetOnlyGT(bool onlyOrNot){onlyGT=onlyOrNot;}
    bool GetOnlyGT(){ return onlyGT;}

	int ExtractHeterSites(int individualToProcess);
	int FillHeterSitesBack(int individualToProcess);
    //Memory management functions
    //virtual bool AllocateMemory(int nIndividuals, int maxHaplos, int nMarkers, float defaultTheta);
    //virtual void EstimateMemoryInfo(int Individuals, int Markers, int States, bool Compact, bool Phased);
    virtual void RetrieveMemoryBlock(int marker);



    virtual bool ForceMemoryAllocation();

    //inline section
	inline float GetTransitionProb(int site, StateIndex from, StateIndex to) {
//		if(Wrapper->transVector.size()<= site) {fprintf(stderr,"%d doesn't exist!\n",site);abort();}
//		if(Wrapper->transVector[site].size()<=from) {fprintf(stderr,"site:%d out of %lu sites, from:%d states too large!\n",site,Wrapper->transVector.size(),from);abort();}
//		if(Wrapper->transVector[site][from].size()<=to) {fprintf(stderr,"site:%d from:%d to:%d states too large!\n",site,from,to);abort();}
//		return Wrapper->transVector[site][from][to];
        uchar allele=GetAllele(site+1,to);
		if((int)Wrapper->Graph.StateNodeMat.size()<= site) {fprintf(stderr,"site %d doesn't exist!\n",site);abort();}
		if((int)Wrapper->Graph.StateNodeMat[site].size()<=from)
        {
//            fprintf(stderr,"site:%d, from:%d states too large!\n",site,from);
            return 0.f;
        }
        if(Wrapper->Graph.StateNodeMat[site][from]->childNodeIndex[allele]== -1||Wrapper->Graph.StateNodeMat[site][from]->childNodeIndex[allele]!=to)
		{
//            fprintf(stderr,"site:%d from:%d to:%d states too large!\n",site,from,to);
            return 0.f;
        }
//		return Wrapper->Graph.StateNodeMat[site][from]->numHapChild[GetAllele(site+1,to)];
        return Wrapper->Graph.GetProbToCurrentNodeConditionalOnParentNode(site, from, allele);//site is for parent
	}
    inline float GetHapProbAt(int site,int index)
    {
        return Wrapper->GetHapProbAt(site,index);
    }
    inline float GetEdgeProbAt(int site, StateIndex from, uchar allele)
    {
        return Wrapper->Graph.GetEdgeProbFromParentNode(site,from,allele);//
    }
	inline uchar GetAllele(int site, StateIndex state)
	{
		return Wrapper->GetAllele(site,state);
	}
    inline StateIndex GetChildNode(int site, StateIndex state, uchar allele)
    {
        return Wrapper->Graph.StateNodeMat[site][state]->childNodeIndex[allele];
    }
	inline int GetStateNumFrom(int site)
	{
		return Wrapper->GetNumStates(site);
	}
	inline int GetCurrentIndividualState(int site, int chrom)
	{
		return Wrapper->haplotypeCluster[site][2 * (individuals - 1) + chrom];
	}
	inline void SetCurrentIndividualState(int site, int chrom, int state)
	{
		Wrapper->haplotypeCluster[site][2 * (individuals - 1) + chrom]=state;
	}
    inline void UpdateStateNum(int num){
         states=num;
    }

    //HMM version two

    struct pairhash {
    public:
        template <typename T, typename U>
        long operator()(const std::pair<T, U> &x) const//Cantor pairing function:
        {
            int a=std::hash<T>()(x.first);
            int b=std::hash<U>()(x.second);
            unsigned long A = (unsigned long)(a >= 0 ? 2 * (long)a : -2 * (long)a - 1);
            unsigned long B = (unsigned long)(b >= 0 ? 2 * (long)b : -2 * (long)b - 1);
            long C = (long)((A >= B ? A * A + A + B : A + B * B) / 2);
            return a < 0 && b < 0 || a >= 0 && b >= 0 ? C : -C - 1;
        }
    };
    typedef std::unordered_map<std::pair<int32_t,int32_t>,float,pairhash> Source;//(nodeA,nodeB)->fwd
    typedef std::unordered_set<std::pair<int32_t,int32_t>,pairhash> NodePair;//(nodeA,nodeB)
    typedef std::unordered_map<StateIndex,std::unordered_map<StateIndex,Source> > ChildToSource;
    std::vector<ChildToSource> genuienParents;

    float* fwdValueSum;
    std::vector<std::unordered_map<int,float> > fwdValueNode1Sum;
    std::vector<std::unordered_map<int,float> > fwdValueNode2Sum;
    float SumFwdValueFromOriginVec(const Source& a)
    {
        float sum(0.f);
        for (auto kv:a) {
            sum+=kv.second;
        }
        return sum;
    }

    int InitialFwdValues();
    int ResetFwdValues();
    double GetGL(int individual, int marker, uchar allele1, uchar allele2);
    float GetRecombRate(int marker);
    StateIndex GetChildNode(int site, int stateIndex, uchar allele)
    {
        if(Wrapper->Graph.StateNodeMat[site][stateIndex]->childNodeIndex[allele]== -1)
            return -1;
        else return Wrapper->Graph.StateNodeMat[site][stateIndex]->childNodeIndex[allele];
    }



    //without recombination
    std::vector<bool> brokenList;
    int ForwardAlgorithm();
    int BackwardSampling(Random *rand, int SampleIndex, char** sampledHaps);
    //with recombination
    class AvailableParentStatePair {
//        std::vector<NodePair > nextAvailableStatePair;//marker/availablePair
        NodePair::iterator nextAvailableStatePairIter;
//        int MarkerIndex;
    public:
        int MarkerIndex;
        std::vector<NodePair > nextAvailableStatePair;//marker/availablePair
        AvailableParentStatePair(int nmarker):nextAvailableStatePair(nmarker,NodePair()),MarkerIndex(0)
        {}
        ~AvailableParentStatePair()
        {
            nextAvailableStatePair.clear();
        }
        void FillNextAvailableStatePair(std::pair<int, int> p) {
            nextAvailableStatePair[MarkerIndex].insert(p);
        }

        void ResetMarkerIndexAt(int index) {
            MarkerIndex = index;
            nextAvailableStatePairIter=nextAvailableStatePair[MarkerIndex].begin();
        }

        std::pair<int, int> GetNextAvailableStatePair() {
            return *(nextAvailableStatePairIter++);
        }

        bool IsEnd() {
            return nextAvailableStatePairIter == nextAvailableStatePair[MarkerIndex].end();
        }
        void NextMarker(){MarkerIndex++;}
        void LastMarker(){
//            nextAvailableStatePair[MarkerIndex].clear();
            MarkerIndex--;
            nextAvailableStatePairIter=nextAvailableStatePair[MarkerIndex].begin();}
    };
    AvailableParentStatePair availablePair;//available parents at each site
    /*!ForwadAlgorithmRec()
     * forward algorithm with recombination model added
     * @return
     */
    int ForwardAlgorithmRec();
    /*!BackwardSamplingRec()
     * sampling haplotype while backward algorithm with recombination model added
     * @param rand random number generator
     * @param SampleIndex index of sample about to sample for
     * @param sampledHaps char[][] array to store sampled haplotype
     * @return
     */
    int BackwardSamplingRec(Random *rand, int SampleIndex, char** sampledHaps);


	char ** sampledHaps;
	int nSampleCopy;


    //KS D value related

    float *** PvalueMatrix;//10k X 10k
    inline int CalculatePvalueMatrix()
    {
        std::cerr<<"Enter CalculatePvalueMatrix() "<<std::endl;
        PvalueMatrix=new float ** [101];
        for (int i = 1; i <=100; ++i) {
            PvalueMatrix[i]=new float* [(int)floor(10000.0/i)+1];
            std::cerr<<"Calculating "<<i<<" thousand"<<std::endl;
            for (int j =i; j <(int)floor(10000.0/i)+1; ++j) {
                PvalueMatrix[i][j]=new float [1000];

                for(int D=1000;D>0;D--)
                {
                    PvalueMatrix[i][j][D-1]=1.-psmirnov2x(double(D)/1000.0, i, j);
//                    std::cerr<<i<<"\t"<<j<<"\t"<<D<<"\t"<<PvalueMatrix[i][j][size_t(D*1000)-1]<<std::endl;
                }
            }
        }
        std::cerr<<"Exit CalculatePvalueMatrix() "<<std::endl;
        return 0;
    }

    inline int WritePvalueMatrix(std::string fileName)
    {
        //std::fstream fout("/Users/fanzhang/Downloads/PlutoTest/PvalueMatrix",std::ios_base::binary|std::ios_base::out);
        GzipFileType fout(fileName.c_str(),"w");
        if(!fout.isOpen())
        {
            std::cerr<<"open file "<<fileName<<" failed!"<<std::endl;
            exit(EXIT_FAILURE);
        }
        for (int i = 1; i <=100 ; ++i) {
            for (int j =i; j <(int)floor(10000.0/i)+1; ++j) {

                for(int D=1000;D>0;D--)
                {
                    fout.write((char*)&(PvalueMatrix[i][j][D-1]),sizeof(float));
//                    std::cerr<<"write:"<<i<<"\t"<<j<<"\t"<<D<<"\t"<<PvalueMatrix[i][j][size_t(D*1000)-1]<<std::endl;

                }
            }
        }
        fout.close();
        return 0;
    }

    inline int ReadPvalueMatrix(std::string fileName)
    {
        //std::fstream fin("/Users/fanzhang/Downloads/PlutoTest/PvalueMatrix",std::ios_base::binary|std::ios_base::in);
        GzipFileType fin(fileName.c_str(),"r");
        if(!fin.isOpen())
        {
            std::cerr<<"open file "<<fileName<<" failed!"<<std::endl;
            std::cerr<<"please specify --calPvalueMatrix to obtain PvalueMatrix file!"<<std::endl;
            exit(EXIT_FAILURE);
        }
        PvalueMatrix=new float ** [101];
        for (int i = 1; i <=100 ; ++i) {
            PvalueMatrix[i]=new float* [(int)floor(10000.0/i)+1];
            for (int j =i; j <(int)floor(10000.0/i)+1; ++j) {
                PvalueMatrix[i][j]=new float [1000];
                for(int D=1000;D>0;D--)
                {
                    fin.read((char*)&(PvalueMatrix[i][j][D-1]),sizeof(float));
//                    std::cerr<<"read:"<<i<<"\t"<<j<<"\t"<<D<<"\t"<<PvalueMatrix[i][j][size_t(D*1000)-1]<<std::endl;

                }
            }
        }
        fin.close();
        return 0;
    }

    inline void DestroyPvalueMatrix()
    {
        //KS table
        for (int i = 1; i <=100; ++i) {
            for (int j =i; j <(int)floor(10000.0/i)+1; ++j) {
                delete [] PvalueMatrix[i][j];
            }
            delete [] PvalueMatrix[i];
        }
        delete [] PvalueMatrix;
    }
protected:

    PBWTWrapper* Wrapper;

    int indexBeingSampled;
	//subset markers related
	char** tmpHaps;
	char ** tmpGeno;
    float* tmpPenetrance;
	int tmpMarkers;
	double max_num;



//    double * phred2prob;

	bool isRev;
    bool onlyGT;

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

	//memory management related
	std::unordered_map<int, std::vector<float *> > memoryBlockList;//size and list of address
	std::unordered_map<int, int> numInUse;//size and list of address
	int totalBlockNum;

	void GetMemoryBlock(int marker);//revise GetMemoryBlock function
	float* GetReuseableBlock();//revise GetReuseableBlock function
	float* GetLargeBlock();
	void ResetMemoryPool();
	void ResetReuseablePool();
	void ReleaseMemoryBlock();
};


#endif //PLUTO_PBWTHAPLOTYPER_H


