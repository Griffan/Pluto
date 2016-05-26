//
// Created by Fan Zhang on 8/6/15.
//

#ifndef PLUTO_PBWTHAPLOTYPER_H
#define PLUTO_PBWTHAPLOTYPER_H


#include "ShotgunHaplotyper.h"
#include "../../pbwtWrapper/PBWTWrapper.h"
class PBWTHaplotyper : public ShotgunHaplotyper{
public:
	bool onlyHeterSite;
    PBWTHaplotyper(int nhaps, int nsnps);
	PBWTHaplotyper();
	void InitAuxillary();
    ~PBWTHaplotyper();

    void RandomSetup(Random * rand);
    int LoopThroughChromosomesViaPBWT();
	int OrderedLoopThroughChromosomesViaPBWT();
	void Transpose(int site, float * source, float * dest);
	void OrderedTranspose(int site, float *source, float *dest);
    virtual void ScoreLeftConditional();
	void OrderedScoreLeftConditional();
    virtual void ConditionOnData(float * matrix, int marker, char phred11, char phred12, char phred22);
	void OrderedConditionOnData(float *matrix, int marker, char phred11, char phred12, char phred22);
    virtual void ImputeAlleles(int marker, int state1, int state2, Random *rand, int currentIndividual, char** haps);
    virtual void ImputeAllele(int haplotype, int marker, int state, char** haps);
    virtual void FillPath(int haplotype, int fromMarker, int toMarker, int state, char** haps);
    virtual void SampleChromosomes(Random * rand);
	void OrderedSampleChromosomes(Random * rand);

    void SetUseRev(bool useOrNot){useRev=useOrNot;}
	bool ReverseInput();

    void SetOnlyGT(bool onlyOrNot){onlyGT=onlyOrNot;}
    bool GetOnlyGT(){ return onlyGT;}

	int ExtractHeterSites(int individualToProcess);
	int FillHeterSitesBack(int individualToProcess);
    //Memory management functions
    //virtual bool AllocateMemory(int nIndividuals, int maxHaplos, int nMarkers, float defaultTheta);
    //virtual void EstimateMemoryInfo(int Individuals, int Markers, int States, bool Compact, bool Phased);
    virtual void RetrieveMemoryBlock(int marker);
	void OrderedRetrieveMemoryBlock(int marker);


    virtual bool ForceMemoryAllocation();

    //inline section
	inline float GetTransitionProb(int site, int from, int to) {
//		if(Wrapper->transVector.size()<= site) {fprintf(stderr,"%d doesn't exist!\n",site);abort();}
//		if(Wrapper->transVector[site].size()<=from) {fprintf(stderr,"site:%d out of %lu sites, from:%d states too large!\n",site,Wrapper->transVector.size(),from);abort();}
//		if(Wrapper->transVector[site][from].size()<=to) {fprintf(stderr,"site:%d from:%d to:%d states too large!\n",site,from,to);abort();}
//		return Wrapper->transVector[site][from][to];
		if(Wrapper->Graph.StateNodeMat.size()<= site) {return 0;/*fprintf(stderr,"%d doesn't exist!\n",site);abort();*/}
		if(Wrapper->Graph.StateNodeMat[site].size()<=from) {return 0;/*fprintf(stderr,"site:%d out of %lu sites, from:%d states too large!\n",site,Wrapper->transVector.size(),from);abort();*/}
		if(Wrapper->Graph.StateNodeMat[site][from].childNodeIndex2NumHap.find(to)==Wrapper->Graph.StateNodeMat[site][from].childNodeIndex2NumHap.end())
		{return 0;/*fprintf(stderr,"site:%d from:%d to:%d states too large!\n",site,from,to);abort();*/}
		return Wrapper->Graph.StateNodeMat[site][from].childNodeIndex2NumHap[to];
	}
	inline uchar GetAllele(int site, int state)
	{
		return Wrapper->clusterAllele[site][state];
	}
	inline int GetStateNumFrom(int site)
	{
		return Wrapper->clusterAllele[site].size();
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


protected:

    PBWTWrapper* Wrapper;

    int indexBeingSampled;
	//subset markers related
	char** tmpHaps;
	char ** tmpGeno;
    float* tmpPenetrance;
	int tmpMarkers;
	double max_num;

	int nSampleCopy;
	char ** sampledHaps;
    double * phred2prob;

	bool useRev;
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


