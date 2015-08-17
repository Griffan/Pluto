//
// Created by Fan Zhang on 8/6/15.
//

#ifndef PLUTO_PBWTHAPLOTYPER_H
#define PLUTO_PBWTHAPLOTYPER_H


#include "ShotgunHaplotyper.h"
#include "../../pbwtWrapper/PBWTWrapper.h"
class PBWTHaplotyper : public ShotgunHaplotyper{
public:
    PBWTHaplotyper(int nhaps, int nsnps);
    ~PBWTHaplotyper();
    int LoopThroughChromosomesViaPBWT();
	void Transpose(int site, float * source, float * dest, float theta);
    virtual void ScoreLeftConditional();
    virtual void ConditionOnData(float * matrix, int marker, char phred11, char phred12, char phred22);
    virtual void ImputeAlleles(int marker, int state1, int state2, Random * rand);
    virtual void SampleChromosomes(Random * rand);

	inline float getTransitionProb(int site, int from, int to);

    //Memory management functions
    //virtual bool AllocateMemory(int nIndividuals, int maxHaplos, int nMarkers, float defaultTheta);
    //virtual void EstimateMemoryInfo(int Individuals, int Markers, int States, bool Compact, bool Phased);
    virtual void RetrieveMemoryBlock(int marker);
    virtual bool ForceMemoryAllocation();

    //inline section
	inline uchar getAllele(int site, int state)
	{
		return Wrapper.clusterAllele[site][state];
	}
	inline int getStateNumFrom(int site)
	{
		return Wrapper.clusterAllele[site].size();
	}
	inline int getCurrentIndividualState(int site, int chrom)
	{
		return Wrapper.haplotypeCluster[site][2 * (individuals - 1) + chrom];
	}
	inline int setCurrentIndividualState(int site, int chrom, int state)
	{
		return Wrapper.haplotypeCluster[site][2 * (individuals - 1) + chrom]=state;
	}
    inline int UpdateStateNum(int num){
        states=num;
    }
	//TODO:
	//1. assign states to current individual
	//2. change transpose function
protected:
    PBWTWrapper Wrapper;

};


#endif //PLUTO_PBWTHAPLOTYPER_H


