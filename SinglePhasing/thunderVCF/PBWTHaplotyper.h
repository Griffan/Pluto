//
// Created by Fan Zhang on 8/6/15.
//

#ifndef PLUTO_PBWTHAPLOTYPER_H
#define PLUTO_PBWTHAPLOTYPER_H




#include <unordered_set>
#include "ShotgunHaplotyper.h"
#include "../../pbwtWrapper/PBWTWrapper.h"
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

	void InitialSampleCopy(Random * rand);
//    void RandomSetup(Random * rand);
    void SwapIndividuals(int a, int b);
    void PrepareRefSetPBWTWrapper();
    void PrepareRefSetPBWTWrapperLeaveOneOut();
    int LoopThroughChromosomesHighPrecision();
    int LoopThroughChromosomesSingleRound();
    int LoopThroughChromosomesLeaveOneOut();

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



    virtual bool ForceMemoryAllocation();

    //inline section
	inline float GetTransitionProb(int site, int from, int to) {
//		if(Wrapper->transVector.size()<= site) {fprintf(stderr,"%d doesn't exist!\n",site);abort();}
//		if(Wrapper->transVector[site].size()<=from) {fprintf(stderr,"site:%d out of %lu sites, from:%d states too large!\n",site,Wrapper->transVector.size(),from);abort();}
//		if(Wrapper->transVector[site][from].size()<=to) {fprintf(stderr,"site:%d from:%d to:%d states too large!\n",site,from,to);abort();}
//		return Wrapper->transVector[site][from][to];
		if(Wrapper->Graph.StateNodeMat.size()<= site) {fprintf(stderr,"site %d doesn't exist!\n",site);abort();}
		if(Wrapper->Graph.StateNodeMat[site].size()<=from) {fprintf(stderr,"site:%d, from:%d states too large!\n",site,from);abort();}
        if(Wrapper->Graph.StateNodeMat[site][from]->childNodeIndex[GetAllele(site+1,to)]== nullptr||*(Wrapper->Graph.StateNodeMat[site][from]->childNodeIndex[GetAllele(site+1,to)])!=to)
		{fprintf(stderr,"site:%d from:%d to:%d states too large!\n",site,from,to);abort();}
		return Wrapper->Graph.StateNodeMat[site][from]->numHapChild[GetAllele(site+1,to)];
	}
	inline uchar GetAllele(int site, int state)
	{
		return Wrapper->GetAllele(site,state);
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
    struct origin
    {
        int firstState;
        int secondState;
		float fwdValue;
        origin(int a,int b, float c)
        {
            firstState=a;
            secondState=b;
            fwdValue=c;
        }
    };
//    typedef std::vector<origin> originVec;


    float* fwdValueSum;
    int InitialFwdValues();
    int ResetFwdValues();
    float GetGL(int individual, int marker, char allele1, char allele2);
    int GetChildNode(int site, int stateIndex, char allele)
    {
        if(Wrapper->Graph.StateNodeMat[site][stateIndex]->childNodeIndex[allele]== nullptr)
            return -1;
        else return *(Wrapper->Graph.StateNodeMat[site][stateIndex]->childNodeIndex[allele]);
    }
//    int SwitchFwdValuePtr()
//    {
//        std::unordered_map<int,std::unordered_map<int,originVec> >* tmp;
//        currentFwdValuePtr->clear();
//        tmp=currentFwdValuePtr;
//        currentFwdValuePtr=nextFwdValuePtr;
//        nextFwdValuePtr=tmp;
//    }
    int ForwardAlgorithm();
    std::vector<bool> brokenList;
    int BackwardSampling(Random *rand, int SampleIndex, char** sampledHaps);
//    std::unordered_map<int,std::unordered_map<int,originVec> > currentFwdValue,nextFwdValue;
//    std::unordered_map<int,std::unordered_map<int,originVec> >* currentFwdValuePtr,*nextFwdValuePtr;

    struct pairhash {
    public:
        template <typename T, typename U>
        long operator()(const std::pair<T, U> &x) const//Cantor pairing function:
        {
            unsigned long a=std::hash<T>()(x.first);
			unsigned long b=std::hash<U>()(x.second);
			unsigned long A = (unsigned long)(a >= 0 ? 2 * a : -2 * a - 1);
			unsigned long B = (unsigned long)(b >= 0 ? 2 * b : -2 * b - 1);
            long C = (long)((A >= B ? A * A + A + B : A + B * B) / 2);
            return a < 0 && b < 0 || a >= 0 && b >= 0 ? C : -C - 1;
        }
    };

    typedef std::unordered_map<std::pair<int32_t,int32_t>,float,pairhash> Source;
	typedef std::unordered_map<int,std::unordered_map<int,Source> > ChildToSource;
    std::vector<ChildToSource> genuienParents;

    float SumFwdValueFromOriginVec(const Source& a)
    {
        float sum(0.f);
        for (auto kv:a) {
            sum+=kv.second;
        }
        return sum;
    }
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
        GzipFileType fout(fileName.c_str(),"r");
        if(!fout.isOpen())
        {
            std::cerr<<"open file /Users/fanzhang/Downloads/PlutoTest/PvalueMatrix failed!"<<std::endl;
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

    inline int DestroyPvalueMatrix()
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

    PBWTWrapper* Wrapper,*fwdWrapper,*backWrapper,*baseWrapper;

    int indexBeingSampled;
	//subset markers related
	char** tmpHaps;
	char ** tmpGeno;
    float* tmpPenetrance;
	int tmpMarkers;
	double max_num;



//    double * phred2prob;

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


