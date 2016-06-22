//
// Created by Fan Zhang on 7/20/15.
//
/*
 * Notice that currently, we only want to forward
 * and backward following thunder's iterating
 * pattern. Hence both forward and backward initiation
 * are not implemented yet.
 */
#ifndef PLUTO_PBWTWRAPPER_H
#define PLUTO_PBWTWRAPPER_H
//#define DEBUG 1


#include "pbwt/pbwt.h"
#include "ks.h"
#include <vector>
#include <unordered_map>
#include <iostream>
#include <algorithm>
#include <queue>
#include <map>
#include <numeric>
#include <fstream>
#include "Rmath.h"


struct max_pair_t
{
    int clusterA;
    size_t sizeA;
    int clusterB;
    size_t sizeB;
    double Dmax;
    bool exact;
    double pval;
    max_pair_t(int a,int b,double c, bool d,double e, size_t A, size_t B)
    {
        clusterA=a;
        clusterB=b;
        Dmax=c;
        exact=d;
        pval=e;
        sizeA=A;
        sizeB=B;
    }
};
bool comparator(const max_pair_t& lhs, const max_pair_t& rhs);


template<typename T>
struct square
{
    T operator()(const T& Left, const T& Right) const
    {
        return (Left + Right*Right);
    }
};

extern const float T_CRITICAL_VALUE[];

inline float GetTCritical(uint32_t df)
{
    if(df<=30) return T_CRITICAL_VALUE[df-1];
    else if(df<=100) return T_CRITICAL_VALUE[29+(df-30)/5];
    else if(df<=200) return T_CRITICAL_VALUE[44];
    else if(df<=500) return T_CRITICAL_VALUE[45];
    else return T_CRITICAL_VALUE[46];
}

struct r_stat_t
{
    float S0x,S1x;
    float Sy;
    float S0xx,S1xx;
    float Syy;
    float S0xy,S1xy;
    int size;
    r_stat_t(const std::vector<float>& yCoordinate)
    {
        size=yCoordinate.size();
        S0x=0;
        S1x=yCoordinate.size();
        Sy=std::accumulate(yCoordinate.begin(), yCoordinate.end(), 0.0);
        S0xx=0;
        S1xx=S1x;
        Syy=std::accumulate(yCoordinate.begin(), yCoordinate.end(), 0.0, square<float>());
        S0xy=0;
        S1xy=Sy;
    }
    r_stat_t()
    {
        size=0;
        S0x=0;
        S1x=0;
        Sy=0;
        S0xx=0;
        S1xx=0;
        Syy=0;
        S0xy=0;
        S1xy=0;
    }
    r_stat_t& operator=(const r_stat_t&a)
    {
        size=a.size;
        S1x=a.S1x;
        Sy=a.Sy;
        S1xx=a.S1xx;
        Syy=a.Syy;
        S0xy=a.S0xy;
        S1xy=a.S1xy;
        return *this;
    }
    r_stat_t Combine(const r_stat_t& a)
    {
        r_stat_t b;
        b.size=size+a.size;
        b.S0x=S0x+a.S1x;
        b.S1x=S1x+a.S0x;
        b.Sy=Sy+a.Sy;
        b.S0xx=S0xx+a.S1xx;
        b.S1xx=S1xx+a.S0xx;
        b.Syy=Syy+a.Syy;
        b.S0xy=S0xy+a.S1xy;
        b.S1xy=S1xy+a.S0xy;
        return b;
    }
    r_stat_t& operator+(const r_stat_t&a)
    {
        size+=a.size;
        S0x+=a.S0x;
        S1x+=a.S1x;
        Sy+=a.Sy;
        S0xx+=a.S0xx;
        S1xx+=a.S1xx;
        Syy+=a.Syy;
        S0xy+=a.S0xy;
        S1xy+=a.S1xy;
        return *this;
    }
    bool IsSignificant()
    {
        double r=(size*S0xy-S0x*Sy)/(sqrt(size*S0xx-S0x*S0x)*sqrt(size*Syy-Sy*Sy));
//        double r2=(size*S1xy-S1x*Sy)/(sqrt(size*S1xx-S1x*S1x)*sqrt(size*Syy-Sy*Sy));
        double t=r*sqrt(size-2)/sqrt(1-r*r);
        return fabs(t)>GetTCritical(size-2);
//        double t2=r2*sqrt(size-2)/sqrt(1-r2*r2);
//        double p=pt(t,size-2,1,0);

//if(p<0.25||p>0.75)        fprintf(stderr,"r:%f,\tbeta hat:%f,\tp value:%f,\t df:%d\n",r,S0xy/S0xx,p,size-2);
//        return p<0.25||p>0.75;

    }
};


//define NODE structure which represents cluster or state in HMM model
class StateNode
{
private:
    StateNode(const StateNode&);
    StateNode & operator=(const StateNode&);
public:
    int nodeIndex;
    float numHap;
//    std::vector<int> parentNodeIndex;
//    std::vector<int> numHapFromParentNode;
//    std::vector<int> childNodeIndex;
//    std::vector<float> numHapToChildNode;
    std::unordered_map<int,StateNode* > parentNodeIndex2NumHap;
//    std::unordered_map<int,float > childNodeIndex2NumHap;
    int* childNodeIndex[2];
    float numHapChild[2];

    StateNode() {
        nodeIndex=0;
        numHap=0;
        childNodeIndex[0]= nullptr;
        childNodeIndex[1]= nullptr;
        numHapChild[0]=0;
        numHapChild[1]=0;
    }

    StateNode(int tIndex, float tNumHap) {
        nodeIndex=tIndex;
        numHap=tNumHap;
        childNodeIndex[0]= nullptr;
        childNodeIndex[1]= nullptr;
        numHapChild[0]=0;
        numHapChild[1]=0;
    }

//    StateNode& operator=(const StateNode& A)
//    {
//        fprintf(stderr,"ever called\n");
//        nodeIndex=A.nodeIndex;
//        numHap=A.numHap;
//        parentNodeIndex2NumHap=A.parentNodeIndex2NumHap;
////        childNodeIndex2NumHap=A.childNodeIndex2NumHap;
//        childNodeIndex[0]=A.childNodeIndex[0];
//        childNodeIndex[1]=A.childNodeIndex[1];
//        numHapChild[0]=A.numHapChild[0];
//        numHapChild[1]=A.numHapChild[1];
//        return *this;
//    }

    StateNode& operator+=(const StateNode& A)
    {
        numHap+=A.numHap;
        for (auto kv:A.parentNodeIndex2NumHap) {
            AddParentNode(kv.first,kv.second);
            if(kv.second->childNodeIndex[0]==&(A.nodeIndex))//A has 0 allele
            {
                kv.second->childNodeIndex[0]=&(this->nodeIndex);
            }
            else if (kv.second->childNodeIndex[1]==&(A.nodeIndex))//A has 1 allele
            {
                kv.second->childNodeIndex[1]=&(this->nodeIndex);
            }
            else{
                fprintf(stderr,"roar from StateNode operator+=!!!!\n%x and %x: %x\n",kv.second->childNodeIndex[0],kv.second->childNodeIndex[1],&(A.nodeIndex));
                exit(EXIT_FAILURE);
            }
        }
//        for (auto kv:A.childNodeIndex2NumHap) {
//            AddChildNode(kv.first,kv.second);
//        }
        return *this;
    }

    void AddParentNode(int index, StateNode* parentAddress) {
            parentNodeIndex2NumHap[index]=parentAddress;
    }

    bool AddChildNode(char allele,int* index, float numHaplotype) {//should only be called once
//        if(index != nullptr) fprintf(stderr,"add %d to allele %d with numHaplotype:%f\n",*index,allele,numHaplotype);
//        else fprintf(stderr,"add nullptr to allele %d with numHaplotype:%f\n",allele,numHaplotype);
        if(childNodeIndex[allele]!=index && childNodeIndex[allele]!= nullptr)
        {
            fprintf(stderr,"roar from StateNode AddChildNode!!!!\n");
            exit(EXIT_FAILURE);
        }
        childNodeIndex[allele]=index;
        numHapChild[allele]+=numHaplotype;
    }

    float GetTransitionProbToChildNode(char allele) {
        return numHapChild[allele];
    }

//    float GetTransitionProbFromParentNode(int index) {
//        return parentNodeIndex2NumHap[index];
//    }
};

class StateNodeContainer
{
public:
    std::vector<std::vector<StateNode*> > StateNodeMat;
    std::vector<StateNode*> tmpNodeVec;
    int nsnps;

    StateNodeContainer();
    StateNodeContainer(int nmarkers):nsnps(nmarkers),StateNodeMat(nmarkers,std::vector<StateNode*>(0,new StateNode()))
    {
    }
    void NormalizeCurrentSiteTransitionProb(int index)
    {
        if(index>=0)
        for (int i = 0; i <StateNodeMat[index].size(); ++i) {

            StateNodeMat[index][i]->numHapChild[0]/= StateNodeMat[index][i]->numHap;
            StateNodeMat[index][i]->numHapChild[1]/= StateNodeMat[index][i]->numHap;
//            fprintf(stderr,"marker:%d\t0:%f\t1:%f\t%f\n",index,StateNodeMat[index][i]->numHapChild[0],StateNodeMat[index][i]->numHapChild[1],StateNodeMat[index][i]->numHap);
        }
    }
};


typedef std::unordered_map<std::string,std::string>  ID2POP;
typedef std::unordered_map<int,std::string>  Index2ID;
typedef std::map<int,std::map<int,bool> > EDGE;

class PBWTWrapper
{
public:
    int N,M;//numSites,numHaps
    int nSamples;
    int nMarkers;

    int prefixLength;

    StateNodeContainer Graph;

    PBWT* pbwtCore;
    PbwtCursor* forwardCursor,*reverseCursor;
    std::vector<std::vector<int> > a,alpha/*reverse*/;
    //std::vector<std::unordered_map<int,int> > aMap,alphaMap;
    std::vector<std::vector<int> > aMap,alphaMap;
    std::vector<std::vector<int> > d,delta;
    std::vector<std::vector<float> > fwdDistance,bkDistance;/*reverse*/;
    std::vector<std::vector<uchar> > sortedY/*only for test*/;
    std::vector<int> c,celta;/*number of zero at each site*/
    std::vector<std::vector<int> > u,ultra;/*relative rank within zeros*/

    std::vector<std::vector<int> > haplotypeCluster;//site, hapID
    std::vector<std::vector<uchar> > clusterAllele;//numCluster;//at each site

    std::vector<std::vector<std::vector<float> > > transVector;	//transition probability: site,from,to
    std::vector<EDGE > inEdges;//valid edges:site, to, from; record indices of states that can arrive at current site and current state
    std::vector<EDGE > outEdges;//valid edges:site, from, to; record number of states that can reach out from current site and current state
    char ** haplotype;//I don't store alleles here, instead I rely on the haplotype storage in libMach
    double * freq1s;

    std::vector<std::vector<int> > clusterMembership;//content is the fwd rank at that site
    std::vector<bool> hasSiblings;

    //mergeSite function variables
    std::vector<std::vector<int> > dist;
    std::priority_queue<max_pair_t,std::vector<max_pair_t>, std::function<bool(const max_pair_t&,const max_pair_t&)> >mergePairList;

    std::unordered_map<int, int> stateOrder;//mapping oldState to newOrder
    int tmpOrder;
    std::vector<uchar> tmpAllele;
    double pval;
    bool EXACT;
    std::unordered_map<int, int> removeMembership;//rankID,state

    //regressionMergeSite function variables
    std::vector<std::vector<float> > leftCoordinate;
    std::vector<std::vector<float> > rightCoordinate;
    std::vector<r_stat_t>  rightCoordinateStat;



    int phased;
    int nSampledCopy;


    PBWTWrapper(){}
    ~PBWTWrapper() {
        if(pbwtCore)
        {
            pbwtDestroy(pbwtCore);

        }
        if(forwardCursor)
        {
            //delete forwardCursor;
            pbwtCursorDestroy(forwardCursor);
        }
        if(reverseCursor)
        {
            //delete forwardCursor;
            pbwtCursorDestroy(reverseCursor);
        }
        //TODO:REVERSE
        //Cannot delete haplotype because haplotype is obtained via SetHap()
//        if(haplotype)
//        {
//            for (int i = 0; i <nSamples ; ++i) {
//                delete [] haplotype[i];
//            }
//            delete [] haplotype;
//        }
        for (int i = 0; i < Graph.StateNodeMat.size(); ++i) {
            for (int j = 0; j <Graph.StateNodeMat[i].size() ; ++j) {
                delete Graph.StateNodeMat[i][j];
            }
            Graph.StateNodeMat[i].clear();
        }

        //KS table
        for (int i = 1; i <=100; ++i) {
            for (int j =i; j <(int)floor(10000.0/i)+1; ++j) {
                delete [] PvalueMatrix[i][j];
            }
            delete [] PvalueMatrix[i];
        }
        delete [] PvalueMatrix;

    }

    void ResetWrapper()
    {

        if(pbwtCore)
        {
            pbwtDestroy(pbwtCore);
            pbwtCore = pbwtCreate(M, N);

        }
        if(forwardCursor)
        {
            //delete forwardCursor;
            pbwtCursorDestroy(forwardCursor);
            forwardCursor = pbwtCursorCreate(pbwtCore, TRUE, TRUE);
        }
        if(reverseCursor)
        {
            pbwtCursorDestroy(reverseCursor);
            reverseCursor = pbwtCursorCreate(pbwtCore, TRUE, TRUE);
        }
        //TODO:REVERSE

        for (int i = 0; i < Graph.StateNodeMat.size(); ++i) {
            for (int j = 0; j <Graph.StateNodeMat[i].size() ; ++j) {
                delete Graph.StateNodeMat[i][j];
            }
            Graph.StateNodeMat[i].clear();
        }
        alphaMap=aMap=a=alpha=d=u=ultra=haplotypeCluster=std::vector<std::vector<int> >(N,std::vector<int>(M,0));
        c=celta=std::vector<int>(N,0);
        clusterAllele=std::vector<std::vector<uchar> >(N,std::vector<uchar>());
    }

    PBWTWrapper(const char ** haps, int nhaps, int nsnps);
    PBWTWrapper(int nhaps,int nsnps);


    int CursorForwards();
    int CursorBackwards();

    int CursorForwardsTo(int k, int T=5);


    int CursorBackwardsTo(int siteBackword, int T=5);

    int CopyHap(int k, PbwtCursor* Cursor);

    bool HasSiblings(int site, int state) {
        if(site ==0) return true;
        for(auto kv:Graph.StateNodeMat[site][state]->parentNodeIndex2NumHap)
        {
//            if(Graph.StateNodeMat[site-1][kv.first].childNodeIndex2NumHap.size()>1) return true;
            if(Graph.StateNodeMat[site-1][kv.first]->childNodeIndex[0]!= nullptr&&Graph.StateNodeMat[site-1][kv.first]->childNodeIndex[1]!= nullptr) return true;
        }
        return false;
    }

    int LabelNoSiblingCluster(int site);
	int UpdateTransVector(int site);

    bool IsEditDistanceOK(int stateA, int stateB, int index, int thresh);
    void MergeSortedArrayToA(std::vector<int> &a, std::vector<int> &b);
    bool IsRecipricalLengthOK(std::vector<int> &a, std::vector<int> &b);
    int CalculateDmax(double & pval, double & Dmax, std::vector<int> & j, std::vector<int>& k);
    int CalculateDmaxBeta(double & pval, double & Dmax, std::vector<int> & j, std::vector<int>& k);

    int MergeAtSite(int site);
    int RegressionMergeAtSite(int site);

    int MoveSegment(const std::unordered_map<int,int>& mergedMembership,int site);

    int SetHaps(char **haps,char **sampledHaps, int nPhase, int nCopy);

    //inline functions
    inline int GetNumStates(int k)
    {
        if(clusterAllele.size()<=k) fprintf(stderr,"site: %d not in clusterAllele\n",k);
        return clusterAllele[k].size();
    }
    inline unsigned long GetNumHaps(int site) const { return haplotypeCluster[site].size(); }

    inline float GetDistanceFromBack(int site, int backRank)const {return bkDistance[site][backRank];}
    inline float GetDistanceFromFwd(int site, int fwdRank)const {return fwdDistance[site][fwdRank];}

    inline int GetHapIDFromBack(int site, int backRank) const { return alpha[site + 1][backRank]; }
    inline int GetHapIDFromFwd(int site, int fwdRank) const { return a[site][fwdRank]; }//you should only use it after a being updated

    inline int GetRankFromBack(int site, int hapID) {return alphaMap[site][hapID];}
    inline int GetRankFromFwd(int site, int hapID) {return aMap[site][hapID];}

    inline int GetHapState(int site, int hapID) { return haplotypeCluster[site][hapID]; }

//    //KS D value related
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
    }

    inline int WritePvalueMatrix()
    {
        std::fstream fout("/Users/fanzhang/Downloads/PlutoTest/PvalueMatrix",std::ios_base::binary|std::ios_base::out);
        if(!fout.is_open())
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
    }

    inline int ReadPvalueMatrix()
    {
        std::fstream fin("/Users/fanzhang/Downloads/PlutoTest/PvalueMatrix",std::ios_base::binary|std::ios_base::in);
        if(!fin.is_open())
        {
            std::cerr<<"open file /Users/fanzhang/Downloads/PlutoTest/PvalueMatrix failed!"<<std::endl;
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
    }
    inline float GetPValue(int n1, int n2, double D)
    {
        if(n1>n2) std::swap(n1,n2);
        return PvalueMatrix[n1][n2][size_t(D*1000)-1];
    }



    //fast update pbwt
    int RemoveIndividualFromPBWT(int individualToProcess);

    int InsertIndividualBackToPBWT(int individualIndex, char** haps);
    /*
    function name: debug functions
    return value: :
    param : :
    author: fanzhang
    time: 8/4/15
    */
    int PrintDistributionAtSite(int state,std::vector<int>& dist);

    template<typename T> void PrintVector(std::vector<T> a,const char* str)
    {
        fprintf(stderr,"debug array %s:\n",str);
        for (int i = 0; i < a.size(); ++i) {
            fprintf(stderr,"%f\t",(float)a[i]);
        }
        fprintf(stderr,"\n");
    }
    template<typename T> inline void PrintMatrix(std::vector<std::vector<T> > a,const char* str)
    {
        fprintf(stderr,"debug matrix size(%d,%d,) %s:\n",a.size(),a[0].size(), str);
        for (int i = 0; i < M; ++i) {
            std::cerr<<i<<"\t";
            for (int j = 0; j < N; ++j) {
                std::cerr<<a[j][i]<<"\t";
            }
            std::cerr<<std::endl;
        }
        std::cerr<<std::endl;
    }
    inline void PrintHap(char** a,std::vector<int>  &b)
    {
        fprintf(stderr,"Hap matrix:\n");
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                std::cerr<<(ushort)a[b[i]][j]<<"\t";
            }
            std::cerr<<std::endl;
        }
        std::cerr<<std::endl;
    }

    int PrintSummary();

    void DoMerge(int site, int retainState, int removeState, std::vector<std::vector<int>> &dist,
                 std::vector<bool, std::allocator<bool>> &removeIndicator,
                 std::vector<bool, std::allocator<bool>> &retainIndicator, std::unordered_map<int, int> &removeMembership);
private:
    PBWTWrapper(const PBWTWrapper&);
    PBWTWrapper & operator=(const PBWTWrapper&);
};

#endif //PLUTO_PBWTWRAPPER_H
