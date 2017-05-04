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

typedef int16_t StateIndex;

struct MaxPair
{
    StateIndex clusterA;
    StateIndex clusterB;
    double Dmax;
    bool exact;
    double pval;
    double f_id;//float id for breaking tie
    MaxPair(int a,int b,double c, bool d,double e,double id)
    {
        clusterA=a;
        clusterB=b;
        Dmax=c;
        exact=d;
        pval=e;
        f_id=id;
    }
};

inline bool comparator(const MaxPair &lhs, const MaxPair &rhs) {
    if(lhs.pval == rhs.pval)
    {
        if(lhs.Dmax==rhs.Dmax) {
            return lhs.f_id<rhs.f_id;
        }
        return lhs.Dmax<rhs.Dmax;
    }
    return lhs.pval < rhs.pval;
}


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

typedef unsigned char uchar;
//define NODE structure which represents cluster or state in HMM model
class StateNode
{
private:
    StateNode(const StateNode&);
    StateNode() {
        numHap[0]=0;
        numHap[1]=0;
        childNodeIndex[0]= -1;
        childNodeIndex[1]= -1;
        allele=0;
    }

public:

    //float numHap;
    uchar allele;
    std::vector<StateIndex> parentNodeIndex;
    StateIndex childNodeIndex[2];
    float numHap[2];


    StateNode(uchar tAllele) {
        numHap[0]=0;
        numHap[1]=0;
        childNodeIndex[0]= -1;
        childNodeIndex[1]= -1;
        allele=tAllele;
    }

    ~StateNode()
    {
        parentNodeIndex.clear();
    }

    StateNode & operator=(const StateNode& A)
    {
        fprintf(stderr,"using operator =\n");
        numHap[0]=A.numHap[0];
        numHap[1]=A.numHap[1];
        allele=A.allele;
        parentNodeIndex=A.parentNodeIndex;
        childNodeIndex[0]=A.childNodeIndex[0];
        childNodeIndex[1]=A.childNodeIndex[1];
        return *this;
    }

    int AddParentNode(StateIndex parentIndex) {
        parentNodeIndex.push_back(parentIndex);
        return 0;
    }

    int AddChildNode(uchar allele, StateIndex index) {//should only be called once
//        if(childNodeIndex[(size_t)allele] != nullptr)
//            fprintf(stderr,"add index: %p to child %p with allele %d(%d) with numHaplotype:%f to\n",index,childNodeIndex[(size_t)allele],allele,(size_t)allele,numHaplotype);
//        else fprintf(stderr,"add nullptr to allele %d with numHaplotype:%f\n",allele,numHaplotype);
        if(childNodeIndex[allele]!=index && childNodeIndex[allele]!= -1)
        {
            fprintf(stderr,"roar from StateNode AddChildNode!!!!allele:%d\t%d\tto\t%d\n",allele,childNodeIndex[allele],index);
            exit(EXIT_FAILURE);
        }
        childNodeIndex[allele]=index;
        numHap[allele]+=1.;
        return 0;
    }

    std::vector<StateIndex>& GetParentIndexVector()
    {
        return parentNodeIndex;
    }


    int writeNode(std::ofstream & fout);
    int readNode(std::ifstream & fin);
};

class StateNodeContainer//actual graph
{
private:
    StateNodeContainer(StateNodeContainer& A)
    {
        nsnps=A.nsnps;
        StateNodeMat=A.StateNodeMat;
        tmpNodeVec=A.tmpNodeVec;
    }
    StateNodeContainer& operator=(StateNodeContainer& A);
public:
    int nsnps;
    std::vector<std::vector<StateNode*> > StateNodeMat;
    std::vector<StateNode*> tmpNodeVec;


    StateNodeContainer()
    {
        nsnps=0;
    }
    ~StateNodeContainer()
    {
        for (auto &vec:StateNodeMat) {
            for (auto & node:vec) {
                delete node;
            }
            vec.clear();
        }
    }

    StateNodeContainer(int nmarkers):nsnps(nmarkers),
                                     StateNodeMat(nmarkers,std::vector<StateNode*>(0, new StateNode(0)))
    {
    }

    int JoinNodes(int marker, StateIndex indexRetain, StateIndex indexRemove)
    {
        StateNodeMat[marker][indexRetain]->numHap[0]+=StateNodeMat[marker][indexRemove]->numHap[0];
        StateNodeMat[marker][indexRetain]->numHap[1]+=StateNodeMat[marker][indexRemove]->numHap[1];
        for (auto kv:StateNodeMat[marker][indexRemove]->parentNodeIndex) {
            StateNodeMat[marker][indexRetain]->AddParentNode(kv);
//            fprintf(stderr,"allele 0:%x, allele 1:%x, address:%x\n",kv.second->childNodeIndex[0],kv.second->childNodeIndex[1],&(A.nodeIndex));
            if(StateNodeMat[marker][indexRemove]->allele==0)//Remove has 0 allele
            {
                StateNodeMat[marker-1][kv]->childNodeIndex[0]=indexRetain;
            }
            else if (StateNodeMat[marker][indexRemove]->allele==1)//Remove has 1 allele
            {
                StateNodeMat[marker-1][kv]->childNodeIndex[1]=indexRetain;
            }
            else{
                fprintf(stderr,"roar from StateNode operator+=!!!!\n%d and %d: %d\n",StateNodeMat[marker-1][kv]->childNodeIndex[0],StateNodeMat[marker-1][kv]->childNodeIndex[1],indexRemove);
                exit(EXIT_FAILURE);
            }
        }
        return 0;
    }

    void UpdateChildNodeIndex(int marker, StateIndex parentIndex, StateIndex newChildIndex, uchar allele)
    {
        StateNodeMat[marker][parentIndex]->childNodeIndex[allele]=newChildIndex;
    }

    float GetProbToCurrentNodeConditionalOnParentNode(int marker, StateIndex parentIndex, uchar allele) {
        return StateNodeMat[marker][parentIndex]->numHap[allele]/(StateNodeMat[marker][parentIndex]->numHap[0]+StateNodeMat[marker][parentIndex]->numHap[1]);
    }

    int writeContainer(const std::string& fileName);
    int readContainer(const std::string& fileName);
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

    float *** PvalueMatrix;//10k X 10k
    inline float GetPValue(int n1, int n2, double D)
    {
        if(n1>n2) std::swap(n1,n2);
        int a=(int)floor(D*1000+0.5);
        return PvalueMatrix[n1][n2][(a<1?1:a)-1];
    }

    int prefixLength;

    StateNodeContainer Graph;

//    PBWT* pbwtCore;
//    PbwtCursor* forwardCursor,*reverseCursor;
    std::vector<int> a,alpha/*reverse*/;//stores array a status after process current column haps;
    std::vector<std::vector<int> > alphaMap;
    std::vector<int> d,delta;
    std::vector<std::vector<int> > allDelta;
    std::vector<std::vector<float> > bkDistance;/*reverse*/;
    std::vector<int> sortedY;
    std::vector<int> c,celta;/*number of zero at each site*/
    //std::vector<std::vector<int> > u,ultra;/*relative rank within zeros*/

    std::vector<std::vector<StateIndex> > haplotypeCluster;//site, hapID
    //std::vector<std::vector<int> > bkHaplotypeCluster;//site, hapID


    char ** haplotype;//I don't store alleles here, instead I rely on the haplotype storage in libMach

    std::vector<std::vector<int> > clusterMembership;//content is the fwd rank at that site


    //mergeSite function variables
    std::vector<float> recomRate;
    std::vector<std::vector<int> > dist;
    std::priority_queue<MaxPair,std::vector<MaxPair>, std::function<bool(const MaxPair&,const MaxPair&)> >mergePairList;

    std::unordered_map<StateIndex, StateIndex> stateOrder;//mapping oldState to newOrder
    StateIndex tmpOrder;
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
//        if(pbwtCore)
//        {
//            pbwtDestroy(pbwtCore);
//
//        }
//        if(forwardCursor)
//        {
//            //delete forwardCursor;
//            pbwtCursorDestroy(forwardCursor);
//        }
//        if(reverseCursor)
//        {
//            //delete forwardCursor;
//            pbwtCursorDestroy(reverseCursor);
//        }
        //TODO:REVERSE
        //Cannot delete haplotype because haplotype is obtained via SetHap()
        if(haplotype)
        {
//            for (int i = 0; i <nSamples ; ++i) {
//                delete [] haplotype[i];
//            }
            delete [] haplotype;
        }
        for (auto &col: Graph.StateNodeMat) {
            for (auto & row:col) {
                delete row;
            }
            col.clear();
        }

    }

    void ResetWrapper()
    {

//        if(pbwtCore)
//        {
//            pbwtDestroy(pbwtCore);
//            pbwtCore = pbwtCreate(M, N);
//
//        }
//        if(forwardCursor)
//        {
//            //delete forwardCursor;
//            pbwtCursorDestroy(forwardCursor);
//            forwardCursor = pbwtCursorCreate(pbwtCore, TRUE, TRUE);
//        }
//        if(reverseCursor)
//        {
//            pbwtCursorDestroy(reverseCursor);
//            reverseCursor = pbwtCursorCreate(pbwtCore, TRUE, TRUE);
//        }
        //TODO:REVERSE

        for (auto & col:Graph.StateNodeMat) {
            for (auto &row:col) {
                delete row;
            }
            col.clear();
        }
        a=std::vector<int>(N,0);
        delta=d=alpha=a;
        alphaMap=std::vector<std::vector<int> >(N,std::vector<int>(M,0));
        haplotypeCluster=std::vector<std::vector<StateIndex > >(N,std::vector<StateIndex>(M,0));
        c=celta=std::vector<int>(N,0);
        sortedY=std::vector<int>(M,0);
//        clusterAllele=std::vector<std::vector<uchar> >(N,std::vector<uchar>());
    }

    void ReleaseWrapperMemory()
    {
        a.clear();
        delta=d=alpha=a;
        alphaMap.clear();
        haplotypeCluster.clear();
        c.clear();
        celta.clear();
        sortedY.clear();
        bkDistance.clear();
        allDelta.clear();
    }

//    PBWTWrapper(const char ** haps, int nhaps, int nsnps);
    PBWTWrapper(int nhaps,int nsnps, float *** _PvalueMatrix);


    int CursorForwards();
    int CursorBackwards();

    int CursorForwardsToDeprecated(int k, int T);
    int CursorForwardsTo(int k, int T);


    int CursorBackwardsTo(int siteBackword, int T=5);

    int CopyHap(int k,std::vector<int>& tmpA);

    bool HasSiblings(int site, StateIndex state) {
        if(site ==0) return true;
        for(auto kv:Graph.StateNodeMat[site][state]->parentNodeIndex)
        {
            if(Graph.StateNodeMat[site-1][kv]->childNodeIndex[0]!= -1 && Graph.StateNodeMat[site-1][kv]->childNodeIndex[0]!=state) return true;
            if(Graph.StateNodeMat[site-1][kv]->childNodeIndex[1]!= -1 && Graph.StateNodeMat[site-1][kv]->childNodeIndex[1]!=state) return true;
        }
        return false;
    }

//    int LabelNoSiblingCluster(int site);
//	int UpdateTransVector(int site);

    bool IsEditDistanceOK(int stateA, int stateB, int index, int thresh);
    bool IsInSameBackCluster(int stateA, int stateB, int site, int error_thresh);

    void MergeSortedArrayToA(std::vector<int> &a, std::vector<int> &b);
    bool IsRecipricalLengthOK(std::vector<int> &a, std::vector<int> &b);
    int CalculateDmax(double & pval, double & Dmax, std::vector<int> & j, std::vector<int>& k);
    int CalculateDmaxBeta(double & pval, double & Dmax, std::vector<int> & j, std::vector<int>& k);

    int MergeAtSite(int site);
    int RegressionMergeAtSite(int site,bool isBaseWrapper);

    int MoveSegment(const std::unordered_map<int,int>& mergedMembership,int site);

    int SetHaps(char **haps, int CopyStart, int CopyEnd, char **sampledHaps, int CopyStart2, int CopyEnd2, float* rate);
    //inline functions
    inline int GetNumStates(int k) const
    {
//        if(clusterAllele.size()<=k) fprintf(stderr,"site: %d not in clusterAllele\n",k);
//        return clusterAllele[k].size();
        return Graph.StateNodeMat[k].size();
    }
    inline uchar GetAllele(int site, StateIndex state)
    {
        return Graph.StateNodeMat[site][state]->allele;
    }
    inline unsigned long GetNumHaps(int site) const { return haplotypeCluster[site].size(); }

    inline float GetDistanceFromBack(int site, int backRank)const {return bkDistance[site][backRank];}
    //inline float GetDistanceFromFwd(int site, int fwdRank)const {return fwdDistance[site][fwdRank];}

    //inline int GetHapIDFromBack(int site, int backRank) const { return alpha[site][backRank]; }
    inline int GetHapIDFromFwd(int fwdRank) const { return a[fwdRank]; }//you should only use it after a being updated

    inline int GetRankFromBack(int site, int hapID) {return alphaMap[site][hapID];}
    //inline int GetRankFromFwd(int site, int hapID) {return aMap[site][hapID];}

    inline StateIndex GetHapStateFromFwd(int site, int hapID) { return haplotypeCluster[site][hapID]; }
    //inline int GetHapStateFromBack(int site, int hapID) { return bkHaplotypeCluster[site][hapID]; }

    inline float GetHapProbAt(int site,int index)
    {
        return (Graph.StateNodeMat[site][index]->numHap[0]+Graph.StateNodeMat[site][index]->numHap[1])/M;
    }

    //fast update pbwt
    //int FastCursorForwards(const PBWTWrapper& motherWrapper);
    //int FastCursorForwardsTo(int k, int T, const PBWTWrapper& baseWrapper);
    /*
    function name: debug functions
    return value: :
    param : :
    author: fanzhang
    time: 8/4/15
    */
    int PrintDistributionAtSite(int state,std::vector<int>& dist);

    template<typename T> void PrintVector(T* a, int size, const char* str)
    {
        fprintf(stderr,"debug array %s:\n",str);
        for (int i = 0; i < size; ++i) {
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

    void DoMerge(int site, StateIndex retainState, StateIndex removeState, std::vector<std::vector<int>> &dist,
                 std::vector<bool, std::allocator<bool>> &removeIndicator,
                 std::vector<bool, std::allocator<bool>> &retainIndicator, std::unordered_map<int, int> &removeMembership);
    int HowManyChildlessState(std::vector<StateNode*> & a)
    {
        int sum=0;
        std::cerr<<sum<<" enter childless state out of "<<a.size()<<" total states"<<std::endl;
        for (auto & node:a) {
            if(node->childNodeIndex[0]== -1 || node->childNodeIndex[1]== -1)
                sum++;
        }
        std::cerr<<sum<<" leave childless state out of "<<a.size()<<" total states"<<std::endl;
        return sum;
    }
    int HowManyChildHapCount(std::vector<StateNode*> & currentNodes, int childSite)
    {
        int sum1(0),sum2(0);
        std::cerr<<sum1<<" enter child haps from sum1 and  "<<sum2<<" child haps from sum2 out of "<<currentNodes.size()<<" total states"<<std::endl;
        int index=0;
        std::unordered_map<int,bool> bag;
        for (auto & node:currentNodes) {
            if(node->childNodeIndex[0]!= -1 )
            {
                index=node->childNodeIndex[0];
                if(bag.find(index) ==bag.end()) bag[index]=true;
                else continue;
                sum1+=Graph.StateNodeMat[childSite][index]->numHap;

            }
            if(node->childNodeIndex[1]!= -1)
            {
                index=node->childNodeIndex[1];
                if(bag.find(index) ==bag.end()) bag[index]=true;
                else continue;
                sum1+=Graph.StateNodeMat[childSite][index]->numHap;

            }

        }
        std::cerr<<sum1<<" leave child haps from sum1 and  "<<sum2<<" child haps from sum2 out of "<<currentNodes.size()<<" total states"<<std::endl;
        return sum1;
    }
private:
    PBWTWrapper(const PBWTWrapper&);
    PBWTWrapper & operator=(const PBWTWrapper&);

    void UpdateAandD(int k);

    int CreateNewCluster(int k, int rank, int i0, StateIndex group);

    void CreateLastSiteCluster(int k, int rank, int i0, StateIndex group);
};

#endif //PLUTO_PBWTWRAPPER_H
