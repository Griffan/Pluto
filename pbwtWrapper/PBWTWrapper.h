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
//#include "../TestUnit/MergingEventSimulator.h"
#include <queue>
#include <map>


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

typedef std::unordered_map<std::string,std::string>  ID2POP;
typedef std::unordered_map<int,std::string>  Index2ID;

typedef std::map<int,std::map<int,bool> > EDGE;

class PBWTWrapper
{
public:
    int N,M;//numSites,numHaps
    int nSamples;
    int nMarkers;
    PBWT* pbwtCore;
    PbwtCursor* forwardCursor,*reverseCursor;
    std::vector<std::vector<int> > a,alpha/*reverse*/;
    //std::vector<std::unordered_map<int,int> > aMap,alphaMap;
    std::vector<std::vector<int> > aMap,alphaMap;
    std::vector<std::vector<int> > d,delta;/*reverse*/;
    std::vector<std::vector<uchar> > sortedY/*only for test*/;
    std::vector<int> c,celta;/*number of zero at each site*/
    std::vector<std::vector<int> > u,ultra;/*relative rank within zeros*/

    std::vector<std::vector<int> > haplotypeCluster;//site, hapID
    std::vector<std::vector<uchar> > clusterAllele;//numCluster;//at each site

    std::vector<std::vector<std::vector<float> > > transVector;	//transition probability: site,from,to
    std::vector<EDGE > inEdges;//valid edges:site, to, from; record indices of states that can arrive at current site and current state
    std::vector<EDGE > outEdges;//valid edges:site, from, to; record number of states that can reach out from current site and current state
    char ** haplotype;//I don't store alleles here, instead I rely on the haplotype storage in libMach

    std::vector<std::vector<int> > clusterMembership;//content is the fwd rank at that site
    std::vector<bool> hasSiblings;

    //mergeSite function variables
    std::vector<std::vector<int> > dist;
    std::priority_queue<max_pair_t,std::vector<max_pair_t>, std::function<bool(const max_pair_t&,const max_pair_t&)> >mergePairList;
    double tmpABS;
    std::unordered_map<int, int> stateOrder;//mapping oldState to newOrder
    int tmpOrder;
    std::vector<uchar> tmpAllele;
    double pval;
    bool EXACT;
    std::unordered_map<int, int> removeMembership;//rankID,state


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
    }

    PBWTWrapper(const char ** haps, int nhaps, int nsnps);
    PBWTWrapper(int nhaps,int nsnps);


    int CursorForwards();
    int CursorBackwards();

    int CursorForwardsTo(int k, int T=5);


    int CursorBackwardsTo(int k, int T=5);

    int CopyHap(int k, PbwtCursor* Cursor);

    int LabelNoSiblingCluster(int site);
	int UpdateTransVector(int site);

    bool IsEditDistanceOK(const std::vector< std::vector<char*> >& backBone,int stateA, int stateB, int index, double thresh);
    void MergeSortedArrayToA(std::vector<int> &a, std::vector<int> &b);
    bool IsRecipricalLengthOK(std::vector<int> &a, std::vector<int> &b);
    int CalculateDmax(double & pval, double & Dmax, std::vector<int> & j, std::vector<int>& k);

    int MergeAtSite(int site);

    int MoveSegment(const std::unordered_map<int,int>& mergedMembership,int site);

    int SetHaps(char **haps);

    //inline functions
    inline int GetNumStates(int k)
    {
        if(clusterAllele.size()<=k) fprintf(stderr,"site: %d not in clusterAllele\n",k);
        return clusterAllele[k].size();
    }
    inline unsigned long GetNumHaps(int site) const { return haplotypeCluster[site].size(); }

    inline int GetHapIDFromBack(int site, int backRank) const { return alpha[site + 1][backRank]; }
    inline int GetHapIDFromFwd(int site, int fwdRank) const { return a[site][fwdRank]; }//you should only use it after a being updated

    inline int GetRankFromBack(int site, int hapID) {return alphaMap[site][hapID];}
    inline int GetRankFromFwd(int site, int hapID) {return aMap[site][hapID];}

    inline int GetHapState(int site, int hapID) { return haplotypeCluster[site][hapID]; }

//    //KS D value related
//    std::vector<std::vector<double> > DvalueMatrix;//10k X 10k
//    float exact_ks_test_p_val;
//    inline int CalculateDvalueMatrix()
//    {
//        DvalueMatrix=std::vector<std::vector<double> >(10000,std::vector<double>(10000,-1.0));
//        for (int i = 1; i <10000 ; ++i) {
//            for (int j = i; j <floor(10000/i+0.5); ++j) {
//
//                double D=1;
//                for(;D>0;D-=0.01)
//                {
//                    if((1 - psmirnov2x(&D, i, j)) > exact_ks_test_p_val) break;
//                }
//                DvalueMatrix[i][j]=D;
//                DvalueMatrix[j][i]=D;
//            }
//        }
//    }
//    inline int GetExactThresh(int n1, int n2)
//    {
//        return DvalueMatrix[n1][n2];
//    }


    inline void resetWrapper()
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
        alphaMap=aMap=a=alpha=d=u=ultra=haplotypeCluster=std::vector<std::vector<int> >(N,std::vector<int>(M,0));
        c=celta=std::vector<int>(N,0);
        clusterAllele=std::vector<std::vector<uchar> >(N,std::vector<uchar>());
        transVector.clear();
        //aMap=alphaMap=std::vector<std::unordered_map<int,int> >(N,std::unordered_map<int,int>());
        inEdges.clear();
        outEdges.clear();
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
            fprintf(stderr,"%d\t",(int)a[i]);
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
