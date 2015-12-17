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
#define DEBUG 1
#include "pbwt/pbwt.h"
#include <vector>
#include <unordered_map>
#include <iostream>
#include <algorithm>

class PBWTWrapper
{
public:
    int N,M;//numSites,numHaps
    PBWT* pbwtCore;
    PbwtCursor* forwardCursor,*reverseCursor;
    std::vector<std::vector<int> > a,alpha/*reverse*/;
    std::vector<std::unordered_map<int,int> > aMap,alphaMap;
    std::vector<std::vector<int> > d;//,delta/*reverse*/;
    std::vector<std::vector<uchar> > sortedY/*only for test*/;
    //std::vector<int> numZero;/*number of zero at each site*/
    std::vector<std::vector<int> > haplotypeCluster;//site, rank
    std::vector<std::vector<uchar> > clusterAllele;//numCluster;//at each site
	//std::unordered_map<int,std::vector<std::vector<float> > > transVector;//transition probability: site,from,to
    std::vector<std::vector<std::vector<float> > > transVector;
    char ** haplotype;//I don't store alleles here, instead I rely on the haplotype storage in libMach
    std::vector<std::vector<int> > clusterMembership;

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
    }

    PBWTWrapper(const char ** haps, int nhaps, int nsnps);
    PBWTWrapper(int nhaps,int nsnps);

    int CursorForwards();
    int CursorBackwards();
    int CursorForwardsTo(int k, int T=5);
    int CursorBackwardsTo(int k, int T=5);
	int ObtainHapFromSinglePhasing(char ** haps);//I implement it here, but not using it for now
    inline int CopyHap(int k, PbwtCursor* Cursor);


	int UpdateTransVector(int site);
    int MergeCluster(int site);

    int MoveSegment(std::vector<std::vector<int> >& MemberShip);
    //bool KStest(std::vector<int>& a,std::vector<int>& b);
    bool KStest(const double& Dmax, const int & sizeA, const int & sizeB);
    int UpdateRankWithinState(std::vector<std::vector<int> > &dist,int stateA, int stateB);

    int SetHaps(char **haps);

    //inline functions
    inline int GetNumStates(int k)
    {
        if(clusterAllele.size()<=k) fprintf(stderr,"site: %d not in clusterAllele\n",k);
        return clusterAllele[k].size();
    }
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
        //TODO:REVERSE
        a=alpha=d=haplotypeCluster=std::vector<std::vector<int> >(pbwtCore->N,std::vector<int>(pbwtCore->M,0));
        clusterAllele=std::vector<std::vector<uchar> >(pbwtCore->N,std::vector<uchar>());
        transVector.clear();
    }
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
            fprintf(stderr,"%d\t",a[i]);
        }
        fprintf(stderr,"\n");
    }
    template<typename T> inline void PrintMatrix(std::vector<std::vector<T> > a,const char* str)
    {
        fprintf(stderr,"debug matrix size(%d,%d,) %s:\n",a.size(),a[0].size(), str);
        for (int i = 0; i < M; ++i) {
            std::cerr<<i<<"\t";
            for (int j = 0; j < N; ++j) {
                std::cerr<<(ushort)a[j][i]<<"\t";
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

    inline float CalDmax(std::vector<int> a, std::vector<int> b)
    {
        int size=a.size()+b.size();
        if(a.size()>b.size())
        {
            std::vector<int> tmp=a;
            a=b;
            b=tmp;
        }
        std::vector<int> A(size,0),B(size,0);
        std::sort(a.begin(),a.end());
        std::sort(b.begin(),b.end());
        int currentA(a[0]),currentB(b[0]);
        for (int totalRank(0),i(0),j(0),cumuNumA(0),cumuNumB(0);totalRank<size;totalRank++) {
            if((currentA<currentB&&i<a.size())||j==b.size())
            {
                cumuNumA=totalRank;
                A[totalRank]=cumuNumA;
                B[totalRank]=cumuNumB;
                i++;
                currentA=a[i];
            }
            else
            {
                cumuNumB=totalRank;
                A[totalRank]=cumuNumA;
                B[totalRank]=cumuNumB;
                j++;
                currentB=b[j];
            }
        }
        float Dmax(0);
        for (int k = 0; k < size; ++k) {
            float tmp=abs(A[k]-B[k]);
            if(tmp>Dmax) Dmax=tmp;
        }

        return Dmax/size;
    }

    int PrintSummary();


};

#endif //PLUTO_PBWTWRAPPER_H
