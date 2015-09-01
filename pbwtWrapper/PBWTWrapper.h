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
#define DEBUG 0
#include "pbwt/pbwt.h"
#include <vector>
#include <unordered_map>
class PBWTWrapper
{
public:
    int N,M;//numSites,numHaps
    PBWT* pbwtCore;
    PbwtCursor* forwardCursor,*reverseCursor;
    std::vector<std::vector<int> > a,alpha/*reverse*/;
    std::vector<std::vector<int> > d;//,delta/*reverse*/;
    //std::vector<int> numZero;/*number of zero at each site*/
    std::vector<std::vector<int> > haplotypeCluster;
    std::vector<std::vector<uchar> > clusterAllele;//numCluster;//at each site
	std::unordered_map<int,std::vector<std::vector<float> > > transVector;//transition probability: site,from,to
    char ** haplotype;//I don't store alleles here, instead I rely on the haplotype storage in libMach

    virtual ~PBWTWrapper() { }

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

    bool KStest(std::vector<int>& a,std::vector<int>& b);

    int setHaps(char ** haps);

    //inline functions
    inline int getNumStates(int k)
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
    template<typename T> inline void PrintVector(T * a,int n, const char* str)
    {
        fprintf(stderr,"debug array %s:\t",str);
        for (int i = 0; i < n; ++i) {
            fprintf(stderr,"%d\t",a[i]);
        }
        fprintf(stderr,"\n");
    }
};

#endif //PLUTO_PBWTWRAPPER_H
