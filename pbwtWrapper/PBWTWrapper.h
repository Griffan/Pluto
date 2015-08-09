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
class PBWTWrapper
{
public:
    PBWT* pbwtCore;
    PbwtCursor* forwardCursor,*reverseCursor;
    std::vector<std::vector<int> > a,alpha/*reverse*/;
    std::vector<std::vector<int> > d,delta/*reverse*/;
    std::vector<int> numZero;/*number of zero at each site*/
    std::vector<std::vector<int> > haplotypeCluster;
    std::vector<int> numCluster;//at each site
    char ** haplotype;//I don't store alleles here, instead I rely on the haplotype storage in libMach

    virtual ~PBWTWrapper() { }

    PBWTWrapper(PBWT* pbwt) { }
    PBWTWrapper(const char ** haps, int nhaps, int nsnps);
    PBWTWrapper(int nhaps,int nsnps);
    int InitializeCursor(BOOL isForwards, BOOL isStart);
	int InitializeReverseCursor(BOOL isForwards, BOOL isStart);
    int CursorForwards();
    int CursorBackwards();
    int CursorForwardsTo(int k, int T=5);
    int CursorBackwardsTo(int k, int T=5);
	int ObtainHapFromSinglePhasing(char ** haps);//I implement it here, but not using it for now
    inline int CopyHap(int k, PbwtCursor* Cursor);

    int MergeCluster(int site);
    bool KStest(std::vector<int>& a,std::vector<int>& b);

    int setHaps(char ** haps);

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
