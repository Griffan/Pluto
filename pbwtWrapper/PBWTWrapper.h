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
#include "../pbwt/pbwt.h"
#include <vector>
class PBWTWrapper
{
    PBWT* pbwtCore;
    PbwtCursor* forwardCursor,backwardCursor;
    std::vector<std::vector<int> > a,alpha/*reverse*/;
    std::vector<std::vector<int> > d,delta/*reverse*/;
    std::vector<int> numZero;/*number of zero at each site*/
    std::vector<std::vector<int> > haplotypeCluster;
    std::vector<int> numCluster;//at each site
    char ** haplotype;//I don't store alleles here, instead I rely on the haplotype storage in libMach
public:
    virtual ~PBWTWrapper() { }

public:
    PBWTWrapper(PBWT* pbwt) { }
    PBWTWrapper(const char ** haps, int nhaps, int nsnps);
    PBWTWrapper(int nhaps,int nsnps);
    int InitializeCursor();
    int InitializeReverseCursor();
    int CursorForwards();
    int CursorBackwards();
    int CursorForwardsTo(int k, int T=5);
    int CursorBackwardsTo(int k, int T=5);
    inline int CopyHap(int k, PbwtCursor* Cursor);
    int PrintDistributionAtSite();
    int IdentifyGroup();//based on content of prefix, acutally this is the function that find set maximal up to length L
    int ObtainRank();//based on content of suffix, this is the function that find rank order of haplotype
};

#endif //PLUTO_PBWTWRAPPER_H
