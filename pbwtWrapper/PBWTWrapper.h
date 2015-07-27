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
#include <vector.h>
class PBWTWrapper
{
    PWBT* pbwtCore,pbwtCoreReverse;
    PbwtCursor* forwardCursor,backwardCursor;
    vector<vector<int> > a,alpha/*reverse*/;
    vector<vector<int> > d,delta/*reverse*/;
    vector<int> numZero;/*number of zero at each site*/
    vector<vector<int> > haplotypeCluster;

public:
    virtual ~PBWTWrapper() { }

public:
    PBWTWrapper(const ::PBWTWrapper::PWBT* pbwt) { }
    PBWTWrapper(const char ** haplotype, int nhaps, int nsnps);
    PBWTWrapper(int nhaps,int nsnps);
    int InitializeCursor();
    int InitializeReverseCursor();
    int CursorForwards();
    int CursorBackwards();
    int CursorForwardsTo();
    int CursorBackwardsTo();
    int PrintDistributionAtSite();
    int IdentifyGroup();//based on content of prefix, acutally this is the function that find set maximal up to length L
    int ObtainRank();//based on content of suffix, this is the function that find rank order of haplotype
};

#endif //PLUTO_PBWTWRAPPER_H
