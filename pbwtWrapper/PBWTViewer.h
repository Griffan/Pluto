//
// Created by Fan Zhang on 1/1/16.
//

#ifndef PLUTO_DEBUGWRAPPER_H
#define PLUTO_DEBUGWRAPPER_H


#include "PBWTWrapper.h"
#include "../TestUnit/MergingEventSimulator.h"

class PBWTViewer: public PBWTWrapper
{

public:

    int numRight;
    int numAltRight;
    std::vector<std::vector<uchar> > clusterAllele;
    PBWTViewer(int a, int b);
//    int CursorForwards(RESULT* result);
//    int CursorForwards(Index2ID& a,ID2POP&b);
//    int CursorForwardsTo(int k, int T = 5, RESULT *result = 0);//for debugging
//    int CursorForwardsTo(int k, Index2ID &a, ID2POP &b,int T=5);
//    int MergeCluster(int site, RESULT* result);
    int MergeCluster(int site, Index2ID &a, ID2POP &b);
    /*The below are implementations for finding maximal length of shared prefix and suffix*/
    int ConfidentOrNot(char** individual, int siteA, int siteB);
    int Phase(char** individual, int siteA, int siteB);
    int MajorVoting(char** individual, char*** voteHapVec, int round, std::vector<int>& heterIndex);
    int FindLengthOfSuffix(char* haplotype, int siteA);
    int FindLengthOfPrefix(char* haplotype, int siteB);
    char** ExtractSubset(int individual);
    int SubsetResume(char** haplotypeX,int individual);
    int Process(int nMarkers, int nSamples, char** haps);
    inline void resetArrays()
    {
        a=alpha=d=u=ultra=haplotypeCluster=std::vector<std::vector<int> >(N,std::vector<int>(M,0));
        c=celta=std::vector<int>(N,0);
        clusterAllele=std::vector<std::vector<uchar> >(N,std::vector<uchar>());
//        transVector.clear();
    }
};


#endif //PLUTO_DEBUGWRAPPER_H
