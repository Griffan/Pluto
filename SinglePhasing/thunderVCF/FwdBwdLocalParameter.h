//
// Created by Fan Zhang on 7/18/18.
//

#ifndef PLUTO_FWDBWDLOCALPARAMETER_H
#define PLUTO_FWDBWDLOCALPARAMETER_H

//    typedef std::unordered_map<std::pair<StateIndex, StateIndex>, float, PairHash> Source;//(nodeA,nodeB)->fwd
//    typedef std::unordered_map<std::pair<StateIndex, StateIndex>, Source, PairHash> DestToSource;
//    typedef std::unordered_map<uint64_t, float> Source;//(nodeA,nodeB)->fwd
//    typedef std::unordered_map<uint64_t, Source> DestToSource;

#include <cstdint>
#include <vector>
#include "../../abseil-cpp/absl/container/flat_hash_map.h"

typedef int16_t StateIndex;
typedef uint64_t NodePair;
typedef std::vector<NodePair> SourceVec;
typedef std::vector<std::vector<NodePair> > Index2SourceVec;

typedef std::vector<float> FwdVec;//associated via index with SourceVec
typedef std::vector<float> BwdVec;//associated via index with SourceVec

typedef std::vector<FwdVec> Index2FwdVec;//associated via index with Index2SourceVec
typedef std::vector<BwdVec> Index2BwdVec;//associated via index with Index2SourceVec


//typedef std::unordered_map<NodePair, int> Dest2SourceVecIndex;//NodePair -> StorageIndex

typedef absl::flat_hash_map<NodePair, int> Dest2SourceVecIndex;//NodePair -> StorageIndex

#define HASH_RESERVE 4096

class FwdBwdLocalParameter {

public:
    int states;
    int markers;
    //from (parentNode1, parentNode2) to (childNode1, childNode2), childNodes are present in conditional graph, but parentNode are not necessarily present
    std::vector<Dest2SourceVecIndex> parentsNodeVec;//edge information stored part 1
    std::vector<Index2SourceVec> megaSourceVec;//edge information stored part 2

    std::vector<Index2FwdVec> megaFwdVec;//StorageIndex -> FwdVec
    std::vector<float> fwdValueSum;
    std::vector<absl::flat_hash_map<int, float> > fwdValueNode1Sum;
    std::vector<absl::flat_hash_map<int, float> > fwdValueNode2Sum;
    std::vector<bool> isRec;

    //backward value
    absl::flat_hash_map<NodePair, float> finalBwdValue;//current marker -> next marker -> bwdValue
    float bwdValueSum;
    absl::flat_hash_map<int, float> bwdValueNode1Sum;
    absl::flat_hash_map<int, float> bwdValueNode2Sum;
    absl::flat_hash_map<NodePair, float> bwdValue;


    explicit FwdBwdLocalParameter(int markers);

    inline NodePair MakePair(StateIndex first, StateIndex second) {
        return static_cast<NodePair>(first) << 32 |  static_cast<NodePair>(second);
    }

    StateIndex GetFirst(NodePair pair);

    StateIndex GetSecond(NodePair pair);

    int FillParentsNodeVec(int i, StateIndex childNode1, StateIndex childNode2, StateIndex parentNode1,
                           StateIndex parentNode2, float tmpFwdValue);

    int ClearParentsNodeVec(int i);

    int ResetParentsNodeVec();

    FwdVec &GetFwdVec(int i, int destIndex);

    FwdVec &GetFwdVec(int i, StateIndex A, StateIndex B);

    float GetFwd(int i, int destIndex, int sourceIndex);

    SourceVec &GetSourceVec(int i, int destIndex);

    SourceVec &GetSourceVec(int i, StateIndex A, StateIndex B);

    int GetDestIndex(int i, StateIndex A, StateIndex B);

    NodePair GetSource(int i, int destIndex, int sourceIndex);

    NodePair GetSource(int i, StateIndex A, StateIndex B, int sourceIndex);

//    template <class valueVec>
    float GetSumValueFromContainer(FwdVec &a);

    float GetSumFwdValueFrom(int i, StateIndex A, StateIndex B);

    float GetSumBwdValueFrom(StateIndex A, StateIndex B);

    std::vector<std::vector<NodePair> > viableStatePair;

    int AddStatePair(int site, StateIndex A, StateIndex B)
    {
        viableStatePair[site].push_back(MakePair(A,B));
        return 0;
    }

    std::vector<NodePair> GetStatePair(int site)
    {
        return viableStatePair[site];
    }


    int ResetBwdValue();

};


#endif //PLUTO_FWDBWDLOCALPARAMETER_H
