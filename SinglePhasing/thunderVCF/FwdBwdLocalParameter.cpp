//
// Created by Fan Zhang on 7/18/18.
//

#include "FwdBwdLocalParameter.h"
#include <iostream>

FwdBwdLocalParameter::FwdBwdLocalParameter(int nmarkers)
{

    states = -1;
    markers = nmarkers;

    Dest2SourceVecIndex dummy;
    dummy.reserve(HASH_RESERVE);
    parentsNodeVec.assign(markers, dummy);

    Index2SourceVec dummy2;
    dummy2.reserve(HASH_RESERVE);
    megaSourceVec.assign(markers, dummy2);

    Index2FwdVec dummy3;
    dummy3.reserve(HASH_RESERVE);
    megaFwdVec.assign(markers, dummy3);

    viableStatePair.assign(markers, std::vector<NodePair>());
}

StateIndex FwdBwdLocalParameter::GetFirst(NodePair pair)
{
    return (StateIndex)(pair >> 32);
}

StateIndex FwdBwdLocalParameter::GetSecond(NodePair pair)
{
    return (StateIndex)(pair & 0xffff);
}

int FwdBwdLocalParameter::FillParentsNodeVec(int i, StateIndex childNode1, StateIndex childNode2, StateIndex parentNode1,
    StateIndex parentNode2, float tmpFwdValue)
{
    int index;
    if (parentsNodeVec[i].find(MakePair(childNode1, childNode2)) != parentsNodeVec[i].end()) //Dest already exists
    {
        index = parentsNodeVec[i][MakePair(childNode1, childNode2)];
        megaFwdVec[i][index].push_back(tmpFwdValue);
        megaSourceVec[i][index].push_back(MakePair(parentNode1, parentNode2));
    } else {
        index = megaFwdVec[i].size();
        parentsNodeVec[i][MakePair(childNode1, childNode2)] = index;

        FwdVec dummyFwdVec(1, tmpFwdValue);
        dummyFwdVec.reserve(HASH_RESERVE);
        megaFwdVec[i].push_back(dummyFwdVec);

        SourceVec dummySourceVec(1, MakePair(parentNode1, parentNode2));
        dummySourceVec.reserve(HASH_RESERVE);
        megaSourceVec[i].push_back(dummySourceVec);
    }
    return 0;
}

int FwdBwdLocalParameter::ClearParentsNodeVec(int i)
{

    parentsNodeVec[i].clear();
    megaFwdVec[i].clear();
    megaSourceVec[i].clear();
    return 0;
}

int FwdBwdLocalParameter::ResetParentsNodeVec()
{
    for (int i = 0; i < markers; ++i) {
        ClearParentsNodeVec(i);
    }
    return 0;
}
FwdVec& FwdBwdLocalParameter::GetFwdVec(int i, int destIndex)
{
    return megaFwdVec[i][destIndex];
}

// std::unordered_map<NodePair , float> &FwdBwdLocalParameter::GetBwdHash(StateIndex A, StateIndex B) {
//    return finalBwdValue[MakePair(A, B)];
//}

FwdVec& FwdBwdLocalParameter::GetFwdVec(int i, StateIndex A, StateIndex B)
{
    return GetFwdVec(i, parentsNodeVec[i][MakePair(A, B)]);
}

float FwdBwdLocalParameter::GetFwd(int i, int destIndex, int sourceIndex)
{
    return megaFwdVec[i][destIndex][sourceIndex];
}

SourceVec& FwdBwdLocalParameter::GetSourceVec(int i, int destIndex)
{
    return megaSourceVec[i][destIndex];
}

SourceVec& FwdBwdLocalParameter::GetSourceVec(int i, StateIndex A, StateIndex B)
{
    return GetSourceVec(i, parentsNodeVec[i][MakePair(A, B)]);
}

int FwdBwdLocalParameter::GetDestIndex(int i, StateIndex A, StateIndex B)
{
    return parentsNodeVec[i][MakePair(A, B)];
}

NodePair FwdBwdLocalParameter::GetSource(int i, int destIndex, int sourceIndex)
{
    return megaSourceVec[i][destIndex][sourceIndex];
}

NodePair FwdBwdLocalParameter::GetSource(int i, StateIndex A, StateIndex B, int sourceIndex)
{
    return megaSourceVec[i][GetDestIndex(i, A, B)][sourceIndex];
}

//template <class valueVec>
float FwdBwdLocalParameter::GetSumValueFromContainer(FwdVec& a)
{
    if (a.size() == 0)
        return 0.f;
    float sum(0.f);
    for (auto kv : a) {
        sum += kv;
    }
    return sum;
}

float FwdBwdLocalParameter::GetSumFwdValueFrom(int i, StateIndex A, StateIndex B)
{
    return GetSumValueFromContainer(GetFwdVec(i, A, B));
}

// float FwdBwdLocalParameter::GetSumBwdValueFrom(StateIndex A, StateIndex B)
//{
//    return GetSumValueFromContainer(GetBwdHash(A, B));//same structure, different alias
//}
int FwdBwdLocalParameter::ResetBwdValue()
{
    bwdValueSum = 0;
    bwdValueNode1Sum.clear();
    bwdValueNode2Sum.clear();
    bwdValue.clear();

    return 0;
}
