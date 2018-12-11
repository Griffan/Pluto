//
// Created by Fan Zhang on 8/10/18.
//

#ifndef PLUTO_DAG_H
#define PLUTO_DAG_H

#include <unordered_set>
#include <string>
#include <vector>
#include <iostream>

typedef int16_t StateIndex;
typedef unsigned char uchar;
typedef std::unordered_set<StateIndex> ParentSet;
//typedef std::vector<StateIndex> ParentSet;

class StateNode
{
private:
    char allele;
    StateIndex childNode[2]={-1,-1};
    float numHap[2]={0.f,0.f};
    ParentSet parentNodeSet;

    StateNode(const StateNode&);
    StateNode() {
        allele=0;
        numHap[0]=0;
        numHap[1]=0;
        childNode[0]= -1;
        childNode[1]= -1;
    }

public:

    explicit StateNode(char tAllele) {
        allele = tAllele;
        numHap[0] = 0;
        numHap[1] = 0;
        childNode[0] = -1;
        childNode[1] = -1;
    }

    ~StateNode()
    {
        parentNodeSet.clear();
    }

    StateNode & operator=(const StateNode& A)
    {
        numHap[0]=A.numHap[0];
        numHap[1]=A.numHap[1];
        allele=A.allele;
        parentNodeSet=A.parentNodeSet;
        childNode[0]=A.childNode[0];
        childNode[1]=A.childNode[1];
        return *this;
    }

    int AddParentNode(StateIndex parentIndex) {
        auto item=parentNodeSet.insert(/*parentNodeSet.end(),*/parentIndex);
//        if(not parentNodeSet.empty())
//        {
//            for(auto k:parentNodeSet)
//            fprintf(stderr,"inserted %d to set of size %d\n",k,parentNodeSet.size());
//        }
        return item.second;
    }

    int AddChildNode(char allele, StateIndex index) {//should only be called once
        if(childNode[allele]!=index && childNode[allele]!= -1)
        {
            fprintf(stderr,"roar from StateNode AddChildNode!!!!allele:%d\t%d\tto\t%d\n",allele,childNode[allele],index);
            exit(EXIT_FAILURE);
        }
        childNode[allele]=index;
        numHap[allele]+=1.;
        return 0;
    }

    inline char GetAllele()const
    {
        return allele;
    }

    inline void SetAllele(char a)
    {
         allele = a;
    }

    inline StateIndex GetChildNodeIndex(char allele)
    {
        return childNode[allele];
    }

    inline void SetChildNodeIndex(char allele, StateIndex value)
    {
        childNode[allele]=value;
    }

    inline float GetNumHap(char allele)
    {
        return numHap[allele];
    }

    inline void SetNumHap(char allele, float value)
    {
        numHap[allele]=value;
    }

    inline void AddNumHap(char allele, float value)
    {
        numHap[allele]+=value;
    }

    ParentSet& GetParentIndexSet()
    {
        return parentNodeSet;
    }

    std::string ToString()
    {
        std::string line;
        line+=std::to_string(allele)+"\t";
        line+=std::to_string(childNode[0])+"\t";
        line+=std::to_string(childNode[1])+"\t";
        line+=std::to_string(numHap[0])+"\t";
        line+=std::to_string(numHap[1])+"\t";
        line+=std::to_string(parentNodeSet.size())+"\t";
        for (auto k:parentNodeSet) {
            line+=std::to_string(k)+"\t";
        }
        return line;
    }

    int WriteNode(std::ofstream &fout);

    int ReadNode(std::ifstream &fin);
};

class DAG//actual graph
{
private:
    DAG(DAG& A);

    DAG& operator=(DAG& A);
public:
    int nsnps=0;
    int nhaps=0;
    std::vector<std::vector<StateNode*> > StateNodeMat;

    DAG();

    ~DAG();

    DAG(int nmarkers, int nHaps);

    void Reset() {
        for (auto &col:StateNodeMat) {
            for (auto &row:col) {
                delete row;
            }
            col.clear();
        }
    }

    void Clear()
    {
        Reset();
    }

    void AddNode(int k, char allele)
    {
        StateNodeMat[k].push_back(new StateNode(allele));
    }

    void DeleteNode(int site, StateIndex stateM)
    {
        delete StateNodeMat[site][stateM];
        StateNodeMat[site][stateM] = nullptr;
    }

    void AddParentNode(int k, StateIndex group, StateIndex prevSiteStateIndex)
    {
        StateNodeMat[k][group]->AddParentNode(prevSiteStateIndex);
    }

    void AddChildNode(int k, StateIndex prevSiteStateIndex, char allele, StateIndex group)
    {
        StateNodeMat[k][prevSiteStateIndex]->AddChildNode(allele, group);
    }

    StateIndex GetChildNode(int site, StateIndex state, char allele) {
        return StateNodeMat[site][state]->GetChildNodeIndex(allele);
    }

    ParentSet GetParentSet(int site, StateIndex stateM)
    {
        return StateNodeMat[site][stateM]->GetParentIndexSet();
    }

    StateNode* GetNode(int site, StateIndex stateM)
    {
        return StateNodeMat[site][stateM];
    }

    void UpdateNodeVec(int site, std::vector<StateNode*> tmpNodeVec)
    {
        StateNodeMat[site] = tmpNodeVec;
    }

    char GetAllele(int site, StateIndex state)const
    {
        return StateNodeMat[site][state]->GetAllele();
    }

    float GetHapProbAt(int site,int index)
    {
        return (StateNodeMat[site][index]->GetNumHap(0)+StateNodeMat[site][index]->GetNumHap(1))/nhaps;
    }

    float GetTransitionProb(int site, StateIndex from, StateIndex to) {
        char allele = GetAllele(site + 1, to);
        if ((int) StateNodeMat.size() <= site) {
            fprintf(stderr, "site %d doesn't exist!\n", site);
            abort();
        }
        if ((int) StateNodeMat[site].size() <= from) {
//            fprintf(stderr, "site:%d, from:%d states too large!\n", site, from);
            return 0.f;
        }
        if (StateNodeMat[site][from]->GetChildNodeIndex(allele) == -1 ||
            StateNodeMat[site][from]->GetChildNodeIndex(allele) != to) {
//            fprintf(stderr, "site:%d from:%d to:%d states too large!\n", site, from, to);
            return 0.f;
        }
        return GetProbToCurrentNodeConditionalOnParentNode(site, from, allele);//site is for parent
    }

    float GetTransitionFreq(int site, StateIndex from, StateIndex to) {
        char allele = GetAllele(site + 1, to);
        if ((int) StateNodeMat.size() <= site) {
            fprintf(stderr, "site %d doesn't exist!\n", site);
            abort();
        }
        if ((int) StateNodeMat[site].size() <= from) {
//            fprintf(stderr, "site:%d, from:%d states too large!\n", site, from);
            return 0.f;
        }
        if (StateNodeMat[site][from]->GetChildNodeIndex(allele) == -1 ||
            StateNodeMat[site][from]->GetChildNodeIndex(allele) != to) {
//            fprintf(stderr, "site:%d from:%d to:%d states too large!\n", site, from, to);
            return 0.f;
        }
        return StateNodeMat[site][from]->GetNumHap(allele);
    }

    int GetNumStates(int k) const
    {
        return StateNodeMat[k].size();
    }

    bool HasSiblings(int site, StateIndex state) {
        if(site ==0) return true;
        for(auto kv:StateNodeMat[site][state]->GetParentIndexSet())
        {
            if(StateNodeMat[site-1][kv]->GetChildNodeIndex(0)!= -1 && StateNodeMat[site-1][kv]->GetChildNodeIndex(0)!=state) return true;
            if(StateNodeMat[site-1][kv]->GetChildNodeIndex(1)!= -1 && StateNodeMat[site-1][kv]->GetChildNodeIndex(1)!=state) return true;
        }
        return false;
    }

    bool Has2ndRelatives(int site, StateIndex state) {
        if(site ==0 or 1) return true;
        for(auto kv:StateNodeMat[site][state]->GetParentIndexSet())
        {
            if(HasSiblings(site - 1, kv)) return true;//parent has sibs
        }
        return false;
    }

    int JoinNodes(int marker, StateIndex indexRetain, StateIndex indexRemove);

    void UpdateChildNodeIndex(int marker, StateIndex parentIndex, StateIndex newChildIndex, char allele);

    float GetProbToCurrentNodeConditionalOnParentNode(int marker, StateIndex parentIndex, char allele);

    float GetEdgeProbFromParentNode(int marker, StateIndex parentIndex, char allele);

    int WriteDAG(const std::string &fileName);

    int ReadDAG(const std::string &fileName);

    int ToJson(const std::string &fileName);

    int FromJson(const std::string &fileName);
};

#endif //PLUTO_DAG_H
