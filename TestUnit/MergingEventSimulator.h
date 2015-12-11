//
// Created by Fan Zhang on 12/9/15.
//

#ifndef PLUTO_MERGINGEVENTSIMULATOR_H
#define PLUTO_MERGINGEVENTSIMULATOR_H

#include <string>
#include <vector>
using namespace std;
typedef vector<int> Branch;
typedef vector<pair<Branch, Branch> > Partition;//Partition of individuals per marker
typedef vector<Partition> DAG;

class MergingEventSimulator {
    int nHaps;
    int markerIndex;
    string originalString;
    vector<string> finalStringArray;
    DAG dag;
    vector<bool> hasMutation;
    vector<pair<int, Branch> > eventRecorder;
public:
    MergingEventSimulator();
    char** ChangeStringArraytoCstringArray();
    void SimulateCoalescentTree();
    int FlipHaps(int start,int end,vector<string>&tmpStringArray, bool op, Branch& lastA);
    void PrintSimulatedEvent();
    void PutSimulatedEventIntoStringArray();
    void PrintTree();
    void PrintDAG();
    inline int GetNumHaps()
    {
        return nHaps;
    }
    inline int GetHapLength()
    {
        return originalString.length();
    }

};


#endif //PLUTO_MERGINGEVENTSIMULATOR_H
