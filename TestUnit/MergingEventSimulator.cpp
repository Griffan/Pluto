//
// Created by Fan Zhang on 12/9/15.
//

#include "MergingEventSimulator.h"
#include "random"
#include <iostream>
#define PREFIX 5
MergingEventSimulator::MergingEventSimulator() {
    originalString="000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    nHaps=600;
    markerIndex=PREFIX;
    finalStringArray=vector<string>(nHaps,originalString);
    Branch a;
    Partition tmp=Partition(0,make_pair(a,a));
    dag=DAG(originalString.length(),tmp);
    hasMutation=vector<bool>(originalString.length(),false);
}

char** MergingEventSimulator::ChangeStringArraytoCstringArray() {
    char** b=new char* [finalStringArray.size()];
    for (int i = 0; i <finalStringArray.size() ; ++i) {
        b[i]= new char [originalString.size()];
        for (int j = 0; j < originalString.size(); ++j) {
            b[i][j]=finalStringArray[i][j];
        }
//        memcpy(b[i],finalStringArray[i].c_str(),finalStringArray[i].length()*sizeof(char));
    }
    return b;
}

void MergingEventSimulator::SimulateCoalescentTree() {
    long long int seed = std::chrono::system_clock::now().time_since_epoch().count();
    srand(seed);
            int numFlipped = rand() % (nHaps-1) + 1;
            vector<string> tmpStringArray=finalStringArray;
            //vector<int> indexFlipped, indexNotFlipped;
           // GenerateRandomnVector(numFlipped, indexFlipped, indexNotFlipped);
    Branch a,b;
    int start=0;
    for (int i = start + 0; i <start + numFlipped ; ++i) {
        a.push_back(i);
    }
    for (int j = start + numFlipped; j <start + tmpStringArray.size() ; ++j) {
        b.push_back(j);
    }
    dag[markerIndex].push_back(make_pair(a,b));
            int oldMarkerIndex=markerIndex;
            FlipHaps(0,numFlipped,tmpStringArray, true,a);
            markerIndex=oldMarkerIndex;
            FlipHaps(numFlipped,nHaps,tmpStringArray, false,b);
    finalStringArray=tmpStringArray;
}

int MergingEventSimulator::FlipHaps(int start,int end, vector<string>&tmpStringArray, bool op, Branch& lastA) {

    if(op) {
        for (int i = start; i < end; ++i) {
            if (i >= tmpStringArray.size()) cerr << "bigger than maximal index" << endl;
            tmpStringArray[i][markerIndex] = '1';
        }
        hasMutation[markerIndex]=true;

    }
    markerIndex++;
    int stringArraySize=end-start;
    if(markerIndex>=originalString.length()) return 0;
    else {

        while(markerIndex<originalString.length()) {
            int randNum = rand();
            if (randNum % 100 > 80) {
                break;
            }
            else
                dag[markerIndex].push_back(make_pair(lastA,lastA));
            markerIndex++;
        }
        if(markerIndex >= originalString.length()) return 0;
        if(stringArraySize>1) {
                    int numFlipped = rand() % (stringArraySize - 1) + 1;
                    int oldMarkerIndex=markerIndex;
                    Branch a,b;
                    for (int i = start + 0; i <start + numFlipped ; ++i) {
                        a.push_back(i);
                    }
                    for (int j = start + numFlipped; j <start + stringArraySize ; ++j) {
                        b.push_back(j);
                    }
                    dag[markerIndex].push_back(make_pair(a,b));
                    FlipHaps(start + 0, start + numFlipped, tmpStringArray, true, a);
                    markerIndex=oldMarkerIndex;
                    FlipHaps(start + numFlipped, start + stringArraySize, tmpStringArray, false,b);
        }
        else
        {
            dag[markerIndex].push_back(make_pair(lastA,lastA));
            if(rand()%2==1)
                FlipHaps(start + 0, start + stringArraySize, tmpStringArray, false, lastA);
            else
                FlipHaps(start + 0, start + stringArraySize, tmpStringArray, true, lastA);
        }
        return 0;
    }
}

void MergingEventSimulator::PrintTree() {
    for (int i = 0; i <nHaps ; ++i) {
        int count=0;
        for (int j = 0; j < finalStringArray[i].length() ; ++j) {
            if(finalStringArray[i][j]=='1') count++;
        }
        std::cerr<<i<<"\t"<<finalStringArray[i]<<"\t"<<count<<endl;
    }
}
static bool IsBranchEqual(Branch&a,Branch&b)
{
    if(a[0]==b[0]&&a.back()==b.back()) return true;
    else
return false;
}
void MergingEventSimulator::PutSimulatedEventIntoStringArray() {
    long long int seed = std::chrono::system_clock::now().time_since_epoch().count();
    srand(seed);
    int indexToPut,partitionToPut;
    Branch a,b;
    do
    {
        do
        {
            indexToPut = rand() % (originalString.length() - 80) + 1;
        }while(hasMutation[indexToPut]||indexToPut<=PREFIX+1);
        partitionToPut=(rand() % dag[indexToPut].size());
        a=dag[indexToPut][partitionToPut].first;
        b=dag[indexToPut][partitionToPut].second;
    }while(a.size()<10);//||!IsBranchEqual(a,b));

    int numFlip=rand()%(a.size())+1;

    for (int i = 0; i <numFlip ; ++i) {
        std::shuffle(a.begin(),a.end(),std::default_random_engine(seed));
        int index=a[i];
        //std::cerr<<"index:"<<index<<" indexToPut:"<<indexToPut<<std::endl;
        finalStringArray[index][indexToPut]='1';
    }
    cerr<<"Simulating event at Marker "<<indexToPut<<", Partition "<<partitionToPut<<", first "<<numFlip<<" individuals"<<endl;
    //eventRecorder.push_back(make_pair(indexToPut,a));
}

void MergingEventSimulator::PrintSimulatedEvent() {
    cerr<<"Final:"<<endl;
    for (int i = 0; i <eventRecorder.size() ; ++i) {
        cerr<<eventRecorder[i].first<<"\t";
        for (int j = 0; j <eventRecorder[i].second.size() ; ++j) {
            cerr<<eventRecorder[i].second.at(j);
        }
        cerr<<endl;
    }
}

void MergingEventSimulator::PrintDAG() {
    cerr<<"DAG:"<<endl;
    for (int i = 0; i <dag.size() ; ++i) {
        cerr << "Marker " << i << ":"<<endl;
        for (int l = 0; l < dag[i].size(); ++l) {

            cerr << "Partition " << l << ":"<<endl;
            cerr << "(left:";
            for (int j = 0; j < dag[i][l].first.size(); ++j) {
                cerr << dag[i][l].first.at(j) << "\t";
            }
            cerr << "right:";
            for (int k = 0; k < dag[i][l].second.size(); ++k) {
                cerr << dag[i][l].second.at(k) << "\t";
            }
            cerr<<")";
            cerr << endl;
        }
    }
}
