//
// Created by Fan Zhang on 10/1/15.
//

#include <iostream>
#include "PBWTWrapper.h"
#include "algorithm"
#include "random"
//#define DEBUG 1
#define DEBUG3
#ifdef DEBUG2
using namespace std;
char**  HapsInit(int M, int N)
{
    char ** haps = new char* [M];
    for(int i=0;i!=M;++i)
        haps[i]=new char [N];
    haps[0]                   ="0000";
    for (int i = 1; i <21 ; ++i) {
        strncpy(haps[i],haps[0],N);
    }
    haps[21]                  ="0001";
    for (int i = 22; i <79+21 ; ++i) {
        strncpy(haps[i],haps[21],N);
    }
    haps[79+21]               ="0011";
    for (int i = 79+21+1; i <95+79+21 ; ++i) {
        strncpy(haps[i],haps[79+21],N);
    }
    haps[95+79+21]            ="0110";
    for (int i = 95+79+21+1; i <116+95+79+21 ; ++i) {
        strncpy(haps[i],haps[95+79+21],N);
    }
    haps[116+95+79+21]        ="1000";
    for (int i = 116+95+79+21+1; i <25+116+95+79+21 ; ++i) {
        strncpy(haps[i],haps[116+95+79+21],N);
    }
    haps[25+116+95+79+21]     ="1001";
    for (int i = 25+116+95+79+21+1; i <112+25+116+95+79+21 ; ++i) {
        strncpy(haps[i],haps[25+116+95+79+21],N);
    }
    haps[112+25+116+95+79+21] ="1011";
    for (int i = 112+25+116+95+79+21+1; i <152+112+25+116+95+79+21 ; ++i) {
        strncpy(haps[i],haps[112+25+116+95+79+21],N);
    }
    return haps;
}
int RandomShuffle(char** &haps,int M)
{
    char* tmp;
    for (int i = 0; i <M*1000 ; ++i) {
        int index=rand() % (M-1);
        int index2=rand() % (M-1);
        tmp=haps[index];
        haps[index]=haps[index2];
        haps[index2]=tmp;
    }

    return 0;
}

char** AddRandomSuffix(char** &haps,int M, unsigned seed)
{
    srand(seed);
    for (int i = 0; i <M ; ++i) {
        int suffix=rand() % 8;
        char* tmp=new char[1+strlen(haps[i])+strlen(string(to_string(suffix)).c_str())];
        strcpy(tmp,haps[i]);
        strcat(tmp,string(to_string(suffix)).c_str());
        haps[i]=tmp;
    }

    return haps;
}
int main(int argc, char ** argv) {


    cerr << "Hello, World!" << endl;
    int M=21+79+95+116+25+112+152;
    int N=8;
    char ** haps=HapsInit(M,N);

    //RandomShuffle(haps,M);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    haps=AddRandomSuffix(haps,M,seed);
    std::shuffle(haps,haps+M,std::default_random_engine(seed));
    std::shuffle(haps,haps+M,std::default_random_engine(seed));
    PBWTWrapper * Graph=new PBWTWrapper(M,N);
    Graph->SetHaps(haps);
    Graph->CursorBackwards();
    Graph->CursorForwards();



        //Graph->PrintHap(haps, Graph->a[0]);
        // Wrapper->PrintHap(tmpHaps,Wrapper->a[6]);
        //Graph->PrintHap(haps, Graph->a[Graph->N - 1]);
        // Wrapper->PrintMatrix(Wrapper->a,"a array matrix");
    //Graph->PrintMatrix(Graph->d, "d array");
        //Wrapper->PrintVector(Wrapper->a[Wrapper->N-7],"last a array");
//TODO:deal with missing data
    return 0;
}
#endif

#ifdef DEBUG3
#include "MergingEventSimulator.h"
int main(int argc, char ** argv) {

    MergingEventSimulator MS;
    cerr << "Hello, World!" << endl;
    char ** haps= nullptr;
    MS.SimulateCoalescentTree();
 //   MS.PrintTree();
    RESULT* result=MS.PutSimulatedEventIntoStringArray();
    //MS.PrintSimulatedEvent();
    MS.PrintTree();
 //   MS.PrintDAG();
    haps=MS.ChangeStringArraytoCstringArray();
    //RandomShuffle(haps,M);
//    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//    std::shuffle(haps,haps+M,std::default_random_engine(seed));
    PBWTWrapper * Graph=new PBWTWrapper(MS.GetNumHaps(),MS.GetHapLength());
    Graph->SetHaps(haps);
    Graph->CursorBackwards();
    Graph->CursorForwards(result);

//TODO:deal with missing data
    return 0;
}
#endif