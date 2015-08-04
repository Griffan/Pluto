#include <iostream>
#include "pbwtWrapper/PBWTWrapper.h"
using namespace std;

int main() {
    cerr << "Hello, World!" << endl;
    int M=10,N=7;
    char ** haps;
    haps = new char* [M];
    //for(int i=0;i!=M;++i)
     //   haps[i]=new char [N];
    haps[0]="0000000";
    haps[1]="0001000";
    haps[2]="0010000";
    haps[3]="0011000";
    haps[4]="0100000";
    haps[5]="0101000";
    haps[6]="0110000";
    haps[7]="0111000";
    haps[8]="1000000";
    haps[9]="1001000";
    PBWTWrapper Graph(M,N);
    Graph.setHaps(haps);
    Graph.CursorForwards();
    Graph.CursorBackwards();

    return 0;
}