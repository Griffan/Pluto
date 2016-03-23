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
#include "../pbwtWrapper/ks.h"
int main(int argc, char ** argv) {

//    /*ks.test test*/
//    vector<int> a={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,18,18,20,21,21,23,24,25,26,27,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,84,85,85,85,88,89,90,91,92,93,94,95,96,97,97,97,97};
//    vector<int> b={0,0,0,0,0,0,0,0,9,10,11,11,11,11,15,16,17,17,17,17,17,22,22,22,22,22,22,22,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,44,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,73,73,73,73,73,73,73,73,73,73,73,73,73,73,73,73,73,73,73,73,73,73,73,73,73,73,73};
//
//    std::cerr<<"P value:"<<ks_test(a,b)<<std::endl;//expected p=1.071e-06
//
//    vector<int> c={4,19,20,25,57,58,59,60,83,84,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,120,122,125,126,127,128,129,130,131,132,133,134,135,165,166,167,168,169,170,171,172,173,174,175,176,178,179,180,181,194};
//    vector<int> d={26,82,85,86,182,183,184,185,186,187,188,189,190,191,192,193,195,196,197,198,199};
//
//    std::cerr<<"P value:"<<ks_test(c,d)<<std::endl;//expected p=3.054e-10
//    exit(0);
    /*MergingEventSimulator test*/
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
    Graph->CursorForwards();

//TODO:deal with missing data
    return 0;
}
#endif

#ifdef DEBUG4

#include <fstream>
#include "../SinglePhasing/libVcf/libVcfVcfFile.h"
#include "DebugWrapper.h"
using namespace libVcf;
typedef std::unordered_map<string,string>  ID2POP;
typedef std::unordered_map<int,string>  Index2ID;

char** readVCF(string & fileName, int & nSamples, int& nMarkers, Index2ID& mapper)
{

    printf("starting LoadHaplotypesFromVCF\n");
    bool warningsPrinted = false;
    try {
        VcfFile* pVcf = new VcfFile;
        pVcf->bSiteOnly = true;
        pVcf->bParseGenotypes = false;
        pVcf->bParseDosages = false;
        pVcf->bParseValues = false;

        VcfMarker* pMarker = new VcfMarker;
        //CalculatePhred2Prob();

        pVcf->openForRead(fileName.c_str());

        nSamples = pVcf->getSampleCount();
        if (nSamples == 0) {
            throw VcfFileException("No individual genotype information exist in the input VCF file %s", fileName.c_str());
        }

        for(int j=0; pVcf->iterateMarker(); ++j) {
            //fprintf(stderr,"j=%d\n",j);
            pMarker = pVcf->getLastMarker();
        }
        nMarkers=pVcf->nNumMarkers;
        delete pVcf;
        //delete pMarker;
    }
    catch ( VcfFileException e ) {
        fprintf(stderr, e.what() );
    }

    char** haplotypes=0;
    try {
        VcfFile* pVcf = new VcfFile;
        pVcf->bSiteOnly = false;
        pVcf->bParseGenotypes = true;
        pVcf->bParseDosages = false;
        pVcf->bParseValues = false;

        VcfMarker* pMarker = new VcfMarker;
        //CalculatePhred2Prob();

        pVcf->openForRead(fileName.c_str());

//        nSamples = pVcf->getSampleCount();
//        if (nSamples == 0) {
//            throw VcfFileException("No individual genotype information exist in the input VCF file %s", fileName.c_str());
//        }
        haplotypes=new char* [2*nSamples];
        for (int k = 0; k <2*nSamples ; ++k) {
            haplotypes[k]= new char [nMarkers];
        }
        for(int j=0; pVcf->iterateMarker(); ++j) {
            //fprintf(stderr,"j=%d\n",j);
            pMarker = pVcf->getLastMarker();
            for(int i=0; i < nSamples; ++i) {
                //fprintf(stderr,"i=%d\n",j);
                mapper[i]=string(pVcf->getSampleID(i).c_str());
                //fprintf(stderr,"now reading %dth individual %s\n", i,mapper[i].c_str());
                unsigned short g = pMarker->vnSampleGenotypes[i];
                char g1, g2;

                // genotype is missing
                if ( g == 0xffff ) {
                    //fprintf(stderr,"ERROR: Observed Missing genotypes");
                    //abort();
                }
                else {
                    // genotype is unphased
                    if ( (g & 0x8000) == 0 ) {
                        if ( !warningsPrinted ) {
                            //fprintf(stderr,"ERROR: Observed unphased genotypes %x",g);
                            //abort();
                        }
                    }
                    g1 = (((g & 0x7f00) >> 8) & 0xff);
                    g2 = (g & 0x7f);
                }

                if ( pMarker->asAlts.Length() > 1 ) {
                    if ( g1 == 0 || g2 == 0 ) {
                        fprintf(stderr,"ERROR: TriAllelic Site, but '0' genotype is observed");
                        abort();
                    }
                    --g1;
                    --g2;
                }
                haplotypes[ i * 2 ][ j ] = g1;
                haplotypes[ i * 2 + 1][ j ] = g2;
            }

        }
        delete pVcf;
        //delete pMarker;
    }
    catch ( VcfFileException e ) {
        fprintf(stderr, e.what() );
    }
    printf("finishing LoadHaplotypesFromVCF\n");
    return haplotypes;

}



int readPed(string & fileName,ID2POP& mapper)
{
    ifstream fin(fileName);
    if(!fin.is_open())
    {
        std::cerr<<"Cannot open file "<<fileName<<std::endl;
    }
    string line;
    string ID,POP;
    getline(fin,line);
    while(getline(fin,line))
    {
        stringstream ss(line);
        ss>>ID>>ID;
        ss>>POP>>POP>>POP>>POP>>POP;
        mapper[ID]=POP;
        //std::cerr<<"reading individuals "<<ID<<" from population:"<<POP<<endl;

    }
    return 0;
}
int main(int argc, char ** argv) {


    cerr << "Hello, World!" << endl;
    char ** haps= nullptr;
    //string inputVcf="/Users/fanzhang/Downloads/PlutoTest/OMNI.merged.chr20.phased_genotypes.20141111.vcf.gz";
    string inputVcf="/Users/fanzhang/Downloads/PlutoTest/test.OMNI.vcf.gz";
    string inputPed="/Users/fanzhang/Downloads/PlutoTest/integrated_call_samples.20130502.ALL.ped";
    int nSamples(0),nMarkers(0);
    Index2ID MAP1;
    haps=readVCF(inputVcf,nSamples,nMarkers, MAP1);
    ID2POP MAP2;
    readPed(inputPed,MAP2);
    DebugWrapper  * Graph=new DebugWrapper(2*nSamples,nMarkers);
    fprintf(stderr,"finished initializing graph\n");
    Graph->SetHaps(haps);
    Graph->CursorBackwards();
    fprintf(stderr,"finished backward procedure\n");
    Graph->CursorForwards(MAP1,MAP2);


//TODO:deal with missing data
    return 0;
}
#endif

#ifdef DEBUG5

#include <fstream>
#include "../SinglePhasing/libVcf/libVcfVcfFile.h"
#include "DebugWrapper.h"
using namespace libVcf;
typedef std::unordered_map<string,string>  ID2POP;
typedef std::unordered_map<int,string>  Index2ID;

char** readVCF(string & fileName, int & nSamples, int& nMarkers)
{

    printf("starting LoadHaplotypesFromVCF\n");
    bool warningsPrinted = false;
    try {
        VcfFile* pVcf = new VcfFile;
        pVcf->bSiteOnly = true;
        pVcf->bParseGenotypes = false;
        pVcf->bParseDosages = false;
        pVcf->bParseValues = false;

        VcfMarker* pMarker = new VcfMarker;
        //CalculatePhred2Prob();

        pVcf->openForRead(fileName.c_str());

        nSamples = pVcf->getSampleCount();
        if (nSamples == 0) {
            throw VcfFileException("No individual genotype information exist in the input VCF file %s", fileName.c_str());
        }

        for(int j=0; pVcf->iterateMarker(); ++j) {
            //fprintf(stderr,"j=%d\n",j);
            pMarker = pVcf->getLastMarker();
        }
        nMarkers=pVcf->nNumMarkers;
        delete pVcf;
        //delete pMarker;
    }
    catch ( VcfFileException e ) {
        fprintf(stderr, e.what() );
    }

    char** haplotypes=0;
    try {
        VcfFile* pVcf = new VcfFile;
        pVcf->bSiteOnly = false;
        pVcf->bParseGenotypes = true;
        pVcf->bParseDosages = false;
        pVcf->bParseValues = false;

        VcfMarker* pMarker = new VcfMarker;
        //CalculatePhred2Prob();

        pVcf->openForRead(fileName.c_str());

//        nSamples = pVcf->getSampleCount();
//        if (nSamples == 0) {
//            throw VcfFileException("No individual genotype information exist in the input VCF file %s", fileName.c_str());
//        }
        haplotypes=new char* [2*nSamples];
        for (int k = 0; k <2*nSamples ; ++k) {
            haplotypes[k]= new char [nMarkers];
        }

        srand(1234567);

        for(int j=0; pVcf->iterateMarker(); ++j) {
            //fprintf(stderr,"j=%d\n",j);
            pMarker = pVcf->getLastMarker();
            for(int i=0; i < nSamples; ++i) {
                //fprintf(stderr,"i=%d\n",j);
                //mapper[i]=string(pVcf->getSampleID(i).c_str());
                //fprintf(stderr,"now reading %dth individual %s\n", i,mapper[i].c_str());
                unsigned short g = pMarker->vnSampleGenotypes[i];
                char g1, g2;

                // genotype is missing
                if ( g == 0xffff ) {
                    //fprintf(stderr,"ERROR: Observed Missing genotypes");
                    //abort();
                }
                else {
                    // genotype is unphased
                    if ( (g & 0x8000) == 0 ) {
                        if ( !warningsPrinted ) {
                            //fprintf(stderr,"ERROR: Observed unphased genotypes %x",g);
                            //abort();
                        }
                    }
                    g1 = (((g & 0x7f00) >> 8) & 0xff);
                    g2 = (g & 0x7f);
                }

                if ( pMarker->asAlts.Length() > 1 ) {
                    if ( g1 == 0 || g2 == 0 ) {
                        fprintf(stderr,"ERROR: TriAllelic Site, but '0' genotype is observed");
                        abort();
                    }
                    --g1;
                    --g2;
                }
//                if(rand()%2==1)
//                {
//                    swap(g1,g2);
//                }
                haplotypes[ i * 2 ][ j ] = g1;
                haplotypes[ i * 2 + 1][ j ] = g2;
            }

        }
        delete pVcf;
        //delete pMarker;
    }
    catch ( VcfFileException e ) {
        fprintf(stderr, e.what() );
    }
    printf("finishing LoadHaplotypesFromVCF\n");
    return haplotypes;

}



int readPed(string & fileName,ID2POP& mapper)
{
    ifstream fin(fileName);
    if(!fin.is_open())
    {
        std::cerr<<"Cannot open file "<<fileName<<std::endl;
    }
    string line;
    string ID,POP;
    getline(fin,line);
    while(getline(fin,line))
    {
        stringstream ss(line);
        ss>>ID>>ID;
        ss>>POP>>POP>>POP>>POP>>POP;
        mapper[ID]=POP;
        //std::cerr<<"reading individuals "<<ID<<" from population:"<<POP<<endl;

    }
    return 0;
}

#include <ctime>
int main(int argc, char ** argv) {


    cerr << "Hello, World!" << endl;
    char ** haps= nullptr;
    //string inputVcf="/Users/fanzhang/Downloads/PlutoTest/OMNI.merged.chr20.phased_genotypes.20141111.vcf.gz";
    string inputVcf="/Users/fanzhang/Downloads/PlutoTest/OMNI_1kg_unrel.recode.vcf";
    string inputPed="/Users/fanzhang/Downloads/PlutoTest/integrated_call_samples.20130502.ALL.ped";
    int nSamples(0),nMarkers(0);

    haps=readVCF(inputVcf,nSamples,nMarkers);


    DebugWrapper  * Graph=new DebugWrapper(2*nSamples,nMarkers);
    Graph->SetHaps(haps);
    std::clock_t    start=std::clock();
    Graph->Process(nMarkers,nSamples,haps);
    std::cout << "Time: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
    delete Graph;

//TODO:deal with missing data
    return 0;
}
#endif