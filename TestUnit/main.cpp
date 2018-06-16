//
// Created by Fan Zhang on 10/1/15.
//

#include <iostream>
#include "PBWTWrapper.h"
#include "random"
#include "gtest/gtest.h"

#define GTEST
#ifdef GTEST
using namespace std;
char**  HapsInit(int M, int N)
{
    char ** haps = new char* [M];
    for(int i=0;i!=M;++i)
        haps[i]=new char [N];
//    haps[0]                   = const_cast<char *>("0000");
    haps[0][0]=0;
    haps[0][1]=0;
    haps[0][2]=0;
    haps[0][3]=0;
    for (int i = 1; i <21 ; ++i) {
        memcpy(haps[i],haps[0],N);
    }
//    haps[21]                  =const_cast<char *>("0001");
    haps[21][0]=0;
    haps[21][1]=0;
    haps[21][2]=0;
    haps[21][3]=1;
    for (int i = 21+1; i <79+21 ; ++i) {
        memcpy(haps[i],haps[21],N);
    }
//    haps[79+21]               =const_cast<char *>("0010");
    haps[79+21][0]=0;
    haps[79+21][1]=0;
    haps[79+21][2]=1;
    haps[79+21][3]=0;
    for (int i = 79+21+1; i <95+79+21 ; ++i) {
        memcpy(haps[i],haps[79+21],N);
    }
//    haps[95+79+21]            =const_cast<char *>("0011");
    haps[95+79+21][0]=0;
    haps[95+79+21][1]=0;
    haps[95+79+21][2]=1;
    haps[95+79+21][3]=1;
    for (int i = 95+79+21+1; i <116+95+79+21 ; ++i) {
        memcpy(haps[i],haps[95+79+21],N);
    }
//    haps[116+95+79+21]        =const_cast<char *>("0100");
    haps[116+95+79+21][0]=0;
    haps[116+95+79+21][1]=1;
    haps[116+95+79+21][2]=0;
    haps[116+95+79+21][3]=0;
    for (int i = 116+95+79+21+1; i <25+116+95+79+21 ; ++i) {
        memcpy(haps[i],haps[116+95+79+21],N);
    }
//    haps[25+116+95+79+21]     =const_cast<char *>("0101");
    haps[25+116+95+79+21][0]=0;
    haps[25+116+95+79+21][1]=1;
    haps[25+116+95+79+21][2]=0;
    haps[25+116+95+79+21][3]=1;
    for (int i = 25+116+95+79+21+1; i <112+25+116+95+79+21 ; ++i) {
        memcpy(haps[i],haps[25+116+95+79+21],N);
    }
//    haps[112+25+116+95+79+21] =const_cast<char *>("0111");
    haps[112+25+116+95+79+21][0]=0;
    haps[112+25+116+95+79+21][1]=1;
    haps[112+25+116+95+79+21][2]=1;
    haps[112+25+116+95+79+21][3]=1;
    for (int i = 112+25+116+95+79+21+1; i <152+112+25+116+95+79+21 ; ++i) {
        memcpy(haps[i],haps[112+25+116+95+79+21],N);
    }

    return haps;
}

char**  HapsInit2(int M, int N)
{
    char ** haps = new char* [M];
    for(int i=0;i!=M;++i)
        haps[i]=new char [N];
//    haps[0]                   = const_cast<char *>("0000");
    haps[0][0]=0;
    haps[0][1]=0;
    haps[0][2]=0;
    haps[0][3]=0;
    haps[0][4]=0;

    for (int i = 1; i <21 ; ++i) {
        memcpy(haps[i],haps[0],N);
    }
//    haps[21]                  =const_cast<char *>("0001");
    haps[21][0]=0;
    haps[21][1]=0;
    haps[21][2]=0;
    haps[21][3]=1;
    haps[21][4]=1;

    for (int i = 21+1; i <79+21 ; ++i) {
        memcpy(haps[i],haps[21],N);
    }
//    haps[79+21]               =const_cast<char *>("0010");
    haps[79+21][0]=0;
    haps[79+21][1]=0;
    haps[79+21][2]=1;
    haps[79+21][3]=1;
    haps[79+21][4]=1;

    for (int i = 79+21+1; i <95+79+21 ; ++i) {
        memcpy(haps[i],haps[79+21],N);
    }
//    haps[95+79+21]            =const_cast<char *>("0011");
    haps[95+79+21][0]=0;
    haps[95+79+21][1]=1;
    haps[95+79+21][2]=1;
    haps[95+79+21][3]=0;
    haps[95+79+21][4]=0;

    for (int i = 95+79+21+1; i <116+95+79+21 ; ++i) {
        memcpy(haps[i],haps[95+79+21],N);
    }
//    haps[116+95+79+21]        =const_cast<char *>("0100");
    haps[116+95+79+21][0]=1;
    haps[116+95+79+21][1]=0;
    haps[116+95+79+21][2]=0;
    haps[116+95+79+21][3]=0;
    haps[116+95+79+21][4]=0;

    for (int i = 116+95+79+21+1; i <25+116+95+79+21 ; ++i) {
        memcpy(haps[i],haps[116+95+79+21],N);
    }
//    haps[25+116+95+79+21]     =const_cast<char *>("0101");
    haps[25+116+95+79+21][0]=1;
    haps[25+116+95+79+21][1]=0;
    haps[25+116+95+79+21][2]=0;
    haps[25+116+95+79+21][3]=1;
    haps[25+116+95+79+21][4]=1;

    for (int i = 25+116+95+79+21+1; i <112+25+116+95+79+21 ; ++i) {
        memcpy(haps[i],haps[25+116+95+79+21],N);
    }
//    haps[112+25+116+95+79+21] =const_cast<char *>("0111");
    haps[112+25+116+95+79+21][0]=1;
    haps[112+25+116+95+79+21][1]=0;
    haps[112+25+116+95+79+21][2]=1;
    haps[112+25+116+95+79+21][3]=1;
    haps[112+25+116+95+79+21][4]=1;

    for (int i = 112+25+116+95+79+21+1; i <152+112+25+116+95+79+21 ; ++i) {
        memcpy(haps[i],haps[112+25+116+95+79+21],N);
    }

    std::random_shuffle(haps, haps + 600);

    return haps;
}
class PBWTWrapperTest : public ::testing::Test {
protected:
    PBWTWrapper* Wrapper;
    int M,N;
    PBWTWrapperTest()
    {
        M=21+79+95+116+25+112+152;
        N=4;
        Wrapper = new PBWTWrapper(M,N,nullptr,0);
    }

    virtual ~PBWTWrapperTest(){}

    virtual void SetUp(){
        cerr << "Hello, World! From PBWTWrapper" << endl;
        char ** haps=HapsInit(M,N);
        Wrapper->SetHaps(haps,0,2*M,0,0,0,0);
    }

    virtual void TearDown(){}
};

/*
   PBWTWrapper(){}
    ~PBWTWrapper(){}

    PBWTWrapper(const char ** haps, int nhaps, int nsnps);
    PBWTWrapper(int nhaps,int nsnps);

    int CursorForwards();
    int CursorBackwards();

    int CursorForwardsTo(int k, int T=5);

    int CursorBackwardsTo(int k, int T=5);

    int CopyHap(int k, PbwtCursor* Cursor);

    int LabelNoSiblingCluster(int site);
	int UpdateTransVector(int site);

    bool IsEditDistanceOK(const std::vector< std::vector<char*> >& backBone,int stateA, int stateB, int index, double thresh);
    int MergeAtSite(int site);

    int MoveSegment(const std::unordered_map<int,int>& mergedMembership,int site);

    int SetHaps(char **haps);

    //inline functions
    inline int GetNumStates(int k)

    inline unsigned long GetNumHaps(int site) const { return haplotypeCluster[site].size(); }

    inline int GetHapIDFromBack(int site, int backRank) const { return alpha[site + 1][backRank]; }
    inline int GetHapIDFromFwd(int site, int fwdRank) const { if(site<0) return fwdRank; else return a[site][fwdRank]; }//you should only use it after a being updated

    inline int GetRankFromBack(int site, int hapID) {return alphaMap[site][hapID];}
    inline int GetRankFromFwd(int site, int hapID) {return aMap[site][hapID];}

    inline int GetHapStateFromFwd(int site, int hapID) { return haplotypeCluster[site][hapID]; }
    inline void ResetWrapper()

 * */
//TEST_F(PBWTWrapperTest, IsEditDistanceOK) {
//    std::vector< std::vector<char*> > backBone(2,std::vector<char*>(0,0));
//    backBone[0].push_back(Wrapper->haplotype[0]);
//    backBone[1].push_back(Wrapper->haplotype[21]);
//    bool result=Wrapper->IsEditDistanceOK(backBone,0,1,0,1);
//    EXPECT_EQ(true,result);
//    result=Wrapper->IsEditDistanceOK(backBone,0,1,0,0);
//    EXPECT_EQ(true,result);
//    result=Wrapper->IsEditDistanceOK(backBone,0,1,0,4);
//    EXPECT_EQ(false,result);
//}

TEST_F(PBWTWrapperTest, IsRecipricalLengthOK)
{
    std::vector<int> a={1,2,3,4,6,7};
    std::vector<int> b={5,8,9};
    bool result=Wrapper->IsRecipricalLengthOK(a,b);
    EXPECT_EQ(false,result);

    std::vector<int> c={1,2,3,4,6,7};
    std::vector<int> d={5};
    result = Wrapper->IsRecipricalLengthOK(c,d);
    EXPECT_EQ(true,result);

    std::vector<int> e={1,2,3,4,6,7};
    std::vector<int> f={9};
    result = Wrapper->IsRecipricalLengthOK(e,f);
    EXPECT_EQ(false,result);

    std::vector<int> g={1,2,3,4,6,7};
    std::vector<int> h={-1,9};
    result = Wrapper->IsRecipricalLengthOK(g,h);
    EXPECT_EQ(true,result);

}

TEST_F(PBWTWrapperTest, MergeSortedArrayToA)
{
    std::vector<int> a={1,2,3,4,6,7};
    std::vector<int> b={5,8,9};
    std::vector<int> c={1,2,3,4,5,6,7,8,9};
    Wrapper->MergeSortedArrayToA(a,b);
    EXPECT_EQ(c,a);

    std::vector<int> a1={1};
    std::vector<int> b1={2};
    std::vector<int> c1={1,2};
    Wrapper->MergeSortedArrayToA(a1,b1);
    EXPECT_EQ(c1,a1);
}

TEST_F(PBWTWrapperTest, CursorBackwards)
{
    std::vector<int> a(M,0);
    for (int k1 = 0; k1 <M; ++k1) {
        a[k1]=k1;
    }
    std::vector<int> d(M,0);
    d[0]=4;//0000
    for (int i = 1; i <21 ; ++i) {
        d[i]=0;
    }
    d[21]=1;//0001
    for (int j = 21+1; j <21+79 ; ++j) {
        d[j]=0;
    }
    d[21+79]=2;//0010
    for (int k = 21+79+1; k <21+79+95; ++k) {
        d[k]=0;
    }
    d[21+79+95]=1;//0011
    for (int l = 21+79+95+1; l <21+79+95+116 ; ++l) {
        d[l]=0;
    }
    d[21+79+95+116]=3;//0100
    for (int m =21+79+95+116+1; m <21+79+95+116+25 ; ++m) {
        d[m]=0;
    }
    d[21+79+95+116+25]=1;//0101
    for (int n = 21+79+95+116+25+1; n <21+79+95+116+25+112 ; ++n) {
        d[n]=0;
    }
    d[21+79+95+116+25+112]=2;//0111
    for (int i1 = 21+79+95+116+25+112+1; i1 <21+79+95+116+25+112+152 ; ++i1) {
        d[i1]=0;
    }

    Wrapper->CursorBackwards();

    EXPECT_EQ(a,Wrapper->alpha);

    EXPECT_EQ(d,Wrapper->delta);
}

TEST_F(PBWTWrapperTest, CursorForwards)
{

    std::vector<int> a(M,0);
    for (int i = 0; i <21 ; ++i) {
        a[i]=i;
    }
    for (int j = 21; j <21+25 ; ++j) {
        a[j]=116+95+79+21+j-21;
    }
    for (int k = 21+25; k <21+25+95; ++k) {
        a[k]=79+21+k-(21+25);
    }
    for (int l = 21+25+95; l <21+25+95+79 ; ++l) {
        a[l]=21 + l-(21+25+95);
    }
    for (int m = 21+25+95+79; m <21+25+95+79+112 ; ++m) {
        a[m]=25+116+95+79+21+m-(21+25+95+79);
    }
    for (int n = 21+25+95+79+112; n <21+25+95+79+112+116 ; ++n) {
        a[n]=95+79+21+n-(21+25+95+79+112);
    }
    for (int i1 = 21+25+95+79+112+116; i1 <21+25+95+79+112+116+152 ; ++i1) {
        a[i1]=112+25+116+95+79+21+i1-(21+25+95+79+112+116);
    }

    std::vector<int> d(M,0);
    d[0]=4;//0000
    for (int i = 1; i <21 ; ++i) {
        d[i]=0;
    }
    d[21]=2;//0100
    for (int j = 21+1; j <21+25 ; ++j) {//46
        d[j]=0;
    }
    d[21+25]=3;//0010
    for (int k = 21+25+1; k <21+25+95; ++k) {//141
        d[k]=0;
    }
    d[21+25+95]=4;//0001
    for (int l = 21+25+95+1; l <21+25+95+79 ; ++l) {//220
        d[l]=0;
    }
    d[21+25+95+79]=2;//0101
    for (int m = 21+25+95+79+1; m <21+25+95+79+112 ; ++m) {//332
        d[m]=0;
    }
    d[21+25+95+79+112]=3;//0011
    for (int n = 21+25+95+79+112+1; n <21+25+95+79+112+116 ; ++n) {//448
        d[n]=0;
    }
    d[21+25+95+79+112+116]=2;//0111
    for (int i1 = 21+25+95+79+112+116+1; i1 <21+25+95+79+112+116+152 ; ++i1) {//600
        d[i1]=0;
    }
    Wrapper->CursorBackwards();
    Wrapper->CursorForwards();
    EXPECT_EQ(a,Wrapper->a);
//    Wrapper->PrintVector(d,"expected");
//    Wrapper->PrintVector(Wrapper->d[Wrapper->N-1],"actual");
    EXPECT_EQ(d,Wrapper->d);
}

TEST_F(PBWTWrapperTest, MoveSegment)
{

    Wrapper->CursorBackwards();

}

class StateNodeTest : public ::testing::Test {
protected:
    PBWTWrapper* Wrapper;
    int M,N;
    StateNodeTest()
    {
        M=21+79+95+116+25+112+152;
        N=4;
        Wrapper = new PBWTWrapper(M,N,0,20);
    }

    virtual ~StateNodeTest(){}

    virtual void SetUp(){
        cerr << "Hello, World! From StateNodeTest" << endl;
        char ** haps=HapsInit(M,N);
        Wrapper->SetHaps(haps,0,2*M,0,0,0,0);
    }

    virtual void TearDown(){}
};

class EdgeMergeTest : public ::testing::Test {
protected:
    PBWTWrapper* Wrapper;
    int M,N;
    EdgeMergeTest()
    {
        M=600;
        N=5;
        Wrapper = new PBWTWrapper(M,N,0,20);
    }

    virtual ~EdgeMergeTest(){}

    virtual void SetUp(){
        cerr << "Hello, World! From StateNodeTest" << endl;
        char ** haps=HapsInit2(M,N);
        Wrapper->SetHaps(haps,0,2*M,0,0,0,0);
    }

    virtual void TearDown(){}
};

TEST_F(EdgeMergeTest, RegressionMergeAtSite)
{
    Wrapper->CursorBackwards();
//    Wrapper->CursorForwardsTo(0,20);
//    Wrapper->CursorForwardsTo(1,20);
//    Wrapper->CursorForwardsTo(2,20);
//    Wrapper->RegressionMergeAtSite(1,true);
    Wrapper->CursorForwards();
}
int main(int argc, char ** argv) {

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
#endif

#ifdef RUN_SIMPLE_SIMULATOR
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
    Graph->SetHaps(haps,0,0,0,0);
    Graph->CursorBackwards();
    Graph->CursorForwards();

//TODO:deal with missing data
    return 0;
}
#endif

#ifdef DRAW_PLOT

#include <fstream>
#include "../SinglePhasing/libVcf/libVcfVcfFile.h"
#include "PBWTViewer.h"
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

/*!
 *
 * @param fileName
 * @param ptr
 * @return
 */
char** readMACS(const string & fileName, int& nSamples, int& nMarkers)
{
    ifstream fin(fileName);
    if(!fin.is_open())
    {
        std::cerr<<"Cannot open file "<<fileName<<std::endl;
    }
    string line,tmp;
    float coord(0.);
    int index(0);
    getline(fin,line);//read cmdline
    stringstream ss(line);
    ss>>tmp>>nSamples;
    getline(fin,line);//read random seed
    getline(fin,line);//read empty line

    getline(fin,line);//read "//"

    while(line.find("segsites",0)==std::string::npos)
        getline(fin,line);//read tree
    //now line contains segsites:totalpos
    ss.clear();
    ss.str(line);

    ss>>tmp>>nMarkers;
    ss.clear();
    float * positions = new float[nMarkers];
    getline(fin,line);//read positions
    ss.str(line);
    ss>>tmp;
    while(ss>>coord) positions[index++] = coord;

    char** hap = new char* [nSamples];
    index = 0;
    while(getline(fin,line))
    {
        hap[index] = new char [nMarkers];
        for (int i = 0; i < nMarkers; ++i) {
            hap[index][i] = line[i]-'0';
        }
        index++;
    }
    nSamples/=2;
    return hap;
}

//#define READVCF 0
int main(int argc, char ** argv) {


    cerr << "Hello, World!" << endl;
    char ** haps= nullptr;
#ifdef READVCF
    string inputVcf="/Users/fanzhang/Downloads/PlutoTest/omni/OMNI.merged.chr20.phased_genotypes.20141111.vcf.gz";
//    string inputVcf="/Users/fanzhang/Downloads/PlutoTest/test.OMNI.vcf.gz";
    string inputPed="/Users/fanzhang/Downloads/PlutoTest/integrated_call_samples.20130502.ALL.ped";
#else
    string inputMACShap="/Users/fanzhang/Downloads/macs/haplotypes.txt";
#endif

    string outputMatrix="/Users/fanzhang/Downloads/PlutoTest/omni/hyun.rank.d.matrix.with.allele.debug.rear";
    string outputIndividual="/Users/fanzhang/Downloads/PlutoTest/omni/hyun.rank.d.individual";
    string outputMarker="/Users/fanzhang/Downloads/PlutoTest/omni/hyun.rank.d.marker";

    int nSamples(0),nMarkers(0);
#ifdef READVCF
    Index2ID MAP1;
    haps=readVCF(inputVcf, nSamples, nMarkers, MAP1);
    ID2POP MAP2;
    readPed(inputPed,MAP2);
#else
    haps=readMACS(inputMACShap, nSamples, nMarkers);
#endif

    PBWTViewer  * Wrapper=new PBWTViewer(2*nSamples,nMarkers);
    fprintf(stderr,"finished initializing graph\n");
    Wrapper->haplotype=haps;
    Wrapper->CursorBackwards();
    fprintf(stderr,"finished backward procedure\n");
    Wrapper->CursorForwards(true);

    ofstream fout(outputMatrix);
    if(!fout.is_open())
    {
        fprintf(stderr,"open file %s failed!",outputMatrix.c_str());
        exit(EXIT_FAILURE);
    }
//    for (int hapID = 0; hapID <Wrapper->M; ++hapID) {
//        for (int marker = 1; marker <Wrapper->N-1; ++marker) {
//            int hapRank=Wrapper->GetRankFromFwd(marker-1,hapID);
//            int backHapRank=Wrapper->GetRankFromBack(marker+1,hapID);
//            int left=(marker-1-Wrapper->d[marker-1][hapRank]);
//            int right =(Wrapper->N-Wrapper->delta[marker+1][backHapRank]-1)-(marker+1);
//            fout<<hapRank<<","<<backHapRank<<","<<(left>0?left:0)<<","<<(right>0?right:0)<<":"<<(int)Wrapper->haplotype[hapID][marker]<<"\t";
//        }
//        fout<<"\n";
//    }
    for (int hapID = 0; hapID <Wrapper->M; ++hapID) {
        for (int marker = Wrapper->N-10; marker <Wrapper->N; ++marker) {
            int hapRank=Wrapper->GetRankFromFwd(marker,hapID);
            int backHapRank=Wrapper->GetRankFromBack(marker,hapID);
            int left=(marker-Wrapper->d[marker][hapRank]);
            int right =(Wrapper->N-Wrapper->delta[marker][backHapRank]-1)-(marker);
            fout<<hapRank<<","<<backHapRank<<","<<Wrapper->d[marker][hapRank]<<","<<Wrapper->delta[marker][backHapRank]<<":"<<(int)Wrapper->haplotype[hapID][marker]<<"\t";
        }
        fout<<"\n";
    }
    fout.close();

#ifdef  READVCF
    fout.open(outputIndividual);
    for (int i = 0; i <Wrapper->M; ++i) {
        fout<<MAP1[i/2]<<endl;
    }
    fout.close();
#endif

    return 0;
}
#endif

#ifdef SIMPLIFIED_PHASING

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