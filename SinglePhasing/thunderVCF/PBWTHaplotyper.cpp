//
// Created by Fan Zhang on 8/6/15.
//

#include <pbwt/pbwt.h>
#include "PBWTHaplotyper.h"
#include "MemoryAllocators.h"
//debug related
static const float UNDERFLOW_MIN = std::numeric_limits<double>::min()*100;
static void printLeftMatrix(float * probability, int numStates)
{
	for (int i = 0; i <numStates; ++i) {
		for (int j = 0; j <= i; ++j, probability++) {
			fprintf(stderr, "(%d,%d):%9.9f\t", i, j, *probability);
		}
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "\n");
}
//initiation
PBWTHaplotyper::PBWTHaplotyper(int nhaps, int nsnps) {

    //Wrapper = new PBWTWrapper(nhaps, nsnps);

}

PBWTHaplotyper::PBWTHaplotyper():genuienParents(10000,ChildToSource()) {//TODO:change the number
    onlyHeterSite=false;
    max_num=1;
    nSampleCopy=0;//additional, the original haps not included
}

void PBWTHaplotyper::InitAuxillary() {

    tmpHaps = AllocateCharMatrix(individuals * 2+(individuals-phased)*nSampleCopy*2, markers);

    tmpGeno = AllocateCharMatrix(individuals, markers*3);

    tmpPenetrance = new float [markers * 9];

    Wrapper=nullptr;
}

PBWTHaplotyper::~PBWTHaplotyper() {

    if(tmpHaps != nullptr) {
        for (int i = 0; i < 2 * individuals+(individuals-phased)*nSampleCopy*2; ++i) {
            delete[] tmpHaps[i];
        }
        delete [] tmpHaps;
    }
    if(tmpGeno != nullptr) {
        for (int i = 0; i < individuals; ++i) {
            delete[] tmpGeno[i];
        }
        delete [] tmpGeno;
    }
    if(tmpPenetrance != nullptr)
        delete [] tmpPenetrance;


    for (int l = 0; l < nSampleCopy*(individuals-phased)*2; ++l) {
        delete [] sampledHaps[l];
    }
    delete [] sampledHaps;

    if(Wrapper != nullptr)
        delete Wrapper;

    if(fwdWrapper != nullptr)
        delete fwdWrapper;
	ReleaseMemoryBlock();
    DestroyPvalueMatrix();
}

//memory management
void PBWTHaplotyper::ReleaseMemoryBlock()
{
	for (std::unordered_map<int, std::vector<float*> >::iterator iter = memoryBlockList.begin(); iter != memoryBlockList.end(); ++iter)
	{
		for (size_t i = 0; i < iter->second.size(); i++)
		{
			if (iter->second[i] != nullptr) delete[] iter->second[i];
		}
	}
}
void PBWTHaplotyper::GetMemoryBlock(int marker)
{
	if (!economyMode || marker == 0 || marker > stack[stackPtr] + gridSize)
	{
		stack[++stackPtr] = marker;
		leftMatrices[marker] = GetLargeBlock();

		ResetReuseablePool();
	}
	else
		leftMatrices[marker] = GetReuseableBlock();
}
float* PBWTHaplotyper::GetLargeBlock()
{
	int blockSize = orderedGenotypes ? states * states : states * (states + 1) / 2;
	if (numInUse.find(blockSize) == numInUse.end())
	{
		numInUse[blockSize] = 0;
		memoryBlockList[blockSize] = std::vector<float*>(0, nullptr);
	}
	if (numInUse[blockSize] < (int)memoryBlockList[blockSize].size())
	{
		numInUse[blockSize]++;
		return memoryBlockList[blockSize][numInUse[blockSize] - 1];
	}
	else
	{
		memoryBlockList[blockSize].push_back(AllocateMemoryBlock());
		numInUse[blockSize]++;
		return memoryBlockList[blockSize][numInUse[blockSize] - 1];
	}
}
float* PBWTHaplotyper::GetReuseableBlock()
{
	int blockSize = orderedGenotypes ? states * states : states * (states + 1) / 2;
	if (numInUse.find(blockSize) == numInUse.end())
	{
		numInUse[blockSize] = 0;
		memoryBlockList[blockSize] = std::vector<float*>(0, nullptr);
	}
	if (numInUse[blockSize] < (int)memoryBlockList[blockSize].size())
	{
		numInUse[blockSize]++;//TODO::reset needed
		return memoryBlockList[blockSize][numInUse[blockSize] - 1];
	}
	else
	{
		memoryBlockList[blockSize].push_back(AllocateMemoryBlock());
		numInUse[blockSize]++;
		return memoryBlockList[blockSize][numInUse[blockSize] - 1];
	}

}
void PBWTHaplotyper::ResetMemoryPool()
{
	nextAvailable = nextSmallAvailable = 0;
	nextReuseable = markers - 1;
	stackPtr = -1;
	for (std::unordered_map<int, int>::iterator iter = numInUse.begin(); iter != numInUse.end();++iter)
	{
		iter->second = 0;
	}
}

void PBWTHaplotyper::ResetReuseablePool()
{
	nextReuseable = markers - 1;
}
void PBWTHaplotyper::RetrieveMemoryBlock(int marker) {
    if (stack[stackPtr] <= marker) {
//        fprintf(stderr, "%d out from RetrieveMemory\n",marker);
        return;
    }
    else {
        ResetReuseablePool();

        float *from = leftMatrices[stack[--stackPtr]];

        for (int i = stack[stackPtr] + 1; i <= marker; i++) {
            int markerindex = i * 3;
            {
				leftMatrices[i] = GetReuseableBlock();

                Transpose(i, from, leftMatrices[i]);
                ConditionOnData(leftMatrices[i], i, genotypes[states / 2][markerindex],
                                genotypes[states / 2][markerindex + 1], genotypes[states / 2][markerindex + 2]);

                from = leftMatrices[i];
            }
        }
    }
}

bool PBWTHaplotyper::ForceMemoryAllocation() {
    // Cycle through individuals, with the exact same steps as the actual
    // haplotyper and request memory ... by requesting all memory upfront,
    // we force crashes to happen early.
    for (int i = 0; i < individuals - phased; i++) {
        ResetMemoryPool();
        GetMemoryBlock(0);

        if (leftMatrices[0] == nullptr)
            return false;

        int skipped = 0;
        for (int j = 1; j < markers; j++)
            //if (genotypes[i][j] != GENOTYPE_MISSING || j == markers - 1)
        {
            GetMemoryBlock(j);

            if (leftMatrices[j] == nullptr)
                return false;
        }
        //else
        //skipped++;

        if (skipped == 0) break;
    }

    if (!phased)
        return true;

    ResetMemoryPool();
    for (int j = 0; j < markers; j++) {
        GetSmallMemoryBlock(j);

        if (leftMatrices[j] == nullptr)
            return false;
    }

    return true;
}
//random setup
void PBWTHaplotyper::InitialSampleCopy(Random * rand)
{

    if (rand == nullptr)
        rand = &globalRandom;
    CalculatePhred2Prob();

    if (nSampleCopy == 0) return;
    sampledHaps = new char* [nSampleCopy*(individuals-phased)*2];
    for (int l = 0; l < nSampleCopy*(individuals-phased)*2; ++l) {
        sampledHaps[l]= new char[markers];
    }

    for (int j = 0; j < markers; j++)
    {
        double mac = 0;
        int markerindex = 3*j;

        double hyperprior11 = freq1s[j] * freq1s[j];
        double hyperprior12 = 2.0 * freq1s[j] * (1.0 - freq1s[j]);
        double hyperprior22 = (1.0 - freq1s[j]) * (1.0 - freq1s[j]);

        for (int i = 0; i < individuals; i++)
        {
            double post11 = hyperprior11 * phred2prob[(size_t)genotypes[i][markerindex]];
            double post12 = hyperprior12 * phred2prob[(size_t)genotypes[i][markerindex+1]];
            double post22 = hyperprior22 * phred2prob[(size_t)genotypes[i][markerindex+2]];
            double sumpost = post11 + post12 + post22;
            post11 /= sumpost;
            post12 /= sumpost;
            post22 /= sumpost;

            // estimated counts of AL2
            mac += post12+ 2*post22;
        }

        //here, each person contributes two alleles
        double freq = 0.5 * mac / (double) individuals;

        double prior_11 = (1.0 - freq) * (1.0 - freq);
        double prior_12 = 2.0 * freq * (1.0 - freq);
        double prior_22 = freq * freq;



        for (int i = 0; i < individuals-phased; i++)
        {
            int observed = (unsigned char) (genotypes[i][j]);

            double posterior_11 = prior_11 * phred2prob[(size_t)genotypes[i][markerindex]];
            double posterior_12 = prior_12 * phred2prob[(size_t)genotypes[i][markerindex+1]];
            double posterior_22 = prior_22 * phred2prob[(size_t)genotypes[i][markerindex+2]];
            double sum = posterior_11 + posterior_12 + posterior_22;

            if (sum == 0)
                printf("Problem!\n");

            posterior_11 /= sum;
            posterior_12 /= sum;

            for (int k = 0; k <nSampleCopy; ++k) {
                double r = rand->Next();
                int index=i*nSampleCopy+k;
                if (r < posterior_11)
                {
                    sampledHaps[index * 2][j] = 0;
                    sampledHaps[index * 2 + 1][j] = 0;
                }
                else if (r < posterior_11 + posterior_12)
                {
                    bool bit = rand->Binary();

                    sampledHaps[index * 2][j] = bit;
                    sampledHaps[index * 2 + 1][j] = bit ^ 1;
                }
                else
                {
                    sampledHaps[index * 2][j] = 1;
                    sampledHaps[index * 2 + 1][j] = 1;
                }
            }
        }
    }
}

//phasing
void PBWTHaplotyper::ConditionOnData(float *matrix, int marker, char phred11, char phred12, char phred22) {
    // We treat missing genotypes as uninformative about the mosaic's
    // underlying state. If we were to allow for deletions and the like,
    // that may no longer be true.

    //if (genotype == GENOTYPE_MISSING)
    //return;
    float* source=matrix;
    double sum=0.;
    double conditional_probs[3];
    uchar ph11 = (unsigned char) phred11;
    uchar ph12 = (unsigned char) phred12;
    uchar ph22 = (unsigned char) phred22;
//    uchar ph11 = 0;
//    uchar ph12 = 0;
//    uchar ph22 = 0;

    CalculatePhred2Prob();
    //if(marker==106) fprintf(stderr,"(ph11,ph12,ph22)=(%d,%d,%d)\n",ph11,ph12,ph22);
    for (int i = 0; i < 3; i++)
        conditional_probs[i] = Penetrance(marker, i, 0) * phred2prob[ph11] +
                               Penetrance(marker, i, 1) * phred2prob[ph12] +
                               Penetrance(marker, i, 2) * phred2prob[ph22];
//    int nan=0;
    for (int i = 0; i < states; i++) {
        double factors[2];

//        factors[0] = conditional_probs[haplotypes[i][marker]];
//        factors[1] = conditional_probs[haplotypes[i][marker] + 1];
        factors[0] = conditional_probs[GetAllele(marker,i)];//Wrapper->clusterAllele[marker][i]];
        factors[1] = conditional_probs[GetAllele(marker,i)+1];//Wrapper->clusterAllele[marker][i] + 1];

        for (int j = 0; j <= i; j++, matrix++) {
            *matrix *= factors[GetAllele(marker,j)];//Wrapper->clusterAllele[marker][j]];
            sum+=*matrix;
//            if(*matrix >0)
//            {
//                fprintf(stderr,"(%d,%d):%g\n",i,j,*matrix);
//                nan=1;
//            }
        }
//        if(*matrix >0)
//        {
//            fprintf(stderr,"(%d,%d):%g\n",i,i,*matrix);
//            nan=1;
//        }
    }

    matrix=source;
    for (int i = 0; i < states; i++) {
        for (int j = 0; j <= i; j++, matrix++) {
            *matrix /=sum;
        }
    }
}

void PBWTHaplotyper::ScoreLeftConditional() {


    ResetMemoryPool();
    UpdateStateNum(GetStateNumFrom(0));
    GetMemoryBlock(0);

    SetupPrior(leftMatrices[0]);
    ConditionOnData(leftMatrices[0], 0, genotypes[individuals - 1][0], genotypes[individuals - 1][1],
                    genotypes[individuals - 1][2]);

    float *from = leftMatrices[0];

    for (int i = 1; i < markers; i++) {
        int markerindex = i * 3;
        //fprintf(stderr, "processing marker %d...\n", i);
        // Cumulative recombination fraction allows us to skip uninformative positions
        // theta = theta + thetas[i - 1] - theta * thetas[i - 1];

        // Skip over uninformative positions to save time
        // maybe check for the difference between min(phred11, phred12, phred22) and the next smallest

        //if (genotypes[states / 2][i] != GENOTYPE_MISSING || i == markers - 1)
        {
            UpdateStateNum(GetStateNumFrom(i));

            GetMemoryBlock(i);

            Transpose(i, from, leftMatrices[i]);
//
//            fprintf(stderr,"marker:%d\tout of (%d) states\n",i,GetStateNumFrom(i));
//            printLeftMatrix(leftMatrices[i],GetStateNumFrom(i));

            ConditionOnData(leftMatrices[i], i, genotypes[individuals - 1][markerindex],genotypes[individuals - 1][markerindex + 1], genotypes[individuals - 1][markerindex + 2]);//based on genotype of individual need to be phased
//            fprintf(stderr,"marker:%d\tout of (%d) states, after condional on data (0,0):%d\t(0,1):%d\t(1,1):%d\n",i,GetStateNumFrom(i),genotypes[individuals - 1][markerindex],genotypes[individuals - 1][markerindex + 1], genotypes[individuals - 1][markerindex + 2]);
////            Wrapper->PrintVector(Wrapper->clusterAllele[i-1],"state allele");
////            Wrapper->PrintVector(Wrapper->haplotypeCluster[i-1],"state");
//            Wrapper->PrintVector(Wrapper->haplotypeCluster[i],"state");
//            Wrapper->PrintVector(Wrapper->clusterAllele[i],"state allele");
//            printLeftMatrix(leftMatrices[i],GetStateNumFrom(i));

            from = leftMatrices[i];

        }

    }

    MarkMemoryPool();
}

void PBWTHaplotyper::SwapIndividuals(int a, int b)
{
    // if (b < 0 || b >= individuals)
    //   printf("Bad Swap!");

    Swap(genotypes[a], genotypes[b]);
    Swap(haplotypes[a * 2], haplotypes[b * 2]);
    Swap(haplotypes[a * 2 + 1], haplotypes[b * 2 + 1]);

    if(nSampleCopy>0) {
        Swap(sampledHaps[a * 2], sampledHaps[(individuals - phased - 1) * 2]);
        Swap(sampledHaps[a * 2 + 1], sampledHaps[(individuals - phased - 1) * 2 + 1]);
    }
    if (diseaseCount)
    {
        Swap(diseaseStatus[a], diseaseStatus[b]);
        Swap(diseaseScores[a], diseaseScores[b]);
    }

    if (weights != nullptr)
    {
        float temp = weights[a];
        weights[a] = weights[b];
        weights[b] = temp;
    }
}

void PBWTHaplotyper::PrepareRefSetPBWTWrapper()
{
    if(Wrapper!= nullptr) {
        delete Wrapper;
        Wrapper = nullptr;
    }
    Wrapper = new PBWTWrapper(2*phased, markers, PvalueMatrix);
    Wrapper->SetHaps(haplotypes,2*(individuals-phased),2*individuals, nullptr,0,0, thetas);
    Wrapper->CursorBackwards();//calculate backwards order of suffix
    Wrapper->CursorForwards(false);
}
int PBWTHaplotyper::LoopThroughChromosomesHighPrecision() {

    ResetCrossovers();

    if(useRev) ReverseInput();
    clock_t t=clock();
    if(Wrapper!= nullptr) {
        delete Wrapper;
        Wrapper = nullptr;
    }
    Wrapper = new PBWTWrapper(2*individuals+(individuals-phased)*nSampleCopy*2, markers, PvalueMatrix);
    Wrapper->SetHaps(haplotypes,0,2*individuals,sampledHaps,0,(individuals-phased)*nSampleCopy*2, thetas);
    Wrapper->CursorBackwards();//calculate backwards order of suffix
    Wrapper->CursorForwards(false);
    clock_t t1=clock();
    printf("[HighPrecision]build model time:%.2f sec\n", (float) (t1 - t) / CLOCKS_PER_SEC);
    for (int i = individuals - 1; i >= 0; i--) {

        if (i < individuals - phased) {
            indexBeingSampled=i;
            SwapIndividuals(i, individuals - 1);


#ifdef DEBUG
             {
                Wrapper->PrintHap(tmpHaps, Wrapper->a[0]);

                // Wrapper->PrintHap(tmpHaps,Wrapper->a[6]);
                Wrapper->PrintHap(tmpHaps, Wrapper->a[Wrapper->N - 1]);
                // Wrapper->PrintMatrix(Wrapper->a,"a array matrix");
                Wrapper->PrintMatrix(Wrapper->d, "d array");
                //Wrapper->PrintVector(Wrapper->a[Wrapper->N-7],"last a array");
            }
#endif

            fprintf(stderr, "[HighPrecision]phasing individual %d...\n", i);
//            ScoreLeftConditional();
            ForwardAlgorithm();
//            t=clock();
//            printf("forward algorithm time:%.2f sec\n", (float) (t-t1) / CLOCKS_PER_SEC);

//            SampleChromosomes(&globalRandom);
            BackwardSampling(&globalRandom,individuals-1,haplotypes);

            for (int j = 0; j <nSampleCopy; ++j) {//n copy per individual
                BackwardSampling(&globalRandom,j+i*nSampleCopy,sampledHaps);
            }
//            t1=clock();
//            printf("sampling time:%.2f sec\n", (float) (t1-t) / CLOCKS_PER_SEC);

#ifdef _DEBUG
            if (!SanityCheck())
               {
               printf("\nProblems above occurred haplotyping individual %d\n\n", i);
               Print();
               }
#endif
            SwapIndividuals(i, individuals - 1);
            ResetFwdValues();
        }
    }
    t=clock();
    printf("[HighPrecision]forward algorithm and sampling time:%.2f sec\n", (float) (t-t1) / CLOCKS_PER_SEC);
    if(useRev) ReverseInput();
    return 0;
}

int PBWTHaplotyper::LoopThroughChromosomesSingleRound() {

    ResetCrossovers();

    if(useRev) ReverseInput();
    clock_t t=clock();
    if(Wrapper!= nullptr) {
        delete Wrapper;
        Wrapper = nullptr;
    }
    Wrapper = new PBWTWrapper(2*phased, markers, PvalueMatrix);
    Wrapper->SetHaps(haplotypes,2*(individuals-phased),2*individuals,nullptr,0,0, thetas);//only copy phased haps into pbwt
    Wrapper->CursorBackwards();//calculate backwards order of suffix
    Wrapper->CursorForwards(false);
    clock_t t1=clock();
    printf("[SingleRound]build model time:%.2f sec\n", (float) (t1 - t) / CLOCKS_PER_SEC);
    for (int i = individuals - 1; i >= 0; i--) {

        if (i < individuals - phased) {
            indexBeingSampled=i;
            SwapIndividuals(i, individuals - 1);
#ifdef DEBUG
            {
                Wrapper->PrintHap(tmpHaps, Wrapper->a[0]);

                // Wrapper->PrintHap(tmpHaps,Wrapper->a[6]);
                Wrapper->PrintHap(tmpHaps, Wrapper->a[Wrapper->N - 1]);
                // Wrapper->PrintMatrix(Wrapper->a,"a array matrix");
                Wrapper->PrintMatrix(Wrapper->d, "d array");
                //Wrapper->PrintVector(Wrapper->a[Wrapper->N-7],"last a array");
            }
#endif

            fprintf(stderr, "[SingleRound]phasing individual %d...\n", i);
//            ScoreLeftConditional();
            ForwardAlgorithm();
//            t=clock();
//            printf("forward algorithm time:%.2f sec\n", (float) (t-t1) / CLOCKS_PER_SEC);

//            SampleChromosomes(&globalRandom);
            BackwardSampling(&globalRandom,individuals-1,haplotypes);

            for (int j = 0; j <nSampleCopy; ++j) {//n copy per individual
                BackwardSampling(&globalRandom,j+i*nSampleCopy,sampledHaps);
            }
//            t1=clock();
//            printf("sampling time:%.2f sec\n", (float) (t1-t) / CLOCKS_PER_SEC);

#ifdef _DEBUG
            if (!SanityCheck())
               {
               printf("\nProblems above occurred haplotyping individual %d\n\n", i);
               Print();
               }
#endif
            SwapIndividuals(i, individuals - 1);
            ResetFwdValues();
        }
    }
    t=clock();
    printf("[SingleRound]forward algorithm and sampling time:%.2f sec\n", (float) (t-t1) / CLOCKS_PER_SEC);
    if(useRev) ReverseInput();
    return 0;
}

int PBWTHaplotyper::LoopThroughChromosomesViaPBWTWithHeterOnly() {

    ResetCrossovers();

    for (int i = individuals - 1; i >= 0; i--) {

        if (i < individuals - phased) {
            indexBeingSampled=i;
            SwapIndividuals(i, individuals - 1);
            if(useRev) ReverseInput();
            clock_t t=clock();
            ExtractHeterSites(individuals-1);
            Wrapper->SetHaps(haplotypes,0,2*individuals,sampledHaps,0,(individuals-phased)*nSampleCopy*2, thetas);
            Wrapper->CursorBackwards();//calculate backwards order of suffix
            Wrapper->CursorForwards(false);

#ifdef DEBUG
            {
                Wrapper->PrintHap(tmpHaps, Wrapper->a[0]);

                // Wrapper->PrintHap(tmpHaps,Wrapper->a[6]);
                Wrapper->PrintHap(tmpHaps, Wrapper->a[Wrapper->N - 1]);
                // Wrapper->PrintMatrix(Wrapper->a,"a array matrix");
                Wrapper->PrintMatrix(Wrapper->d, "d array");
                //Wrapper->PrintVector(Wrapper->a[Wrapper->N-7],"last a array");
            }
#endif
            printf("%d markers used for individual %d\n",markers,i);
            clock_t t1=clock();
            printf("build model time:%.2f sec\n", (float) (t1 - t) / CLOCKS_PER_SEC);
            //exit(EXIT_SUCCESS);
//            if (weights != NULL)
//                ScaleWeights();
//
//            if (updateDiseaseScores)
//                ScoreNPL();

            //if (i < individuals - phased) {
            fprintf(stderr, "phasing individual %d...\n", i);
//            ScoreLeftConditional();
            ForwardAlgorithm();
            t=clock();
            printf("forward algorithm time:%.2f sec\n", (float) (t-t1) / CLOCKS_PER_SEC);

//            SampleChromosomes(&globalRandom);
            BackwardSampling(&globalRandom,individuals-1,haplotypes);

            for (int j = 0; j <nSampleCopy; ++j) {//n copy per individual
                BackwardSampling(&globalRandom,j+i*nSampleCopy,sampledHaps);
            }
            t1=clock();
            printf("sampling time:%.2f sec\n", (float) (t1-t) / CLOCKS_PER_SEC);
            //exit(EXIT_SUCCESS);

//            if (updateDiseaseScores && diseaseCount)
//                IntegrateNPL();

#ifdef _DEBUG
            if (!SanityCheck())
               {
               printf("\nProblems above occurred haplotyping individual %d\n\n", i);
               Print();
               }
#endif
            FillHeterSitesBack(individuals - 1);
            SwapIndividuals(i, individuals - 1);
            //Wrapper->ResetWrapper();
            if(useRev) ReverseInput();
        }


    }

    return 0;
}



static float * GetProbability(float* source, int stateA, int stateB)
{
    //00,10,11,20,21,22,30,31,32,33...
    if(stateB>stateA) std::swap(stateA,stateB);
    return &source[stateA*(stateA+1)/2+stateB];
}

static float * GetOutput(float *dest, int stateA, int stateB)
{
    //00,10,11,20,21,22,30,31,32,33...
    if(stateB>stateA) std::swap(stateA,stateB);
    return &dest[stateA*(stateA+1)/2+stateB];
}

bool PBWTHaplotyper::ReverseInput()
{
    int begin=0;
    int end=markers-1;
    for (; begin <end; ++begin,--end) {
        for (int i = 0; i <individuals; ++i) {
            //haplotypes
            std::swap(haplotypes[i*2][begin],haplotypes[i*2][end]);
            std::swap(haplotypes[i*2+1][begin],haplotypes[i*2+1][end]);
            //genotypes
            std::swap(genotypes[i][begin*3],genotypes[i][end*3]);
            std::swap(genotypes[i][begin*3+1],genotypes[i][end*3+1]);
            std::swap(genotypes[i][begin*3+2],genotypes[i][end*3+2]);
        }

        for (int i = 0; i <(individuals-phased)*nSampleCopy; ++i) {
            //haplotypes
            std::swap(sampledHaps[i*2][begin],sampledHaps[i*2][end]);
            std::swap(sampledHaps[i*2+1][begin],sampledHaps[i*2+1][end]);
        }
        //penetrance
        for (int j = 0; j <3 ; ++j) {
            for (int k = 0; k < 3; ++k) {
                std::swap(Penetrance(begin, j, k),Penetrance(end, j, k));
            }
        }
    }

    begin=0;
    end=markers-2;
    for (; begin <end; ++begin,--end) {
        std::swap(thetas[begin],thetas[end]);
    }
    return true;
}

#define HMM_PRUNE 1
#define HMM_PRUNE2 1
void PBWTHaplotyper::Transpose(int site, float *source, float *dest)//site indicate dest marker index
{
    bool passOnce=false;
    float sum = 0.0;
    float factor=1.0;
    float *probability = source;
    float *output = dest;
    *output=0.0;
    int fromWhere = site - 1;//used by transvector, because vector start from 0
    int toWhere = site;
    int numFromStates = GetStateNumFrom(fromWhere);
    int numToStates = GetStateNumFrom(toWhere);//number of cluster
//    fprintf(stderr,"site:%d\tstates:%d\tfromwhere:%d\ttowhere:%d\n",site,states,GetStateNumFrom(fromWhere),GetStateNumFrom(toWhere));
    //int currentIndividualOriginalState1 = GetCurrentIndividualState(fromWhere, 0);
    //int currentIndividualOriginalState2 = GetCurrentIndividualState(fromWhere, 1);
    // This final loop actually transposes the probabilities for each state
    if(max_num < 1e-2) {
       factor = 10e6;
    }
    max_num=0;
//    if (weights == NULL)

        for (int k = 0; k < numToStates; ++k) {

            for (int l = 0; l < k; ++l) {
                *output = 0.;
#ifdef HMM_PRUNE2
//                if(site==8356)fprintf(stderr,"marker:%d,(%d,%d),geno:%d and GL:%d,%d,%d\n",toWhere,k,l,GetAllele(toWhere, k) + GetAllele(toWhere, l),genotypes[individuals - 1][toWhere * 3],genotypes[individuals - 1][toWhere * 3 +1],genotypes[individuals - 1][toWhere * 3 +2]);

                if ( genotypes[individuals - 1][toWhere * 3 + GetAllele(toWhere, k) + GetAllele(toWhere, l)] > 10) {
                    output++;
//                    fprintf(stderr,"marker:%d,(%d,%d),geno:%d and GL:%d,%d,%d\n",toWhere,k,l,GetAllele(toWhere, k) + GetAllele(toWhere, l),genotypes[individuals - 1][toWhere * 3],genotypes[individuals - 1][toWhere * 3 +1],genotypes[individuals - 1][toWhere * 3 +2]);

                    continue;
                }
#endif
//                probability=source;
//                for (int i = 0; i < numFromStates; i++) {
//                    for (int j = 0; j < i; j++) {
//                        *output +=(*probability) *
//                                  (GetTransitionProb(fromWhere, i, k) * GetTransitionProb(fromWhere, j, l)+
//                                  GetTransitionProb(fromWhere, i, l) * GetTransitionProb(fromWhere, j, k));
//                        probability++;
//                    }//end of inner loop, then deal with homo haplotype
//                    *output +=
//                            (*probability) * GetTransitionProb(fromWhere, i, k) * GetTransitionProb(fromWhere, i, l)*2;
//                    probability++;
//                }
//                for (auto kvi:Wrapper->inEdges[toWhere][k]) {
//                    for(auto kvj:Wrapper->inEdges[toWhere][l]){
                for (auto kvi:Wrapper->Graph.StateNodeMat[toWhere][k]->parentNodeIndex2NumHap) {
                    for(auto kvj:Wrapper->Graph.StateNodeMat[toWhere][l]->parentNodeIndex2NumHap){
                        sum=(*GetProbability(source,kvi.first,kvj.first))
                            * GetTransitionProb(fromWhere, kvi.first, k)
                            * GetTransitionProb(fromWhere, kvj.first, l);

                        if(kvi==kvj) sum*=2;//compensate for 2 combinations
//                        if((site==8355||site==8356)&&*GetProbability(source,kvi.first,kvj.first)>0)
//                        {
//                            fprintf(stderr,"(%d,%d)to (%d,%d):probability:%g\tk-first:%g\tl-second:%g\tk-second:%g\tl-first:%g\t\n",k,l,kvi.first,kvj.first,*GetProbability(source,kvi.first,kvj.first),
//                                    GetTransitionProb(fromWhere, k, kvi.first), GetTransitionProb(fromWhere, l, kvj.first),GetTransitionProb(fromWhere, k, kvj.first) , GetTransitionProb(fromWhere, l, kvi.first));
//                        }
                        if((sum < std::numeric_limits<float>::min()) &&
                                    (*GetProbability(source,kvi.first,kvj.first) > 0.)) sum=std::numeric_limits<float>::min();
                        *output += sum;
                    }
                }

                *output *=factor;
                if(*output>max_num) {
                    max_num=*output;
                }
                output++;

            }
            //end of inner loop, then deal with homo haplotype
            *output=0.;

#ifdef HMM_PRUNE2
//            if(site==8356)fprintf(stderr,"marker:%d,(%d,%d),geno:%d and GL:%d,%d,%d\n",toWhere,k,k,GetAllele(toWhere, k) + GetAllele(toWhere, k),genotypes[individuals - 1][toWhere * 3],genotypes[individuals - 1][toWhere * 3 +1],genotypes[individuals - 1][toWhere * 3 +2]);
            if (genotypes[individuals - 1][toWhere * 3 + GetAllele(toWhere, k) + GetAllele(toWhere, k)] > 10) {
                output++;
//                fprintf(stderr,"(%d,%d),geno:%d and GL:%d,%d,%d\n",k,k,GetAllele(toWhere, k) + GetAllele(toWhere, k),genotypes[individuals - 1][toWhere * 3],genotypes[individuals - 1][toWhere * 3 +1],genotypes[individuals - 1][toWhere * 3 +2]);
                continue;
            }
#endif
//            probability=source;
//            for (int i = 0; i < numFromStates; i++) {
//                for (int j = 0; j < i; j++) {
//                    *output +=
//                            *probability * (GetTransitionProb(fromWhere, i, k) * GetTransitionProb(fromWhere, j, k));
//                    probability++;
//                }//end of inner loop, then deal with homo haplotype
//                *output +=
//                        *probability * GetTransitionProb(fromWhere, i, k) * GetTransitionProb(fromWhere, i, k);
//                probability++;
//            }

//            for (int m = 0; m <numFromStates ; ++m) {
//                fprintf(stderr,"site:%d:(from state %d to state %d) %f\t\n",site,m,k,GetTransitionProb(fromWhere, m, k));
//            }

            for (auto iter=Wrapper->Graph.StateNodeMat[toWhere][k]->parentNodeIndex2NumHap.begin();iter!=Wrapper->Graph.StateNodeMat[toWhere][k]->parentNodeIndex2NumHap.end();++iter) {
                for(auto iter2=Wrapper->Graph.StateNodeMat[toWhere][k]->parentNodeIndex2NumHap.begin();iter2!=iter;++iter2){
                    sum=(*GetProbability(source,iter->first,iter2->first))
                        * GetTransitionProb(fromWhere, iter->first, k)
                        * GetTransitionProb(fromWhere, iter2->first, k);

                    if((sum < std::numeric_limits<float>::min()) &&
                       (*GetProbability(source,iter->first,iter2->first) > 0.)) sum=std::numeric_limits<float>::min();
                    *output += sum;
                }
                sum=(*GetProbability(source,iter->first,iter->first)) *  GetTransitionProb(fromWhere, iter->first, k)* GetTransitionProb(fromWhere, iter->first, k);
//                if((site==8355||site==8356)&&*GetProbability(source,iter->first,iter->first)>0)
//                {
//                    fprintf(stderr,"(%d,%d)to (%d,%d):probability:%g\tk-first:%g\tl-second:%g\tk-second:%g\tl-first:%g\t\n",k,k,iter->first,iter->first,*GetProbability(source,iter->first,iter->first),
//                            GetTransitionProb(fromWhere, k, iter->first), GetTransitionProb(fromWhere, k, iter->first),GetTransitionProb(fromWhere, k, iter->first) , GetTransitionProb(fromWhere, k, iter->first));
//                }
                if((sum < std::numeric_limits<float>::min()) &&
                   (*GetProbability(source,iter->first,iter->first) > 0.)) sum=std::numeric_limits<float>::min();
                *output += sum;
            }
            *output *= factor;
            if(*output>max_num) {
                max_num=*output;
            }
            output++;

        }
}

void  PBWTHaplotyper::ImputeAlleles(int marker, int state1, int state2, Random *rand, int currentIndividual, char** haps) {


    int currentHap1 = 2 * currentIndividual;
    int currentHap2 = currentHap1 + 1;

    int copied1 =GetAllele(marker,state1);//Wrapper->clusterAllele[marker][state1];//haplotypes[state1][marker];
    int copied2 =GetAllele(marker,state2); //Wrapper->clusterAllele[marker][state2];//haplotypes[state2][marker];
//    fprintf(stdout,"marker %d copied genotype: %d|%d\n",marker,copied1,copied2);

    int markerindex = marker * 3;
//    int ph11 = (unsigned char) genotypes[states / 2][markerindex];
//    int ph12 = (unsigned char) genotypes[states / 2][markerindex + 1];
//    int ph22 = (unsigned char) genotypes[states / 2][markerindex + 2];
    int ph11 = (unsigned char) genotypes[individuals - 1][markerindex];
    int ph12 = (unsigned char) genotypes[individuals - 1][markerindex + 1];
    int ph22 = (unsigned char) genotypes[individuals - 1][markerindex + 2];

    CalculatePhred2Prob();

    double posterior_11 = Penetrance(marker, copied1 + copied2, 0) * phred2prob[ph11];
    double posterior_12 = Penetrance(marker, copied1 + copied2, 1) * phred2prob[ph12];
    double posterior_22 = Penetrance(marker, copied1 + copied2, 2) * phred2prob[ph22];
    double sum = posterior_11 + posterior_12 + posterior_22;

    posterior_11 /= sum;
    posterior_22 /= sum;

    double r = rand->Next();

    if (r < posterior_11)//homo ref alleles
    {
        if(copied1!=copied2||copied1!=0)
            fprintf(stdout,"individidual: %d,Homo ref Marker:%d from else (%d,%d) and orginal gl:(%d,%d,%d)\n",currentIndividual,marker,copied1,copied2,ph11,ph12,ph22);

        haps[currentHap1][marker] = 0;
        haps[currentHap2][marker] = 0;
    }
    else if (r < posterior_11 + posterior_22)//home alt alleles
    {
        if(copied1!=copied2||copied1!=1)
            fprintf(stdout,"individual: %d,Homo alt Marker:%d from else (%d,%d) and orginal gl:(%d,%d,%d)\n",currentIndividual,marker,copied1,copied2,ph11,ph12,ph22);
        haps[currentHap1][marker] = 1;
        haps[currentHap2][marker] = 1;
    }
    else if (copied1 != copied2)//heter states and heter alleles
    {
//        double rate = GetErrorRate(marker);

//        if (rand->Next() < rate * rate / ((rate * rate) + (1.0 - rate) * (1.0 - rate)))//if both alleles mutated
//        {
//            copied1 = !copied1;
//            copied2 = !copied2;
//        }

        haps[currentHap1][marker] = copied1;
        haps[currentHap2][marker] = copied2;
    }
    else//homo states but heter alleles
    {
        fprintf(stdout,"individual: %d, Heter Marker:%d from else (%d,%d) and orginal gl:(%d,%d,%d)\n",currentIndividual,marker,copied1,copied2,ph11,ph12,ph22);
        bool bit = rand->Binary();

        haps[currentHap1][marker] = bit;
        haps[currentHap2][marker] = bit ^ 1;
    }

    int imputed1 = haps[currentHap1][marker];
    int imputed2 = haps[currentHap2][marker];
    //fprintf(stdout,"imputed genotype: %d|%d\n",imputed1,imputed2);
    //int differences = abs(copied1 - imputed1) + abs(copied2 - imputed2);
    int differences = abs(copied1 + copied2 - imputed1 - imputed2);

    error_models[marker].matches += 2 - differences;
    error_models[marker].mismatches += differences;
}

void  PBWTHaplotyper::ImputeAllelesRaw(int marker, int state1, int state2, Random *rand, int currentIndividual, char** haps) {


    int currentHap1 = 2 * currentIndividual;
    int currentHap2 = currentHap1 + 1;

    int copied1 =GetAllele(marker,state1);
    int copied2 =GetAllele(marker,state2);

    haps[currentHap1][marker] = copied1;
    haps[currentHap2][marker] = copied2;

    int imputed1 = haps[currentHap1][marker];
    int imputed2 = haps[currentHap2][marker];

    int differences = abs(copied1 + copied2 - imputed1 - imputed2);

    error_models[marker].matches += 2 - differences;
    error_models[marker].mismatches += differences;
}

void PBWTHaplotyper::ImputeAllele(int haplotype, int marker, int state, char** haps) {
    // if (updateDiseaseScores) UpdateDiseaseScores(marker, state);

    haps[haplotype][marker] = GetAllele(marker, state);
}

void PBWTHaplotyper::FillPath(int haplotype, int fromMarker, int toMarker, int state, char** haps) {
    fromMarker++;

    while (fromMarker < toMarker)
        ImputeAllele(haplotype, fromMarker++, state, haps);
}

void PBWTHaplotyper::SampleChromosomes(Random *rand) {
    int fromWhere = markers - 1;

    //Print(markers - 1);
//    RewindMemoryPool();
//    UpdateStateNum(GetStateNumFrom(fromWhere));
//    RetrieveMemoryBlock(fromWhere);

    float *probability = leftMatrices[fromWhere];
    float sum = 0.0;
    float sampleSum = 0.0;
    float sumPhase = 0.0;//initial phase could be either side

    int currentIndividual;
    int first = 0, second = 0;
    float choice;

    for (int copy = 0; copy < nSampleCopy + 1; ++copy) {
        if (copy == nSampleCopy) currentIndividual = individuals - 1;
        else currentIndividual = indexBeingSampled * nSampleCopy + copy;
        fromWhere = markers - 1;

        RewindMemoryPool();
        UpdateStateNum(GetStateNumFrom(fromWhere));
        RetrieveMemoryBlock(fromWhere);

        probability = leftMatrices[fromWhere];
        sum=0.0;
        // Calculate sum over all states
        for (int i = 0; i < states; i++) {
            for (int j = 0; j <= i; j++) {
#ifdef HMM_PRUNE//this is for sync problem
                if (genotypes[individuals - 1][fromWhere * 3 + GetAllele(fromWhere, i) + GetAllele(fromWhere, j)] >
                    10) {
                    probability++;
//                fprintf(stderr, "(%d,%d):geno:%d\t", i, j, GetAllele(fromWhere, i) + GetAllele(fromWhere, j));
                    continue;
                }
#endif
                sum += *probability;
//            fprintf(stderr, "(%d,%d):%f\t", i, j, *probability);
                probability++;
            }
//        fprintf(stderr,"\n");
        }
        // Sample number and select state
        choice = rand->Uniform(0, sum);
        sampleSum = 0.0;

        first = 0;
        second = 0;
        for (probability = leftMatrices[fromWhere]; first < states; first++) {
            for (second = 0; second <= first; second++) {
#ifdef HMM_PRUNE
                if (genotypes[individuals - 1][fromWhere * 3 + GetAllele(fromWhere, first) +
                                               GetAllele(fromWhere, second)] > 10) {
                    probability++;
                    continue;
                }
#endif
                sampleSum += *probability;
                probability++;

                if (sampleSum >= choice) break;
            }
            if (second <= first) break;
        }


    float max_prob(0),tmp_sum(0);
        float last_sum(sum);
        int j0;
        for (int j = markers - 2; j >= 0; j--) {

//            fprintf(stderr,
//                    "marker:%d\tTotalSum:%g\tsampleSum: %9.9g, choice:%9.9g,Chose (%d,%d) out of (%d) states and geno:%d and GL:%d,%d,%d\n",
//                    j + 1, last_sum, sampleSum, choice, first, second, GetStateNumFrom(j + 1),
//                    GetAllele(fromWhere, first) + GetAllele(fromWhere, second),
//                    genotypes[individuals - 1][fromWhere * 3], genotypes[individuals - 1][fromWhere * 3 + 1],
//                    genotypes[individuals - 1][fromWhere * 3 + 2]);
//        max_prob=0;
            last_sum = 0.;

            if (copy == nSampleCopy) {
                ImputeAlleles(j + 1, first, second, rand, currentIndividual, haplotypes);
            }
            else {
                ImputeAlleles(j + 1, first, second, rand, currentIndividual, sampledHaps);
            }

            j0 = j;
            fromWhere = j;//TODO: remain question
            sum = 0.0;
            UpdateStateNum(GetStateNumFrom(fromWhere));
            RetrieveMemoryBlock(fromWhere);
            probability = leftMatrices[fromWhere];


            for (int k = 0; k < states; k++) {
                for (int l = 0; l < k; l++, probability++) {

#ifdef HMM_PRUNE
                    if (genotypes[individuals - 1][fromWhere * 3 + GetAllele(fromWhere, k) + GetAllele(fromWhere, l)] >
                        10) {
//                            fprintf(stderr,"(%d,%d):geno:%d\t",k,l,GetAllele(fromWhere, k) + GetAllele(fromWhere, l));
                        continue;
                    }
#endif

//                        tmp_sum=sum;
                    sum += *probability *
                           (GetTransitionProb(fromWhere, k, first) *
                            GetTransitionProb(fromWhere, l, second)
                            + GetTransitionProb(fromWhere, k, second) *
                              GetTransitionProb(fromWhere, l, first));
//                    if(sum-tmp_sum>0)fprintf(stderr,"(%d,%d)to (%d,%d):probability:%g\tk-first:%g\tl-second:%g\tk-second:%g\tl-first:%g\t\n",k,l,first,second,*probability,GetTransitionProb(fromWhere, k, first), GetTransitionProb(fromWhere, l, second),GetTransitionProb(fromWhere, k, second) , GetTransitionProb(fromWhere, l, first));

                }
#ifdef HMM_PRUNE
                if (genotypes[individuals - 1][fromWhere * 3 + GetAllele(fromWhere, k) + GetAllele(fromWhere, k)] >
                    10) {
                    probability++;
//                    fprintf(stderr,"(%d,%d):geno:%d\n",k,k,GetAllele(fromWhere, k) + GetAllele(fromWhere, k));
                    continue;
                }
#endif

//                tmp_sum=sum;
                sum += *probability *
                       GetTransitionProb(fromWhere, k, first) *
                       GetTransitionProb(fromWhere, k, second);

//                if(sum-tmp_sum>0) fprintf(stderr,"(%d,%d) to (%d,%d):probability:%g\tmid:%g\n",k,k,first,second,*probability ,GetTransitionProb(fromWhere, k, first)*GetTransitionProb(fromWhere, k, second));
                probability++;
            }
            last_sum = sum;
            // Sample number and decide how many state changes occurred between the
            // two positions
            choice = rand->Uniform(0, sum);


            // But perhaps the first or second haplotype recombined
            probability = leftMatrices[j];

            // Try to select any other state
            sampleSum = 0.0;
            sumPhase = 0.0;//resolve phasing
            // Save the original states
            int first0 = first;
            int second0 = second;

            for (first = 0; first < states; first++) {
                for (second = 0; second < first; second++, probability++) {
#ifdef HMM_PRUNE
                    if (genotypes[individuals - 1][fromWhere * 3 + GetAllele(fromWhere, first) +
                                                   GetAllele(fromWhere, second)] > 10) {
                        continue;
                    }
#endif
                    sumPhase = sampleSum + *probability *
                                     GetTransitionProb(fromWhere, first, first0) *
                                     GetTransitionProb(fromWhere, second, second0);
                    sampleSum += *probability *
                           (GetTransitionProb(fromWhere, first, first0) *
                            GetTransitionProb(fromWhere, second, second0)
                            + GetTransitionProb(fromWhere, first, second0) *
                              GetTransitionProb(fromWhere, second, first0));

                    if (sampleSum >= choice) {
                        goto JUMPOUT;
                        //break;
                    }
                }
#ifdef HMM_PRUNE
                if (genotypes[individuals - 1][fromWhere * 3 + GetAllele(fromWhere, first) +
                                               GetAllele(fromWhere, second)] > 10) {
                    probability++;
                    continue;
                }
#endif
                sampleSum += *probability * (GetTransitionProb(fromWhere, first, first0) *
                                       GetTransitionProb(fromWhere, first, second0));
                probability++;

                if (sampleSum >= choice) {
                    //second=first;
                    break;
                }
                // if (second <= first) break;
            }
            JUMPOUT:
            if (sumPhase < choice) {//TODO:To confirm if it's the former phase, greater or less than
                int temp = first;
                first = second;
                second = temp;
            }
            if (first >= states or second >= states) {
                std::cerr << "Sampled states out of state space!Abort" << std::endl;
                exit(EXIT_FAILURE);
            }

            // Record outcomes for intermediate, uninformative, positions
            if(copy == nSampleCopy) {
                FillPath(currentIndividual * 2, j, j0 + 1, first,haplotypes);
                FillPath(currentIndividual * 2 + 1, j, j0 + 1, second,haplotypes);
            }
            else
            {
                FillPath(currentIndividual * 2, j, j0 + 1, first,sampledHaps);
                FillPath(currentIndividual * 2 + 1, j, j0 + 1, second, sampledHaps);
            }
        }


        if (copy == nSampleCopy) {
            ImputeAlleles(0, first, second, rand, currentIndividual, haplotypes);
        }
        else {
            ImputeAlleles(0, first, second, rand, currentIndividual, sampledHaps);
        }
    }


}

//HMM version two
int PBWTHaplotyper::InitialFwdValues()
{
    brokenList.assign(markers,false);
    fwdValueSum=new float [markers];
    for (int k = 0; k <markers ; ++k) {
        fwdValueSum[k]=0.f;
    }
    int SampleIndex=individuals-1;
//    currentFwdValuePtr=&currentFwdValue;
//    nextFwdValuePtr=&nextFwdValue;
    float prior=1./(states*states);
    float tmpFwdValue(0.),gl(0.);

    for (int i = 0; i < states; ++i) {
        for (int j = 0; j <states ; ++j) {
            int allele1=GetAllele(0,i);
            int allele2=GetAllele(0,j);
            gl=GetGL(SampleIndex,0,allele1,allele2);
//            fprintf(stderr,"site:%d,prev(0,0) to current(%d,%d) : allell1:%d\tallele2:%d\tgl:%g\n",
//                    0,i,j,allele1,allele2,gl);
            if(gl>1e-1)
            {
                tmpFwdValue = prior * gl;
                if (tmpFwdValue < UNDERFLOW_MIN) {
                    tmpFwdValue = UNDERFLOW_MIN;
                }
//                (*currentFwdValuePtr)[i][j].push_back(origin(0, 0, tmpFwdValue));
                fwdValueSum[0] += tmpFwdValue;
                genuienParents[0][i][j][std::make_pair(0,0)]=tmpFwdValue;
            }
        }
    }

    for(auto iter= genuienParents[0].begin();iter!=genuienParents[0].end();++iter)
        for(auto iter2=iter->second.begin();iter2!=iter->second.end();++iter2)
            for (auto iter3 = iter2->second.begin(); iter3 != iter2->second.end(); ++iter3) {
                iter3->second/=fwdValueSum[0];
            }
    return 0;
}

int PBWTHaplotyper::ResetFwdValues()
{
    delete [] fwdValueSum;

    for (int k = 0; k <markers; ++k) {
        genuienParents[k].clear();
    }
    return 0;
}

double PBWTHaplotyper::GetGL(int individual, int marker, char allele1, char allele2)
{
//    std::cerr<<individual<<"\t"<<marker<<"\t"<<(int)allele1<<"\t"<<(int)allele2<<std::endl;

    return phred2prob[(size_t)genotypes[individual][3*marker+allele1+allele2]];
}

void PBWTHaplotyper::RandomSetup(Random * rand)
{
    if (rand == NULL)
        rand = &globalRandom;

    CalculatePhred2Prob();

    for (int j = 0; j < markers; j++)
    {
        int markerindex = 3*j;
        for (int i = 0; i < individuals-phased; i++)
        {

            int posterior_11 =  genotypes[i][markerindex];
            int posterior_12 =  genotypes[i][markerindex+1];
            int posterior_22 =  genotypes[i][markerindex+2];
            int min =std::min(std::min(posterior_11,posterior_12),posterior_22);//phred score
            if (min == posterior_11)
            {
                haplotypes[i * 2][j] = 0;
                haplotypes[i * 2 + 1][j] = 0;
            }
            else if (min== posterior_12)
            {
                bool bit = rand->Binary();

                haplotypes[i * 2][j] = bit;
                haplotypes[i * 2 + 1][j] = bit ^ 1;
            }
            else
            {
                haplotypes[i * 2][j] = 1;
                haplotypes[i * 2 + 1][j] = 1;
            }

        }
    }
}



int PBWTHaplotyper::ForwardAlgorithm()
{
    int SampleIndex=individuals-1;
    UpdateStateNum(GetStateNumFrom(0));
    InitialFwdValues();
    float prevFwdValue(0.f);
    float tmpFwdValue(0.f);
    float gl(0.f);

    int fitPair=0;
    int totalPair=0;
    int noChildPair=0;

    char allele1,allele2;
    int childNode1,childNode2;
    for (int i = 1; i < markers; i++) {
        fitPair=0;
        totalPair=0;
        noChildPair=0;
        fwdValueSum[i]=0;

        for(auto iter=genuienParents[i-1].begin();iter!=genuienParents[i-1].end();++iter)//all states at site i-1, iter->first: hap1
            for (auto iter2 = iter->second.begin(); iter2 !=iter->second.end(); ++iter2)// iter2->first:hap2; iter2->second:all the source states to current state
            {
                prevFwdValue=SumFwdValueFromOriginVec(iter2->second);//all the parents that lead to current child pair
                totalPair++;
                for ( allele1 = 0; allele1 <2 ; ++allele1) {
                    childNode1=GetChildNode(i-1,iter->first,allele1);//i child of i-1 child state
                    if(childNode1==-1) continue;
                    for ( allele2 = 0; allele2 <2 ; ++allele2) {
                        childNode2=GetChildNode(i-1,iter2->first,allele2);//i child of i-1 child state
                        if(childNode2==-1) continue;
                        gl = GetGL(SampleIndex, i, allele1, allele2);//i gl
//                        fprintf(stderr,"site:%d,prev(%d,%d) to current(%d,%d) : allell1:%d\tallele2:%d\tedgeNumHap1:%g\tedgeNumHap2:%g\tgl:%g\n",
//                                    i,iter->first,iter2->first,childNode1,childNode2,allele1,allele2,GetTransitionProb(i-1,iter->first,childNode1),GetTransitionProb(i-1,iter2->first,childNode2),gl);
                        if (gl > 1e-1 or i < 10)
                        {
                            fitPair++;
                            tmpFwdValue=prevFwdValue*
                                        GetTransitionProb(i-1,iter->first,childNode1)*
                                        GetTransitionProb(i-1,iter2->first,childNode2)*
                                        gl;//i fwdValue
                            if(tmpFwdValue < UNDERFLOW_MIN && prevFwdValue > 0)
                            {
                                tmpFwdValue= UNDERFLOW_MIN;
                            }
//                            fprintf(stderr,"site:%d,prev(%d,%d) to current(%d,%d) : %g and prevFwd:%g\ttp1:%g\ttp2:%g\tgl:%g fwdValueSum:%g\n",
//                                    i,iter->first,iter2->first,childNode1,childNode2,tmpFwdValue,prevFwdValue,
//                                    GetTransitionProb(i-1,iter->first,childNode1),
//                                    GetTransitionProb(i-1,iter2->first,childNode2),gl,fwdValueSum[i]);
//                            (*nextFwdValuePtr)[childNode1][childNode2].push_back(origin(iter->first,iter2->first,tmpFwdValue));
                            fwdValueSum[i]+=tmpFwdValue;
                            genuienParents[i][childNode1][childNode2][std::make_pair(iter->first,iter2->first)]=tmpFwdValue;
                        }
//                        else
//                            fprintf(stderr,"site:%d,prev(%d,%d) to current(%d,%d) pass: 0 and prevFwd:%g\tgl:%g\n",
//                                    i,iter->first,iter2->first,childNode1,childNode2,prevFwdValue,gl);
                    }
                }
            }
        if (fitPair == 0)
        {

            fprintf(stderr,"[Warning]marker %d broken! %d totalPair, %d noChildPair\n",i,totalPair,noChildPair);
            exit(EXIT_FAILURE);
            //brokenList[i] = true;

            UpdateStateNum(GetStateNumFrom(i));
            prevFwdValue = 1.f / (states * states);
            for (int j = 0; j < states; ++j)
                for (int k = 0; k < states; ++k) {
                    allele1 = GetAllele(i, j);
                    allele2 = GetAllele(i, k);
                    gl = GetGL(SampleIndex, i, allele1, allele2);//i gl
                    if (gl > 1e-1) {
                        tmpFwdValue = prevFwdValue * gl;//i fwdValue
                        if (tmpFwdValue < UNDERFLOW_MIN && prevFwdValue > 0) {
                            tmpFwdValue = UNDERFLOW_MIN;
                        }
                        fwdValueSum[i] += tmpFwdValue;
                        genuienParents[i][j][k][std::make_pair(0, 0)] = tmpFwdValue;
                    }
                }
        }
//        fprintf(stderr,"site:%d report overall fwd:%g\n",i,fwdValueSum[i]);
        for(auto iter=genuienParents[i].begin();iter!=genuienParents[i].end();++iter)
            for(auto iter2=iter->second.begin();iter2!=iter->second.end();++iter2)
                for (auto iter3 = iter2->second.begin(); iter3 != iter2->second.end(); ++iter3) {
                    iter3->second/=fwdValueSum[i];
                }

//        SwitchFwdValuePtr();

//        if(fitPair > 0) continue;
//        else
//        {
//            fprintf(stderr,"at marker %d, not viable path!!!!\n",i);
//            exit(EXIT_FAILURE);
//        }

    }
    return 0;
}
int PBWTHaplotyper::BackwardSampling(Random *rand, int SampleIndex, char** sampledHaps)
{
    double choice(0.);
    double sum(0.);
    float gl0(0.f);
    float transProb0(0.f),transProb1(0.f);
    float prevFwdValue(0.f);
    int sampledFirst(0),sampledSecond(0),first0(0),second0(0);
    origin chosenOrigin{0,0,0.f};

    choice=rand->Uniform(0,1);
//    for(auto iter=genuienParents[markers-1].begin();iter!=genuienParents[markers-1].end();++iter)
//        for(auto iter2=iter->second.begin();iter2!=iter->second.end();++iter2)
//        {
//            prevFwdValue = SumFwdValueFromOriginVec(iter2->second);
//            fprintf(stderr,"(%d,%d)\tvalue:%g\n",iter->first,iter2->first,prevFwdValue);
//
//        }
    for(auto iter=genuienParents[markers-1].begin();iter!=genuienParents[markers-1].end();++iter)
        for(auto iter2=iter->second.begin();iter2!=iter->second.end();++iter2)
        {
                prevFwdValue = SumFwdValueFromOriginVec(iter2->second);
                sum += prevFwdValue;
                if (sum > choice) {
                    sampledFirst = iter->first;
                    sampledSecond = iter2->first;
                    goto SAMPLE_BREAK;
                }
        }
    SAMPLE_BREAK:
    ImputeAlleles(markers-1, sampledFirst, sampledSecond, rand, SampleIndex, sampledHaps);
//    fprintf(stderr,"marker:%d\tsum:%g\tchoice:%g\t(%d,%d)\tvalue:%g\toverallFwd:%g\n",markers-1,sum,choice,sampledFirst,sampledSecond,prevFwdValue,fwdValueSum[markers-1]);
    sum=0.;
        choice=rand->Uniform(0,prevFwdValue);
        for (auto kv: genuienParents[markers-1][sampledFirst][sampledSecond])
        {
//            transProb0 = GetTransitionProb(markers - 2, kv.first.first, sampledFirst);
//            transProb1 = GetTransitionProb(markers - 2, kv.first.second, sampledSecond);
//            prevFwdValue = kv.second * fwdValueSum[markers - 1] / (transProb0 * transProb1 * gl0);
            sum += kv.second;
            if(sum>choice)
            {
                chosenOrigin.firstState=kv.first.first;
                chosenOrigin.secondState=kv.first.second;
                gl0 = GetGL(individuals-1, markers-1, GetAllele(markers-1, sampledFirst), GetAllele(markers-1, sampledSecond));
                transProb0 = GetTransitionProb(markers - 2, chosenOrigin.firstState, sampledFirst);
                transProb1 = GetTransitionProb(markers - 2, chosenOrigin.secondState, sampledSecond);
                chosenOrigin.fwdValue=kv.second*fwdValueSum[markers-1]/(transProb0 * transProb1 * gl0);
                break;
            }
        }


    for (int i = markers-2; i >0; --i) {
//        fprintf(stderr,"marker:%d\tsum:%g\tchoice:%g\t(%d,%d)\tvalue:%g\n",i,sum,choice,chosenOrigin.firstState,chosenOrigin.secondState,chosenOrigin.fwdValue);
//        fprintf(stderr,"marker:%d\tsum:%g\n",i,sum);

        sum=0.;
        first0=chosenOrigin.firstState;
        second0=chosenOrigin.secondState;

        ImputeAlleles(i, first0, second0, rand, SampleIndex, sampledHaps);

        {
            if (fabs(chosenOrigin.fwdValue - UNDERFLOW_MIN) < std::numeric_limits<double>::epsilon())
                chosenOrigin.fwdValue = SumFwdValueFromOriginVec(genuienParents[i][first0][second0]);

            choice = rand->Uniform(0, chosenOrigin.fwdValue);
//        fprintf(stderr,"fwd:%g and summation:%g\n",chosenOrigin.fwdValue,SumFwdValueFromOriginVec(genuienParents[i][first0][second0]));
            if (genuienParents[i].find(first0) == genuienParents[i].end() ||
                genuienParents[i][first0].find(second0) == genuienParents[i][first0].end()) {
                fprintf(stderr, "marker %d does not have state %d or %d\n", i, first0, second0);
                exit(EXIT_FAILURE);
            }

            for (auto kv:genuienParents[i][first0][second0])//we are actually sampling i-1's states
            {
                sum += kv.second;//kv.second is normalized
//                fprintf(stderr, "from (%d,%d) to (%d,%d) : %g\n", first0, second0, kv.first.first, kv.first.second,
//                        kv.second);
                if (sum > choice) {

                    chosenOrigin.firstState = kv.first.first;
                    chosenOrigin.secondState = kv.first.second;

                    gl0 = GetGL(individuals - 1, i, GetAllele(i, first0), GetAllele(i, second0));

                    transProb0 = GetTransitionProb(i - 1, chosenOrigin.firstState, first0);
                    transProb1 = GetTransitionProb(i - 1, chosenOrigin.secondState, second0);
                    chosenOrigin.fwdValue = kv.second * fwdValueSum[i] / (transProb0 * transProb1 * gl0);
//                fprintf(stderr, "try parents:(%d,%d) at marker %d with fwdvalue:%g\tgl:%g\ttp1:%g\ttp2:%g\ttotalFwd:%g\n",
//                        kv.first.first,kv.first.second,
//                        i, kv.second, gl0, transProb0, transProb1,fwdValueSum[i]);

                    break;
                }
//            else
//            {
//                fprintf(stderr, "try parents:(%d,%d) failed at marker %d with fwdvalue:%g\tgl:%g\ttp1:%g\ttp2:%g\n",
//                        kv.first.first,kv.first.second, i, kv.second, gl0,
//                        transProb0, transProb1);
//            }
            }
        }
//        FillPath(SampleIndex * 2, i, i0 + 1, chosenOrigin.firstState,haplotypes);
//        FillPath(SampleIndex * 2 + 1, i, i0 + 1, chosenOrigin.secondState,haplotypes);
    }

    ImputeAlleles(0, chosenOrigin.firstState, chosenOrigin.secondState, rand, SampleIndex, sampledHaps);
    return 0;
}

int PBWTHaplotyper::ExtractHeterSites(int individualToProcess) {//apply after swap individualToProcess to the back

    if(Wrapper!= nullptr) {
        delete Wrapper;
        Wrapper = nullptr;
    }
    absoluteIndexToRelative.clear();
    relativeIndexToAbsolute.clear();
    tmpMarkers=0;

    if(onlyHeterSite){

        std::vector<bool> HeterIndex(markers,false);
        for (int i = 0; i < markers ; ++i) {
//            fprintf(stderr,"through heter partat marker %d: %d\t%d\t%d\n",i,genotypes[individualToProcess][i*3],genotypes[individualToProcess][i*3+1],genotypes[individualToProcess][i*3+2]);
            if(genotypes[individualToProcess][i*3+1] < genotypes[individualToProcess][i*3]&&
                    genotypes[individualToProcess][i*3+1] < genotypes[individualToProcess][i*3+2])
            {
                HeterIndex[i] = true;
                absoluteIndexToRelative[i] = tmpMarkers++;
                relativeIndexToAbsolute.push_back(i);
            }
        }

        for (int i = 0; i < tmpMarkers ; ++i) {
            int markerAbsoluteNow=relativeIndexToAbsolute[i];
            if(!HeterIndex[markerAbsoluteNow])
            {
                continue;
            }//homo
            for (int j = 0; j <individuals; ++j) {

                tmpHaps[2*j][i]=haplotypes[2*j][markerAbsoluteNow];
                tmpHaps[2*j+1][i]=haplotypes[2*j+1][markerAbsoluteNow];
                tmpGeno[j][i*3]=genotypes[j][markerAbsoluteNow*3];
                tmpGeno[j][i*3+1]=genotypes[j][markerAbsoluteNow*3+1];
                tmpGeno[j][i*3+2]=genotypes[j][markerAbsoluteNow*3+2];
            }
            for (int j = 0; j <(individuals-phased)*nSampleCopy; ++j) {
                tmpHaps[individuals*2+2*j][i]=sampledHaps[2*j][markerAbsoluteNow];
                tmpHaps[individuals*2+2*j+1][i]=sampledHaps[2*j+1][markerAbsoluteNow];
            }
            for (int j = 0; j <3 ; ++j) {
                for (int k = 0; k < 3; ++k) {
                    tmpPenetrance[i * 9 + j * 3 + k]=Penetrance(markerAbsoluteNow, j, k);
                }
            }
        }
        SwapTempHaps();
    }
    else
    {
        tmpMarkers=markers;
    }

    if(tmpMarkers==0)
    {
        fprintf(stderr,"found 0 markers available...abort!\n");
        abort();
    }
    Wrapper = new PBWTWrapper(2*individuals+(individuals-phased)*nSampleCopy*2, tmpMarkers, PvalueMatrix);

    return 0;
}

int PBWTHaplotyper::FillHeterSitesBack(int individualToProcess) {

    if(onlyHeterSite) {
        SwapTempHaps();
        int markerAbsoluteNow = 0;
        for (int i = 0; i < tmpMarkers; ++i) {
            markerAbsoluteNow = relativeIndexToAbsolute[i];
            haplotypes[2 * individualToProcess][markerAbsoluteNow] = tmpHaps[2 * individualToProcess][i];
            haplotypes[2 * individualToProcess + 1][markerAbsoluteNow] = tmpHaps[2 * individualToProcess+1][i];
        }
    }
    for (auto & parents: genuienParents) {
        parents.clear();
    }
    return 0;
}
