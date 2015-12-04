//
// Created by Fan Zhang on 8/6/15.
//

#include <pbwt/pbwt.h>
#include "PBWTHaplotyper.h"
#include "MemoryAllocators.h"
//debug related
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

PBWTHaplotyper::PBWTHaplotyper() {
    onlyHeterSite=false;
    max_num=0;
}

void PBWTHaplotyper::InitAuxillary() {

    tmpHaps = AllocateCharMatrix(individuals * 2, markers);

    tmpGeno = AllocateCharMatrix(individuals, markers*3);

    Wrapper=NULL;
}

PBWTHaplotyper::~PBWTHaplotyper() {

    if(tmpHaps != NULL) {
        for (int i = 0; i < 2 * individuals; ++i) {
            delete[] tmpHaps[i];
        }
        delete [] tmpHaps;
    }
    if(tmpGeno != NULL) {
        for (int i = 0; i < individuals; ++i) {
            delete[] tmpGeno[i];
        }
        delete [] tmpGeno;
    }
    if(Wrapper != NULL)
        delete Wrapper;
	ReleaseMemoryBlock();
}

//memory management
void PBWTHaplotyper::ReleaseMemoryBlock()
{
	for (std::unordered_map<int, std::vector<float*> >::iterator iter = memoryBlockList.begin(); iter != memoryBlockList.end(); ++iter)
	{
		for (size_t i = 0; i < iter->second.size(); i++)
		{
			if (iter->second[i] != NULL) delete[] iter->second[i];
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
		memoryBlockList[blockSize] = std::vector<float*>(0, NULL);
	}
	if (numInUse[blockSize] < memoryBlockList[blockSize].size())
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
		memoryBlockList[blockSize] = std::vector<float*>(0, NULL);
	}
	if (numInUse[blockSize] < memoryBlockList[blockSize].size())
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
       // fprintf(stderr, "out from RetrieveMemory\n");
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

        if (leftMatrices[0] == NULL)
            return false;

        int skipped = 0;
        for (int j = 1; j < markers; j++)
            //if (genotypes[i][j] != GENOTYPE_MISSING || j == markers - 1)
        {
            GetMemoryBlock(j);

            if (leftMatrices[j] == NULL)
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

        if (leftMatrices[j] == NULL)
            return false;
    }

    return true;
}

//phasing
void PBWTHaplotyper::ConditionOnData(float *matrix, int marker, char phred11, char phred12, char phred22) {
    // We treat missing genotypes as uninformative about the mosaic's
    // underlying state. If we were to allow for deletions and the like,
    // that may no longer be true.

    //if (genotype == GENOTYPE_MISSING)
    //return;

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

    for (int i = 0; i < states; i++) {
        double factors[2];

//        factors[0] = conditional_probs[haplotypes[i][marker]];
//        factors[1] = conditional_probs[haplotypes[i][marker] + 1];
        factors[0] = conditional_probs[Wrapper->clusterAllele[marker][i]];
        factors[1] = conditional_probs[Wrapper->clusterAllele[marker][i] + 1];

        for (int j = 0; j <= i; j++, matrix++)
            *matrix *= factors[Wrapper->clusterAllele[marker][j]];
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

            //fprintf(stderr,"marker:%d\tout of (%d) states\n",i,GetStateNumFrom(i-1));
            //printLeftMatrix(from,GetStateNumFrom(i-1));

            ConditionOnData(leftMatrices[i], i, genotypes[individuals - 1][markerindex],genotypes[individuals - 1][markerindex + 1], genotypes[individuals - 1][markerindex + 2]);//based on genotype of individual need to be phased

            from = leftMatrices[i];

        }

    }

    MarkMemoryPool();
}

int PBWTHaplotyper::LoopThroughChromosomesViaPBWT() {

    ResetCrossovers();

	//char* tmphap1 = new char[markers];
	//char* tmphap2 = new char[markers];

    for (int i = individuals - 1; i >= 0; i--) {

        if (i < individuals - phased) {
            SwapIndividuals(i, individuals - 1);
            clock_t t=clock();
            ExtractHeterSites(individuals-1);
            Wrapper->SetHaps(haplotypes);
            Wrapper->CursorBackwards();//calculate backwards order of suffix
            Wrapper->CursorForwards();

            if(DEBUG) {
                Wrapper->PrintHap(tmpHaps, Wrapper->a[0]);

                // Wrapper->PrintHap(tmpHaps,Wrapper->a[6]);
                Wrapper->PrintHap(tmpHaps, Wrapper->a[Wrapper->N - 1]);
                // Wrapper->PrintMatrix(Wrapper->a,"a array matrix");
                Wrapper->PrintMatrix(Wrapper->d, "d array");
                //Wrapper->PrintVector(Wrapper->a[Wrapper->N-7],"last a array");
            }
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
           // fprintf(stderr, "phasing individual %d...\n", i);
            ScoreLeftConditional();
            t=clock();
            printf("forward algorithm time:%.2f sec\n", (float) (t-t1) / CLOCKS_PER_SEC);
            SampleChromosomes(&globalRandom);
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
            //Wrapper->resetWrapper();
        }

    }

    return 0;
}


void PBWTHaplotyper::Transpose(int site, float *source, float *dest)//site indicate dest marker index
{


    float sum = 0.0;
    float factor=1.0;
    float *probability = source;
    float *output = dest;
    *output=0.0;
    int fromWhere = site - 1;//used by transvector, because vector start from 0
    int toWhere = site;
    int numFromStates = GetStateNumFrom(fromWhere);
    int numToStates = GetStateNumFrom(toWhere);//number of cluster
    //fprintf(stderr,"site:%d\tstates:%d\tfromwhere:%d\ttowhere:%d\n",site,states,getStateNumFrom(fromWhere),GetStateNumFrom(toWhere));
    //int currentIndividualOriginalState1 = GetCurrentIndividualState(fromWhere, 0);
    //int currentIndividualOriginalState2 = GetCurrentIndividualState(fromWhere, 1);
    // This final loop actually transposes the probabilities for each state
   if(max_num < 10e-10) {
       factor = 10e5;
       max_num=0;
   }
//    if (weights == NULL)
        for (int k = 0; k < numToStates; ++k) {

            for (int l = 0; l < k; ++l) {

                *output=0.0;
                probability=source;
                for (int i = 0; i < numFromStates; i++) {
                    for (int j = 0; j < i; j++) {
//                        fprintf(stderr,"a:%f\t%d\t%d\n",GetTransitionProb(fromWhere, i, k),i,k);
//                        fprintf(stderr,"b:%f\t%d\t%d\n",GetTransitionProb(fromWhere, j, l),j,l);
//                        fprintf(stderr,"c:%f\t%d\t%d\n",GetTransitionProb(fromWhere, i, l),i,l);
//                        fprintf(stderr,"d:%f\t%d\t%d\n",GetTransitionProb(fromWhere, j, k),j,k);
//                        fprintf(stderr,"e:%f\t\n",*probability);
                        *output +=
                                *probability * (GetTransitionProb(fromWhere, i, k) * GetTransitionProb(fromWhere, j, l)
                                                + GetTransitionProb(fromWhere, i, l) * GetTransitionProb(fromWhere, j,
                                                                                                         k)
                                               );

                        probability++;
                    }//end of inner loop, then deal with homo haplotype
                    *output +=
                            *probability * GetTransitionProb(fromWhere, i, k) * GetTransitionProb(fromWhere, i, l);
                    probability++;
                }
                //fprintf(stderr,"site:%d:(%d,%d) %f\t\n",site,k,l,*output);
                if(*output>max_num) max_num=*output;
                *output *=factor;
                output++;

            }
            //end of inner loop, then deal with homo haplotype
            *output=0.0;
            probability=source;
            for (int i = 0; i < numFromStates; i++) {
                for (int j = 0; j < i; j++) {
                    *output +=
                            *probability * (GetTransitionProb(fromWhere, i, k) * GetTransitionProb(fromWhere, j, k));
                    probability++;
                }//end of inner loop, then deal with homo haplotype
                *output +=
                        *probability * GetTransitionProb(fromWhere, i, k) * GetTransitionProb(fromWhere, i, k);
                probability++;
            }
            //fprintf(stderr,"site:%d:(%d,%d) %f\t\n",site,k,k,*output);
            if(*output>max_num) max_num=*output;
            *output *= factor;
            output++;

        }

}

void  PBWTHaplotyper::ImputeAlleles(int marker, int state1, int state2, Random *rand) {


    int currentHap1 = 2 * (individuals - 1);
    int currentHap2 = currentHap1 + 1;

    int copied1 = Wrapper->clusterAllele[marker][state1];//haplotypes[state1][marker];
    int copied2 = Wrapper->clusterAllele[marker][state2];//haplotypes[state2][marker];
    //if(marker==106)fprintf(stdout,"marker %d copied genotype: %d|%d\t",marker,copied1,copied2);

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
    { //fprintf(stdout,"from 00\t");
        haplotypes[currentHap1][marker] = 0;
        haplotypes[currentHap2][marker] = 0;
    }
    else if (r < posterior_11 + posterior_22)//home alt alleles
    {//fprintf(stdout,"from 11\t");
        haplotypes[currentHap1][marker] = 1;
        haplotypes[currentHap2][marker] = 1;
    }
    else if (copied1 != copied2)//heter states and heter alleles
    {
        //fprintf(stdout,"from 01\t");
        double rate = GetErrorRate(marker);

        if (rand->Next() < rate * rate / ((rate * rate) + (1.0 - rate) * (1.0 - rate)))//if both alleles mutated
        {
            copied1 = !copied1;
            copied2 = !copied2;
        }

        haplotypes[currentHap1][marker] = copied1;
        haplotypes[currentHap2][marker] = copied2;
    }
    else//homo states but heter alleles
    {
        //fprintf(stdout,"from else\t");
        bool bit = rand->Binary();

        haplotypes[currentHap1][marker] = bit;
        haplotypes[currentHap2][marker] = bit ^ 1;
    }

    int imputed1 = haplotypes[currentHap1][marker];
    int imputed2 = haplotypes[currentHap2][marker];
    //fprintf(stdout,"imputed genotype: %d|%d\n",imputed1,imputed2);
    //int differences = abs(copied1 - imputed1) + abs(copied2 - imputed2);
    int differences = abs(copied1 + copied2 - imputed1 - imputed2);

    error_models[marker].matches += 2 - differences;
    error_models[marker].mismatches += differences;
}


void PBWTHaplotyper::ImputeAllele(int haplotype, int marker, int state) {
    // if (updateDiseaseScores) UpdateDiseaseScores(marker, state);

    haplotypes[haplotype][marker] = GetAllele(marker, state);
}

void PBWTHaplotyper::FillPath(int haplotype, int fromMarker, int toMarker, int state) {
    fromMarker++;

    while (fromMarker < toMarker)
        ImputeAllele(haplotype, fromMarker++, state);
}

void PBWTHaplotyper::SampleChromosomes(Random *rand) {
    int fromWhere = markers-1;
    int currentHap1 = 2 * (individuals - 1);
    int currentHap2 = currentHap1 + 1;
    //Print(markers - 1);
    RewindMemoryPool();
    UpdateStateNum(GetStateNumFrom(markers - 1));
    RetrieveMemoryBlock(markers - 1);

    float *probability = leftMatrices[markers - 1];
    float sum = 0.0;

    float lastSum=sum;

    // Calculate sum over all states
    for (int i = 0; i < states; i++)
        for (int j = 0; j <= i; j++) {
            sum += *probability;
            probability++;
        }

    // Sample number and select state
    float choice = rand->Uniform(0, sum);
    lastSum=sum;
    sum = 0.0;

    int first = 0, second = 0;
    for (probability = leftMatrices[markers - 1]; first < states; first++) {
        for (second = 0; second <= first; second++) {
            sum += *probability;
            //fprintf(stderr,"(%d,%d):%f\t",first,second,*probability);
            probability++;

            if (sum >= choice) break;
        }

        if (second <= first) break;
    }

    // printf("Cumulative probability: %g\n", sum);
    // printf("           Random draw: %g\n", choice);
    // printf("        Selected state: %g\n", *(probability - 1));

    for (int j = markers - 2; j >= 0; j--) {
 //       fprintf(stderr,"marker:%d\tlastSum:%9.9f\tSum: %9.9f, choice:%9.9f,Chose (%d,%d) out of (%d) states\n",j+1, lastSum,sum, choice,first, second,GetStateNumFrom(j+1));
//        //Wrapper->PrintVector(Wrapper->haplotypeCluster[j],"states");
//        //if(j==106)
//        printLeftMatrix(leftMatrices[j+1],GetStateNumFrom(j+1));

        ImputeAlleles(j + 1, first, second, rand);//TODO:modify needed, like the part before fillpath
        int j0 = j;

        sum=0.0;
        UpdateStateNum(GetStateNumFrom(j));
        RetrieveMemoryBlock(j);
        probability = leftMatrices[j];


        fromWhere=j;//TODO: remain question
        if(first!=second) {
            for (int k = 0; k < states; k++) {
                for (int l = 0; l < k; l++, probability++) {
                    sum += *probability *
                           (GetTransitionProb(fromWhere, k, first) * GetTransitionProb(fromWhere, l, second)
                            + GetTransitionProb(fromWhere, k, second) *
                              GetTransitionProb(fromWhere, l, first));
                    //fprintf(stderr,"(%d,%d):sum:%f\tleft:%f\tright:%f\t",k,l,*probability * (getTransitionProb(fromWhere, k, first) * getTransitionProb(fromWhere, l, second) + getTransitionProb(fromWhere, k, second) * getTransitionProb(fromWhere, l, first)),getTransitionProb(fromWhere, k, first) * getTransitionProb(fromWhere, l, second),GetTransitionProb(fromWhere, k, second) * getTransitionProb(fromWhere, l, first));
                }
                sum += *probability * GetTransitionProb(fromWhere, k, first) * GetTransitionProb(fromWhere, k, second);
                //fprintf(stderr,"(%d,%d):sum:%f\tmid1:%f\tmid2:%f\n",k,k,*probability * getTransitionProb(fromWhere, k, first) * GetTransitionProb(fromWhere, k, second),getTransitionProb(fromWhere, k, first),getTransitionProb(fromWhere, k, second));
                probability++;
            }
        }
        else
        {
            for (int k = 0; k < states; k++) {
                for (int l = 0; l < k; l++, probability++) {
                    sum += *probability * GetTransitionProb(fromWhere, k, first) * GetTransitionProb(fromWhere, l, second);
                    //fprintf(stderr,"(%d,%d):sum:%f\tleft:%f\tright:%f\t",k,l,*probability * (getTransitionProb(fromWhere, k, first) * getTransitionProb(fromWhere, l, second) + getTransitionProb(fromWhere, k, second) * getTransitionProb(fromWhere, l, first)),getTransitionProb(fromWhere, k, first) * getTransitionProb(fromWhere, l, second),GetTransitionProb(fromWhere, k, second) * getTransitionProb(fromWhere, l, first));
                }
                sum += *probability * GetTransitionProb(fromWhere, k, first) * GetTransitionProb(fromWhere, k, second);
                //fprintf(stderr,"(%d,%d):sum:%f\tmid1:%f\tmid2:%f\n",k,k,*probability * getTransitionProb(fromWhere, k, first) * GetTransitionProb(fromWhere, k, second),getTransitionProb(fromWhere, k, first),getTransitionProb(fromWhere, k, second));
                probability++;
            }
        }

        //fprintf(stderr,"site:%d finished\n",fromWhere);
        // Sample number and decide how many state changes occurred between the
        // two positions
        choice = rand->Uniform(0, sum);


        // But perhaps the first or second haplotype recombined
        probability = leftMatrices[j];
        lastSum=sum;
        // Try to select any other state
        sum = 0.0;
        float sumPerPair = 0.0;//resolve phasing
        // Save the original states
        int first0 = first;
        int second0 = second;
        if(first0 != second0) {
            for (first = 0; first < states; first++) {
                for (second = 0; second < first; second++, probability++) {
                    sumPerPair = sum + *probability * GetTransitionProb(fromWhere, first, first0) *
                                       GetTransitionProb(fromWhere, second, second0);
                    sum += *probability *
                           (GetTransitionProb(fromWhere, first, first0) * GetTransitionProb(fromWhere, second, second0)
                            + GetTransitionProb(fromWhere, first, second0) *
                              GetTransitionProb(fromWhere, second, first0));

                    if (sum >= choice) {
                        goto JUMPOUT;
                        //break;
                    }
                }

                sumPerPair = sum + *probability * GetTransitionProb(fromWhere, first, first0) *
                                   GetTransitionProb(fromWhere, first, second0);
                sum += *probability * (GetTransitionProb(fromWhere, first, first0) *
                                       GetTransitionProb(fromWhere, first, second0));
                probability++;

                if (sum >= choice) {
                    //second=first;
                    break;
                }

                // if (second <= first) break;
            }
            JUMPOUT:
            if (sumPerPair < choice) {//TODO:To confirm if it's the former phase, greater or less than
                int temp = first;
                first = second;
                second = temp;
            }
        }
        else {
            for (first = 0; first < states; first++) {
                for (second = 0; second < first; second++, probability++) {

                    sum += *probability *
                           (GetTransitionProb(fromWhere, first, first0) * GetTransitionProb(fromWhere, second, second0));

                    if (sum >= choice) {
                        goto JUMPOUT2;
                        //break;
                    }
                }

                sum += *probability * (GetTransitionProb(fromWhere, first, first0) *
                                       GetTransitionProb(fromWhere, first, second0));
                probability++;

                if (sum >= choice) {
                    //second=first;
                    break;
                }

                // if (second <= first) break;
            }

        }
        JUMPOUT2:
        // Record outcomes for intermediate, uninformative, positions
        FillPath(currentHap1, j, j0 + 1, first);
        FillPath(currentHap2, j, j0 + 1, second);

    }

    ImputeAlleles(0, first, second, rand);

}


int PBWTHaplotyper::ExtractHeterSites(int individualToProcess) {//apply after swap individualToProcess to the back

    if(Wrapper!=NULL) {
        delete Wrapper;

    }
    absoluteIndexToRelative.clear();
    relativeIndexToAbsolute.clear();
    tmpMarkers=0;

    if(onlyHeterSite){
        //fprintf(stderr,"through heeter part\n");
        std::vector<int> HeterIndex(markers,0);
        for (int i = 0; i < markers ; ++i) {
            if(haplotypes[2*individualToProcess][i]!=haplotypes[2*individualToProcess+1][i])
            {
                HeterIndex[i] = 1;
                absoluteIndexToRelative[i] = tmpMarkers++;
                relativeIndexToAbsolute.push_back(i);
            }
        }

        for (int j = 0; j <individuals-DEBUG; ++j) {
            for (int i = 0; i < tmpMarkers ; ++i) {
                int markerAbsoluteNow=relativeIndexToAbsolute[i];
                if(!HeterIndex[markerAbsoluteNow])
                {
                    continue;
                }//homo
                tmpHaps[2*j][i]=haplotypes[2*j][markerAbsoluteNow];
                tmpHaps[2*j+1][i]=haplotypes[2*j+1][markerAbsoluteNow];
                tmpGeno[j][i*3]=genotypes[j][markerAbsoluteNow*3];
                tmpGeno[j][i*3+1]=genotypes[j][markerAbsoluteNow*3+1];
                tmpGeno[j][i*3+2]=genotypes[j][markerAbsoluteNow*3+2];
            }
        }
        SwapTempHaps();
    }
    else
    {
//        for (int j = 0; j <individuals-DEBUG; ++j) {
//            for (int i = 0; i < markers ; ++i) {
//                tmpHaps[2*j][i]=haplotypes[2*j][i];
//                tmpHaps[2*j+1][i]=haplotypes[2*j+1][i];
//                tmpGeno[j][i*3]=genotypes[j][i*3];
//                tmpGeno[j][i*3+1]=genotypes[j][i*3+1];
//                tmpGeno[j][i*3+2]=genotypes[j][i*3+2];
//                if(j==0) {
//                    absoluteIndexToRelative[i] = tmpMarkers++;
//                    relativeIndexToAbsolute.push_back(i);
//                }
//
//            }
//        }
        tmpMarkers=markers;
    }

    if(tmpMarkers==0)
    {
        fprintf(stderr,"found 0 markers available...abort!\n");
        abort();
    }
    Wrapper = new PBWTWrapper(2*(individuals-1-DEBUG)/*not include individualToProcess*/, tmpMarkers);

    return 0;
}

int PBWTHaplotyper::FillHeterSitesBack(int individualToProcess) {

    if(onlyHeterSite) {
        SwapTempHaps();
        int markerAbsoluteNow = 0;
        for (int i = 0; i < tmpMarkers; ++i) {
            markerAbsoluteNow = relativeIndexToAbsolute[i];
            haplotypes[2 * (individualToProcess - 1)][markerAbsoluteNow] = tmpHaps[2 * (individualToProcess - 1)][i];
            haplotypes[2 * (individualToProcess - 1) + 1][markerAbsoluteNow] = tmpHaps[2 * (individualToProcess - 1) +
                                                                                       1][i];
        }
    }

//    if(Wrapper!=NULL) {
//        delete Wrapper;
//        Wrapper=NULL;
//    }
//    absoluteIndexToRelative.clear();
//    relativeIndexToAbsolute.clear();
//    tmpMarkers=0;

    return 0;
}
