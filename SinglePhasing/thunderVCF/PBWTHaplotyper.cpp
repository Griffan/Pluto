//
// Created by Fan Zhang on 8/6/15.
//

#include "PBWTHaplotyper.h"

PBWTHaplotyper::PBWTHaplotyper(int nhaps, int nsnps) : Wrapper(nhaps, nsnps) {

}

PBWTHaplotyper::~PBWTHaplotyper() {

}


bool PBWTHaplotyper::ForceMemoryAllocation() {
    // Cycle through individuals, with the exact same steps as the actual
    // haplotyper and request memory ... by requesting all memory upfront,
    // we force crashes to happen early.
    for (int i = 0; i < individuals - phased; i++) {
        ResetMemoryPool();
        UpdateStateNum(Wrapper.clusterAllele[0].size());
        GetMemoryBlock(0);

        if (leftMatrices[0] == NULL)
            return false;

        int skipped = 0;
        for (int j = 1; j < markers; j++)
            //if (genotypes[i][j] != GENOTYPE_MISSING || j == markers - 1)
        {
            UpdateStateNum(Wrapper.clusterAllele[j].size());
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
        UpdateStateNum(Wrapper.clusterAllele[j].size());
        GetSmallMemoryBlock(j);

        if (leftMatrices[j] == NULL)
            return false;
    }

    return true;
}

void PBWTHaplotyper::ConditionOnData(float *matrix, int marker, char phred11, char phred12, char phred22) {
    // We treat missing genotypes as uninformative about the mosaic's
    // underlying state. If we were to allow for deletions and the like,
    // that may no longer be true.

    //if (genotype == GENOTYPE_MISSING)
    //return;

    double conditional_probs[3];
    int ph11 = (unsigned char) phred11;
    int ph12 = (unsigned char) phred12;
    int ph22 = (unsigned char) phred22;

    CalculatePhred2Prob();

    for (int i = 0; i < 3; i++)
        conditional_probs[i] = Penetrance(marker, i, 0) * phred2prob[ph11] +
                               Penetrance(marker, i, 1) * phred2prob[ph12] +
                               Penetrance(marker, i, 2) * phred2prob[ph22];

    for (int i = 0; i < states; i++) {
        double factors[2];

//        factors[0] = conditional_probs[haplotypes[i][marker]];
//        factors[1] = conditional_probs[haplotypes[i][marker] + 1];
        factors[0] = conditional_probs[Wrapper.clusterAllele[marker][i]];
        factors[1] = conditional_probs[Wrapper.clusterAllele[marker][i] + 1];

        for (int j = 0; j <= i; j++, matrix++)
            *matrix *= factors[Wrapper.clusterAllele[marker][j]];
    }
}

void PBWTHaplotyper::ScoreLeftConditional() {
    ResetMemoryPool();
    GetMemoryBlock(0);

    UpdateStateNum(Wrapper.clusterAllele[0].size());

    SetupPrior(leftMatrices[0]);
    ConditionOnData(leftMatrices[0], 0, genotypes[individuals - 1][0], genotypes[individuals - 1][1],
                    genotypes[individuals - 1][2]);

    double theta = 0.0;
    float *from = leftMatrices[0];
    for (int i = 1; i < markers; i++) {
        int markerindex = i * 3;

        // Cumulative recombination fraction allows us to skip uninformative positions
        theta = theta + thetas[i - 1] - theta * thetas[i - 1];//TODO:need to be updated accordingly

        // Skip over uninformative positions to save time
        // maybe check for the difference between min(phred11, phred12, phred22) and the next smallest

        //if (genotypes[states / 2][i] != GENOTYPE_MISSING || i == markers - 1)
        {
            UpdateStateNum(Wrapper.clusterAllele[i].size());
            GetMemoryBlock(i);

            Transpose(i,from, leftMatrices[i], theta);
            ConditionOnData(leftMatrices[i], i, genotypes[individuals - 1][markerindex],
                            genotypes[individuals - 1][markerindex + 1], genotypes[individuals - 1][markerindex +
                                                                                                    2]);//based on genotype of individual need to be phased

            theta = 0;
            from = leftMatrices[i];
        }
    }

    MarkMemoryPool();
}

int PBWTHaplotyper::LoopThroughChromosomesViaPBWT() {
    ResetCrossovers();


    for (int i = individuals - 1; i >= 0; i--) {
        SwapIndividuals(i, individuals - 1);

        if (weights != NULL)
            ScaleWeights();

        if (updateDiseaseScores)
            ScoreNPL();

        if (i < individuals - phased) {

            ScoreLeftConditional();
            SampleChromosomes(&globalRandom);

            if (updateDiseaseScores && diseaseCount)
                IntegrateNPL();

#ifdef _DEBUG
         if (!SanityCheck())
            {
            printf("\nProblems above occurred haplotyping individual %d\n\n", i);
            Print();
            }
#endif
        }
        else {
            ScoreLeftConditionalForHaplotype();
            SampleHaplotypeSource(&globalRandom);
            SwapHaplotypes(states, states + 1);
            ScoreLeftConditionalForHaplotype();
            SampleHaplotypeSource(&globalRandom);
            SwapHaplotypes(states, states + 1);
        }

        SwapIndividuals(i, individuals - 1);
    }

    return 0;
}

float PBWTHaplotyper::getTransitionProb(int site, int from, int to) {
    return Wrapper.transVector[site][from][to];
}

void PBWTHaplotyper::Transpose(int site, float *source, float *dest,
                                      float theta)//site indicate dest marker index
{
    if (theta == 0.0) {
        for (int i = 0; i < states; i++)
            for (int j = 0; j <= i; j++, dest++, source++)
                *dest = *source;

        return;
    }

    float sum = 0.0;
    float *probability = source;
    float *output = dest;

    int fromWhere = site - 1;//used by transvector, because vector start from 0
    int toWhere = site;
    int fromStates = states;
    int toStates = Wrapper.clusterAllele[toWhere].size();

    // This final loop actually transposes the probabilities for each state
    if (weights == NULL)
        for (int k = 0; k < toStates; ++k)
        {
            for (int l = 0; l < k; ++l)
            {
                for (int i = 0; i < fromStates; i++)
                {
                    for (int j = 0; j < i; j++)
                    {
                        *output += *probability * (getTransitionProb(fromWhere, i, k) * getTransitionProb(fromWhere, j, l)+getTransitionProb(fromWhere, i, l) * getTransitionProb(fromWhere, j, k));
                        probability++;
                    }//end of inner loop, then deal with homo haplotype
                    //TODO: homo haplotype
                    *output += *probability * getTransitionProb(fromWhere, i, k) * getTransitionProb(fromWhere, i, l) * 2;
                    probability++;
                }
                output++;
            }//end of inner loop, then deal with homo haplotype
        }
}

void PBWTHaplotyper::SampleChromosomes(Random *rand) {
    ShotgunHaplotyper::SampleChromosomes(rand);
}
