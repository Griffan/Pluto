//
// Created by Fan Zhang on 8/6/15.
//

#include "PBWTHaplotyper.h"
#include "../../pbwtWrapper/PBWTWrapper.h"

PBWTHaplotyper::PBWTHaplotyper() {

}

PBWTHaplotyper::~PBWTHaplotyper() {

}

void PBWTHaplotyper::ScoreLeftConditional() {
    ResetMemoryPool();
    GetMemoryBlock(0);

    UpdateStateNum(Wrapper.numCluster[0]);

    SetupPrior(leftMatrices[0]);
    ConditionOnData(leftMatrices[0], 0, genotypes[states / 2][0], genotypes[states / 2][1], genotypes[states / 2][2]);

    double theta = 0.0;
    float *from = leftMatrices[0];
    for (int i = 1; i < markers; i++)
    {
        int markerindex = i*3;

        // Cumulative recombination fraction allows us to skip uninformative positions
        theta = theta + thetas[i - 1] - theta * thetas[i - 1];

        // Skip over uninformative positions to save time
        // maybe check for the difference between min(phred11, phred12, phred22) and the next smallest

        //if (genotypes[states / 2][i] != GENOTYPE_MISSING || i == markers - 1)
        {
            GetMemoryBlock(i);

            Transpose(from, leftMatrices[i], theta);
            ConditionOnData(leftMatrices[i], i, genotypes[states / 2][markerindex],
                            genotypes[states / 2][markerindex+1], genotypes[states / 2][markerindex+2]);

            theta = 0;
            from = leftMatrices[i];
        }
    }

    MarkMemoryPool();
}

int PBWTHaplotyper::LoopThroughChromosomesViaPBWT() {
    ResetCrossovers();


    for (int i = individuals - 1; i >= 0; i--)
    {
        SwapIndividuals(i, individuals - 1);

        if (weights != NULL)
            ScaleWeights();

        if (updateDiseaseScores)
            ScoreNPL();

        if (i < individuals - phased)
        {

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
        else
        {
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


