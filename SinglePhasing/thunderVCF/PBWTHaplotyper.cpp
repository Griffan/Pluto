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
		UpdateStateNum(getStateNumFrom(0));
        GetMemoryBlock(0);

        if (leftMatrices[0] == NULL)
            return false;

        int skipped = 0;
        for (int j = 1; j < markers; j++)
            //if (genotypes[i][j] != GENOTYPE_MISSING || j == markers - 1)
        {
			UpdateStateNum(getStateNumFrom(j));
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
		UpdateStateNum(getStateNumFrom(j));
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

	UpdateStateNum(getStateNumFrom(0));

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
			UpdateStateNum(getStateNumFrom(i));
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
    int numFromStates = states;
	int numToStates = getStateNumFrom(toWhere);//number of cluster

	int currentIndividualOriginalState1 = getCurrentIndividualState(fromWhere,0);
	int currentIndividualOriginalState2 = getCurrentIndividualState(fromWhere,1);
    // This final loop actually transposes the probabilities for each state
    if (weights == NULL)
		for (int k = 0; k < numToStates; ++k)
        {
            for (int l = 0; l < k; ++l)
            {
				for (int i = 0; i < numFromStates; i++)
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

void  PBWTHaplotyper::ImputeAlleles(int marker, int state1, int state2, Random * rand)
{
	// if (updateDiseaseScores)
	//   {
	//   UpdateDiseaseScores(marker, state1);
	//   UpdateDiseaseScores(marker, state2);
	//   }

	int imputed1 = getAllele(marker,state1);
	int imputed2 = getAllele(marker,state2);

	int genotype = genotypes[individuals-1][marker];

	int currentHap1 = 2 * (individuals - 1);
	int currentHap2 = currentHap1 + 1;

	if (genotype != GENOTYPE_HOMOZYGOUS_FOR_ONE &&
		genotype != GENOTYPE_HOMOZYGOUS_FOR_TWO)
	{
		haplotypes[currentHap1][marker] = imputed1;
		haplotypes[currentHap2][marker] = imputed2;
	}

	if (genotype == GENOTYPE_MISSING) return;

	int differences = abs(genotype - imputed1 - imputed2 - 1);

	if (genotype == GENOTYPE_HETEROZYGOUS && differences == 0)
		error_models[marker].uncertain_pairs++;
	else
	{
		error_models[marker].matches += 2 - differences;
		error_models[marker].mismatches += differences;
	}

	if (genotype != GENOTYPE_HETEROZYGOUS) return;

	if (imputed1 == imputed2)
		if (rand->Binary())
			haplotypes[currentHap1][marker] = !imputed2;
		else
			haplotypes[currentHap2][marker] = !imputed1;
}
void PBWTHaplotyper::SampleChromosomes(Random *rand) {
      // Print(markers - 1);
   RewindMemoryPool();
   UpdateStateNum(getStateNumFrom(markers - 1));
   RetrieveMemoryBlock(markers - 1);

   float * probability = leftMatrices[markers - 1];
   float sum = 0.0;

   // Calculate sum over all states
   for (int i = 0; i < states; i++)
      for (int j = 0; j <= i; j++)
         {
         sum += *probability;
         probability++;
         }

   // Sample number and select state
   float choice = rand->Uniform(0, sum);

   sum = 0.0;

   int first = 0, second = 0;
   for (probability = leftMatrices[markers - 1]; first < states; first++)
      {
      for (second = 0; second <= first; second++)
         {
         sum += *probability;
         probability++;

         if (sum >= choice) break;
         }

      if (second <= first) break;
      }

   // printf("Cumulative probability: %g\n", sum);
   // printf("           Random draw: %g\n", choice);
   // printf("        Selected state: %g\n", *(probability - 1));

   for (int j = markers - 2; j >= 0; j--)
      {
      //printf("Sum: %f, Chose (%d,%d)\n", sum, first, second);
      ImputeAlleles(j + 1, first, second, rand);//TODO:modify needed

      // Starting marker for this iteration
      int   j0 = j;

      // Cumulative recombination fraction, skipping over uninformative
      // positions
      float theta = thetas[j];
      //while (genotypes[states / 2][j] == GENOTYPE_MISSING && j > 0)
         //{
         //--j;
         //theta = theta + thetas[j] - theta * thetas[j];
         //}

      // When examining the previous location we consider three alternatives:
      // states that could be reached when both haplotypes recombine (11),
      // states that can be reached when the first (10) or second (01) haplotype recombines,
      // and the states that can be reached without recombination.

      float sum00 = 0.0, sum01 = 0.0, sum10 = 0.0, sum11 = 0.0;
	  UpdateStateNum(getStateNumFrom(j));
      RetrieveMemoryBlock(j);
      probability = leftMatrices[j];

      for (int k = 0; k < states; k++)
         for (int l = 0; l <= k; l++, probability++)
            {
            sum11 += *probability;
            if (first == k || first == l) sum01 += *probability;
            if (second == k || second == l) sum10 += *probability;
            if (first == k && second == l || first == l && second == k) sum00 += *probability;
            }

      if (weights != NULL)
         {
         sum01 *= weights[second / 2];
         sum10 *= weights[first / 2];
         sum11 *= weights[second / 2] * weights[first / 2];
         }

      sum = sum11 * theta * theta / (states * states) +
            (sum10 + sum01) * theta * (1.0 - theta) / states +
            sum00 * (1.0 - theta) * (1.0 - theta);

      // Sample number and decide how many state changes occurred between the
      // two positions
      choice = rand->Uniform(0, sum);

      // The most likely outcome is that no changes occur ...
      choice -= sum00 * (1.0 - theta) * (1.0 - theta);
      if (choice <= 0.0)
         {
         // Record outcomes for intermediate, uninformative, positions
         FillPath(states, j, j0 + 1, first);
         FillPath(states + 1, j, j0 + 1, second);

         continue;
         }

      // But perhaps the first or second haplotype recombined
      probability = leftMatrices[j];

      choice -= sum10 * theta * (1.0 - theta) / states;
      if (choice <= 0.0)
         {
         // The first haplotype changed ...
         choice = choice * states / (theta * (1.0 - theta));

         // Record the original state
         int first0 = first;

         if (weights != NULL) choice /= weights[first / 2];

         for (first = 0; first < states; first++)
            {
            if (first >= second)
               choice += probability[first * (first + 1) / 2 + second];
            else
               choice += probability[second * (second + 1) / 2 + first];

            if (choice >= 0.0) break;
            }

         // Record outcomes for intermediate, uninformative, positions
         SamplePath(states, j, j0 + 1, first, first0, rand);
         FillPath(states + 1, j, j0 + 1, second);

         continue;
         }

      choice -= sum01 * theta * (1.0 - theta) / states;
      if (choice <= 0.0)
         {
         // The second haplotype changed ...
         choice = choice * states / (theta * (1.0 - theta));

         // Save the original state
         int second0 = second;

         if (weights != NULL) choice /= weights[second / 2];

         for (second = 0; second < states; second++)
            {
            if (first >= second)
               choice += probability[first * (first + 1) / 2 + second];
            else
               choice += probability[second * (second + 1) / 2 + first];

            if (choice >= 0.0) break;
            }

         // Record outcomes for intermediate, uninformative, positions
         FillPath(states, j, j0 + 1, first);
         SamplePath(states + 1, j, j0 + 1, second, second0, rand);

         continue;
         }

      // Try to select any other state
      choice *= states * states / (theta * theta);
      sum = 0.0;

      // Save the original states
      int first0 = first;
      int second0 = second;

      if (weights != NULL) choice /= weights[first / 2] * weights[second / 2];

      for (first = 0; first < states; first++)
         {
         for (second = 0; second <= first; second++, probability++)
            {
            sum += *probability;

            if (sum > choice) break;
            }

         if (second <= first) break;
         }

      if (rand->Binary())
         {
         int temp = first;
         first = second;
         second = temp;
         }

      // Record outcomes for intermediate, uninformative, positions
      SamplePath(states, j, j0 + 1, first, first0, rand);
      SamplePath(states + 1, j, j0 + 1, second, second0, rand);
      }

   ImputeAlleles(0, first, second, rand);
}
