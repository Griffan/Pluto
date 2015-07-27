////////////////////////////////////////////////////////////////////// 
// thunder/Main.cpp 
// (c) 2000-2008 Goncalo Abecasis
// 
// This file is distributed as part of the MaCH source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile MaCH.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Saturday April 12, 2008
// 

//static int NUM_NON_GLF_FIELDS = 9;

#include <vector>
#include <map>
#include <unordered_map>
#include <fstream>

#include "ShotgunHaplotyper.h"
#include "ShotgunManners.h"
#include "OutputHandlers.h"
#include "DosageCalculator.h"
#include "MergeHaplotypes.h"
#include "HaplotypeLoader.h"
#include "Parameters.h"
#include "InputFile.h"
#include "Error.h"
#include "libVcfVcfFile.h"
using namespace libVcf;

float * thetas = NULL;
int     nthetas = 0;

float * error_rates = NULL;
int     nerror_rates = 0;
std::unordered_map<std::string, int> unphaseMarkerIdx;//record marker name and relative index in unphased vcf
std::unordered_map<std::string, int> unphaseMarkerUIdx;//record marker name and relative index in phased ref vcf
std::unordered_map<std::string, bool> unphaseMarkerFlag;//record if this marker only shows in unphased vcf

std::unordered_map<std::string, bool> pidIncludedInUnphasedVcf;
std::unordered_map<std::string, bool> pidIncludedInPhasedVcf;
std::unordered_map<std::string, bool> pidExcludedInUnphasedVcf;
std::unordered_map<std::string, bool> pidExcludedInPhasedVcf;

std::unordered_map<std::string, std::pair<int, int> > DuplicatedIndividualPair;

// print output files directly in VCF format
// inVcf contains skeleton of VCF information to copy from
// consensus contains the haplotype information to replace GT field
// dosage contains dosage information to be added as DS field
// thetas contains recombination rate information between markers
// error-rates contains per-marker error rates
// rsqs contains rsq_hat estimates ??
void OutputVCFConsensus(const String& inVcf, Pedigree & ped, ConsensusBuilder & consensus, DosageCalculator &doses, const String& filename, float* thetas, float* error_rates, ShotgunHaplotyper& engine)
{
	consensus.Merge(); // calculate consensus sequence

	if (consensus.stored)
		printf("Merged sampled haplotypes to generate consensus\n"
		"   Consensus estimated to have %.1f errors in missing genotypes and %.1f flips in haplotypes\n\n",
		consensus.errors, consensus.flips);

	// read and write VCF inputs
	try {

		fprintf(stderr, "Outputing VCF file %s\n", filename.c_str());

		VcfFile* pVcf = new VcfFile;
		IFILE outVCF = ifopen(filename.c_str(), "wb");
		if (outVCF == NULL) {
			error("Cannot open output file %s for writing", filename.c_str());
			exit(-1);
		}

		pVcf->bSiteOnly = false;
		pVcf->bParseGenotypes = false;
		pVcf->bParseDosages = false;
		pVcf->bParseValues = true;
		pVcf->openForRead(inVcf.c_str());

		// add proper header information
		//for (int i = 1; i < pVcf->asMetaKeys.Length(); ++i) {
		//	if ((pVcf->asMetaKeys[i - 1].SubStr(0, 4).Compare("INFO") == 0) && (pVcf->asMetaKeys[i].SubStr(0, 4).Compare("INFO") != 0)) {
		//		pVcf->asMetaKeys.InsertAt(i, "INFO");
		//		pVcf->asMetaValues.InsertAt(i, "<ID=AVGPOST,Number=1,Type=Float,Description=\"Average Posterior Probability from thunderVCF\">");
		//		++i;

		//		pVcf->asMetaKeys.InsertAt(i, "INFO");
		//		pVcf->asMetaValues.InsertAt(i, "<ID=RSQ,Number=1,Type=Float,Description=\"Imputation Quality from thunderVCF\">");
		//		++i;

		//		pVcf->asMetaKeys.InsertAt(i, "INFO");
		//		pVcf->asMetaValues.InsertAt(i, "<ID=ERATE,Number=1,Type=Float,Description=\"Per-marker error rate from thunderVCF\">");
		//		++i;

		//		pVcf->asMetaKeys.InsertAt(i, "INFO");
		//		pVcf->asMetaValues.InsertAt(i, "<ID=THETA,Number=1,Type=Float,Description=\"Recombination parameter with next marker from thunderVCF\">");
		//		++i;
		//	}
		//	if ((pVcf->asMetaKeys[i - 1].SubStr(0, 6).Compare("FORMAT") != 0) && (pVcf->asMetaKeys[i].SubStr(0, 6).Compare("FORMAT") == 0)) {
		//		pVcf->asMetaKeys.InsertAt(i + 1, "FORMAT");
		//		pVcf->asMetaValues.InsertAt(i + 1, "<ID=DS,Number=1,Type=Integer,Description=\"Genotype dosage from thunderVCF\">");
		//		++i;
		//	}
		//}

		

		char** haplotypes = consensus.consensus;

		// check the sanity of data
		if (pVcf->getSampleCount() == 0) {
			throw VcfFileException("No individual genotype information exist in the input VCF file %s", filename.c_str());
		}

		int nSamples = pVcf->getSampleCount();

		// build map of personID -> sampleIndex
		std::map<std::string, int> pedMap;
		for (int i = 0; i < (engine.individuals-engine.phased)/*ped.count*/; ++i) {
			pedMap[ped[i].pid.c_str()] = i;
			//fprintf(stderr,"Adding (%s,%d)\n", ped[i].pid.c_str(), i);
		}

		std::vector<int> vcf2ped;
		std::vector<int> outputSubset;
		for (int i = 0; i < nSamples; ++i) {
			//if (pidIncludedInUnphasedVcf.size() == 0 || pidIncludedInUnphasedVcf.find(std::string(pVcf->vpVcfInds[i]->sIndID.c_str())) != pidIncludedInUnphasedVcf.end()){//in list
				std::map<std::string, int>::iterator found = pedMap.find(pVcf->vpVcfInds[i]->sIndID.c_str());

				if (found == pedMap.end()) {
					//error("Cannot find individual ID %s", pVcf->vpVcfInds[i]->sIndID.c_str());
					//exit(-1);
					continue;
				}
				else {
					fprintf(stderr, "Found (%s,%d)\n", pVcf->vpVcfInds[i]->sIndID.c_str(), found->second);
					vcf2ped.push_back(found->second);
					outputSubset.push_back(i);
				}
				
			//}
		}
		pVcf->printVCFHeaderSubset(outVCF, outputSubset); // print header file

		// read VCF lines
		VcfMarker* pMarker = new VcfMarker;

		char sDose[255];
		double freq, maf, avgPost, rsq;
        String markerName;
		for (int m = 0; pVcf->iterateMarker(); ++m) {
			//fprintf(stderr,"m=%d\n",m);

			pMarker = pVcf->getLastMarker();
			markerName.printf("%s:%d", pMarker->sChrom.c_str(), pMarker->nPos);

			if(unphaseMarkerFlag[std::string(markerName.c_str())]==false) {//if not unphase unique marker, update content otherwise remains as before
				doses.CalculateMarkerInfo(m, freq, maf, avgPost, rsq);

				////fprintf(stderr,"foo1\n");
				int nInfo = pMarker->asInfoKeys.Find("LDAF");
				if (nInfo < -1) {
					sprintf(sDose, "%.4lf", 1. - freq);
					pMarker->asInfoKeys.Add("LDAF");
					pMarker->asInfoValues.Add(sDose);
				}
				else
					pMarker->asInfoValues[nInfo].printf("%.4lf", 1. - freq);

				nInfo = pMarker->asInfoKeys.Find("AVGPOST");
				if (nInfo < -1) {
					sprintf(sDose, "%.4lf", avgPost);
					pMarker->asInfoKeys.Add("AVGPOST");
					pMarker->asInfoValues.Add(sDose);
				}
				else
					pMarker->asInfoValues[nInfo].printf("%.4lf", avgPost);

				nInfo = pMarker->asInfoKeys.Find("RSQ");
				if (nInfo < -1) {
					sprintf(sDose, "%.4lf", rsq);
					pMarker->asInfoKeys.Add("RSQ");
					pMarker->asInfoValues.Add(sDose);
				}
				else
					pMarker->asInfoValues[nInfo].printf("%.4lf", rsq);

				nInfo = pMarker->asInfoKeys.Find("ERATE");
				if (nInfo < -1) {
					sprintf(sDose, "%.4lf", nerror_rates ? error_rates[m] / nerror_rates : 0);
					pMarker->asInfoKeys.Add("ERATE");
					pMarker->asInfoValues.Add(sDose);
				}
				else
					pMarker->asInfoValues[nInfo].printf("%.4lf", nerror_rates ? error_rates[m] / nerror_rates : 0);

				nInfo = pMarker->asInfoKeys.Find("THETA");
				if (nInfo < -1) {
					if (m != engine.markers - 1)
						sprintf(sDose, "%.4lf", nthetas ? thetas[m] / nthetas : 0);
					else
						sprintf(sDose, "%.4lf", nthetas ? thetas[m - 1] / nthetas : 0);
					pMarker->asInfoKeys.Add("THETA");
					pMarker->asInfoValues.Add(sDose);
				}
				else {
					if (m != engine.markers - 1)
						pMarker->asInfoValues[nInfo].printf("%.4lf", nthetas ? thetas[m] / nthetas : 0);
					else
						pMarker->asInfoValues[nInfo].printf("%.4lf", nthetas ? thetas[m - 1] / nthetas : 0);
				}

				int GTidx = pMarker->asFormatKeys.Find("GT");
				if (GTidx < -1) {
					throw VcfFileException("Cannot recognize GT key in FORMAT field");
				}

				int DSidx = GTidx + 1;
				int PLidx = GTidx + 2;

				DSidx = pMarker->asFormatKeys.Find("DS");
				if (DSidx < -1) {
					//throw VcfFileException("Cannot recognize DS key in FORMAT field");
					pMarker->asFormatKeys.InsertAt(GTidx + 1, "DS");
				}

				int nFormats = pMarker->asFormatKeys.Length();

				pMarker->setSampleSize(vcf2ped.size(), pVcf->bParseGenotypes, pVcf->bParseDosages, pVcf->bParseValues);

				//fprintf(stderr,"nFormats=%d\tGTidx=%d\tDSidx=%d\n",nFormats,GTidx,DSidx);

				for (int i = 0; i < vcf2ped.size(); ++i) {
					int pi = vcf2ped[i];
					// modify GT values;
					//fprintf(stderr,"i=%d, pi=%d, GTidx = %d, nFormats = %d, asSampleValues.Length() = %d, haplotypes = %x\n",i,pi,GTidx,nFormats,pMarker->asSampleValues.Length(), haplotypes);
					if (pMarker->asAlts.Length() == 1) {
						pMarker->asSampleValues[nFormats*i + GTidx].printf("%d|%d", haplotypes[pi * 2][m], haplotypes[pi * 2 + 1][m]);
					}
					else {
						pMarker->asSampleValues[nFormats*i + GTidx].printf("%d|%d", haplotypes[pi * 2][m] + 1, haplotypes[pi * 2 + 1][m] + 1);
					}
					// add DS values
					sprintf(sDose, "%.3lf", 2 - doses.GetDosage(pi, m));
					if (DSidx < -1)
					{
					    pMarker->asSampleValues.InsertAt(nFormats*i + DSidx, sDose);
					}
					else
					    pMarker->asSampleValues[nFormats * i + DSidx].printf("%.3lf", sDose);
				}
			}
			pMarker->printVCFMarker(outVCF,false); // print marker to output file
		}
		delete pVcf;
		//delete pMarker;
		ifclose(outVCF);
	}
	catch (VcfFileException e) {
		error(e.what());
	}
}
void UnphasedSamplesOutputVCF(const String& inVcf, Pedigree & ped, DosageCalculator &doses, const String& filename, float* thetas, float* error_rates, ShotgunHaplotyper &engine)
{

	// read and write VCF inputs
	try {

		fprintf(stderr, "Outputing VCF file %s\n", filename.c_str());

		VcfFile* pVcf = new VcfFile;
		IFILE outVCF = ifopen(filename.c_str(), "wb");
		if (outVCF == NULL) {
			error("Cannot open output file %s for writing", filename.c_str());
			exit(-1);
		}

		pVcf->bSiteOnly = false;
		pVcf->bParseGenotypes = false;
		pVcf->bParseDosages = false;
		pVcf->bParseValues = true;
		pVcf->openForRead(inVcf.c_str());

		// add proper header information
		//for (int i = 1; i < pVcf->asMetaKeys.Length(); ++i) {
		//	if ((pVcf->asMetaKeys[i - 1].SubStr(0, 4).Compare("INFO") == 0) && (pVcf->asMetaKeys[i].SubStr(0, 4).Compare("INFO") != 0)) {
		//		pVcf->asMetaKeys.InsertAt(i, "INFO");
		//		pVcf->asMetaValues.InsertAt(i, "<ID=AVGPOST,Number=1,Type=Float,Description=\"Average Posterior Probability from thunderVCF\">");
		//		++i;

		//		pVcf->asMetaKeys.InsertAt(i, "INFO");
		//		pVcf->asMetaValues.InsertAt(i, "<ID=RSQ,Number=1,Type=Float,Description=\"Imputation Quality from thunderVCF\">");
		//		++i;

		//		pVcf->asMetaKeys.InsertAt(i, "INFO");
		//		pVcf->asMetaValues.InsertAt(i, "<ID=ERATE,Number=1,Type=Float,Description=\"Per-marker error rate from thunderVCF\">");
		//		++i;

		//		pVcf->asMetaKeys.InsertAt(i, "INFO");
		//		pVcf->asMetaValues.InsertAt(i, "<ID=THETA,Number=1,Type=Float,Description=\"Recombination parameter with next marker from thunderVCF\">");
		//		++i;
		//	}
		//	if ((pVcf->asMetaKeys[i - 1].SubStr(0, 6).Compare("FORMAT") != 0) && (pVcf->asMetaKeys[i].SubStr(0, 6).Compare("FORMAT") == 0)) {
		//		pVcf->asMetaKeys.InsertAt(i + 1, "FORMAT");
		//		pVcf->asMetaValues.InsertAt(i + 1, "<ID=DS,Number=1,Type=Integer,Description=\"Genotype dosage from thunderVCF\">");
		//		++i;
		//	}
		//}

		// check the sanity of data
		if (pVcf->getSampleCount() == 0) {
			throw VcfFileException("No individual genotype information exist in the input VCF file %s", filename.c_str());
		}

		int nSamples = pVcf->getSampleCount();//num of samples need to be phased
		//fprintf(stderr, "we got %d samples\n", nSamples);
		// build map of personID -> sampleIndex
		std::map<std::string, int> pedMap;
		for (int i = 0; i < (engine.individuals-engine.phased)/*ped.count*/; ++i) {
			pedMap[ped[i].pid.c_str()] = i;
			//fprintf(stderr,"Adding (%s,%d)\n", ped[i].pid.c_str(), i);
		}

		std::vector<int> vcf2ped;
		std::vector<int> outputSubset;
		for (int i = 0; i < nSamples; ++i) {
				std::map<std::string, int>::iterator found = pedMap.find(pVcf->vpVcfInds[i]->sIndID.c_str());
				if (found == pedMap.end()) {
					//error("Cannot find individual ID %s", pVcf->vpVcfInds[i]->sIndID.c_str());
					//exit(-1);
					continue;
				}
				else {
					fprintf(stderr, "Found (%s,%d)\n", pVcf->vpVcfInds[i]->sIndID.c_str(), found->second);
					vcf2ped.push_back(found->second);
					outputSubset.push_back(i);
				}
		}
		pVcf->printVCFHeaderSubset(outVCF, outputSubset); // print header file
		// read VCF lines
		VcfMarker* pMarker = new VcfMarker;
		char sDose[255];
		double freq, maf, avgPost, rsq;
		String markerName;

		for (int m = 0; pVcf->iterateMarker(); ++m) {
			pMarker = pVcf->getLastMarker();
			markerName.printf("%s:%d", pMarker->sChrom.c_str(), pMarker->nPos);

			if(unphaseMarkerFlag[std::string(markerName.c_str())]==false) {//if not unphase unique marker, update content otherwise remains as before
				doses.CalculateMarkerInfo(m, freq, maf, avgPost, rsq);

				////fprintf(stderr,"foo1\n");
				int nInfo = pMarker->asInfoKeys.Find("LDAF");
				if (nInfo < -1) {
					sprintf(sDose, "%.4lf", 1. - freq);
					pMarker->asInfoKeys.Add("LDAF");
					pMarker->asInfoValues.Add(sDose);
				}
				else
					pMarker->asInfoValues[nInfo].printf("%.4lf", 1. - freq);

				nInfo = pMarker->asInfoKeys.Find("AVGPOST");
				if (nInfo < -1) {
					sprintf(sDose, "%.4lf", avgPost);
					pMarker->asInfoKeys.Add("AVGPOST");
					pMarker->asInfoValues.Add(sDose);
				}
				else
					pMarker->asInfoValues[nInfo].printf("%.4lf", avgPost);

				nInfo = pMarker->asInfoKeys.Find("RSQ");
				if (nInfo < -1) {
					sprintf(sDose, "%.4lf", rsq);
					pMarker->asInfoKeys.Add("RSQ");
					pMarker->asInfoValues.Add(sDose);
				}
				else
					pMarker->asInfoValues[nInfo].printf("%.4lf", rsq);

				nInfo = pMarker->asInfoKeys.Find("ERATE");
				if (nInfo < -1) {
					sprintf(sDose, "%.4lf", nerror_rates ? error_rates[m] / nerror_rates : 0);
					pMarker->asInfoKeys.Add("ERATE");
					pMarker->asInfoValues.Add(sDose);
				}
				else
					pMarker->asInfoValues[nInfo].printf("%.4lf", nerror_rates ? error_rates[m] / nerror_rates : 0);

				nInfo = pMarker->asInfoKeys.Find("THETA");
				if (nInfo < -1) {
					if (m != engine.markers - 1)
						sprintf(sDose, "%.4lf", nthetas ? thetas[m] / nthetas : 0);
					else
						sprintf(sDose, "%.4lf", nthetas ? thetas[m - 1] / nthetas : 0);
					pMarker->asInfoKeys.Add("THETA");
					pMarker->asInfoValues.Add(sDose);
				}
				else {
					if (m != engine.markers - 1)
						pMarker->asInfoValues[nInfo].printf("%.4lf", nthetas ? thetas[m] / nthetas : 0);
					else
						pMarker->asInfoValues[nInfo].printf("%.4lf", nthetas ? thetas[m - 1] / nthetas : 0);
				}

				int GTidx = pMarker->asFormatKeys.Find("GT");
				if (GTidx < -1) {
					throw VcfFileException("Cannot recognize GT key in FORMAT field");
				}

				int DSidx = GTidx + 1;
				int PLidx = GTidx + 2;

				DSidx = pMarker->asFormatKeys.Find("DS");
				if (DSidx < -1) {
					//throw VcfFileException("Cannot recognize DS key in FORMAT field");
					pMarker->asFormatKeys.InsertAt(GTidx + 1, "DS");
				}

				int nFormats = pMarker->asFormatKeys.Length();

				pMarker->setSampleSize(vcf2ped.size(), pVcf->bParseGenotypes, pVcf->bParseDosages, pVcf->bParseValues);

				//fprintf(stderr,"nFormats=%d\tGTidx=%d\tDSidx=%d\n",nFormats,GTidx,DSidx);

				for (int i = 0; i < vcf2ped.size(); ++i) {
					int pi = vcf2ped[i];
					// modify GT values;
					//fprintf(stderr,"i=%d, tok=%d, pi=%d, GTidx = %d, nFormats = %d, asSampleValues.Length() = %d, haplotypes = %x\n",i, tok,pi,GTidx,nFormats,pMarker->asSampleValues.Length(), engine.haplotypes);
					if (pMarker->asAlts.Length() == 1) {
						pMarker->asSampleValues[nFormats * i + GTidx].printf("%d|%d", engine.haplotypes[pi * 2][m],
																			 engine.haplotypes[pi * 2 + 1][m]);
					}
					else {
						pMarker->asSampleValues[nFormats * i + GTidx].printf("%d|%d", engine.haplotypes[pi * 2][m] + 1,
																			 engine.haplotypes[pi * 2 + 1][m] + 1);
					}
					// add DS values
					sprintf(sDose, "%.3lf", 2 - doses.GetDosage(pi, m));
					if (DSidx < -1)
					{
					    pMarker->asSampleValues.InsertAt(nFormats*i + DSidx, sDose);
					}
					else
					    pMarker->asSampleValues[nFormats * i + DSidx].printf("%.3lf", sDose);
				}
			}
			pMarker->printVCFMarker(outVCF,false); // print marker to output file
		}
		delete pVcf;
		ifclose(outVCF);
	}
	catch (VcfFileException e) {
		error(e.what());
	}
}

void UpdateVector(float * current, float * & vector, int & n, int length)
{
	if (n++ == 0)
	{
		vector = new float[length];

		for (int i = 0; i < length; i++)
			vector[i] = current[i];
	}
	else
		for (int i = 0; i < length; i++)
			vector[i] += current[i];
}

void UpdateErrorRates(Errors * current, float * & vector, int & n, int length)
{
	if (n++ == 0)
	{
		vector = new float[length];

		for (int i = 0; i < length; i++)
			vector[i] = current[i].rate;
	}
	else
		for (int i = 0; i < length; i++)
			vector[i] += current[i].rate;
}

int MemoryAllocationFailure()
{
	printf("FATAL ERROR - Memory allocation failed\n");
	return -1;
}
void LoadPidToBeIncluded(String & filename, String & filename2)
{
	if (filename != "")
	{
		std::ifstream fin(filename.c_str());
		if (!fin.is_open()) {
			std::cerr << "Open file " << filename << " failed!" << std::endl; exit(EXIT_FAILURE);
		}
		std::string tmpLine;
		while (getline(fin, tmpLine))
		{
			pidIncludedInUnphasedVcf[tmpLine] = true;
		}
		fin.close();
	}
	if (filename2 != "")
	{
		std::ifstream fin2(filename2.c_str());
		if (!fin2.is_open()) {
			std::cerr << "Open file " << filename2 << " failed!" << std::endl; exit(EXIT_FAILURE);
		}
		std::string tmpLine;
		while (getline(fin2, tmpLine))
		{
			pidIncludedInPhasedVcf[tmpLine] = true;
		}
		fin2.close();
	}
}

void LoadPidToBeExcluded(String & filename, String & filename2)
{
	if (filename != "")
	{
		std::ifstream fin(filename.c_str());
		if (!fin.is_open()) {
			std::cerr << "Open file " << filename << " failed!" << std::endl; exit(EXIT_FAILURE);
		}
		std::string tmpLine;
		while (getline(fin, tmpLine))
		{
			pidExcludedInUnphasedVcf[tmpLine] = true;
		}
		fin.close();
	}
	if (filename2 != "")
	{
		std::ifstream fin2(filename2.c_str());
		if (!fin2.is_open()) {
			std::cerr << "Open file " << filename2 << " failed!" << std::endl; exit(EXIT_FAILURE);
		}
		std::string tmpLine;
		while (getline(fin2, tmpLine))
		{
			pidExcludedInPhasedVcf[tmpLine] = true;
		}
		fin2.close();
	}
}
void LoadShotgunSamples(Pedigree &ped, const String& filename, std::unordered_map<std::string, bool>& pidIncluded, std::unordered_map<std::string, bool>& pidExcluded, int& num) {
	//printf("starting LoadShotgunSamples\n\n");

	try {

		VcfFile* pVcf = new VcfFile;
		pVcf->bSiteOnly = true;
		pVcf->bParseGenotypes = false;
		pVcf->bParseDosages = false;
		pVcf->bParseValues = false;
		pVcf->openForRead(filename.c_str());

		int PreviousCount = ped.count;
		// check the sanity of data
		if (pVcf->getSampleCount() == 0) {
			throw VcfFileException("No individual genotype information exist in the input VCF file %s", filename.c_str());
		}
		//std::cerr << "Sample Size:" << pVcf->getSampleCount() << std::endl;
		//std::cerr << "Include Size:" << pidIncluded.size() << std::endl;
		//std::cerr << "Exlude Size:" << pidExcluded.size() << std::endl;
		for (int i = 0; i < pVcf->getSampleCount(); ++i) {
			//std::cerr << "input:" << std::string(pVcf->vpVcfInds[i]->sIndID.c_str()) << "\t" << pidIncluded[std::string(pVcf->vpVcfInds[i]->sIndID.c_str())] << std::endl;
			if ((pidIncluded.size() == 0 || pidIncluded.find(std::string(pVcf->vpVcfInds[i]->sIndID.c_str())) != pidIncluded.end()) && (pidExcluded.size()==0||pidExcluded.find(std::string(pVcf->vpVcfInds[i]->sIndID.c_str())) == pidExcluded.end()))
			{
				if (DuplicatedIndividualPair.find(std::string(pVcf->vpVcfInds[i]->sIndID.c_str())) == DuplicatedIndividualPair.end())//never showed before
				{
					ped.AddPerson(pVcf->vpVcfInds[i]->sIndID, pVcf->vpVcfInds[i]->sIndID, "0", "0", 1, 1);
					DuplicatedIndividualPair[std::string(pVcf->vpVcfInds[i]->sIndID.c_str())] = std::make_pair(i,-1);
					//std::cerr << "Never saw" << pVcf->vpVcfInds[i]->sIndID << std::endl;
				}
				else
				{
					warning("Individual %s duplicated!", pVcf->vpVcfInds[i]->sIndID.c_str());
					ped.AddPerson(pVcf->vpVcfInds[i]->sIndID, pVcf->vpVcfInds[i]->sIndID+String("_dup"), "0", "0", 1, 1);
					DuplicatedIndividualPair[std::string(pVcf->vpVcfInds[i]->sIndID.c_str())].second = i + PreviousCount;
				}
				num++;
			}

		}

		delete pVcf;
	}
	catch (VcfFileException e) {
		error(e.what());
	}

	//ped.Sort();
	printf("Loaded %d individuals from file %s\n\n", ped.count, filename.c_str());
}

void LoadPolymorphicSites(const String& filename) {
	try {
		VcfFile* pVcf = new VcfFile;
		pVcf->bSiteOnly = true;
		pVcf->bParseGenotypes = false;
		pVcf->bParseDosages = false;
		pVcf->bParseValues = false;
		pVcf->openForRead(filename.c_str());

		VcfMarker* pMarker = new VcfMarker;

		StringArray altalleles;
		String markerName;

		while (pVcf->iterateMarker()) {
			int markers = Pedigree::markerCount;
			pMarker = pVcf->getLastMarker();

			markerName.printf("%s:%d", pMarker->sChrom.c_str(), pMarker->nPos);
			int marker = Pedigree::GetMarkerID(markerName);// that's where we fill up marker numbers for ped file
			//initialize flag map to allocate memory equals to marker number in ref set
			unphaseMarkerFlag[std::string(markerName.c_str())] = false;
			unphaseMarkerUIdx[std::string(markerName.c_str())] = marker;
			unphaseMarkerIdx[std::string(markerName.c_str())] = -1;

			int al1, al2;

			//printf("Re-opening VCF file\n");

			if (pMarker->asAlts.Length() == 2) {
				al1 = Pedigree::LoadAllele(marker, pMarker->asAlts[0]);
				al2 = Pedigree::LoadAllele(marker, pMarker->asAlts[1]);
			}
			else {
				al1 = Pedigree::LoadAllele(marker, pMarker->sRef);
				al2 = Pedigree::LoadAllele(marker, pMarker->asAlts[0]);
			}


			if (markers != marker) {
				error("Each polymorphic site should only occur once, but site %s is duplicated\n", markerName.c_str());
			}

			if (al1 != 1 || al2 != 2) {
				error("Allele labels '%s' and '%s' for polymorphic site '%s' are not valid\n", (const char *)altalleles[0], (const char *)altalleles[1], markerName.c_str());
			}
		}
		delete pVcf;
		//delete pMarker;
	}
	catch (VcfFileException e) {
		error(e.what());
	}
}
void LoadUnphasedPolymorphicSites(const String& filename) {
	try {
		VcfFile* pVcf = new VcfFile;
		pVcf->bSiteOnly = true;
		pVcf->bParseGenotypes = false;
		pVcf->bParseDosages = false;
		pVcf->bParseValues = false;
		pVcf->openForRead(filename.c_str());

		VcfMarker* pMarker = new VcfMarker;
		StringArray altalleles;
		String markerName;
		int localIdx = 0;


		while (pVcf->iterateMarker()) {
			
			pMarker = pVcf->getLastMarker();
			markerName.printf("%s:%d", pMarker->sChrom.c_str(), pMarker->nPos);
			int idx =  Pedigree::markerLookup.Integer(markerName);//only look up, no add in
			if (idx != -1)//shown in ref panel marker set
			{
				unphaseMarkerFlag[std::string(markerName.c_str())] = true;
				unphaseMarkerIdx[std::string(markerName.c_str())] = localIdx;
			}
			else
			{
				unphaseMarkerIdx[std::string(markerName.c_str())] = localIdx;
			}
			++localIdx;
			//int marker = Pedigree::GetMarkerID(markerName);
			//int al1, al2;

			////printf("Re-opening VCF file\n");

			//if (pMarker->asAlts.Length() == 2) {
			//	al1 = Pedigree::LoadAllele(marker, pMarker->asAlts[0]);
			//	al2 = Pedigree::LoadAllele(marker, pMarker->asAlts[1]);
			//}
			//else {
			//	al1 = Pedigree::LoadAllele(marker, pMarker->sRef);
			//	al2 = Pedigree::LoadAllele(marker, pMarker->asAlts[0]);
			//}

			////printf("Re-opening VCF file\n");

			//if (markers != marker) {
			//	error("Each polymorphic site should only occur once, but site %s is duplicated\n", markerName.c_str());
			//}


			//if (al1 != 1 || al2 != 2) {
			//	error("Allele labels '%s' and '%s' for polymorphic site '%s' are not valid\n", (const char *)altalleles[0], (const char *)altalleles[1], markerName.c_str());
			//}
		}
		delete pVcf;
		//delete pMarker;
	}
	catch (VcfFileException e) {
		error(e.what());
	}
}
void LoadShotgunResults(Pedigree &ped, char** genotypes, /*char* refalleles, double* freq1s,*/ const String & filename, int maxPhred, ShotgunHaplotyper&engine) {
	//printf("starting LoadShotgunResults\n\n");
	try {
		VcfFile* pVcf = new VcfFile;
		pVcf->bSiteOnly = false;
		pVcf->bParseGenotypes = false;
		pVcf->bParseDosages = false;
		pVcf->bParseValues = true;
		pVcf->openForRead(filename.c_str());

		// check the sanity of data
		if (pVcf->getSampleCount() == 0) {
			throw VcfFileException("No individual genotype information exist in the input VCF file %s", filename.c_str());
		}

		int nSamples = pVcf->getSampleCount();
		//vector<int> personIndices(ped.count,-1);
		std::unordered_map<int, int> personIndices;
		StringIntHash originalPeople; // key: famid+subID, value: original order (0 based);
		int person = 0;
		for (int i = 0; i < nSamples; i++) {
				originalPeople.Add(pVcf->vpVcfInds[i]->sIndID + "." + pVcf->vpVcfInds[i]->sIndID, person);
				person++;
		}

		for (int i = 0; i < (engine.individuals-engine.phased);/* ped.count;*/ i++) {
			if (originalPeople.Integer(ped[i].famid + "." + ped[i].pid) != -1)
			{
				personIndices[originalPeople.Integer(ped[i].famid + "." + ped[i].pid)] = i;
			}

		}

		int markerindex = 0;
		VcfMarker* pMarker = new VcfMarker;
		String markerName;
		while (pVcf->iterateMarker()) {//for each marker

			pMarker = pVcf->getLastMarker();
			markerName.printf("%s:%d", pMarker->sChrom.c_str(), pMarker->nPos);
			//printf("now for marker %s:%d\t", pMarker->sChrom.c_str(), pMarker->nPos);
			if (unphaseMarkerFlag.find(std::string(markerName.c_str()))!=unphaseMarkerFlag.end()&&unphaseMarkerFlag[std::string(markerName.c_str())] == true)
				markerindex = unphaseMarkerUIdx[std::string(markerName.c_str())];
			else
				continue;
			//int AFidx = pMarker->asInfoKeys.Find("AF");
			int PLidx = pMarker->asFormatKeys.Find("PL");
			int GLflag = 0;
			if (PLidx < 0) {
				PLidx = pMarker->asFormatKeys.Find("GL");
				if (PLidx >= 0) GLflag = 1;
				else throw VcfFileException("Cannot recognize GL or PL key in FORMAT field");
			}
			//printf("reading vcf 1\n\n");
			int formatLength = pMarker->asFormatKeys.Length();
			int idx11 = 0, idx12 = 1, idx22 = 2;

			StringArray phred;
			int genoindex = markerindex * 3;

			for (int i = 0; i < nSamples; i++)//for each individual
			{
				//printf("%s\t", pMarker->asSampleValues[PLidx + i*formatLength].c_str());
				if (personIndices.find(i)!= personIndices.end()){
					phred.ReplaceTokens(pMarker->asSampleValues[PLidx + i*formatLength], ",");

					int phred11 = GLflag ? static_cast<int>(-10. * phred[idx11].AsDouble()) : phred[idx11].AsInteger();
					int phred12 = GLflag ? static_cast<int>(-10. * phred[idx12].AsDouble()) : phred[idx12].AsInteger();
					int phred22 = GLflag ? static_cast<int>(-10. * phred[idx22].AsDouble()) : phred[idx22].AsInteger();

					if ((phred11 < 0) || (phred11 < 0) || (phred12 < 0)) {
						error("Negative PL or Positive GL observed");
					}

					//printf("phred scores are %f, %f, %f\n", phred[idx11].AsDouble(), phred[idx12].AsDouble(), phred[idx22].AsDouble());

					if (phred11 > maxPhred) phred11 = maxPhred;
					if (phred12 > maxPhred) phred12 = maxPhred;
					if (phred22 > maxPhred) phred22 = maxPhred;

					genotypes[personIndices[i]][genoindex] = phred11;
					genotypes[personIndices[i]][genoindex + 1] = phred12;
					genotypes[personIndices[i]][genoindex + 2] = phred22;
				}
			}
		}

		delete pVcf;
		//delete pMarker;
	}
	catch (VcfFileException e) {
		error(e.what());
	}
}
void LoadGenotypeFromPhasedVcf(Pedigree &ped, char** genotypes, char* refalleles, double* freq1s, const String & filename, int maxPhred, int phased, ShotgunHaplotyper& engine, double defaultErrorRate, double defaultTransRate) {
	//printf("starting LoadPhasedVcf\n\n");

	try {
		VcfFile* pVcf = new VcfFile;
		pVcf->bSiteOnly = false;
		pVcf->bParseGenotypes = false;
		pVcf->bParseDosages = false;
		pVcf->bParseValues = true;
		pVcf->openForRead(filename.c_str());

		// check the sanity of data
		if (pVcf->getSampleCount() == 0) {
			throw VcfFileException("No individual genotype information exist in the input VCF file %s", filename.c_str());
		}
		vector<int> phaseIdx(ped.count, -1);//initially assume all the individuals are phased
		int nSamples = pVcf->getSampleCount();

		//vector<int> personIndices(ped.count, -1);
		std::unordered_map<int, int> personIndices;
		StringIntHash originalPeople; // key: famid+subID, value: original order (0 based);
		int person = 0;
		for (int i = 0; i < nSamples; i++) {
			{
				//std::cerr << "if ordered:" << pVcf->vpVcfInds[i]->sIndID << std::endl;
				originalPeople.Add(pVcf->vpVcfInds[i]->sIndID + "." + pVcf->vpVcfInds[i]->sIndID, person);
				person++;
			}
		}

		for (int i = (engine.individuals-engine.phased); i < engine.individuals; i++) {
			int idx = originalPeople.Integer(ped[i].famid + "." + ped[i].pid);
			if (idx != -1)//phased in this vcf
			{
				personIndices[originalPeople.Integer(ped[i].famid + "." + ped[i].pid)] = i;

			}
		}

		int markerindex = 0;
		VcfMarker* pMarker = new VcfMarker;
		String markerName;
		while (pVcf->iterateMarker()) {//for each marker
			refalleles[markerindex] = 1;
			//printf("reading vcf 1\n\n");
			pMarker = pVcf->getLastMarker();
			markerName.printf("%s:%d", pMarker->sChrom.c_str(), pMarker->nPos);
			//fprintf(stderr,"now adding marker:%s\n",markerName.c_str());
			int AFidx = pMarker->asInfoKeys.Find("AF");
			if (AFidx == -1)
			{
				pMarker->asInfoKeys.PrintLine();
				pMarker->asInfoValues.PrintLine();
			}
			int PLidx = pMarker->asFormatKeys.Find("PL");
			/*setting error rate*/
			int ERidx = pMarker->asInfoKeys.Find("ERATE");
			int THidx = pMarker->asInfoKeys.Find("THETA");
			if (ERidx != -1 && THidx != -1)// ERATE, THETA exist
			{
				//error_rates[markerindex] = pMarker->asInfoValues[THidx].AsDouble();
				engine.SetErrorRate(markerindex, pMarker->asInfoValues[THidx].AsDouble());
				if (markerindex!=engine.markers-1)
				engine.thetas[markerindex] = pMarker->asInfoValues[THidx].AsDouble();
			}
			else
			{
				fprintf(stderr, "No ERATE or THETA tag found in input vcf, now using command line(--errorRate and --transRate)settings:\n Error Rate:%f\tTrans Rate(Theta):%f\n",defaultErrorRate,defaultTransRate);
				engine.SetErrorRate(defaultErrorRate);
				engine.thetas[markerindex] = defaultTransRate;
			}
			int GLflag = 0;
			if (PLidx < 0) {
				PLidx = pMarker->asFormatKeys.Find("GL");
				if (PLidx >= 0) GLflag = 1;
				else throw VcfFileException("Cannot recognize GL or PL key in FORMAT field");
			}
			//printf("reading vcf 2\n\n");
			int formatLength = pMarker->asFormatKeys.Length();
			int idx11 = 0, idx12 = 1, idx22 = 2;
			if (AFidx == -1) {
				int ANidx = pMarker->asInfoKeys.Find("AN");
				int ACidx = pMarker->asInfoKeys.Find("AC");
				if ((ANidx < 0) || (ACidx < 0)) {
					throw VcfFileException("Cannot recognize AF key in FORMAT field");
				}
				else {
					freq1s[markerindex] = 1. - (pMarker->asInfoValues[ACidx].AsDouble() + .5) / (pMarker->asInfoValues[ANidx].AsDouble() + 1.);
				}
			}
			else if (pMarker->asAlts.Length() == 1) {
				freq1s[markerindex] = (1. - pMarker->asInfoValues[AFidx].AsDouble());

				if (PLidx < 0) {
					error("Missing PL key in FORMAT field");
				}
			}
			else {
				// AF1,AF2 -- freq1s is AF1
				freq1s[markerindex] = pMarker->asInfoValues[AFidx].AsDouble();

				if (PLidx < 0) {
					PLidx = pMarker->asFormatKeys.Find("PL3");
					if (PLidx < 0) {
						PLidx = pMarker->asFormatKeys.Find("GL3");
						if (PLidx >= 0) GLflag = 1;
					}
					idx11 = 2;
					idx12 = 4;
					idx22 = 5;
				}
				if (PLidx < 0) {
					error("Missing PL key in FORMAT field");
				}
			}
			//printf("reading vcf 3\n\n");
			StringArray phred;
			int genoindex = markerindex * 3;

			for (int i = 0; i < nSamples; i++)//for each phased individual
			{
				//printf("phred scores are   %d, %d\n", i, ped.count);
				if (personIndices.find(i) != personIndices.end()){//pidIncludedInPhasedVcf.size()==0||pidIncludedInPhasedVcf.find(std::string(pVcf->vpVcfInds[i]->sIndID.c_str())) != pidIncludedInPhasedVcf.end()){
					phred.ReplaceTokens(pMarker->asSampleValues[PLidx + i*formatLength], ",");
					int phred11 = GLflag ? static_cast<int>(-10. * phred[idx11].AsDouble()) : phred[idx11].AsInteger();
					int phred12 = GLflag ? static_cast<int>(-10. * phred[idx12].AsDouble()) : phred[idx12].AsInteger();
					int phred22 = GLflag ? static_cast<int>(-10. * phred[idx22].AsDouble()) : phred[idx22].AsInteger();

					if ((phred11 < 0) || (phred11 < 0) || (phred12 < 0)) {
						error("Negative PL or Positive GL observed");
					}

					//printf("phred scores are %d, %d, %d, %d, %d\n", phred11, phred12, phred22,i, ped.count );

					if (phred11 > maxPhred) phred11 = maxPhred;
					if (phred12 > maxPhred) phred12 = maxPhred;
					if (phred22 > maxPhred) phred22 = maxPhred;

					genotypes[personIndices[i]][genoindex] = phred11;
					genotypes[personIndices[i]][genoindex + 1] = phred12;
					genotypes[personIndices[i]][genoindex + 2] = phred22;
				}
			}
			//printf("reading vcf 4\n\n");

			++markerindex;
		}
		//initialize missing value in unphased file as 0
		for (std::unordered_map<std::string, bool>::const_iterator iter = unphaseMarkerFlag.begin(); iter != unphaseMarkerFlag.end(); ++iter)
		{
			if (iter->second == false)// if this marker not shown in unphased set but in reference set
			{
				markerindex = unphaseMarkerUIdx[iter->first];
				int genoindex = markerindex * 3;
				for (int i = 0; i < (engine.individuals-engine.phased); i++)//for each unphased individual
				{
					//if(phaseIdx[i]==-1) continue;
					//phred.ReplaceTokens(pMarker->asSampleValues[PLidx + i*formatLength], ",");
					int phred11 = 0;// GLflag ? static_cast<int>(-10. * phred[idx11].AsDouble()) : phred[idx11].AsInteger();
					int phred12 = 0;// GLflag ? static_cast<int>(-10. * phred[idx12].AsDouble()) : phred[idx12].AsInteger();
					int phred22 = 0;// GLflag ? static_cast<int>(-10. * phred[idx22].AsDouble()) : phred[idx22].AsInteger();

					if ((phred11 < 0) || (phred11 < 0) || (phred12 < 0)) {
						error("Negative PL or Positive GL observed");
					}

					//printf("phred scores are %d, %d, %d\n", phred11, phred12, phred22);

					if (phred11 > maxPhred) phred11 = maxPhred;
					if (phred12 > maxPhred) phred12 = maxPhred;
					if (phred22 > maxPhred) phred22 = maxPhred;

					genotypes[i][genoindex] = phred11;
					genotypes[i][genoindex + 1] = phred12;
					genotypes[i][genoindex + 2] = phred22;
				}
			}
		}

		delete pVcf;
		//delete pMarker;
	}
	catch (VcfFileException e) {
		error(e.what());
	}
}
int main(int argc, char ** argv)
{
	String shotgunfile, mapfile, outfile("mach1.out"), phasedfile, pidIncludeFromUnphased(""), pidIncludeFromPhased(""), pidExcludeFromUnphased(""), pidExcludeFromPhased("");
	String crossFile, errorFile;
	clock_t t;
	t = clock();
	double errorRate = 0.01;
	double transRate = 0.01;
	int seed = 123456, warmup = 0, states = 0, weightedStates = 0;
	int burnin = 0, rounds = 0, polling = 0, samples = 0, SamplingRounds=0;
	int maxPhred = 255;
	bool compact = false;
	bool mle = false, mledetails = false, uncompressed = false;
	bool inputPhased = false;
	bool phaseByRef = false;
	bool randomPhase = true;
	bool fixTrans = false;

	SetupCrashHandlers();
	SetCrashExplanation("reading command line options");

	printf("Thunder_Glf 1.0.9 -- Markov Chain Haplotyping for Shotgun Sequence Data\n"
		"(c) 2005-2007 Goncalo Abecasis, Yun Li, with thanks to Paul Scheet\n\n");

	ParameterList pl;

	BEGIN_LONG_PARAMETERS(longParameters)
		LONG_PARAMETER_GROUP("Shotgun Sequences")
		LONG_STRINGPARAMETER("unphasedVcf", &shotgunfile)
		LONG_STRINGPARAMETER("refVcf", &phasedfile)
		LONG_INTPARAMETER("maxPhred", &maxPhred)
		LONG_PARAMETER_GROUP("Optional Files")
		LONG_STRINGPARAMETER("includeUnphasedIDs", &pidIncludeFromUnphased)
		LONG_STRINGPARAMETER("includePhasedIDs", &pidIncludeFromPhased)
		LONG_STRINGPARAMETER("excludeUnphasedIDs", &pidExcludeFromUnphased)
		LONG_STRINGPARAMETER("excludePhasedIDs", &pidExcludeFromPhased)
		LONG_STRINGPARAMETER("crossoverMap", &crossFile)
		LONG_STRINGPARAMETER("errorMap", &errorFile)
		LONG_STRINGPARAMETER("physicalMap", &mapfile)
		LONG_PARAMETER_GROUP("Markov Sampler")
		LONG_INTPARAMETER("seed", &seed)
		LONG_INTPARAMETER("burnin", &burnin)
		LONG_INTPARAMETER("rounds", &rounds)
		LONG_INTPARAMETER("SamplingRounds", &SamplingRounds)
		LONG_PARAMETER_GROUP("Haplotyper")
		LONG_INTPARAMETER("states", &states)
		LONG_DOUBLEPARAMETER("errorRate", &errorRate)
		LONG_DOUBLEPARAMETER("transRate", &transRate)
		LONG_INTPARAMETER("weightedStates", &weightedStates)
		LONG_PARAMETER("compact", &compact)
		LONG_PARAMETER("fixTrans", &fixTrans)
		LONG_PARAMETER_GROUP("Phasing")
		EXCLUSIVE_PARAMETER("randomPhase", &randomPhase)
		EXCLUSIVE_PARAMETER("inputPhased", &inputPhased)
		EXCLUSIVE_PARAMETER("refPhased", &phaseByRef)
		LONG_PARAMETER_GROUP("Imputation")
		LONG_PARAMETER("geno", &OutputManager::outputGenotypes)
		LONG_PARAMETER("quality", &OutputManager::outputQuality)
		LONG_PARAMETER("dosage", &OutputManager::outputDosage)
		LONG_PARAMETER("probs", &OutputManager::outputProbabilities)
		LONG_PARAMETER("mle", &mle)
		LONG_PARAMETER_GROUP("Output Files")
		LONG_STRINGPARAMETER("prefix", &outfile)
		LONG_PARAMETER("phase", &OutputManager::outputHaplotypes)
		LONG_PARAMETER("uncompressed", &OutputManager::uncompressed)
		LONG_PARAMETER("mldetails", &mledetails)
		LONG_PARAMETER_GROUP("Interim Output")
		LONG_INTPARAMETER("sampleInterval", &samples)
		LONG_INTPARAMETER("interimInterval", &polling)
		END_LONG_PARAMETERS();

	pl.Add(new LongParameters("Available Options", longParameters));

	pl.Add(new HiddenString('m', "Map File", mapfile));
	pl.Add(new HiddenString('o', "Output File", outfile));
	pl.Add(new HiddenInteger('r', "Haplotyping Rounds", rounds));
	pl.Add(new HiddenDouble('e', "Error Rate", errorRate));

	pl.Read(argc, argv);
	pl.Status();

	if (OutputManager::outputDosage == false) { // hmkang 
		error("--dosage flag must be set in this implementation");
	}

	// Setup random seed ...
	globalRandom.Reset(seed);

	SetCrashExplanation("loading information on polymorphic sites");


	ShotgunHaplotyper engine;//declaration of engine, also will call default constructor

	// Setup and load a list of polymorphic sites, each with two allele labels ...
	Pedigree ped; 

	SetCrashExplanation("loading shotgun data - first pass");
	LoadPidToBeIncluded(pidIncludeFromUnphased, pidIncludeFromPhased);
	LoadPidToBeExcluded(pidExcludeFromUnphased, pidExcludeFromPhased);
	/*We add unphased samples first*/
	int numUnphased(0);
	LoadShotgunSamples(ped, shotgunfile, pidIncludedInUnphasedVcf, pidExcludedInUnphasedVcf, numUnphased);// here shotgunfile is the vcf file, here vcf is used for filling up the first five column of PED file(check the PED format).
	std::cerr << "Load unphased individuals:"<<numUnphased<<std::endl;
	if(ped.count < 1)
    {
       error("SinglePhasing requires more than 0 sample.");
    }

	/*now loading phased individuals*/
	LoadShotgunSamples(ped, phasedfile, pidIncludedInPhasedVcf, pidExcludedInPhasedVcf, engine.phased);// here shotgunfile is the vcf file, here vcf is used for filling up the first five column of PED file(check the PED format).
	std::cerr << "Load phased individuals:" << engine.phased << std::endl;
	/*Notice that now we adding markers as subset of phased markers*/
	LoadPolymorphicSites(phasedfile);// here only extracted site information only, used for site check
	LoadUnphasedPolymorphicSites(shotgunfile);
	SetCrashExplanation("loading map information for polymorphic sites");

	printf("Loaded information on %d polymorphic sites\n\n", Pedigree::markerCount);

	Pedigree::LoadMarkerMap(mapfile);//the format of mapfiles is:	chrome\tmarker_name\tposition
	

	// Check if physical map is available
	bool positionsAvailable = true;

	for (int i = 0; i < ped.markerCount; i++)
		if (Pedigree::GetMarkerInfo(i)->chromosome < 0)
		{
		positionsAvailable = false;//no physical map available
		break;
		}

	if (positionsAvailable)
	{
		printf("    Physical map will be used to improve crossover rate estimates.\n");

		for (int i = 1; i < ped.markerCount; i++)
			if (ped.GetMarkerInfo(i)->position <= ped.GetMarkerInfo(i - 1)->position ||
				ped.GetMarkerInfo(i)->chromosome != ped.GetMarkerInfo(i - 1)->chromosome)
			{
			printf("    FATAL ERROR -- Problems with physical map ...\n\n"
				"    Before continuing, check the following:\n"
				"    * All markers are on the same chromosome\n"
				"    * All marker positions are unique\n"
				"    * Markers in pedigree and haplotype files are ordered by physical position\n\n");
			return -1;
			}
	}

	printf("\n");

	printf("Processing input files and allocating memory for haplotyping\n");

	SetCrashExplanation("allocating memory for haplotype engine and consensus builder");



	engine.economyMode = compact;//

	engine.EstimateMemoryInfo(ped.count, ped.markerCount, states, compact, false);
	engine.AllocateMemory(ped.count, states, ped.markerCount, (float)transRate);

	printf("Copy unphased genotypes into haplotyping engine\n");
	// Copy genotypes into haplotyping engine
	if (engine.readyForUse)
		LoadShotgunResults(ped, engine.genotypes, /*engine.refalleles, engine.freq1s, */shotgunfile, maxPhred, engine);//this is where thunder copy GL into genotype arrays.helpful for understand genotype datastructre.

	printf("Done loading shotgun file\n\n");
	// Copy phased haplotypes into haplotyping engine, but we put phased haps in the end
	
	printf("Copy phased genotypes into haplotyping engine\n");
	if (engine.readyForUse)
		LoadGenotypeFromPhasedVcf(ped, engine.genotypes, engine.refalleles, engine.freq1s, phasedfile, maxPhred, engine.phased, engine, errorRate, transRate);//this is where thunder copy GL into genotype arrays.helpful for understand genotype datastructre.



	if (engine.readyForUse == false || engine.ForceMemoryAllocation() == false)
		return MemoryAllocationFailure();//check error

	if (positionsAvailable && engine.AllocateDistances())//it is interesting to notice that there are two position information sources, one is from VCF the other is from markerMap
	{
		for (int i = 1; i < ped.markerCount; i++)//here the distance is based on markerMap file 
			engine.distances[i - 1] = ped.GetMarkerInfo(i)->position -
			ped.GetMarkerInfo(i - 1)->position;
	}

	engine.ShowMemoryInfo();

	if (mle)//not sure 
	{
		engine.ShowMLEMemoryInfo();
		if (!engine.AllocateMLEMemory())
			return MemoryAllocationFailure();
	}

	ConsensusBuilder::EstimateMemoryInfo(SamplingRounds, ped.count * 2, ped.markerCount);
	ConsensusBuilder consensus(SamplingRounds, ped.count * 2, ped.markerCount);

	if (consensus.readyForUse == false)
		return MemoryAllocationFailure();
	
	DosageCalculator::storeDistribution = OutputManager::outputDosage ||
		OutputManager::outputQuality ||
		OutputManager::outputGenotypes;


	DosageCalculator::EstimateMemoryInfo(1, ped.count, ped.markerCount);
	DosageCalculator doses(1, ped.count, ped.markerCount);

	if (doses.readyForUse == false)
		return MemoryAllocationFailure();

	if (states < weightedStates) {
		error("Total number of states (--states) must be equal or greater than the total number of weighted states (--weightStates)");
	}
	engine.weightedStates = weightedStates;

	printf("Memory allocated successfully\n\n");

	SetCrashExplanation("loading error rate and cross over maps");

	bool newline = engine.LoadCrossoverRates(crossFile);
	newline |= engine.LoadErrorRates(errorFile);
	if (newline) printf("\n");

	//engine.SetErrorRate(errorRate);
	UpdateVector(engine.thetas, thetas, nthetas, engine.markers - 1);
	UpdateErrorRates(engine.error_models, error_rates, nerror_rates, engine.markers);

	SetCrashExplanation("searching for initial haplotype set");

	if (inputPhased) {
		printf("Loading phased information from the input VCF file\n\n");
		engine.LoadHaplotypesFromVCF(shotgunfile);
	}
	else if (phaseByRef) {
		printf("Assigning haplotypes based on reference genome\n\n");
		engine.PhaseByReferenceSetup();
	}
	else {
		printf("Assigning random set of haplotypes\n\n");
		engine.RandomSetup();
	}
	printf("Found initial haplotype set\n\n");
	//OutputManager::WriteHaplotypes(outfile, ped, engine.haplotypes);
	//return 0;
	engine.LoadHaplotypesFromPhasedVCF(ped, phasedfile);// , pidIncludedInPhasedVcf);

	SetCrashExplanation("revving up haplotyping engine");

	SetCrashExplanation("interating through markov chain haplotyping procedure");

	//for (int i = 0; i < rounds; i++)
	//{
	engine.LoopThroughChromosomes(consensus,SamplingRounds,ped, DuplicatedIndividualPair);
	if (!fixTrans) engine.UpdateThetas();
	errorRate = engine.UpdateErrorRate();

	//printf("Markov Chain iteration %d [%d mosaic crossovers]\n",
	//	i + 1, engine.TotalCrossovers());

	//if (i < burnin)
	//	continue;

	//if (OutputManager::outputHaplotypes)
	//	consensus.Store(engine.haplotypes);

	if (doses.storeDosage || doses.storeDistribution)
		doses.Update(engine.haplotypes);

	UpdateVector(engine.thetas, thetas, nthetas, engine.markers - 1);
	UpdateErrorRates(engine.error_models, error_rates, nerror_rates, engine.markers);

	//	//if (polling > 0 && ((i - burnin) % polling) == 0) {
	//	int i = 0;// adjust for following code
	//		OutputVCFConsensus(shotgunfile, ped, consensus, doses, outfile + ".prelim" + (i + 1) + ".vcf.gz", thetas, error_rates);
	//		OutputManager::OutputConsensus(ped, consensus, doses, outfile + ".prelim" + (i + 1));
	//	//}

	//	//if (samples > 0 && ((i - burnin) % samples) == 0)
	//		OutputManager::WriteHaplotypes(outfile + ".sample" + (i + 1) + ".gz", ped, engine.haplotypes);

	////}

	if (rounds) printf("\n");

	SetCrashExplanation("estimating maximum likelihood solution, conditional on current state");

	if (mle)
	{
		// Use best available error and crossover rates for MLE
		if (nerror_rates)
			for (int i = 0; i < engine.markers; i++)
				engine.SetErrorRate(i, error_rates[i] / nerror_rates);

		if (nthetas)
			for (int i = 0; i < engine.markers - 1; i++)
				engine.thetas[i] = thetas[i] / nthetas;

		engine.OutputMLEs(ped, outfile, mledetails);
	}

	//   ParseHaplotypes(engine.haplotypes, engine.individuals * 2 - 2, engine.markers, 32);

	SetCrashExplanation("outputing solution");
	fprintf(stderr, "%d %d\n", ped.count, ped.markerCount);
	// If we did multiple rounds of haplotyping, then generate consensus
	//if (rounds > 1)
	//	OutputVCFConsensus(shotgunfile, ped, consensus, doses, outfile + ".vcf.gz", thetas, error_rates);
	//OutputManager::OutputConsensus(ped, consensus, doses, outfile);
	//else
	//if (OutputManager::outputHaplotypes)
	//	OutputManager::WriteHaplotypes(outfile, ped, engine.haplotypes);
	//else
	{
		UnphasedSamplesOutputVCF(shotgunfile, ped, doses, outfile + ".vcf.gz", thetas, error_rates, engine);
		OutputVCFConsensus(shotgunfile, ped, consensus, doses, outfile + ".consensus.vcf.gz", thetas, error_rates,engine);
	}

	printf("Estimated mismatch rate in Markov model is: %.5f\n", errorRate);
	printf("Total time:%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
}
