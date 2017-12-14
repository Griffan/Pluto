/*The MIT License (MIT)

Copyright (c) 2017 Fan Zhang, Hyun Min Kang

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */


#include <vector>
#include <map>
#include <unordered_map>
#include <fstream>

#include "ShotgunManners.h"
#include "OutputHandlers.h"
#include "DosageCalculator.h"
#include "MergeHaplotypes.h"
#include "HaplotypeLoader.h"
#include "Parameters.h"
#include "InputFile.h"
#include "Error.h"
#include "libVcfVcfFile.h"
#include "PBWTHaplotyper.h"

using namespace libVcf;

float *thetas = NULL;
int nthetas = 0;

float *error_rates = NULL;
int nerror_rates = 0;
//record marker name and relative index in unphased vcf
std::unordered_map<std::string, int> unphaseMarkerIdx;
//record marker name and relative index in phased ref vcf
std::unordered_map<std::string, int> refMarkerIdx;

std::unordered_map<std::string, bool> unifiedMarkerSet;//false:only shows in ref vcf; true:showed in both phased and unphased vcf

std::unordered_map<std::string, bool> pidIncludedInUnphasedVcf;
std::unordered_map<std::string, bool> pidIncludedInPhasedVcf;
std::unordered_map<std::string, bool> pidExcludedInUnphasedVcf;
std::unordered_map<std::string, bool> pidExcludedInPhasedVcf;

std::unordered_map<std::string, std::pair<int, int> > DuplicatedIndividualPair;

void SetErrorRateAndTransRate(PBWTHaplotyper &engine, double defaultErrorRate, double defaultTransRate, int markerindex,
                              const VcfMarker *pMarker, float& prevGeneticDistance, float& currentGeneticDistance) {
    if(engine.geneticMapAvailable)
    {
        engine.SetErrorRate(markerindex, defaultErrorRate);
        currentGeneticDistance=engine.GDMap.InferGeneticDistance(pMarker->sChrom.c_str(), pMarker->nPos);
        if (markerindex != 0) {//start from 2nd marker
            engine.thetas[markerindex-1] = engine.GDMap.CalculateRecombinationRate(prevGeneticDistance,currentGeneticDistance);
        }
        prevGeneticDistance=currentGeneticDistance;
    }
    else {
        //fprintf(stderr,"No ERATE or THETA tag found in input vcf, now using command line(--errorRate and --transRate)settings:\n Error Rate:%f\tTrans Rate(Theta):%f\n",defaultErrorRate, defaultTransRate);
        engine.SetErrorRate(markerindex, defaultErrorRate);
        if (markerindex != 0)
            engine.thetas[markerindex-1] = defaultTransRate;
    }
}
// print output files directly in VCF format
// inVcf contains skeleton of VCF information to copy from
// consensus contains the haplotype information to replace GT field
// dosage contains dosage information to be added as DS field
// thetas contains recombination rate information between markers
// error-rates contains per-marker error rates
// rsqs contains rsq_hat estimates ??
void OutputVCFConsensus(const String &inVcf, Pedigree &ped, ConsensusBuilder &consensus, const String &filename,
                        float *thetas, float *error_rates, PBWTHaplotyper &engine) {
    consensus.Merge(); // calculate consensus sequence

    if (consensus.stored)
        printf("Merged sampled haplotypes to generate consensus\n"
                       "   Consensus estimated to have %.1f errors in missing genotypes and %.1f flips in haplotypes\n\n",
               consensus.errors, consensus.flips);

    // read and write VCF inputs
    try {

        fprintf(stderr, "Outputing VCF file %s\n", filename.c_str());

        VcfFile *pVcf = new VcfFile;
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
        for (int i = 1; i < pVcf->asMetaKeys.Length(); ++i) {
        	if ((pVcf->asMetaKeys[i - 1].SubStr(0, 4).Compare("INFO") == 0) && (pVcf->asMetaKeys[i].SubStr(0, 4).Compare("INFO") != 0)) {
        		pVcf->asMetaKeys.InsertAt(i, "INFO");
        		pVcf->asMetaValues.InsertAt(i, "<ID=AVGPOST,Number=1,Type=Float,Description=\"Average Posterior Probability from thunderVCF\">");
        		++i;

        		pVcf->asMetaKeys.InsertAt(i, "INFO");
        		pVcf->asMetaValues.InsertAt(i, "<ID=RSQ,Number=1,Type=Float,Description=\"Imputation Quality from thunderVCF\">");
        		++i;

        		pVcf->asMetaKeys.InsertAt(i, "INFO");
        		pVcf->asMetaValues.InsertAt(i, "<ID=ERATE,Number=1,Type=Float,Description=\"Per-marker error rate from thunderVCF\">");
        		++i;

        		pVcf->asMetaKeys.InsertAt(i, "INFO");
        		pVcf->asMetaValues.InsertAt(i, "<ID=THETA,Number=1,Type=Float,Description=\"Recombination parameter with next marker from thunderVCF\">");
        		++i;
        	}
        	if ((pVcf->asMetaKeys[i - 1].SubStr(0, 6).Compare("FORMAT") != 0) && (pVcf->asMetaKeys[i].SubStr(0, 6).Compare("FORMAT") == 0)) {
        		pVcf->asMetaKeys.InsertAt(i + 1, "FORMAT");
        		pVcf->asMetaValues.InsertAt(i + 1, "<ID=DS,Number=1,Type=Integer,Description=\"Genotype dosage from thunderVCF\">");
        		++i;
        	}
        }



        char **haplotypes = consensus.consensus;

        // check the sanity of data
        if (pVcf->getSampleCount() == 0) {
            throw VcfFileException("No individual genotype information exist in the input VCF file %s",
                                   filename.c_str());
        }

        int nSamples = pVcf->getSampleCount();

        // build map of personID -> sampleIndex
        std::map<std::string, int> pedMap;
        for (int i = 0; i < (engine.individuals - engine.phased)/*ped.count*/; ++i) {
//            fprintf(stderr,"Adding (%s,%d)\n", ped[i].pid.c_str(), i);
            pedMap[ped[i].pid.c_str()] = i;

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
        VcfMarker *pMarker= nullptr;// = new VcfMarker;

        char sDose[255];
        double freq(0.), maf(0.), avgPost(0.), rsq(0.);
        String markerName;
        for (int m = 0; pVcf->iterateMarker(); ++m) {
            //fprintf(stderr,"m=%d\n",m);

            pMarker = pVcf->getLastMarker();
            markerName.printf("%s:%d:%s", pMarker->sChrom.c_str(), pMarker->nPos,pMarker->asAlts[0].c_str());//assume tri-alleles was divided into two lines

            if (unifiedMarkerSet[std::string(markerName.c_str())] ==
                true) {//if unphase also has this marker, update content otherwise remains as before
//                doses.CalculateMarkerInfo(m, freq, maf, avgPost, rsq);

                ////fprintf(stderr,"foo1\n");
                int nInfo = pMarker->asInfoKeys.Find("LDAF");
                if (nInfo < 0) {
                    sprintf(sDose, "%.4f", 1. - freq);
                    pMarker->asInfoKeys.Add("LDAF");
                    pMarker->asInfoValues.Add(sDose);
                }
                else
                    pMarker->asInfoValues[nInfo].printf("%.4f", 1. - freq);

                nInfo = pMarker->asInfoKeys.Find("AVGPOST");
                if (nInfo < 0) {
                    sprintf(sDose, "%.4f", avgPost);
                    pMarker->asInfoKeys.Add("AVGPOST");
                    pMarker->asInfoValues.Add(sDose);
                }
                else
                    pMarker->asInfoValues[nInfo].printf("%.4f", avgPost);

                nInfo = pMarker->asInfoKeys.Find("RSQ");
                if (nInfo < 0) {
                    sprintf(sDose, "%.4f", rsq);
                    pMarker->asInfoKeys.Add("RSQ");
                    pMarker->asInfoValues.Add(sDose);
                }
                else
                    pMarker->asInfoValues[nInfo].printf("%.4f", rsq);

                nInfo = pMarker->asInfoKeys.Find("ERATE");
                if (nInfo < 0) {
                    sprintf(sDose, "%.4f", nerror_rates ? error_rates[m] / nerror_rates : 0);
                    pMarker->asInfoKeys.Add("ERATE");
                    pMarker->asInfoValues.Add(sDose);
                }
                else
                    pMarker->asInfoValues[nInfo].printf("%.4f", nerror_rates ? error_rates[m] / nerror_rates : 0);

                nInfo = pMarker->asInfoKeys.Find("THETA");
                if (nInfo < 0) {
                    if (m != engine.markers - 1)
                        sprintf(sDose, "%.4f", nthetas ? thetas[m] / nthetas : 0);
                    else
                        sprintf(sDose, "%.4f", nthetas ? thetas[m - 1] / nthetas : 0);
                    pMarker->asInfoKeys.Add("THETA");
                    pMarker->asInfoValues.Add(sDose);
                }
                else {
                    if (m != engine.markers - 1)
                        pMarker->asInfoValues[nInfo].printf("%.4f", nthetas ? thetas[m] / nthetas : 0);
                    else
                        pMarker->asInfoValues[nInfo].printf("%.4f", nthetas ? thetas[m - 1] / nthetas : 0);
                }

                int GTidx = pMarker->asFormatKeys.Find("GT");
                if (GTidx < 0) {
                    throw VcfFileException("Cannot recognize GT key in FORMAT field");
                }

                int nFormats = pMarker->asFormatKeys.Length();

                pMarker->setSampleSize(static_cast<int>(vcf2ped.size()), pVcf->bParseGenotypes, pVcf->bParseDosages, pVcf->bParseValues);

                //fprintf(stderr,"nFormats=%d\tGTidx=%d\tDSidx=%d\n",nFormats,GTidx,DSidx);

                for (int i = 0; i < (int)vcf2ped.size(); ++i) {
                    int pi = vcf2ped[i];
                    // modify GT values;
                    if (pMarker->asAlts.Length() == 1) {
                        pMarker->asSampleValues[nFormats * i + GTidx].printf("%d|%d", haplotypes[pi * 2][m],
                                                                             haplotypes[pi * 2 + 1][m]);
                    }
                    else {
                        pMarker->asSampleValues[nFormats * i + GTidx].printf("%d|%d", haplotypes[pi * 2][m] + 1,
                                                                             haplotypes[pi * 2 + 1][m] + 1);
                    }
                }
            }
            pMarker->printVCFMarker(outVCF, false); // print marker to output file
        }
        delete pVcf;
        //delete pMarker;
        ifclose(outVCF);
    }
    catch (VcfFileException e) {
        error(e.what());
    }
}

void UnphasedSamplesOutputVCF(const String &inVcf, Pedigree &ped, const String &filename, float *thetas, float *error_rates,
                              PBWTHaplotyper &engine) {

    // read and write VCF inputs
    try {

        fprintf(stderr, "Outputing VCF file %s\n", filename.c_str());

        VcfFile *pVcf = new VcfFile;
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
        for (int i = 1; i < pVcf->asMetaKeys.Length(); ++i) {
        	if ((pVcf->asMetaKeys[i - 1].SubStr(0, 4).Compare("INFO") == 0) && (pVcf->asMetaKeys[i].SubStr(0, 4).Compare("INFO") != 0)) {
        		pVcf->asMetaKeys.InsertAt(i, "INFO");
        		pVcf->asMetaValues.InsertAt(i, "<ID=AVGPOST,Number=1,Type=Float,Description=\"Average Posterior Probability from thunderVCF\">");
        		++i;

        		pVcf->asMetaKeys.InsertAt(i, "INFO");
        		pVcf->asMetaValues.InsertAt(i, "<ID=RSQ,Number=1,Type=Float,Description=\"Imputation Quality from thunderVCF\">");
        		++i;

        		pVcf->asMetaKeys.InsertAt(i, "INFO");
        		pVcf->asMetaValues.InsertAt(i, "<ID=ERATE,Number=1,Type=Float,Description=\"Per-marker error rate from thunderVCF\">");
        		++i;

        		pVcf->asMetaKeys.InsertAt(i, "INFO");
        		pVcf->asMetaValues.InsertAt(i, "<ID=THETA,Number=1,Type=Float,Description=\"Recombination parameter with next marker from thunderVCF\">");
        		++i;
        	}
        	if ((pVcf->asMetaKeys[i - 1].SubStr(0, 6).Compare("FORMAT") != 0) && (pVcf->asMetaKeys[i].SubStr(0, 6).Compare("FORMAT") == 0)) {
        		pVcf->asMetaKeys.InsertAt(i + 1, "FORMAT");
        		pVcf->asMetaValues.InsertAt(i + 1, "<ID=DS,Number=1,Type=Integer,Description=\"Genotype dosage from thunderVCF\">");
        		++i;
        	}
        }

        // check the sanity of data
        if (pVcf->getSampleCount() == 0) {
            throw VcfFileException("No individual genotype information exist in the input VCF file %s",
                                   filename.c_str());
        }

        int nSamples = pVcf->getSampleCount();//num of samples need to be phased
        //fprintf(stderr, "we got %d samples\n", nSamples);
        // build map of personID -> sampleIndex
        std::map<std::string, int> pedMap;
        for (int i = 0; i < (engine.individuals - engine.phased)/*ped.count*/; ++i) {
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
        VcfMarker *pMarker = nullptr;// = new VcfMarker;
        char sDose[255];
        double freq(0.), maf(0.), avgPost(0.), rsq(0.);
        String markerName;
        int markerIndex(0);
        for (int m = 0; pVcf->iterateMarker(); ++m) {
            pMarker = pVcf->getLastMarker();
            markerName.printf("%s:%d:%s", pMarker->sChrom.c_str(), pMarker->nPos,pMarker->asAlts[0].c_str());//assume tri-alleles was divided into two lines
            if (unifiedMarkerSet.find(std::string(markerName.c_str())) != unifiedMarkerSet.end() &&
                unifiedMarkerSet[std::string(markerName.c_str())] == true)
            {//if unphase also has this marker i.e. shared marker, update content otherwise remains as before

                markerIndex = refMarkerIdx[std::string(markerName.c_str())];

//                doses.CalculateMarkerInfo(m, freq, maf, avgPost, rsq);

//                fprintf(stderr,"foo1,marker:%d\t%f\t%f\t%f\t%f\n",m,freq,maf,avgPost,rsq);

                int nInfo = pMarker->asInfoKeys.Find("LDAF");
                if (nInfo < 0) {
                    sprintf(sDose, "%.4f", 1. - freq);
                    pMarker->asInfoKeys.Add("LDAF");
                    pMarker->asInfoValues.Add(sDose);
                }
                else {
                    //fprintf(stderr,"pMarker->asInfoValues[nInfo]:%s\tfreq:%f\n",pMarker->asInfoValues[nInfo].c_str(),freq);
                    pMarker->asInfoValues[nInfo].printf("%.4f", 1. - freq);
                }

                nInfo = pMarker->asInfoKeys.Find("AVGPOST");
                if (nInfo < 0) {
                    sprintf(sDose, "%.4f", avgPost);
                    pMarker->asInfoKeys.Add("AVGPOST");
                    pMarker->asInfoValues.Add(sDose);
                }
                else
                    pMarker->asInfoValues[nInfo].printf("%.4f", avgPost);

                nInfo = pMarker->asInfoKeys.Find("RSQ");
                if (nInfo < 0) {
                    sprintf(sDose, "%.4f", rsq);
                    pMarker->asInfoKeys.Add("RSQ");
                    pMarker->asInfoValues.Add(sDose);
                }
                else
                    pMarker->asInfoValues[nInfo].printf("%.4f", rsq);

                nInfo = pMarker->asInfoKeys.Find("ERATE");
                if (nInfo < 0) {
                    sprintf(sDose, "%.4f", nerror_rates ? error_rates[m] / nerror_rates : 0);
                    pMarker->asInfoKeys.Add("ERATE");
                    pMarker->asInfoValues.Add(sDose);
                }
                else
                    pMarker->asInfoValues[nInfo].printf("%.4f", nerror_rates ? error_rates[m] / nerror_rates : 0);

                nInfo = pMarker->asInfoKeys.Find("THETA");
                if (nInfo < 0) {
                    if (m != engine.markers - 1)
                        sprintf(sDose, "%.4f", nthetas ? thetas[m] / nthetas : 0);
                    else
                        sprintf(sDose, "%.4f", nthetas ? thetas[m - 1] / nthetas : 0);
                    pMarker->asInfoKeys.Add("THETA");
                    pMarker->asInfoValues.Add(sDose);
                }
                else {
                    if (m != engine.markers - 1)
                        pMarker->asInfoValues[nInfo].printf("%.4f", nthetas ? thetas[m] / nthetas : 0);
                    else
                        pMarker->asInfoValues[nInfo].printf("%.4f", nthetas ? thetas[m - 1] / nthetas : 0);
                }

                int GTidx = pMarker->asFormatKeys.Find("GT");
                if (GTidx < 0) {
                    throw VcfFileException("Cannot recognize GT key in FORMAT field");
                }

                int nFormats = pMarker->asFormatKeys.Length();

                pMarker->setSampleSize(static_cast<int>(vcf2ped.size()), pVcf->bParseGenotypes, pVcf->bParseDosages, pVcf->bParseValues);

                for (int i = 0; i < (int)vcf2ped.size(); ++i) {
                    int pi = vcf2ped[i];
                    if (pMarker->asAlts.Length() == 1) {
                        pMarker->asSampleValues[nFormats * i + GTidx].printf("%d|%d", engine.haplotypes[pi * 2][markerIndex],
                                                                             engine.haplotypes[pi * 2 + 1][markerIndex]);
                    }
                    else {
                        pMarker->asSampleValues[nFormats * i + GTidx].printf("%d|%d", engine.haplotypes[pi * 2][markerIndex] + 1,
                                                                             engine.haplotypes[pi * 2 + 1][markerIndex] + 1);
                    }
                }
            }
            pMarker->printVCFMarker(outVCF, false); // print marker to output file
        }
        delete pVcf;
        //delete pMarker;
        ifclose(outVCF);
    }
    catch (VcfFileException e) {
        error(e.what());
    }
}

void UpdateVector(float *current, float *&vector, int &n, int length) {
    if (n++ == 0) {
        vector = new float[length];

        for (int i = 0; i < length; i++)
            vector[i] = current[i];
    }
    else
        for (int i = 0; i < length; i++)
            vector[i] += current[i];
}

void UpdateErrorRates(Errors *current, float *&vector, int &n, int length) {
    if (n++ == 0) {
        vector = new float[length];

        for (int i = 0; i < length; i++)
            vector[i] = current[i].rate;
    }
    else
        for (int i = 0; i < length; i++)
            vector[i] += current[i].rate;
}

int MemoryAllocationFailure() {
    printf("FATAL ERROR - Memory allocation failed\n");
    return -1;
}

void LoadPidToBeIncluded(String filename, String filename2) {
    if (filename != "") {
        std::ifstream fin(filename.c_str());
        if (!fin.is_open()) {
            std::cerr << "Open file " << filename << " failed!" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string tmpLine;
        while (getline(fin, tmpLine)) {
            pidIncludedInUnphasedVcf[tmpLine] = true;
        }
        fin.close();
    }
    if (filename2 != "") {
        std::ifstream fin2(filename2.c_str());
        if (!fin2.is_open()) {
            std::cerr << "Open file " << filename2 << " failed!" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string tmpLine;
        while (getline(fin2, tmpLine)) {
            pidIncludedInPhasedVcf[tmpLine] = true;
        }
        fin2.close();
    }
}

void LoadPidToBeExcluded(String filename, String filename2) {
    if (filename != "") {
        std::ifstream fin(filename.c_str());
        if (!fin.is_open()) {
            std::cerr << "Open file " << filename << " failed!" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string tmpLine;
        while (getline(fin, tmpLine)) {
            pidExcludedInUnphasedVcf[tmpLine] = true;
        }
        fin.close();
    }
    if (filename2 != "") {
        std::ifstream fin2(filename2.c_str());
        if (!fin2.is_open()) {
            std::cerr << "Open file " << filename2 << " failed!" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string tmpLine;
        while (getline(fin2, tmpLine)) {
            pidExcludedInPhasedVcf[tmpLine] = true;
        }
        fin2.close();
    }
}

void LoadSamples(Pedigree &ped, const String &filename, std::unordered_map<std::string, bool> &pidIncluded,
                 std::unordered_map<std::string, bool> &pidExcluded, int &num) {
    //printf("starting LoadSamples\n\n");

    try {

        VcfFile *pVcf = new VcfFile;
        pVcf->bSiteOnly = true;
        pVcf->bParseGenotypes = false;
        pVcf->bParseDosages = false;
        pVcf->bParseValues = false;
        pVcf->openForRead(filename.c_str());

        int PreviousCount = ped.count;
        // check the sanity of data
        if (pVcf->getSampleCount() == 0) {
            throw VcfFileException("No individual genotype information exist in the input VCF file %s",
                                   filename.c_str());
        }
        //std::cerr << "Sample Size:" << pVcf->getSampleCount() << std::endl;
        //std::cerr << "Include Size:" << pidIncluded.size() << std::endl;
        //std::cerr << "Exlude Size:" << pidExcluded.size() << std::endl;
        for (int i = 0; i < pVcf->getSampleCount(); ++i) {
            //std::cerr << "input:" << std::string(pVcf->vpVcfInds[i]->sIndID.c_str()) << "\t" << pidIncluded[std::string(pVcf->vpVcfInds[i]->sIndID.c_str())] << std::endl;
            if ((pidIncluded.size() == 0 ||
                 pidIncluded.find(std::string(pVcf->vpVcfInds[i]->sIndID.c_str())) != pidIncluded.end()) &&
                (pidExcluded.size() == 0 ||
                 pidExcluded.find(std::string(pVcf->vpVcfInds[i]->sIndID.c_str())) == pidExcluded.end()))
            {
                if (DuplicatedIndividualPair.find(std::string(pVcf->vpVcfInds[i]->sIndID.c_str())) ==
                    DuplicatedIndividualPair.end())//never showed before
                {
                    ped.AddPerson(pVcf->vpVcfInds[i]->sIndID, pVcf->vpVcfInds[i]->sIndID, "0", "0", 1, 1);
                    DuplicatedIndividualPair[std::string(pVcf->vpVcfInds[i]->sIndID.c_str())] = std::make_pair(i, -1);
                    //std::cerr << "Never saw" << pVcf->vpVcfInds[i]->sIndID << std::endl;
                }
                else {
                    warning("Individual %s duplicated!", pVcf->vpVcfInds[i]->sIndID.c_str());
                    ped.AddPerson(pVcf->vpVcfInds[i]->sIndID, pVcf->vpVcfInds[i]->sIndID + String("_dup"), "0", "0", 1,
                                  1);
                    DuplicatedIndividualPair[std::string(pVcf->vpVcfInds[i]->sIndID.c_str())].second =
                            i + PreviousCount;
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
//    printf("Loaded %d individuals from file %s\n\n", ped.count, filename.c_str());
}

void LoadRefPanelPolymorphicSites(const String &filename) {
    try {
        VcfFile *pVcf = new VcfFile;
        pVcf->bSiteOnly = true;
        pVcf->bParseGenotypes = false;
        pVcf->bParseDosages = false;
        pVcf->bParseValues = false;
        pVcf->openForRead(filename.c_str());

        VcfMarker *pMarker = nullptr;// = new VcfMarker;

        StringArray altalleles;
        String markerName;

        while (pVcf->iterateMarker()) {
            int markers = Pedigree::markerCount;
            pMarker = pVcf->getLastMarker();

            markerName.printf("%s:%d:%s", pMarker->sChrom.c_str(), pMarker->nPos,pMarker->asAlts[0].c_str());//assume tri-alleles was divided into two lines
            int marker = Pedigree::GetMarkerID(markerName);// that's where we fill up marker numbers for ped file
            //initialize flag map to allocate memory equals to marker number in ref set
            unifiedMarkerSet[std::string(markerName.c_str())] = false;
            refMarkerIdx[std::string(markerName.c_str())] = marker;
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
                warning("Each polymorphic site should only occur once, but site %s is duplicated\n", markerName.c_str());
            }

            if (al1 != 1 || al2 != 2) {
                error("Allele labels '%s' and '%s' for polymorphic site '%s' are not valid\n",
                      (const char *) altalleles[0], (const char *) altalleles[1], markerName.c_str());
            }
        }
        delete pVcf;
        //delete pMarker;
    }
    catch (VcfFileException e) {
        error(e.what());
    }
}
// use ref vcf markers as backbone, ignore sites that are not shown in ref vcf
void LoadUnphasedPolymorphicSites(const String &filename) {
    try {
        VcfFile *pVcf = new VcfFile;
        pVcf->bSiteOnly = true;
        pVcf->bParseGenotypes = false;
        pVcf->bParseDosages = false;
        pVcf->bParseValues = false;
        pVcf->openForRead(filename.c_str());

        VcfMarker *pMarker = nullptr;//new VcfMarker;
        StringArray altalleles;
        String markerName;
        int localIdx = 0;


        while (pVcf->iterateMarker()) {

            pMarker = pVcf->getLastMarker();
            markerName.printf("%s:%d:%s", pMarker->sChrom.c_str(), pMarker->nPos,pMarker->asAlts[0].c_str());//assume tri-alleles was divided into two lines
            int idx = Pedigree::markerLookup.Integer(markerName);//only look up, no add in
            if (idx != -1)//shown in ref panel marker set
            {
                unifiedMarkerSet[std::string(markerName.c_str())] = true;//also shown in unphased vcf
                unphaseMarkerIdx[std::string(markerName.c_str())] = localIdx;//index in unphased vcf
            }
            else {//not shown in ref
                unphaseMarkerIdx[std::string(markerName.c_str())] = localIdx;
            }
            ++localIdx;
        }
        delete pVcf;
        //delete pMarker;
    }
    catch (VcfFileException e) {
        error(e.what());
    }
}

void LoadGenotypeFromUnphasedVCF(Pedigree &ped, const String &filename, int maxPhred, PBWTHaplotyper &engine) {
//printf("starting LoadGenotypeFromUnphasedVCF\n\n");
    try {
        VcfFile *pVcf = new VcfFile;
        pVcf->bSiteOnly = false;
        pVcf->bParseGenotypes = false;
        pVcf->bParseDosages = false;
        pVcf->bParseValues = true;
        pVcf->openForRead(filename.c_str());

        // check the sanity of data
        if (pVcf->getSampleCount() == 0) {
            throw VcfFileException("No individual genotype information exist in the input VCF file %s",
                                   filename.c_str());
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

        for (int i = 0; i < (engine.individuals - engine.phased);/* ped.count;*/ i++) {//first N samples in ped, not necessary the unphased sample
            if (originalPeople.Integer(ped[i].famid + "." + ped[i].pid) != -1) {//find index of this sample in current vcf
                personIndices[originalPeople.Integer(ped[i].famid + "." + ped[i].pid)] = i;//map index in current vcf to index in ped file
            }

        }

        int markerindex = 0;
        VcfMarker *pMarker = nullptr;//new VcfMarker;
        String markerName;
        while (pVcf->iterateMarker()) {//for each marker

            pMarker = pVcf->getLastMarker();
            markerName.printf("%s:%d:%s", pMarker->sChrom.c_str(), pMarker->nPos,pMarker->asAlts[0].c_str());//assume tri-alleles was divided into two lines
//            printf("now for marker %s:%d\t", pMarker->sChrom.c_str(), pMarker->nPos);
            if (unifiedMarkerSet.find(std::string(markerName.c_str())) != unifiedMarkerSet.end() &&
                unifiedMarkerSet[std::string(markerName.c_str())] == true)
                markerindex = refMarkerIdx[std::string(markerName.c_str())];
            else
                continue;
            //int AFidx = pMarker->asInfoKeys.Find("AF");

            int PLidx = pMarker->asFormatKeys.Find("PL");
            int GLidx = pMarker->asFormatKeys.Find("GL");
            int GTidx = pMarker->asFormatKeys.Find("GT");

            if (PLidx < 0 && GLidx < 0 && GTidx <0) {
                    throw VcfFileException("Cannot recognize GT, GL or PL key in FORMAT field");
            }
            //printf("reading vcf 1\n\n");
            int formatLength = pMarker->asFormatKeys.Length();
            int idx11 = 0, idx12 = 1, idx22 = 2;

            StringArray phred;
            int genoindex = markerindex * 3;

            long phred11(-1), phred12(-1), phred22(-1);
            for (int i = 0; i < nSamples; i++)//for each individual in current vcf
            {
                if (personIndices.find(i) != personIndices.end()) {//not in index mapping relations, which is impossible

                    if(PLidx >=0)//found PL
                    {
                        phred.ReplaceTokens(pMarker->asSampleValues[PLidx + i * formatLength], ",");
                        phred11=phred[idx11].AsInteger();
                        phred12=phred[idx12].AsInteger();
                        phred22=phred[idx22].AsInteger();
                    }
                    else if(GLidx >=0)//found GL
                    {
                        phred.ReplaceTokens(pMarker->asSampleValues[GLidx + i * formatLength], ",");
                        phred11=static_cast<int>(-10. * phred[idx11].AsDouble());
                        phred12=static_cast<int>(-10. * phred[idx12].AsDouble());
                        phred22=static_cast<int>(-10. * phred[idx22].AsDouble());
                    }
                    if(GTidx >=0)//found GT
                    {
                        phred.ReplaceTokens(pMarker->asSampleValues[GTidx + i * formatLength], "|/");
                        long geno=phred[0].AsInteger()+phred[1].AsInteger();
                        if(geno==0)
                        {
                            phred11=0;
                            phred12=30;
                            phred22=50;
                        }
                        else if(geno==1)
                        {
                            phred11=50;
                            phred12=0;
                            phred22=50;
                        }
                        else
                        {
                            phred11=50;
                            phred12=30;
                            phred22=0;
                        }
                    }

                    if ((phred11 < 0) || (phred12 < 0) || (phred22 < 0)) {
                        error("Negative PL or Positive GL observed");
                    }

                    if (phred11 > maxPhred) phred11 = maxPhred;
                    if (phred12 > maxPhred) phred12 = maxPhred;
                    if (phred22 > maxPhred) phred22 = maxPhred;
//
//                    printf("phred scores are %f, %f, %f;\tphred11/12/22 %d, %d, %d\n", phred[idx11].AsDouble(), phred[idx12].AsDouble(), phred[idx22].AsDouble(),phred11,phred12,phred22);

                    engine.genotypes[personIndices[i]][genoindex] = static_cast<char>(phred11);
                    engine.genotypes[personIndices[i]][genoindex + 1] = static_cast<char>(phred12);
                    engine.genotypes[personIndices[i]][genoindex + 2] = static_cast<char>(phred22);
//                    fprintf(stderr,"marker:%d\t%d\t%d\t%d\n",markerindex,engine.genotypes[personIndices[i]][genoindex],engine.genotypes[personIndices[i]][genoindex + 1],engine.genotypes[personIndices[i]][genoindex + 2] );
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

void LoadGenotypeAndHaplotypeFromPhasedVCF(Pedigree &ped, const String &filename, int maxPhred, int phased,
                               PBWTHaplotyper &engine, double defaultErrorRate, double defaultTransRate) {
    //printf("starting LoadPhasedVcf\n\n");

    try {
        VcfFile *pVcf = new VcfFile;
        pVcf->bSiteOnly = false;
        pVcf->bParseGenotypes = false;
        pVcf->bParseDosages = false;
        pVcf->bParseValues = true;
        pVcf->openForRead(filename.c_str());

        // check the sanity of data
        if (pVcf->getSampleCount() == 0) {
            throw VcfFileException("No individual genotype information exist in the input VCF file %s",
                                   filename.c_str());
        }
//        std::vector<int> phaseIdx(ped.count, -1);//initially assume all the individuals are phased
        int nSamples = pVcf->getSampleCount();

        //vector<int> personIndices(ped.count, -1);
        std::unordered_map<int, int> personIndices;
        StringIntHash sampleOrderInCurrentVcf; // key: famid+subID, value: original order (0 based); in phased file
        int person = 0;
        for (int i = 0; i < nSamples; i++) {//add samples in current phased file into sampleOrderInCurrentVcf
            {
                //std::cerr << "if ordered:" << pVcf->vpVcfInds[i]->sIndID << std::endl;
                sampleOrderInCurrentVcf.Add(pVcf->vpVcfInds[i]->sIndID + "." + pVcf->vpVcfInds[i]->sIndID, person);
                person++;
            }
        }
//        for (int j = 0; j <engine.individuals; ++j) {
//            std::cerr<<ped[j].famid + "." + ped[j].pid<<std::endl;
//        }

        for (int i = (engine.individuals - engine.phased); i < engine.individuals; i++) {
            int idx = sampleOrderInCurrentVcf.Integer(ped[i].famid + "." + ped[i].pid);//this requires assumption indivisuals in ped stored as unphased individuals + phased individuals
            if (idx != -1)//phased in this vcf
            {
                personIndices[idx] = i;//put idx sample in this vcf into i th position in engine

            }
        }

        int markerindex = 0;
        VcfMarker *pMarker = nullptr;
        String markerName;
        float prevGeneticDistance=0.f;
        float currentGeneticDistance=0.f;
        while (pVcf->iterateMarker()) {//for each marker

            pMarker = pVcf->getLastMarker();
            //markerName.printf("%s:%d:%s", pMarker->sChrom.c_str(), pMarker->nPos,pMarker->asAlts[0].c_str());//assume tri-alleles was divided into two lines

            engine.refalleles[markerindex] = pMarker->sRef[0];

            int AFidx = pMarker->asInfoKeys.Find("AF");
            if (AFidx == -1) {
                pMarker->asInfoKeys.PrintLine();
                pMarker->asInfoValues.PrintLine();
            }
            /*setting error rate*/
//            int ERidx = pMarker->asInfoKeys.Find("ERATE");
//            int THidx = pMarker->asInfoKeys.Find("THETA");
//            if (ERidx != -1 && THidx != -1)// ERATE, THETA exist
//            {
//                //error_rates[markerindex] = pMarker->asInfoValues[THidx].AsDouble();
//                engine.SetErrorRate(markerindex, pMarker->asInfoValues[ERidx].AsDouble());
//                if (markerindex != engine.markers - 1)
//                    engine.thetas[markerindex] = pMarker->asInfoValues[THidx].AsDouble();
//            }
//            else
            SetErrorRateAndTransRate(engine, defaultErrorRate, defaultTransRate, markerindex, pMarker,
                                     prevGeneticDistance, currentGeneticDistance);

            int PLidx = pMarker->asFormatKeys.Find("PL");
            int GLidx = pMarker->asFormatKeys.Find("GL");
            int GTidx = pMarker->asFormatKeys.Find("GT");
            if (GTidx <0) {
                throw VcfFileException("Cannot recognize GT key in FORMAT field");
            }
            //printf("reading vcf 2\n\n");
            int formatLength = pMarker->asFormatKeys.Length();
            int idx11 = 0, idx12 = 1, idx22 = 2;
            if (AFidx == -1) {
                int ANidx = pMarker->asInfoKeys.Find("AN");
                int ACidx = pMarker->asInfoKeys.Find("AC");
                if ((ANidx < 0) || (ACidx < 0)) {
                    //throw VcfFileException("Cannot recognize AF key in FORMAT field");
                    engine.freq1s[markerindex]=1e-6;
                }
                else {
                    engine.freq1s[markerindex] = 1. - (pMarker->asInfoValues[ACidx].AsDouble() + .5) /
                                               (pMarker->asInfoValues[ANidx].AsDouble() + 1.);
                }
            }
            else if (pMarker->asAlts.Length() == 1) {
                engine.freq1s[markerindex] = (1. - pMarker->asInfoValues[AFidx].AsDouble());
            }
            else {
                // AF1,AF2 -- freq1s is AF1
                engine.freq1s[markerindex] = pMarker->asInfoValues[AFidx].AsDouble();
            }

            StringArray phred;
            int genoindex = markerindex * 3;
            long phred11(-1),phred12(-1),phred22(-1);
            for (int i = 0; i < nSamples; i++)//for each phased individual
            {
                //printf("phred scores are   %d, %d\n", i, ped.count);
                if (personIndices.find(i) != personIndices.end())
                {
                    if(PLidx>=0)//found PL
                    {
                        phred.ReplaceTokens(pMarker->asSampleValues[PLidx + i * formatLength], ",");
                        phred11=phred[idx11].AsInteger();
                        phred12=phred[idx12].AsInteger();
                        phred22=phred[idx22].AsInteger();
                    }
                    else if(GLidx >= 0)//found GL
                    {
                        phred.ReplaceTokens(pMarker->asSampleValues[GLidx + i * formatLength], ",");
                        phred11=static_cast<int>(-10. * phred[idx11].AsDouble());
                        phred12=static_cast<int>(-10. * phred[idx12].AsDouble());
                        phred22=static_cast<int>(-10. * phred[idx22].AsDouble());
                    }
                    if(GTidx >= 0)//found GT
                    {
                        phred.ReplaceTokens(pMarker->asSampleValues[GTidx + i * formatLength], "|/");
                        long geno=phred[0].AsInteger()+phred[1].AsInteger();
                        if(geno==0)
                        {
                            phred11=0;
                            phred12=30;
                            phred22=50;
                        }
                        else if(geno==1)
                        {
                            phred11=50;
                            phred12=0;
                            phred22=50;
                        }
                        else
                        {
                            phred11=50;
                            phred12=30;
                            phred22=0;
                        }
                        engine.haplotypes[personIndices[i] * 2][markerindex] = static_cast<char>(phred[0].AsInteger());
                        engine.haplotypes[personIndices[i] * 2 + 1][markerindex] = static_cast<char>(phred[1].AsInteger());
                    }
                    if ((phred11 < 0) || (phred12 < 0) || (phred22 < 0)) {
                        error("Negative PL or Positive GL observed");
                    }

                    //printf("phred scores are %d, %d, %d, %d, %d\n", phred11, phred12, phred22,i, ped.count );

                    if (phred11 > maxPhred) phred11 = maxPhred;
                    if (phred12 > maxPhred) phred12 = maxPhred;
                    if (phred22 > maxPhred) phred22 = maxPhred;

                    engine.genotypes[personIndices[i]][genoindex] = static_cast<char>(phred11);
                    engine.genotypes[personIndices[i]][genoindex + 1] = static_cast<char>(phred12);
                    engine.genotypes[personIndices[i]][genoindex + 2] = static_cast<char>(phred22);
                }
            }
            //printf("reading vcf 4\n\n");

            ++markerindex;
        }
        //initialize missing value in unphased file as 0
        for (std::unordered_map<std::string, bool>::const_iterator iter = unifiedMarkerSet.begin();
             iter != unifiedMarkerSet.end(); ++iter) {
            if (iter->second == false)// if this marker not shown in unphased set but in reference set
            {
                markerindex = refMarkerIdx[iter->first];
                int genoindex = markerindex * 3;
                for (int i = 0; i < (engine.individuals - engine.phased); i++)//for each unphased individual
                {
                    //if(phaseIdx[i]==-1) continue;
                    //phred.ReplaceTokens(pMarker->asSampleValues[PLidx + i*formatLength], ",");
                    int phred11 = 0;// GLflag ? static_cast<int>(-10. * phred[idx11].AsDouble()) : phred[idx11].AsInteger();
                    int phred12 = 0;// GLflag ? static_cast<int>(-10. * phred[idx12].AsDouble()) : phred[idx12].AsInteger();
                    int phred22 = 0;// GLflag ? static_cast<int>(-10. * phred[idx22].AsDouble()) : phred[idx22].AsInteger();

                    if ((phred11 < 0) || (phred12 < 0) || (phred22 < 0)) {
                        error("Negative PL or Positive GL observed");
                    }

                    //printf("phred scores are %d, %d, %d\n", phred11, phred12, phred22);

                    if (phred11 > maxPhred) phred11 = maxPhred;
                    if (phred12 > maxPhred) phred12 = maxPhred;
                    if (phred22 > maxPhred) phred22 = maxPhred;

                    engine.genotypes[i][genoindex] = static_cast<char>(phred11);
                    engine.genotypes[i][genoindex + 1] = static_cast<char>(phred12);
                    engine.genotypes[i][genoindex + 2] = static_cast<char>(phred22);
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


//graph reading version
std::vector<float> errorRateHolder;
std::vector<float> transRateHolder;
void LoadRefPanelPolymorphicSites(const std::string &filename) {
    try {
        std::ifstream fin(filename+".site");
        if(!fin.is_open())
        {
            fprintf(stderr,"Cannot open file %s.site! Make sure run index first",filename.c_str());
        }
        std::string line,chrom,nPos,asAlt,errorRate,transRate;
        StringArray altalleles;
        String markerName;

        while (std::getline(fin,line)) {
            std::stringstream ss(line);
            ss>>chrom>>nPos>>asAlt>>errorRate>>transRate;
            markerName.printf("%s:%s:%s", chrom.c_str(), nPos.c_str(),asAlt.c_str());//assume tri-alleles was divided into two lines
            int marker = Pedigree::GetMarkerID(markerName);// that's where we fill up marker numbers for ped file
            //initialize flag map to allocate memory equals to marker number in ref set
            unifiedMarkerSet[std::string(markerName.c_str())] = false;
            refMarkerIdx[std::string(markerName.c_str())] = marker;
            unphaseMarkerIdx[std::string(markerName.c_str())] = -1;
            errorRateHolder.push_back(atof(errorRate.c_str()));
            transRateHolder.push_back(atof(transRate.c_str()));
        }
    }
    catch (VcfFileException e) {
        error(e.what());
    }
}
void LoadGenotypeFromUnphasedVCF(Pedigree &ped, const String &filename, PBWTHaplotyper &engine,double defaultErrorRate, double defaultTransRate) {
    //printf("starting LoadGenotypeFromUnphasedVCF\n\n");
    try {
        VcfFile *pVcf = new VcfFile;
        pVcf->bSiteOnly = false;
        pVcf->bParseGenotypes = false;
        pVcf->bParseDosages = false;
        pVcf->bParseValues = true;
        pVcf->openForRead(filename.c_str());

        // check the sanity of data
        if (pVcf->getSampleCount() == 0) {
            throw VcfFileException("No individual genotype information exist in the input VCF file %s",
                                   filename.c_str());
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

        for (int i = 0; i < engine.individuals;i++) {//first N samples in ped, not necessary the unphased sample
            int idx=originalPeople.Integer(ped[i].famid + "." + ped[i].pid);
            if ( idx != -1) {//find index of this sample in current vcf
                personIndices[idx] = i;//map index in current vcf to index in ped file
            }

        }

        int markerindex = 0;
        VcfMarker *pMarker = nullptr;//new VcfMarker;
        String markerName;
        float prevGeneticDistance=0.f;
        float currentGeneticDistance=0.f;
        while (pVcf->iterateMarker()) {//for each marker

            pMarker = pVcf->getLastMarker();
            markerName.printf("%s:%d:%s", pMarker->sChrom.c_str(), pMarker->nPos,pMarker->asAlts[0].c_str());//assume tri-alleles was divided into two lines

            if (unifiedMarkerSet.find(std::string(markerName.c_str())) != unifiedMarkerSet.end()) {
                unifiedMarkerSet[std::string(markerName.c_str())] = true;
                markerindex = refMarkerIdx[std::string(markerName.c_str())];
            }
            else {
                fprintf(stderr,"skip marker %s which not shown in reference panel!\n",markerName.c_str());
                continue;//leave out markers that not in ref panel
            }
            //int AFidx = pMarker->asInfoKeys.Find("AF");

            int PLidx = pMarker->asFormatKeys.Find("PL");
            int GLidx = pMarker->asFormatKeys.Find("GL");
            int GTidx = pMarker->asFormatKeys.Find("GT");

            if (PLidx < 0 && GLidx < 0 && GTidx <0) {
                throw VcfFileException("Cannot recognize GT, GL or PL key in FORMAT field");
            }
            //printf("reading vcf 1\n\n");
            int formatLength = pMarker->asFormatKeys.Length();
            int idx11 = 0, idx12 = 1, idx22 = 2;

            StringArray phred;
            int genoindex = markerindex * 3;

            long phred11(-1), phred12(-1), phred22(-1);
            for (int i = 0; i < nSamples; i++)//for each individual in current vcf
            {
                if (personIndices.find(i) != personIndices.end()) {//not in index mapping relations, which is impossible

                    if(PLidx >=0)//found PL
                    {
                        phred.ReplaceTokens(pMarker->asSampleValues[PLidx + i * formatLength], ",");
                        phred11=phred[idx11].AsInteger();
                        phred12=phred[idx12].AsInteger();
                        phred22=phred[idx22].AsInteger();
                    }
                    else if(GLidx >=0)//found GL
                    {
                        phred.ReplaceTokens(pMarker->asSampleValues[GLidx + i * formatLength], ",");
                        phred11=static_cast<int>(-10. * phred[idx11].AsDouble());
                        phred12=static_cast<int>(-10. * phred[idx12].AsDouble());
                        phred22=static_cast<int>(-10. * phred[idx22].AsDouble());
                    }
                    if(GTidx >=0)//found GT
                    {
                        phred.ReplaceTokens(pMarker->asSampleValues[GTidx + i * formatLength], "|/");
                        long geno=phred[0].AsInteger()+phred[1].AsInteger();
                        if(geno==0)
                        {
                            phred11=0;
                            phred12=30;
                            phred22=50;
                        }
                        else if(geno==1)
                        {
                            phred11=50;
                            phred12=0;
                            phred22=50;
                        }
                        else
                        {
                            phred11=50;
                            phred12=30;
                            phred22=0;
                        }
                    }

                    if ((phred11 < 0) || (phred12 < 0) || (phred22 < 0)) {
                        error("Negative PL or Positive GL observed");
                    }

                    if (phred11 > 255) phred11 = 255;
                    if (phred12 > 255) phred12 = 255;
                    if (phred22 > 255) phred22 = 255;
//
//                    printf("phred scores are %f, %f, %f;\tphred11/12/22 %d, %d, %d\n", phred[idx11].AsDouble(), phred[idx12].AsDouble(), phred[idx22].AsDouble(),phred11,phred12,phred22);

                    engine.genotypes[personIndices[i]][genoindex] = static_cast<char>(phred11);
                    engine.genotypes[personIndices[i]][genoindex + 1] = static_cast<char>(phred12);
                    engine.genotypes[personIndices[i]][genoindex + 2] = static_cast<char>(phred22);
//                    fprintf(stderr,"marker:%d\t%d\t%d\t%d\n",markerindex,engine.genotypes[personIndices[i]][genoindex],engine.genotypes[personIndices[i]][genoindex + 1],engine.genotypes[personIndices[i]][genoindex + 2] );
                }
            }
        }

        for (std::unordered_map<std::string, bool>::const_iterator iter = unifiedMarkerSet.begin();
             iter != unifiedMarkerSet.end(); ++iter) {
            if (iter->second == false)// if this marker not shown in unphased set but in reference set
            {
                markerindex = refMarkerIdx[iter->first];
                int genoindex = markerindex * 3;
                for (int i = 0; i < (engine.individuals - engine.phased); i++)//for each unphased individual
                {
                    engine.genotypes[i][genoindex] = 0;
                    engine.genotypes[i][genoindex + 1] = 0;
                    engine.genotypes[i][genoindex + 2] = 0;
                }
            }
        }

        delete pVcf;
    }
    catch (VcfFileException e) {
        error(e.what());
    }
}
//graph building version
void LoadGenotypeAndHaplotypeFromPhasedVCF(Pedigree &ped, const String &filename, PBWTHaplotyper &engine, double defaultErrorRate, double defaultTransRate) {
    try {
        VcfFile *pVcf = new VcfFile;
        pVcf->bSiteOnly = false;
        pVcf->bParseGenotypes = false;
        pVcf->bParseDosages = false;
        pVcf->bParseValues = true;
        pVcf->openForRead(filename.c_str());
        std::ofstream refVCFSite(std::string(filename.c_str())+".site");
        // check the sanity of data
        if (pVcf->getSampleCount() == 0) {
            throw VcfFileException("No individual genotype information exist in the input VCF file %s",
                                   filename.c_str());
        }
        int nSamples = pVcf->getSampleCount();

        std::unordered_map<int, int> personIndices;//order in VCF -> order in haplotype memory
        StringIntHash sampleOrderInCurrentVcf; // key: famid+subID, value: original order (0 based); in phased file
        int person = 0;
        for (int i = 0; i < nSamples; i++) {//add samples in current phased file into sampleOrderInCurrentVcf
            {
                //std::cerr << "if ordered:" << pVcf->vpVcfInds[i]->sIndID << std::endl;
                sampleOrderInCurrentVcf.Add(pVcf->vpVcfInds[i]->sIndID + "." + pVcf->vpVcfInds[i]->sIndID, person);
                person++;
            }
        }

        for (int i = 0; i < engine.individuals; i++) {
            int idx = sampleOrderInCurrentVcf.Integer(ped[i].famid + "." + ped[i].pid);//this requires assumption indivisuals in ped stored as unphased individuals + phased individuals
            if (idx != -1)//phased in this vcf
            {
                personIndices[idx] = i;//put idx sample in this vcf into i th position in engine
            }
        }

        int markerindex = 0;
        VcfMarker *pMarker = nullptr;
        String markerName;

        double prevGeneticDistance(0.),currentGeneticDistance(0.);
        double tmpTransRate(0.);
        while (pVcf->iterateMarker()) {//for each marker
            pMarker = pVcf->getLastMarker();

            markerName.printf("%s:%d:%s", pMarker->sChrom.c_str(), pMarker->nPos,pMarker->asAlts[0].c_str());//assume tri-alleles was divided into two lines
            //initialize flag map to allocate memory equals to marker number in ref set
            unifiedMarkerSet[std::string(markerName.c_str())] = false;
            refMarkerIdx[std::string(markerName.c_str())] = markerindex;

            if(engine.geneticMapAvailable)
            {
                    currentGeneticDistance=engine.GDMap.InferGeneticDistance(pMarker->sChrom.c_str(), pMarker->nPos);
                    if (markerindex != 0) {//start from 2nd marker
                        tmpTransRate = engine.GDMap.CalculateRecombinationRate(prevGeneticDistance,currentGeneticDistance);
                    }
                    prevGeneticDistance=currentGeneticDistance;
            }
            else {
                    //fprintf(stderr,"No ERATE or THETA tag found in input vcf, now using command line(--errorRate and --transRate)settings:\n Error Rate:%f\tTrans Rate(Theta):%f\n",defaultErrorRate, defaultTransRate);
                    if (markerindex != 0)
                        tmpTransRate = defaultTransRate;
            }
            refVCFSite<<pMarker->sChrom.c_str()<<"\t"<<pMarker->nPos<<"\t"<<pMarker->asAlts[0].c_str()<<"\t"<<defaultErrorRate<<"\t"<<tmpTransRate<<std::endl;

            int PLidx = pMarker->asFormatKeys.Find("PL");
            int GLidx = pMarker->asFormatKeys.Find("GL");
            int GTidx = pMarker->asFormatKeys.Find("GT");
            if (GTidx <0) {
                throw VcfFileException("Cannot recognize GT key in FORMAT field");
            }
            int formatLength = pMarker->asFormatKeys.Length();

            StringArray phred;
//            int genoindex = markerindex * 3;
            long phred11(-1),phred12(-1),phred22(-1);
            for (int i = 0; i < nSamples; i++)//for each phased individual
            {
                //printf("phred scores are   %d, %d\n", i, ped.count);
                if (personIndices.find(i) != personIndices.end())
                {
                    if(PLidx>=0)//found PL
                    {
                        phred.ReplaceTokens(pMarker->asSampleValues[PLidx + i * formatLength], ",");
                        phred11=phred[0].AsInteger();
                        phred12=phred[1].AsInteger();
                        phred22=phred[2].AsInteger();
                    }
                    else if(GLidx >= 0)//found GL
                    {
                        phred.ReplaceTokens(pMarker->asSampleValues[GLidx + i * formatLength], ",");
                        phred11=static_cast<int>(-10. * phred[0].AsDouble());
                        phred12=static_cast<int>(-10. * phred[1].AsDouble());
                        phred22=static_cast<int>(-10. * phred[2].AsDouble());
                    }
                    if(GTidx >= 0)//found GT
                    {
                        phred.ReplaceTokens(pMarker->asSampleValues[GTidx + i * formatLength], "|/");
//                        long geno=phred[0].AsInteger()+phred[1].AsInteger();
//                        if(geno==0)
//                        {
//                            phred11=0;
//                            phred12=30;
//                            phred22=50;
//                        }
//                        else if(geno==1)
//                        {
//                            phred11=50;
//                            phred12=0;
//                            phred22=50;
//                        }
//                        else
//                        {
//                            phred11=50;
//                            phred12=30;
//                            phred22=0;
//                        }
                        engine.haplotypes[personIndices[i] * 2][markerindex] = static_cast<char>(phred[0].AsInteger());
                        engine.haplotypes[personIndices[i] * 2 + 1][markerindex] = static_cast<char>(phred[1].AsInteger());
                    } else
                    {
                        fprintf(stderr,"no GT field found for phased VCF!");
                        exit(EXIT_FAILURE);
                    }
//                    if ((phred11 < 0) || (phred12 < 0) || (phred22 < 0)) {
//                        error("Negative PL or Positive GL observed");
//                    }

                    //printf("phred scores are %d, %d, %d, %d, %d\n", phred11, phred12, phred22,i, ped.count );

//                    if (phred11 > 255) phred11 = 255;
//                    if (phred12 > 255) phred12 = 255;
//                    if (phred22 > 255) phred22 = 255;
//
//                    engine.genotypes[personIndices[i]][genoindex] = static_cast<char>(phred11);
//                    engine.genotypes[personIndices[i]][genoindex + 1] = static_cast<char>(phred12);
//                    engine.genotypes[personIndices[i]][genoindex + 2] = static_cast<char>(phred22);
                }
            }
            ++markerindex;
        }
        refVCFSite.close();
        delete pVcf;
    }
    catch (VcfFileException e) {
        error(e.what());
    }
}
int BuildGraph(int argc, char **argv) {

    String outfile("mach1.out"), phasedfile("Empty"), pidIncludeFromUnphased(
            ""), pidIncludeFromPhased(
            ""), pidExcludeFromUnphased(""), pidExcludeFromPhased(""), PMatrix(""),calPMatrix("");
    String GDFile;

    clock_t t;
    t = clock();
    int seed = 123456,samplingRounds = 1;

    int prefixLength = 120;
    double errorRate = 0.01;
    double transRate = 0.01;
    bool onlyHeterSite = false;

    SetupCrashHandlers();
    SetCrashExplanation("reading command line options");

    printf("Pluto 0.01 -- Markov Chain Haplotyping for Shotgun Sequence Data\n"
                   "(c) 2015-2017 Fan Zhang, Goncalo Abecasis, and Hyun Min Kang\n\n");

    ParameterList pl;

    BEGIN_LONG_PARAMETERS(longParameters)
                    LONG_PARAMETER_GROUP("Shotgun Sequences")
                    LONG_STRINGPARAMETER("refVCF", &phasedfile)
                    LONG_PARAMETER_GROUP("Optional Files")
                    LONG_STRINGPARAMETER("includeUnphasedIDs", &pidIncludeFromUnphased)
                    LONG_STRINGPARAMETER("includePhasedIDs", &pidIncludeFromPhased)
                    LONG_STRINGPARAMETER("excludeUnphasedIDs", &pidExcludeFromUnphased)
                    LONG_STRINGPARAMETER("excludePhasedIDs", &pidExcludeFromPhased)
                    LONG_PARAMETER_GROUP("Graph Builder")
                    LONG_INTPARAMETER("graphComplexity", &prefixLength)
                    LONG_STRINGPARAMETER("PvalueMatrix", &PMatrix)
                    LONG_STRINGPARAMETER("calPvalueMatrix", &calPMatrix)
                    LONG_STRINGPARAMETER("geneticDistance", &GDFile)
                    LONG_INTPARAMETER("seed", &seed)
                    LONG_PARAMETER("onlyHeterSite", &onlyHeterSite)
                    LONG_PARAMETER_GROUP("Output Files")
                    LONG_STRINGPARAMETER("outPrefix", &outfile)
    END_LONG_PARAMETERS();

    pl.Add(new LongParameters("Available Options", longParameters));

    pl.Read(argc, argv);
    pl.Status();

    // Setup random seed ...
    globalRandom.Reset(seed);

    PBWTHaplotyper engine;//declaration of engine, also will call default constructor
    engine.nSampleCopy = samplingRounds;
    engine.onlyHeterSite = onlyHeterSite;
    engine.prefixLength = prefixLength;
    engine.outputPrefix = std::string(phasedfile.c_str())+".DAG";

    SetCrashExplanation("loading Pvalue Matrix");

    if(PMatrix.IsEmpty() and calPMatrix.IsEmpty()) {
        std::cerr<<"parameter --PvalueMatrix [PATH] or --calPvalueMatrix [PATH] required!"<<std::endl;
        exit(EXIT_FAILURE);
    }
    else if(!PMatrix.IsEmpty()) {
        engine.ReadPvalueMatrix(std::string(PMatrix.c_str()));
        fprintf(stderr, "Done reading P-Value Matrix.\n\n");
    }
    else if(!calPMatrix.IsEmpty()) {
        fprintf(stderr, "Calculate P-Value Matrix...\n");
        engine.CalculatePvalueMatrix();
        engine.WritePvalueMatrix(std::string(calPMatrix.c_str()));
        std::cerr<<"Pvalue Matrix calculated, next time you can specify parameter --PvalueMatrix [PATH] to skip calculation stage!"<<std::endl;
    }

    SetCrashExplanation("loading information of Genetic Map");

    if (!GDFile.IsEmpty()) {
        engine.GDMap.InputGeneticDistanceMap(std::string(GDFile.c_str()));
        engine.geneticMapAvailable = true;
	fprintf(stderr, "Done reading genetic map file.\n\n");
    }

    SetCrashExplanation("loading information of individuals");
    // Setup and load a list of individuals
    Pedigree ped;
    LoadPidToBeIncluded(String(""), pidIncludeFromPhased);
    LoadPidToBeExcluded(String(""), pidExcludeFromPhased);

    /*now loading phased individuals*/
    // here unphasedfile is the vcf file and is used for filling up the first five column of PED file(check the PED format).
    if (phasedfile != "Empty")
        LoadSamples(ped, phasedfile, pidIncludedInPhasedVcf, pidExcludedInPhasedVcf, engine.phased);
    fprintf(stderr,"Done loading %d phased individuals.\n\n",engine.phased);

    /*Notice that now we adding markers as subset of phased markers*/
    // here only extracted site information only, used for site check

    SetCrashExplanation("loading information for polymorphic sites");

    if (phasedfile != "Empty") {
        LoadRefPanelPolymorphicSites(phasedfile);
    } else
    {
        std::cerr<<"--refVCF is required" <<std::endl;
        exit(EXIT_FAILURE);
    }

    fprintf(stderr,"Done loading information on %d polymorphic sites.\n\n", ped.markerCount);

//    fprintf(stderr,"Processing input files and allocating memory for haplotyping\n");

    SetCrashExplanation("allocating memory for haplotype engine graph builder");

    engine.AllocateMemory(ped.count, ped.markerCount);

    SetCrashExplanation("loading genotype");

    LoadGenotypeAndHaplotypeFromPhasedVCF(ped, phasedfile, engine, errorRate, transRate);//this is where we copy GL into genotype arrays
    fprintf(stderr, "Done loading phased genotype file\n\n");


    SetCrashExplanation("building graph...");
    engine.ConstructGraph();
    fprintf(stderr, "Done building graph.\n\n");

    fprintf(stderr,"Total time:%.2f sec\n", (float) (clock() - t) / CLOCKS_PER_SEC);
    return 0;
}

int PhaseByRefGraph(int argc, char **argv) {

    String unphasedfile, mapfile, outfile("mach1.out"), phasedfile("Empty"), pidIncludeFromUnphased(
            ""), pidIncludeFromPhased(
            ""), pidExcludeFromUnphased(""), pidExcludeFromPhased(""), PMatrix(""),calPMatrix("");
    String crossFile, errorFile;

    clock_t t;
    t = clock();
    double errorRate = 0.01;
    double transRate = 0.01;
    int seed = 123456, warmup = 0, states = 0, weightedStates = 0;
    int burnin = 5, rounds = 10, polling = 0, samples = 0, samplingRounds = 0;
    int maxPhred = 255;

    int prefixLength = 120;

    bool compact = false;
    bool mle = false, mledetails = false, uncompressed = false;

    bool inputPhased = false;
    bool phaseByRef = false;
    bool randomPhase = false;
    bool fixTrans = true;

    bool onlyHeterSite = false;

    SetupCrashHandlers();
    SetCrashExplanation("reading command line options");

    printf("Pluto 0.01 -- Markov Chain Haplotyping for Shotgun Sequence Data\n"
                   "(c) 2015-2017 Fan Zhang, Goncalo Abecasis, and Hyun Min Kang\n\n");

    ParameterList pl;

    BEGIN_LONG_PARAMETERS(longParameters)
                    LONG_PARAMETER_GROUP("Input Files")
                    LONG_STRINGPARAMETER("unphasedVCF", &unphasedfile)
                    LONG_STRINGPARAMETER("refVCF", &phasedfile)
                    LONG_PARAMETER_GROUP("Optional Files")
                    LONG_STRINGPARAMETER("includeUnphasedIDs", &pidIncludeFromUnphased)
                    LONG_STRINGPARAMETER("includePhasedIDs", &pidIncludeFromPhased)
                    LONG_STRINGPARAMETER("excludeUnphasedIDs", &pidExcludeFromUnphased)
                    LONG_STRINGPARAMETER("excludePhasedIDs", &pidExcludeFromPhased)
                    LONG_STRINGPARAMETER("crossoverMap", &crossFile)
                    LONG_STRINGPARAMETER("errorMap", &errorFile)
                    LONG_STRINGPARAMETER("physicalMap", &mapfile)//TODO:decide which of these two, GD and physicalMap, to use
                    LONG_PARAMETER_GROUP("Markov Sampler")
                    LONG_INTPARAMETER("seed", &seed)
                    LONG_PARAMETER_GROUP("Haplotyper")
                    LONG_DOUBLEPARAMETER("errorRate", &errorRate)
                    LONG_DOUBLEPARAMETER("transRate", &transRate)
                    LONG_PARAMETER("fixTrans", &fixTrans)
                    LONG_PARAMETER("onlyHeterSite", &onlyHeterSite)
                    LONG_PARAMETER_GROUP("Phasing")
                    EXCLUSIVE_PARAMETER("randomPhase", &randomPhase)
                    EXCLUSIVE_PARAMETER("inputPhased", &inputPhased)
                    EXCLUSIVE_PARAMETER("refPhased", &phaseByRef)
                    LONG_PARAMETER_GROUP("Output Files")
                    LONG_STRINGPARAMETER("outPrefix", &outfile)
    END_LONG_PARAMETERS();

    pl.Add(new LongParameters("Available Options", longParameters));


    pl.Read(argc, argv);
    pl.Status();


    // Setup random seed ...
    globalRandom.Reset(seed);

    if (rounds < burnin) burnin = 0;

    PBWTHaplotyper engine;//declaration of engine, also will call default constructor
    engine.nSampleCopy = samplingRounds;
    engine.onlyHeterSite = onlyHeterSite;
    engine.geneticMapAvailable = false;
    engine.prefixLength = prefixLength;

    SetCrashExplanation("loading information of individuals");
    // Setup and load a list of individuals
    Pedigree ped;//ped file determined the dimension of haplotype and genotype matrix
    LoadPidToBeIncluded(pidIncludeFromUnphased, String(""));
    LoadPidToBeExcluded(pidExcludeFromUnphased, String(""));

    /*We add unphased individuals first*/
    int numUnphased(0);
    LoadSamples(ped, unphasedfile, pidIncludedInUnphasedVcf, pidExcludedInUnphasedVcf, numUnphased);
    std::cerr << "Load unphased individuals:" << numUnphased << std::endl;
    if (ped.count < 1) {
        error("SinglePhasing requires more than 0 sample.");
    }
    /*now loading phased individuals*/
    // here unphasedfile is the vcf file, here vcf is used for filling up the first five column of PED file(check the PED format).
    if (phasedfile != "Empty") {
        LoadSamples(ped, phasedfile, pidIncludedInPhasedVcf, pidExcludedInPhasedVcf,
                    engine.phasedForByRef);//keep engine.phased==0 while knowing the ref size
    }
    std::cerr << "Detected phased individuals in reference panel:" << engine.phasedForByRef << std::endl;

    /*Notice that now we adding markers as subset of phased markers*/
    // here only extracted site information only, used for site check

    SetCrashExplanation("loading information for polymorphic sites");

    if (phasedfile != "Empty")
    {
        LoadRefPanelPolymorphicSites(std::string(phasedfile.c_str()));
    } else
    {
        fprintf(stderr,"--refVCF is still required, please use the one that you built graph from\n");
        exit(EXIT_FAILURE);
    }

    fprintf(stderr,"Load information on %d polymorphic sites\n\n", Pedigree::markerCount);

    SetCrashExplanation("allocating memory for haplotype engine");

//    engine.EstimateMemoryInfo(ped.count, ped.markerCount, states, compact, false);
    engine.AllocateMemory(ped.count,ped.markerCount);
    engine.SetErrorAndTheta(errorRateHolder,transRateHolder);
//    engine.InitAuxillary();

    SetCrashExplanation("loading genotype");
    fprintf(stderr,"Copy unphased genotypes into haplotyping engine\n");
    // Copy genotypes into haplotyping engine
    if (engine.readyForUse)
        LoadGenotypeFromUnphasedVCF(ped, unphasedfile, engine, errorRate, transRate);//this is where we copy GL into genotype arrays

    fprintf(stderr,"Done loading unphased genotype file\n\n");
    // Copy phased haplotypes into haplotyping engine, but we put phased haps in the end

    SetCrashExplanation("searching for initial haplotype set");

    if (inputPhased) {
        fprintf(stderr,"Loading phased information from the input VCF file\n\n");
        engine.LoadHaplotypesFromVCF(unphasedfile);
//        engine.InitialSampleCopy(NULL);
    }
    else if (phaseByRef) {
        fprintf(stderr,"Assigning haplotypes based on reference genome\n\n");
        engine.PhaseByReferenceSetup();
//        engine.InitialSampleCopy(NULL);
    }
    else {
        fprintf(stderr,"Assigning random set of haplotypes\n\n");
        engine.RandomSetup(NULL);
//        engine.InitialSampleCopy(NULL);
    }
    fprintf(stderr,"Found initial haplotype set\n\n");

    SetCrashExplanation("haplotyping procedure");

    engine.loadGraph=std::string(phasedfile.c_str())+".DAG";

//            engine.LoopThroughChromosomesSingleRound();
    engine.LoopThroughChromosomesRecomb();


    SetCrashExplanation("outputing solution");
    fprintf(stderr, "In total we phased %d individuals at %d markers.\n", ped.count-engine.phasedForByRef, ped.markerCount);
    // If we did multiple rounds of haplotyping, then generate consensus
    {
        UnphasedSamplesOutputVCF(unphasedfile, ped, outfile + ".vcf.gz", thetas, error_rates, engine);
    }
    fprintf(stderr,"Total time:%.2f sec\n", (float) (clock() - t) / CLOCKS_PER_SEC);
    return 0;
}

int PhasingMain(int argc, char **argv) {


    String unphasedfile, mapfile, outfile("mach1.out"), phasedfile("Empty"), pidIncludeFromUnphased(
            ""), pidIncludeFromPhased(
            ""), pidExcludeFromUnphased(""), pidExcludeFromPhased(""), PMatrix(""),calPMatrix("");
    String crossFile, errorFile;
    String GDFile;

    clock_t t;
    t = clock();
    double errorRate = 0.01;
    double transRate = 0.01;
    int seed = 123456, warmup = 0, states = 0, weightedStates = 0;
    int burnin = 5, rounds = 10, polling = 0, samples = 0, samplingRounds = 1;
    int maxPhred = 255;

    int prefixLength = 120;

    bool compact = false;
    bool mle = false, mledetails = false, uncompressed = false;

    bool inputPhased = false;
    bool phaseByRef = false;
    bool randomPhase = false;
    bool fixTrans = true;

    bool isSingleRound = false;
    bool onlyHeterSite = false;

    SetupCrashHandlers();
    SetCrashExplanation("reading command line options");

    printf("Pluto 0.01 -- Markov Chain Haplotyping for Shotgun Sequence Data\n"
                   "(c) 2015-2017 Fan Zhang, Goncalo Abecasis, and Hyun Min Kang\n\n");

    ParameterList pl;

    BEGIN_LONG_PARAMETERS(longParameters)
                    LONG_PARAMETER_GROUP("Shotgun Sequences")
                    LONG_STRINGPARAMETER("unphasedVCF", &unphasedfile)
                    LONG_STRINGPARAMETER("refVCF", &phasedfile)
                    LONG_INTPARAMETER("maxPhred", &maxPhred)
                    LONG_PARAMETER_GROUP("Optional Files")
                    LONG_STRINGPARAMETER("includeUnphasedIDs", &pidIncludeFromUnphased)
                    LONG_STRINGPARAMETER("includePhasedIDs", &pidIncludeFromPhased)
                    LONG_STRINGPARAMETER("excludeUnphasedIDs", &pidExcludeFromUnphased)
                    LONG_STRINGPARAMETER("excludePhasedIDs", &pidExcludeFromPhased)
                    LONG_STRINGPARAMETER("crossoverMap", &crossFile)
                    LONG_STRINGPARAMETER("errorMap", &errorFile)
                    LONG_STRINGPARAMETER("geneticDistance", &GDFile)
                    LONG_STRINGPARAMETER("physicalMap", &mapfile)//decide which of these two, GD and physicalMap, to use
                    LONG_PARAMETER_GROUP("Graph Builder")
                    LONG_INTPARAMETER("graphComplexity", &prefixLength)
                    LONG_PARAMETER_GROUP("Markov Sampler")
                    LONG_INTPARAMETER("seed", &seed)
                    LONG_INTPARAMETER("burnin", &burnin)
                    LONG_INTPARAMETER("rounds", &rounds)
                    LONG_INTPARAMETER("samplingRounds", &samplingRounds)
                    LONG_PARAMETER_GROUP("Haplotyper")
                    LONG_INTPARAMETER("states", &states)
                    LONG_DOUBLEPARAMETER("errorRate", &errorRate)
                    LONG_DOUBLEPARAMETER("transRate", &transRate)
                    LONG_INTPARAMETER("weightedStates", &weightedStates)
                    LONG_PARAMETER("compact", &compact)
                    LONG_PARAMETER("fixTrans", &fixTrans)
                    LONG_PARAMETER("onlyHeterSite", &onlyHeterSite)
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
                    LONG_STRINGPARAMETER("PvalueMatrix", &PMatrix)
                    LONG_STRINGPARAMETER("calPvalueMatrix", &calPMatrix)
    END_LONG_PARAMETERS();

    pl.Add(new LongParameters("Available Options", longParameters));

    pl.Add(new HiddenString('m', "Map File", mapfile));
    pl.Add(new HiddenString('o', "Output File", outfile));
    pl.Add(new HiddenInteger('r', "Haplotyping Rounds", rounds));
    pl.Add(new HiddenDouble('e', "Error Rate", errorRate));

    pl.Read(argc, argv);
    pl.Status();


    // Setup random seed ...
    globalRandom.Reset(seed);




    if (rounds < burnin) burnin = 0;

    PBWTHaplotyper engine;//declaration of engine, also will call default constructor
    engine.nSampleCopy = samplingRounds;
    engine.onlyHeterSite = onlyHeterSite;
    engine.geneticMapAvailable = false;
    engine.prefixLength = prefixLength;

    SetCrashExplanation("loading Pvalue Matrix");

    if(PMatrix.IsEmpty() and calPMatrix.IsEmpty()) {
        std::cerr<<"parameter --PvalueMatrix [PATH] or --calPvalueMatrix [PATH] required!"<<std::endl;
        exit(EXIT_FAILURE);
    }
    else if(!PMatrix.IsEmpty())
        engine.ReadPvalueMatrix(std::string(PMatrix.c_str()));
    else if(!calPMatrix.IsEmpty()) {
        engine.CalculatePvalueMatrix();
        engine.WritePvalueMatrix(std::string(calPMatrix.c_str()));
        std::cerr<<"Pvalue Matrix calculated, next time you can specify parameter --PvalueMatrix [PATH] to skip calculation stage!"<<std::endl;
    }

    SetCrashExplanation("loading information of individuals");
    // Setup and load a list of individuals
    Pedigree ped;
    LoadPidToBeIncluded(pidIncludeFromUnphased, pidIncludeFromPhased);
    LoadPidToBeExcluded(pidExcludeFromUnphased, pidExcludeFromPhased);

    /*We add unphased individuals first*/
    int numUnphased(0);
    LoadSamples(ped, unphasedfile, pidIncludedInUnphasedVcf, pidExcludedInUnphasedVcf, numUnphased);
    std::cerr << "Load unphased individuals:" << numUnphased << std::endl;
    if (ped.count < 1) {
        error("SinglePhasing requires more than 0 sample.");
    }
    /*now loading phased individuals*/
    // here unphasedfile is the vcf file, here vcf is used for filling up the first five column of PED file(check the PED format).
    if (phasedfile != "Empty")
        LoadSamples(ped, phasedfile, pidIncludedInPhasedVcf, pidExcludedInPhasedVcf, engine.phased);
    std::cerr << "Load phased individuals:" << engine.phased << std::endl;

    /*Notice that now we adding markers as subset of phased markers*/
    // here only extracted site information only, used for site check

    SetCrashExplanation("loading information for polymorphic sites");

    if (phasedfile != "Empty")
    {
        LoadRefPanelPolymorphicSites(phasedfile);
        LoadUnphasedPolymorphicSites(unphasedfile);
    }
    else
    {
        LoadRefPanelPolymorphicSites(unphasedfile);
        LoadUnphasedPolymorphicSites(unphasedfile);
    }

    fprintf(stderr,"Load information on %d polymorphic sites\n\n", Pedigree::markerCount);

    SetCrashExplanation("loading information of Genetic Map");

    if (!GDFile.IsEmpty()) {
        engine.GDMap.InputGeneticDistanceMap(std::string(GDFile.c_str()));
        engine.geneticMapAvailable = true;
    }

    Pedigree::LoadMarkerMap(mapfile);//the format of mapfiles is:	chrome\tmarker_name\tposition


    // Check if physical map is available
    bool positionsAvailable = true;

    for (int i = 0; i < ped.markerCount; i++)
        if (Pedigree::GetMarkerInfo(i)->chromosome < 0) {
            positionsAvailable = false;//no physical map available
            break;
        }

    if (positionsAvailable) {
        printf("    Physical map will be used to improve crossover rate estimates.\n");

        for (int i = 1; i < ped.markerCount; i++)
            if (ped.GetMarkerInfo(i)->position <= ped.GetMarkerInfo(i - 1)->position ||
                ped.GetMarkerInfo(i)->chromosome != ped.GetMarkerInfo(i - 1)->chromosome) {
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
    engine.ShotgunHaplotyper::AllocateMemory(ped.count, states, ped.markerCount, (float) transRate);
    engine.InitAuxillary();

    SetCrashExplanation("loading genotype");
    fprintf(stderr,"Copy unphased genotypes into haplotyping engine\n");
    // Copy genotypes into haplotyping engine
    if (engine.readyForUse)
        LoadGenotypeFromUnphasedVCF(ped, unphasedfile, maxPhred, engine);//this is where we copy GL into genotype arrays

    fprintf(stderr,"Done loading unphased genotype file\n\n");
    // Copy phased haplotypes into haplotyping engine, but we put phased haps in the end


    if (phasedfile!="Empty") {
        fprintf(stderr, "Copy phased genotypes into haplotyping engine\n");
        LoadGenotypeAndHaplotypeFromPhasedVCF(ped, phasedfile, maxPhred, engine.phased, engine, errorRate,
                                              transRate);//this is where we copy GL into genotype arrays
        fprintf(stderr, "Done loading phased genotype file\n\n");
    }

    //TODO:decide which version to use
    if (positionsAvailable &&
        engine.AllocateDistances())//notice that there are two position information sources, one is from VCF the other is from markerMap
    {
        for (int i = 1; i < ped.markerCount; i++)//here the distance is based on markerMap file
            engine.distances[i - 1] = ped.GetMarkerInfo(i)->position -
                                      ped.GetMarkerInfo(i - 1)->position;
    }

    engine.ShowMemoryInfo();


    int ConsensusBuilderRounds=0;
    if(rounds-burnin==1)
    {
        isSingleRound=true;
        ConsensusBuilderRounds=samplingRounds;
    }
    else ConsensusBuilderRounds=rounds-burnin;
    ConsensusBuilder::EstimateMemoryInfo(ConsensusBuilderRounds, (ped.count-engine.phased) * 2, ped.markerCount);
    ConsensusBuilder consensus(ConsensusBuilderRounds, (ped.count-engine.phased) * 2, ped.markerCount);

    if (consensus.readyForUse == false)
        return MemoryAllocationFailure();

    printf("Memory allocated successfully\n\n");

    SetCrashExplanation("loading error rate and cross over maps");

    //TODO:decide which version to use
    bool newline = engine.LoadCrossoverRates(crossFile);
    newline |= engine.LoadErrorRates(errorFile);
    if (newline) printf("\n");

    SetCrashExplanation("searching for initial haplotype set");

    if (inputPhased) {
        printf("Loading phased information from the input VCF file\n\n");
        engine.LoadHaplotypesFromVCF(unphasedfile);
        engine.InitialSampleCopy(NULL);
    }
    else if (phaseByRef) {
        printf("Assigning haplotypes based on reference genome\n\n");
        engine.PhaseByReferenceSetup();
        engine.InitialSampleCopy(NULL);
    }
    else {
        printf("Assigning random set of haplotypes\n\n");
        engine.RandomSetup(NULL);
        engine.InitialSampleCopy(NULL);
    }
    printf("Found initial haplotype set\n\n");


    SetCrashExplanation("phasing procedure");

//    engine.loadGraph="HG00535.phased.DAG";

        if (isSingleRound)
//            engine.LoopThroughChromosomesSingleRound();
            engine.LoopThroughChromosomesRecomb();
        else
            engine.LoopThroughChromosomesHighPrecision();

        if (!fixTrans) engine.UpdateThetas();
        errorRate = engine.UpdateErrorRate();

        if (OutputManager::outputHaplotypes) {
            if (isSingleRound)
                consensus.StoreForSingleRound(engine.sampledHaps, engine.nSampleCopy);
            else
                consensus.Store(engine.haplotypes);
        }


    SetCrashExplanation("outputing solution");
    // If we did multiple rounds of haplotyping, then generate consensus
    {
        UnphasedSamplesOutputVCF(unphasedfile, ped, outfile + ".vcf.gz", thetas, error_rates, engine);
        if (OutputManager::outputHaplotypes)
            OutputVCFConsensus(unphasedfile, ped, consensus, outfile + ".consensus.vcf.gz", thetas, error_rates,
                               engine);
    }
    printf("Total time:%.2f sec\n", (float) (clock() - t) / CLOCKS_PER_SEC);
    return 0;
}
