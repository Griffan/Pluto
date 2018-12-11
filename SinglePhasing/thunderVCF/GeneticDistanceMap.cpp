//
// Created by Fan Zhang on 7/29/16.
//

#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include "GeneticDistanceMap.h"
#include <algorithm>

void GeneticDistanceMap::InputGeneticDistanceMap(const std::string &inputFile) {
    std::ifstream FIN(inputFile);
    if(!FIN.is_open())
    {
        std::cerr<<"Open file "<<inputFile<< " falied! Exit."<<std::endl;
        exit(EXIT_FAILURE);
    }
    std::string line,chr;
    int start;
    float gen;
    getline(FIN,line);//discard header
    while(getline(FIN,line))
    {
        std::stringstream ss(line);
        ss>>chr>>start>>gen;
        chrBpGeneticDistance[std::string("chr")+chr][start]=gen;
    }
}

float GeneticDistanceMap::InferGeneticDistance(std::string chr, int start) {

    if(chr.find("chr")==chr.npos)
    {
        chr=std::string("chr")+chr;
    }
    if(chrBpGeneticDistance.find(chr)==chrBpGeneticDistance.end())
    {
        std::cerr<<"Unknown chromsome ID:"<<chr<<" for genetic distance map file!"<<std::endl;
        exit(EXIT_FAILURE);
    }

    std::map<int, float>::iterator before;
    before = chrBpGeneticDistance[chr].upper_bound(start);
    if(chrBpGeneticDistance[chr].begin() == before) return before->second;
    else if(chrBpGeneticDistance[chr].end()==before ) return (--before)->second;
    else {
        auto after = before;
        before--;
//        std::cerr<<"start:"<<start<<"\tbefore.first:"<<before->first
//                 <<"\tafter.second:"<<after->second<<"\tbefore.second:"<<before->second
//                 <<"\tafter.first:"<<after->first <<"\tbefore.first:"<<before->first<<std::endl;

        return (start - before->first) * (after->second - before->second) / (after->first - before->first)+before->second;
    }
}

float GeneticDistanceMap::CalculateRecombinationRate(float prev, float current) {
    const float factor=2.0f;
    prev/=100;
    current/=100;
    return std::max(1.f/powf((1.f+prev/factor),2.f)-1.f/powf((1.f+current/factor),2.f),10e-6f);
}


