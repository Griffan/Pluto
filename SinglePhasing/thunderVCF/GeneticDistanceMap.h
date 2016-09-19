//
// Created by Fan Zhang on 7/29/16.
//

#ifndef PLUTO_GENETICDISTANCEMAP_H
#define PLUTO_GENETICDISTANCEMAP_H

#include <map>

class GeneticDistanceMap {
    std::map<std::string,std::map<int, float> > chrBpGeneticDistance;
public:
    void InputGeneticDistanceMap(const std::string & inputFile);
    float InferGeneticDistance(std::string chr, int start);
    float CalculateRecombinationRate(float prev, float current);
};


#endif //PLUTO_GENETICDISTANCEMAP_H
