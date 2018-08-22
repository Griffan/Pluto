//
// Created by Fan Zhang on 8/10/18.
//

#include "DAG.h"
#include <fstream>

int StateNode::WriteNode(std::ofstream &fout) {
    int totalSize = 0;
    fout.write((char *) &allele, sizeof(uchar));
    totalSize += 2 * sizeof(uchar);
    fout.write((char *) &childNode, 2 * sizeof(StateIndex));
    totalSize += 2 * sizeof(StateIndex);
    fout.write((char *) &numHap, 2 * sizeof(float));
    totalSize += 2 * sizeof(float);
    unsigned long sizeOfSet = parentNodeSet.size();
    fout.write((char *) &sizeOfSet, sizeof(unsigned long));
    totalSize += sizeof(unsigned long);
    for (auto k:parentNodeSet) {
        fout.write((char *) &k, sizeof(StateIndex));
        totalSize += sizeof(StateIndex);
    }
    return totalSize;
}

int StateNode::ReadNode(std::ifstream &fin) {
    fin.read((char *) &allele, sizeof(uchar));
    fin.read((char *) &childNode, 2 * sizeof(StateIndex));
    fin.read((char *) &numHap, 2 * sizeof(float));
    unsigned long sizeOfSet;
    fin.read((char *) &sizeOfSet, sizeof(unsigned long));
    StateIndex k;
    for (int i = 0; i != sizeOfSet; i++) {
        fin.read((char *) &k, sizeof(StateIndex));
        AddParentNode(k);
    }
    return 0;
}

DAG::DAG(DAG &A) {
    nsnps = A.nsnps;
    nhaps = A.nhaps;
    StateNodeMat = A.StateNodeMat;
//        tmpNodeVec=A.tmpNodeVec;
}

DAG::DAG() {
    nsnps = 0;
    nhaps = 0;
}

DAG::~DAG() {
    for (auto &vec:StateNodeMat) {
        for (auto &node:vec) {
            delete node;
        }
        vec.clear();
    }
}

DAG::DAG(int nmarkers, int nHaps) : nsnps(nmarkers), nhaps(nHaps),
                                    StateNodeMat(nmarkers, std::vector<StateNode *>(0, new StateNode(0))) {
}

int DAG::JoinNodes(int marker, StateIndex indexRetain, StateIndex indexRemove) {
    StateNodeMat[marker][indexRetain]->SetNumHap(0, StateNodeMat[marker][indexRetain]->GetNumHap(0) +
                                                    StateNodeMat[marker][indexRemove]->GetNumHap(0));
    StateNodeMat[marker][indexRetain]->SetNumHap(1, StateNodeMat[marker][indexRetain]->GetNumHap(1) +
                                                    StateNodeMat[marker][indexRemove]->GetNumHap(1));

    for (auto kv:StateNodeMat[marker][indexRemove]->GetParentIndexSet()) {
        StateNodeMat[marker][indexRetain]->AddParentNode(kv);
        if (StateNodeMat[marker][indexRemove]->GetAllele() == 0)//Remove has 0 allele
        {
            StateNodeMat[marker - 1][kv]->SetChildNodeIndex(0, indexRetain);
        } else if (StateNodeMat[marker][indexRemove]->GetAllele() == 1)//Remove has 1 allele
        {
            StateNodeMat[marker - 1][kv]->SetChildNodeIndex(1, indexRetain);
        } else {
            fprintf(stderr, "roar from StateNode operator+=!!!!\n%d and %d: %d\n",
                    StateNodeMat[marker - 1][kv]->GetChildNodeIndex(0),
                    StateNodeMat[marker - 1][kv]->GetChildNodeIndex(1), indexRemove);
            exit(EXIT_FAILURE);
        }
    }
    return 0;
}

void DAG::UpdateChildNodeIndex(int marker, StateIndex parentIndex, StateIndex newChildIndex, char allele) {
    StateNodeMat[marker][parentIndex]->SetChildNodeIndex(allele, newChildIndex);
}

float DAG::GetProbToCurrentNodeConditionalOnParentNode(int marker, StateIndex parentIndex, char allele) {
    return StateNodeMat[marker][parentIndex]->GetNumHap(allele) /
           (StateNodeMat[marker][parentIndex]->GetNumHap(0) + StateNodeMat[marker][parentIndex]->GetNumHap(1));
}

float DAG::GetEdgeProbFromParentNode(int marker, StateIndex parentIndex, char allele) {
    return StateNodeMat[marker][parentIndex]->GetNumHap(allele) / nhaps;
}

int DAG::WriteDAG(const std::string &fileName) {
    int totalSize = 0;
    std::ofstream fout(fileName, std::ifstream::binary);

    if (!fout.is_open()) {
        std::cerr << "open file " << fileName << " failed!" << std::endl;
        exit(EXIT_FAILURE);
    }
    fout.write((char *) &nhaps, sizeof(int));
    totalSize += sizeof(int);
    fout.write((char *) &nsnps, sizeof(int));
    totalSize += sizeof(int);
    for (int i = 0; i < nsnps; ++i) {
        unsigned long currentSiteNodeNum = StateNodeMat[i].size();
        fout.write((char *) &currentSiteNodeNum, sizeof(unsigned long));
        totalSize += sizeof(unsigned long);
        for (int j = 0; j < currentSiteNodeNum; ++j) {
            totalSize += StateNodeMat[i][j]->WriteNode(fout);
        }
    }
    fout.close();

//        std::ofstream fout2(fileName+".txt");
//        fout2<<"nhaps:"<<nhaps<<std::endl;
//        fout2<<"nsnps"<<nsnps<<std::endl;
//        for (int i = 0; i <nsnps ; ++i) {
//            fout2<<i<<"th snp nodeVec size:"<<StateNodeMat[i].size()<<std::endl;
//            for (int j = 0; j <StateNodeMat[i].size(); ++j) {
//                fout2<<j<<"th node:"<<StateNodeMat[i][j]->ToString()<<std::endl;
//            }
//        }
//        fout2<<std::endl;
//        fout2<<std::endl;
//        fout2.close();

    return totalSize;
}

int DAG::ReadDAG(const std::string &fileName) {
    std::cerr << "Reading graph from [ReadDAG]..." << std::endl;
    std::ifstream fin(fileName, std::ifstream::binary);

    if (!fin.is_open()) {
        std::cerr << "open file " << fileName << " failed!" << std::endl;
        exit(EXIT_FAILURE);
    }
    fin.read((char *) &nhaps, sizeof(int));
    fin.read((char *) &nsnps, sizeof(int));
    for (int i = 0; i < nsnps; ++i) {
        unsigned long currentSiteNodeNum(0);
        fin.read((char *) &currentSiteNodeNum, sizeof(unsigned long));
        for (int j = 0; j < currentSiteNodeNum; ++j) {
            StateNode *tmpNode = new StateNode(255);
            tmpNode->ReadNode(fin);
            StateNodeMat[i].push_back(tmpNode);
        }
    }
    fin.close();

//        std::ofstream fout2(fileName+".txt");
//        fout2<<"nhaps:"<<nhaps<<std::endl;
//        fout2<<"nsnps"<<nsnps<<std::endl;
//        for (int i = 0; i <nsnps ; ++i) {
//            fout2<<i<<"th snp node size:"<<StateNodeMat[i].size()<<std::endl;
//            for (int j = 0; j <StateNodeMat[i].size(); ++j) {
//                fout2<<j<<"th node:"<<StateNodeMat[i][j]->ToString()<<std::endl;
//            }
//        }
//        fout2<<std::endl;
//        fout2<<std::endl;
//        fout2.close();
    return 0;
}

/*
 *
{
"nodes": [
{"blahbalbal": "d3","id":5},
{"name": "d3.svg"},
{"name": "d3.svg.area"},
{"name": "d3.svg.line"},
{"name": "d3.scale"},
{"name": "d3.scale.linear"},
{"name": "d3.scale.ordinal"}
],
"links": [
{"source": 0, "target": 1, "weight":2},
{"source": 1, "target": 2},
{"source": 1, "target": 3},
{"source": 0, "target": 4},
{"source": 4, "target": 5},
{"source": 4, "target": 6}
]
}
 */
int DAG::ToJson(const std::string &fileName) {
    int minorIndex = 0;
    int majorIndex = 0;
    std::ofstream fout(fileName);

    if (!fout.is_open()) {
        std::cerr << "open file " << fileName << " failed!" << std::endl;
        exit(EXIT_FAILURE);
    }
    int nMarkers = nsnps;
    fout << "{\n"
         << "  \"nodes\": [" << std::endl;
    for (int i = 0; i </*nsnps*/nMarkers; ++i) {
        unsigned long currentSiteNodeNum = StateNodeMat[i].size();
        for (int j = 0; j < currentSiteNodeNum; ++j) {
            if (i == nMarkers - 1 and j == currentSiteNodeNum - 1)
                fout << "{\"name\": \"" << i << "." << j << "\"}"
                     << std::endl;
            else
                fout << "{\"name\": \"" << i << "." << j << "\"},"
                     << std::endl;
        }
    }
    fout << "  ]," << std::endl;
    fout << "  \"links\": [" << std::endl;
    for (int i = 0; i </*nsnps-1*/nMarkers - 1; ++i) {
        unsigned long currentSiteNodeNum = StateNodeMat[i].size();
        majorIndex += currentSiteNodeNum;
        for (int j = 0; j < currentSiteNodeNum; ++j) {
            if (StateNodeMat[i][j]->GetChildNodeIndex(0) != -1)
                fout << "{\"source\": " << minorIndex
                     << ", \"target\": " << majorIndex + StateNodeMat[i][j]->GetChildNodeIndex(0)
                     << ", \"weight\": " << StateNodeMat[i][j]->GetNumHap(0)
                     << ", \"allele\": " << 0
                     << "}," << std::endl;
            if (i == nMarkers - 2 and j == currentSiteNodeNum - 1)
                fout << "{\"source\": " << minorIndex
                     << ", \"target\": " << majorIndex + StateNodeMat[i][j]->GetChildNodeIndex(1)
                     << ", \"weight\": " << StateNodeMat[i][j]->GetNumHap(1)
                     << ", \"allele\": " << 1
                     << "}" << std::endl;
            else if (StateNodeMat[i][j]->GetChildNodeIndex(1) != -1)
                fout << "{\"source\": " << minorIndex
                     << ", \"target\": " << majorIndex + StateNodeMat[i][j]->GetChildNodeIndex(1)
                     << ", \"weight\": " << StateNodeMat[i][j]->GetNumHap(1)
                     << ", \"allele\": " << 1
                     << "}," << std::endl;
            minorIndex++;
        }
    }
    fout << "  ]\n}" << std::endl;

    fout.close();
    return minorIndex;
}

int DAG::FromJson(const std::string &fileName) {
    std::cerr << "Reading graph from [FromJson]..." << std::endl;
    std::ifstream fin(fileName);

    if (!fin.is_open()) {
        std::cerr << "open file " << fileName << " failed!" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string line;
    std::getline(fin, line);//fin<<"{\n"
    std::getline(fin, line);//<<"  \"nodes\": ["<<std::endl;
    /*
    for (int i = 0; i <nMarkers; ++i) {
        unsigned long currentSiteNodeNum = StateNodeMat[i].size();
        for (int j = 0; j <currentSiteNodeNum; ++j) {
            if(i == nMarkers - 1 and j == currentSiteNodeNum - 1)
                fin<<"{\"name\": \""<<i<<"."<<j<<":"<<(int)StateNodeMat[i][j]->GetAllele()<<"\"}"<<std::endl;
            else
                fin<<"{\"name\": \""<<i<<"."<<j<<":"<<(int)StateNodeMat[i][j]->GetAllele()<<"\"},"<<std::endl;
        }
    }*/
    std::vector<std::pair<StateIndex, StateNode *> > nodeVec;
    long curMarker(0);
    StateIndex curIndex(0);
    unsigned long pivot(0);
    unsigned long len(0);
    char allele(-1);
    while (getline(fin, line) and line.find("name") != std::string::npos) {
        pivot = line.find(':');
        len = line.find('.', pivot + 2) - (pivot + 3);
        curMarker = std::stol(line.substr(pivot + 3, len));
        pivot = line.find('.', pivot + 1);
        len = line.find(':', pivot + 1) - (pivot + 1);
        curIndex = static_cast<StateIndex >(std::stoi(line.substr(pivot + 1, len)));

        StateNode *tmpNode = new StateNode(-1);
        nodeVec.emplace_back(curIndex, tmpNode);
        assert(curIndex == StateNodeMat[curMarker].size());
        StateNodeMat[curMarker].push_back(tmpNode);

    }
//    std::getline(fin, line);//fin<<"  ],"<<std::endl;
    std::getline(fin, line);//fin<<"  \"links\": ["<<std::endl;

    long sourceIndex(0);
    long targetIndex(0);
    float weight(0);
    while (getline(fin, line) and line.find("source") != std::string::npos) {
        pivot = line.find(':');
        len = line.find(',', pivot + 1) - (pivot + 1);
        sourceIndex = std::stol(line.substr(pivot + 1, len));
        pivot = line.find(':', pivot + 1);
        len = line.find(',', pivot + 1) - (pivot + 1);
        targetIndex = std::stol(line.substr(pivot + 1, len));
        pivot = line.find(':', pivot + 1);
        len = line.find(',', pivot + 1) - (pivot + 1);
        weight = std::stof(line.substr(pivot + 1, len));
        pivot = line.find(':', pivot + 1);
        len = line.find('}', pivot + 1) - (pivot + 1);
        allele = std::stoi(line.substr(pivot + 1, len));

        if (sourceIndex >= 0) {
            if (nodeVec[sourceIndex].second->GetNumHap(allele) < weight/*work around for graph inconsistency*/) {
                nodeVec[sourceIndex].second->SetChildNodeIndex(allele, nodeVec[targetIndex].first);
                nodeVec[sourceIndex].second->SetNumHap(allele, weight);
                nodeVec[targetIndex].second->AddParentNode(nodeVec[sourceIndex].first);
                nodeVec[targetIndex].second->SetAllele(allele);
            }
        }
        else
            nodeVec[targetIndex].second->SetAllele(allele);
    }
    std::getline(fin, line);//fin<<"  ]\n}"<<std::endl;
    fin.close();
    return 0;
}