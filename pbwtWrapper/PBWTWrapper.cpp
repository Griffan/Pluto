//
// Created by Fan Zhang on 7/20/15.
//


#include <cmath>
#include "PBWTWrapper.h"
#include "iostream"
#include <fstream>
#include <functional>
#include <cstring>

const float T_CRITICAL_VALUE[] =
        {12.71, 4.3, 3.18, 2.78, 2.57,
         2.45, 2.37, 2.31, 2.26, 2.23,
         2.2, 2.18, 2.16, 2.15, 2.13,
         2.12, 2.11, 2.1, 2.09, 2.09,
         2.08, 2.07, 2.07, 2.06, 2.06,
         2.06, 2.05, 2.05, 2.05, 2.04,/*first 30*/
         2.03, 2.02, 2.01, 2.01, 2,
         2, 2, 1.99, 1.99, 1.99, 1.99,
         1.99, 1.99, 1.98,/*every 5 degrees untill 100*/
         1.97,/*200*/ 1.97,/*500*/ 1.96/*infinity*/
        };

float P_thresh = 0.01;

PBWTWrapper::PBWTWrapper(int nhaps, int nsnps) : Graph(nsnps, nhaps) {
    nSamples = nhaps / 2;
    nMarkers = nsnps;
    N = nsnps;
    M = nhaps;
    haplotype = nullptr;
}

PBWTWrapper::PBWTWrapper(int nhaps, int nsnps, float ***t_PvalueMatrix, int prefixLen)
        : prefixLength(prefixLen), Graph(nsnps, nhaps),
          a(nhaps, 0), alpha(a),
          alphaMap(nsnps, a),
          d(nhaps, 0),
          delta(d),
          allDelta(alphaMap),
          bkDistance(nsnps, std::vector<float>(nhaps,0.f)),
          sortedY(a),
          c(nsnps, 0), celta(c),
          haplotypeCluster(nsnps,std::vector<StateIndex>(nhaps, 0)),
          mergePairList(std::function<bool(const MaxPair &,
                                     const MaxPair &)>(comparator)) {
    nSamples = nhaps / 2;
    nMarkers = nsnps;
    N = nsnps;
    M = nhaps;//last two haps are slots for current individual need to be phased
    PvalueMatrix = t_PvalueMatrix;
    //cerr<<"Inside PBWTWrapper M:"<<M<<endl;
//    pbwtCore = pbwtCreate(nhaps, nsnps);
    //pbwtCore->CompressedAllele = arrayCreate(4096 * 32, uchar);

//    forwardCursor = pbwtCursorCreate(pbwtCore, TRUE, TRUE);

//    reverseCursor = pbwtCursorCreate(pbwtCore, FALSE, TRUE);
    for (int i = 0; i < M; ++i) {
        a[i] = i;
    }
    alpha = a;

    haplotype = nullptr;

}


int PBWTWrapper::CursorForwards() {//so far only implemented for test purpose


    //PrintVector(forwardCursor->a,M,"end arrary aFend check 0");
    for (int k = 0; k != N; ++k) {
//        fprintf(stderr,"at site %d\n",k);
        CursorForwardsTo(k, prefixLength);
//        PrintNonZeroVector(d, M, "end arrary aFend check 1");

    }
    //copy end of a to PBWT
    //PrintVector(forwardCursor->a,M,"end arrary aFend check 1");

    PrintSummary();
    //update crossover rate?
    return 0;
}

/*
int PBWTWrapper::CursorForwardsToDeprecated(int k, int T) {

    int rank, i0 = 0, ia;
    int group = 0;

    int hapID(0), prevSiteStateIndex(0);

    clusterMembership.clear();
    dist.clear();
    rightCoordinate.clear();
    rightCoordinateStat.clear();


    int tmpNumHap(0);
    char allele(-1);

    int tmpT = k > T ? T : k;



    //cluster of the previous site k-1
    for (rank = 0; rank < forwardCursor->M; ++rank) {

        if (forwardCursor->d[rank] >
            (k - tmpT) &&
            (k != 0)) {//new cluster if current sequence and last sequence have common sequence less than T
            if (rank != 0) {
                tmpNumHap = rank - i0;
//                if(k-1>1&&k-1<12&&rank>3220&&rank<3242){printf("now site %d rank:%d d:%d tmpT:%d allle:%d a:%d new cluster detected\n",
//                                                               k,rank,forwardCursor->d[rank],tmpT,forwardCursor->sortedY[rank],forwardCursor->a[rank]);}

                allele = haplotype[GetHapIDFromFwd(k - 1, i0)][k - 1];

                Graph.StateNodeMat[k - 1].push_back(new StateNode(group, tmpNumHap, allele));
                Graph.StateNodeMat[k - 1][group]->nodeIndex = group;

                dist.push_back(std::vector<int>(tmpNumHap, 0));//state->rank_occupied
                rightCoordinate.push_back(std::vector<float>(tmpNumHap, (float) 0.));//state->rank_occupied
                std::vector<int> tmpMem;
                for (ia = i0; ia < rank; ++ia) {
                    hapID = GetHapIDFromFwd(k - 1, ia);//original ID, by treating rank as backward ranking
                    haplotypeCluster[k - 1][hapID] = group;
                    tmpMem.push_back(ia);
                    dist[group][ia - i0] = GetRankFromBack(k - 1, hapID);
                    rightCoordinate[group][ia - i0] = GetDistanceFromBack(k - 1, dist[group][ia - i0]);
                    if (k >= 2) {
                        prevSiteStateIndex = GetHapStateFromFwd(k - 2, hapID);
                        Graph.StateNodeMat[k - 1][group]->AddParentNode(prevSiteStateIndex, Graph.StateNodeMat[k - 2][prevSiteStateIndex]);//site k-1
                        Graph.StateNodeMat[k - 2][prevSiteStateIndex]->AddChildNode(allele, &(Graph.StateNodeMat[k - 1][group]->nodeIndex), hapID <= M ? 1.f : 1e-6f);
                    }
                }

                if (k == 1) {
                    Graph.StateNodeMat[0][group]->AddParentNode(0, nullptr);//site 1
                    Graph.StateNodeMat[0][group]->SetID(0,0);
                }
                else
                    Graph.StateNodeMat[k - 1][group]->SetID(prevSiteStateIndex,Graph.StateNodeMat[k - 2][prevSiteStateIndex]->ID);

                Graph.RegisterState(k-1, Graph.StateNodeMat[k - 1][group]->ID, &(Graph.StateNodeMat[k - 1][group]->nodeIndex));

                rightCoordinateStat.push_back(rightCoordinate[group]);
                std::sort(dist[group].begin(),
                          dist[group].end());//TODO:only sort nodes need to be compared, skip those otherwise
//                clusterAllele[k - 1].push_back(allele);
                clusterMembership.push_back(tmpMem);

                i0 = rank;
                group++;
            }
        }
    }
    //finish the last segment if i0 didn't reach the end
    if (i0 < forwardCursor->M) {
        if (k != 0) {
            tmpNumHap = forwardCursor->M - i0;

            allele = haplotype[GetHapIDFromFwd(k - 1, i0)][k - 1];
//            if (allele != 0 and allele != 1) {
//                fprintf(stderr, "last state marker:%d\ti0:%d\tgroup:%d\tID:%d\n", k - 1, i0, group,
//                        GetHapIDFromFwd(k - 1, i0));
//                exit(EXIT_FAILURE);
//            }
            Graph.StateNodeMat[k - 1].push_back(new StateNode(group, tmpNumHap, allele));
            Graph.StateNodeMat[k - 1][group]->nodeIndex = group;

            dist.push_back(std::vector<int>(tmpNumHap, 0));//state->rank_occupied
            rightCoordinate.push_back(std::vector<float>(tmpNumHap, (float) 0.));//state->rank_occupied
            std::vector<int> tmpMem;
            for (ia = i0; ia < forwardCursor->M; ++ia) {
                hapID = GetHapIDFromFwd(k - 1, ia);//original ID, by treating rank as backward ranking
                haplotypeCluster[k - 1][hapID] = group;
                tmpMem.push_back(ia);
                dist[group][ia - i0] = GetRankFromBack(k - 1, hapID);
                rightCoordinate[group][ia - i0] = GetDistanceFromBack(k - 1, dist[group][ia - i0]);
                if (k >= 2) {
                    prevSiteStateIndex = GetHapStateFromFwd(k - 2, hapID);
                    Graph.StateNodeMat[k - 1][group]->AddParentNode(prevSiteStateIndex,
                                                                    Graph.StateNodeMat[k - 2][prevSiteStateIndex]);
                    Graph.StateNodeMat[k - 2][prevSiteStateIndex]->AddChildNode(allele, &(Graph.StateNodeMat[k - 1][group]->nodeIndex),
                                                                                hapID <= M ? 1.f : 1e-6f);
                }
            }

            if (k == 1) {
                Graph.StateNodeMat[0][group]->AddParentNode(0, nullptr);//site 1
                Graph.StateNodeMat[0][group]->SetID(0,0);
            }
            else
                Graph.StateNodeMat[k - 1][group]->SetID(prevSiteStateIndex,Graph.StateNodeMat[k - 2][prevSiteStateIndex]->ID);

            Graph.RegisterState(k-1, Graph.StateNodeMat[k - 1][group]->ID, &(Graph.StateNodeMat[k - 1][group]->nodeIndex));

            rightCoordinateStat.push_back(rightCoordinate[group]);
            std::sort(dist[group].begin(), dist[group].end());
//            clusterAllele[k - 1].push_back(allele);
            clusterMembership.push_back(tmpMem);
        }
    }
//    if(k>=2) {
//        std::cerr << "site:" << k - 2 << ":\t ";
//        HowManyChildlessState(Graph.StateNodeMat[k - 2]);
//        HowManyChildHapCount(Graph.StateNodeMat[k - 2], k - 1);
//    }
    //merge cluster based on KS test
    int merged = 0;
#ifdef DEBUG
    if(k!=0){
        fprintf(stderr, "at site:%d\n", k-1);
        PrintVector(haplotypeCluster[k-1],"haplotype state before merge state");
        PrintVector(clusterAllele[k-1],"state allele before merge allele");
    }
#endif
    if (k >= 1) {


        if (GetNumStates(k-1) != 1) merged = RegressionMergeAtSite(k - 1, true);//TODO:implement this function
//        if (merged) std::cerr<<"merged at site:"<<k-1<<std::endl;
        Graph.NormalizeCurrentSiteTransitionProb(k - 2);

        //if(k!=0&&clusterAllele[k-1].size()!=1) test = MergeAtSiteExperiment(k-1);
#ifdef DEBUG
        if(test) {
            fprintf(stderr, "at site:%d\n", k-1);
            PrintVector(haplotypeCluster[k-1], "haplotype state after merge state");
            PrintVector(clusterAllele[k-1],"state allele after merge allele");
            fprintf(stderr,"\n");
        }
#endif
//        UpdateTransVector(k - 1);


    }
//    std::cerr<<"site:"<<k-2<<"After merge:";
//    int tmpNum=HowManyChildlessState(Graph.StateNodeMat[k - 2]);
//    if(tmpNum>0) abort();


    //copy haplotypes into forwardCursor->y, here is the start of new site, before this
    //section, we are grouping haplotypes based on previous site's alleles
    CopyHap(k, forwardCursor);
    //now use haplotype alleles on current site k, to update array a and array d

    int u = 0, v = 0;
    int p = k + 1;
    int q = k + 1;

    for (rank = 0; rank < forwardCursor->M; ++rank) {

        if (forwardCursor->d[rank] > p) p = forwardCursor->d[rank];
        if (forwardCursor->d[rank] > q) q = forwardCursor->d[rank];

        if (forwardCursor->sortedY[rank] == 0) {
            forwardCursor->a[u] = forwardCursor->a[rank];
            forwardCursor->d[u] = p;
            ++u;
            p = 0;
            forwardCursor->c++;
            this->u[k][rank] = u;
        }
        else {
            forwardCursor->b[v] = forwardCursor->a[rank];
            forwardCursor->e[v] = q;
            ++v;
            q = 0;
            this->u[k][rank] = u;
        }
    }

    memcpy(forwardCursor->a + u, forwardCursor->b, v * sizeof(int));
    memcpy(forwardCursor->d + u, forwardCursor->e, v * sizeof(int));
    c[k] = u;

    forwardCursor->d[forwardCursor->M] = k + 2;
    a[k].assign(forwardCursor->a, forwardCursor->a + forwardCursor->M);//stores array a status after process current column haps

    for (int j = 0; j < (int)a[k].size(); ++j) {
        aMap[k][a[k][j]] = j;//at site k, hapID a[k] is stored at jth row/rank
    }
    d[k].assign(forwardCursor->d, forwardCursor->d + forwardCursor->M);

    //section below are dealing with last site
    if (k == N - 1)
    {
        i0 = 0;
        group = 0;

        clusterMembership.clear();
//        hasSiblings.clear();
        dist.clear();
        for (rank = 0; rank < forwardCursor->M; ++rank) {

            if (forwardCursor->d[rank] >
                (k - tmpT)) {//new cluster if current sequence and last sequence have common sequence less than T
                if (rank != 0) {
//                    fprintf(stderr,"d:%d\tk-tmpT:%d\n",forwardCursor->d[rank],k-tmpT);
                    tmpNumHap = rank - i0;
                    allele = haplotype[GetHapIDFromFwd(k, i0)][k];
                    Graph.StateNodeMat[k].push_back(new StateNode(group, tmpNumHap, allele));
                    Graph.StateNodeMat[k][group]->nodeIndex = group;
                    dist.push_back(std::vector<int>(tmpNumHap, 0));//state->rank_occupied
                    std::vector<int> tmpMem;
                    for (ia = i0; ia < rank; ++ia) {
                        hapID = GetHapIDFromFwd(k, ia);//original ID, by treating rank as backward ranking
                        haplotypeCluster[k][hapID] = group;
                        tmpMem.push_back(ia);
                        dist[group][ia - i0] = GetRankFromBack(k, hapID);
                        if (k >= 1) {
                            prevSiteStateIndex = GetHapStateFromFwd(k - 1, hapID);
                            Graph.StateNodeMat[k][group]->AddParentNode(prevSiteStateIndex,
                                                                        Graph.StateNodeMat[k - 1][prevSiteStateIndex]);
                            Graph.StateNodeMat[k - 1][prevSiteStateIndex]->AddChildNode(allele,
                                                                                        &(Graph.StateNodeMat[k][group]->nodeIndex),
                                                                                        hapID <= M ? 1.f : 1e-6f);
                        }
                    }
                    Graph.StateNodeMat[k][group]->SetID(prevSiteStateIndex,Graph.StateNodeMat[k - 1][prevSiteStateIndex]->ID);
                    Graph.RegisterState(k, Graph.StateNodeMat[k][group]->ID, &(Graph.StateNodeMat[k][group]->nodeIndex));
                    Graph.StateNodeMat[k][group]->AddChildNode(0, nullptr, 0);//end of the chain
                    std::sort(dist[group].begin(), dist[group].end());
//                    clusterAllele[k].push_back(allele);
                    clusterMembership.push_back(tmpMem);

                    i0 = rank;
                    group++;
                }
            }
        }
        //finish the last segment if i0 didn't reach the end
        if (i0 < forwardCursor->M) {
            tmpNumHap = forwardCursor->M - i0;
            allele = haplotype[GetHapIDFromFwd(k, i0)][k];

            dist.push_back(std::vector<int>(tmpNumHap, 0));//state->rank_occupied
            Graph.StateNodeMat[k].push_back(new StateNode(group, tmpNumHap, allele));
            Graph.StateNodeMat[k][group]->nodeIndex = group;
            std::vector<int> tmpMem;
            for (ia = i0; ia < forwardCursor->M; ++ia) {
                hapID = GetHapIDFromFwd(k, ia);//original ID, by treating rank as backward ranking
                haplotypeCluster[k][hapID] = group;
                tmpMem.push_back(ia);
                dist[group][ia - i0] = GetRankFromBack(k, hapID);
                if (k >= 1) {
                    prevSiteStateIndex = GetHapStateFromFwd(k - 1, hapID);
                    Graph.StateNodeMat[k][group]->AddParentNode(prevSiteStateIndex,
                                                                Graph.StateNodeMat[k - 1][prevSiteStateIndex]);
                    Graph.StateNodeMat[k - 1][prevSiteStateIndex]->AddChildNode(allele,
                                                                                &(Graph.StateNodeMat[k][group]->nodeIndex),
                                                                                hapID <= M ? 1.f : 1e-6f);
                }
            }
            Graph.StateNodeMat[k][group]->SetID(prevSiteStateIndex,Graph.StateNodeMat[k - 1][prevSiteStateIndex]->ID);
            Graph.RegisterState(k, Graph.StateNodeMat[k][group]->ID, &(Graph.StateNodeMat[k][group]->nodeIndex));
            Graph.StateNodeMat[k][group]->AddChildNode(0, nullptr, 0);
            std::sort(dist[group].begin(), dist[group].end());
//            clusterAllele[k].push_back(allele);
            clusterMembership.push_back(tmpMem);
        }

//            LabelNoSiblingCluster(k);
//            UpdateTransVector(k);//for the last second site
        Graph.NormalizeCurrentSiteTransitionProb(k - 1);
    }
    return 0;
}
*/
int PBWTWrapper::CursorForwardsTo(int k, int T) {
/*T is the length that how far you look back
 *This function must be called along the sites, no skip permitted;
 *Mask the site you want to skip at the begining if you have to.
 */
    int rank = 0;
    int i0 = 0;
    StateIndex group = 0;

    clusterMembership.clear();
    dist.clear();
    rightCoordinate.clear();
    rightCoordinateStat.clear();

    int tmpT = k > T ? T : k;
    if (k >= 1) {
        for (rank = 0; rank < M; ++rank) {
            if (d[rank] >
                (k - tmpT)) {//new cluster if current sequence and last sequence have common sequence less than T
                if (rank != 0) {
                    CreateNewCluster(k, rank, i0, group);
                    i0 = rank;
                    group++;
                }
            }
        }
        //finish the last segment if i0 didn't reach the end
        if (i0 < M) {
            CreateNewCluster(k, M, i0, group);
        }
        int numState = GetNumStates(k - 1);
        if (numState > 1 && k < (N - 20)) RegressionMergeAtSite(k - 1, true);//TODO:implement this function

    }

    UpdateAandD(k);
    //section below are dealing with last site
    if (k == N - 1) {
        i0 = 0;
        group = 0;
        clusterMembership.clear();
        dist.clear();
        for (rank = 0; rank < M; ++rank) {
            /*assign states of last column based on previous d and sortedY*/
            if (d[rank] >
                (k - tmpT)) {//new cluster if current sequence and last sequence have common sequence less than T
                if (rank != 0) {
                    CreateLastSiteCluster(k, rank, i0, group);
                    i0 = rank;
                    group++;
                }
            }
        }
        //finish the last segment if i0 didn't reach the end
        if (i0 < M) {
            CreateLastSiteCluster(k, M, i0, group);
        }
//        Graph.NormalizeCurrentSiteTransitionProb(k - 1);
    }
    return 0;
}

void PBWTWrapper::CreateLastSiteCluster(int k, int rank, int i0, StateIndex group) {

    int tmpNumHap = rank - i0;
    uchar allele = haplotype[GetHapIDFromFwd(i0)][k];
    StateIndex prevSiteStateIndex = 0;
    int hapID = 0;

    Graph.StateNodeMat[k].push_back(new StateNode(allele));
    dist.push_back(std::vector<int>(tmpNumHap, 0));//state->rank_occupied
    std::vector<int> tmpMem;
    for (int ia = i0; ia < rank; ++ia) {
        hapID = GetHapIDFromFwd(ia);//original ID, by treating rank as backward ranking
        haplotypeCluster[k][hapID] = group;
        tmpMem.push_back(ia);
        dist[group][ia - i0] = GetRankFromBack(k, hapID);
        prevSiteStateIndex = GetHapStateFromFwd(k - 1, hapID);
        Graph.StateNodeMat[k][group]->AddParentNode(prevSiteStateIndex);
        Graph.StateNodeMat[k - 1][prevSiteStateIndex]->AddChildNode(allele, group);
    }
//    Graph.StateNodeMat[k][group]->AddChildNode(0, -1);//end of the chain
    sort(dist[group].begin(), dist[group].end());
    clusterMembership.push_back(tmpMem);
}

int PBWTWrapper::CreateNewCluster(int k, int rank, int i0, StateIndex group) {

    int tmpNumHap = rank - i0;
    int hapID = GetHapIDFromFwd(i0);
    uchar allele = haplotype[hapID][k - 1];
    StateIndex prevSiteStateIndex = 0;


    Graph.StateNodeMat[k - 1].push_back(new StateNode(allele));

    dist.push_back(std::vector<int>(tmpNumHap, 0));//state->rank_occupied
    rightCoordinate.push_back(std::vector<float>(tmpNumHap, (float) 0.));//state->rank_occupied
    std::vector<int> tmpMem;
    for (int ia = i0; ia < rank; ++ia) {
        hapID = GetHapIDFromFwd(ia);//original ID, by treating rank as backward ranking
        haplotypeCluster[k - 1][hapID] = group;
        tmpMem.push_back(ia);
        dist[group][ia - i0] = GetRankFromBack(k - 1, hapID);
        rightCoordinate[group][ia - i0] = GetDistanceFromBack(k - 1, dist[group][ia - i0]);
        if (k >= 2) {
            prevSiteStateIndex = GetHapStateFromFwd(k - 2, hapID);
            Graph.StateNodeMat[k - 1][group]->AddParentNode(prevSiteStateIndex);//site k-1
            Graph.StateNodeMat[k - 2][prevSiteStateIndex]->AddChildNode(allele, group);
        }
    }
//    if (k == 1) {
//        Graph.StateNodeMat[0][group]->AddParentNode(0);//site 1
//    }

    rightCoordinateStat.push_back(rightCoordinate[group]);
    sort(dist[group].begin(),
         dist[group].end());//TODO:only sort nodes need to be compared, skip those otherwise

    clusterMembership.push_back(tmpMem);
    return prevSiteStateIndex;
}

void PBWTWrapper::UpdateAandD(int k) {
    //copy haplotypes into forwardCursor->y, here is the start of new site, before this
    //section, we are grouping haplotypes based on previous site's alleles
    CopyHap(k, a);
    //now use haplotype alleles on current site k, to update array a and array d
    std::vector<int> b(M, 0), e(M, 0);
    int c = 0;
    int u = 0, v = 0;
    int p = k + 1;
    int q = k + 1;

    for (int rank = 0; rank < M; ++rank) {
        if (d[rank] > p) p = d[rank];
        if (d[rank] > q) q = d[rank];
        if (sortedY[rank] == 0) {
            a[u] = a[rank];
            d[u] = p;
            ++u;
            p = 0;
            c++;
//                PBWTWrapper::u[k][rank] = u;
        } else {
            b[v] = a[rank];
            e[v] = q;
            ++v;
            q = 0;
//                PBWTWrapper::u[k][rank] = u;
        }
    }

//    memcpy(a + u, b, v * sizeof(int));
    std::copy(b.begin(), b.begin() + v, a.begin() + u);

//    memcpy(d + u, e, v * sizeof(int));
    std::copy(e.begin(), e.begin() + v, d.begin() + u);

//    c[k] = u;
//    d[M] = k + 2; /* sentinels */
}

/*
int PBWTWrapper::FastCursorForwardsTo(int k, int T, const PBWTWrapper& baseWrapper) {

    int rank, i0 = 0, ia;
    int group = 0;

    int hapID(0), prevSiteStateIndex(0);

    clusterMembership.clear();
    dist.clear();
    rightCoordinate.clear();
    rightCoordinateStat.clear();


    int tmpNumHap(0);
    char allele(-1);

    int tmpT = k > T ? T : k;


    //cluster of the previous site k-1
    for (rank = 0; rank < forwardCursor->M; ++rank) {

            (k - tmpT) &&
            (k >=1 )) {//new cluster if current sequence and last sequence have common sequence less than T
            if (rank != 0) {
                tmpNumHap = rank - i0;
                allele = haplotype[GetHapIDFromFwd(i0)][k - 1];

                Graph.StateNodeMat[k - 1].push_back(new StateNode(group, tmpNumHap, allele));
                Graph.StateNodeMat[k - 1][group]->nodeIndex = group;

                dist.push_back(std::vector<int>(tmpNumHap, 0));//state->rank_occupied
                rightCoordinate.push_back(std::vector<float>(tmpNumHap, (float) 0.));//state->rank_occupied
                std::vector<int> tmpMem;
                for (ia = i0; ia < rank; ++ia) {
                    hapID = GetHapIDFromFwd(ia);//original ID, by treating rank as backward ranking
                    haplotypeCluster[k - 1][hapID] = group;
                    tmpMem.push_back(ia);
                    dist[group][ia - i0] = GetRankFromBack(k - 1, hapID);
                    rightCoordinate[group][ia - i0] = GetDistanceFromBack(k - 1, dist[group][ia - i0]);
                    if (k >= 2) {
                        prevSiteStateIndex = GetHapStateFromFwd(k - 2, hapID);
                        Graph.StateNodeMat[k - 1][group]->AddParentNode(prevSiteStateIndex, Graph.StateNodeMat[k - 2][prevSiteStateIndex]);//site k-1
                        Graph.StateNodeMat[k - 2][prevSiteStateIndex]->AddChildNode(allele, &(Graph.StateNodeMat[k - 1][group]->nodeIndex), 1);
                    }
                }

                if (k == 1) {
                    Graph.StateNodeMat[0][group]->AddParentNode(0, nullptr);//site 1
                    Graph.StateNodeMat[0][group]->SetID(0,0);
                }
                else
                    Graph.StateNodeMat[k - 1][group]->SetID(prevSiteStateIndex,Graph.StateNodeMat[k - 2][prevSiteStateIndex]->ID);

                //update old state, add in new state
                auto iter = baseWrapper.Graph.StateNodeID2IndexPtr[k - 1].equal_range(Graph.StateNodeMat[k - 1][group]->ID);
                if (iter.first == baseWrapper.Graph.StateNodeID2IndexPtr[k - 1].end())
                {
                    Graph.StateNodeMat[k - 1][group]->needMergeUpdate = true;
//                    std::cerr<<"appear here"<<std::endl;
                }
                else{
                    Graph.StateNodeMat[k - 1][group]->needMergeUpdate = true;
                    for (auto i = iter.first; i !=iter.second; ++i) {
                        if(Graph.StateNodeMat[k - 1][group]->IsStateIdentical(
                                *(baseWrapper.Graph.StateNodeMat[k - 1][*(i->second)])))
                        {
                            Graph.StateNodeMat[k - 1][group]->needMergeUpdate = false;
                            break;
                        }
                    }
                }

                rightCoordinateStat.push_back(rightCoordinate[group]);
                std::sort(dist[group].begin(),dist[group].end());
                clusterMembership.push_back(tmpMem);

                i0 = rank;
                group++;


            }
        }
    }
    //finish the last segment if i0 didn't reach the end
    if (i0 < forwardCursor->M) {
        if (k != 0) {
            tmpNumHap = forwardCursor->M - i0;
            allele = haplotype[GetHapIDFromFwd(i0)][k - 1];

            Graph.StateNodeMat[k - 1].push_back(new StateNode(group, tmpNumHap, allele));
            Graph.StateNodeMat[k - 1][group]->nodeIndex = group;

            dist.push_back(std::vector<int>(tmpNumHap, 0));//state->rank_occupied
            rightCoordinate.push_back(std::vector<float>(tmpNumHap, (float) 0.));//state->rank_occupied
            std::vector<int> tmpMem;
            for (ia = i0; ia < forwardCursor->M; ++ia) {
                hapID = GetHapIDFromFwd(ia);//original ID, by treating rank as backward ranking
                haplotypeCluster[k - 1][hapID] = group;
                tmpMem.push_back(ia);
                dist[group][ia - i0] = GetRankFromBack(k - 1, hapID);
                rightCoordinate[group][ia - i0] = GetDistanceFromBack(k - 1, dist[group][ia - i0]);
                if (k >= 2) {
                    prevSiteStateIndex = GetHapStateFromFwd(k - 2, hapID);
                    Graph.StateNodeMat[k - 1][group]->AddParentNode(prevSiteStateIndex,
                                                                    Graph.StateNodeMat[k - 2][prevSiteStateIndex]);
                    Graph.StateNodeMat[k - 2][prevSiteStateIndex]->AddChildNode(allele, &(Graph.StateNodeMat[k -
                                                                                                             1][group]->nodeIndex),
                                                                                1);
                }
            }

            if (k == 1) {
                Graph.StateNodeMat[0][group]->AddParentNode(0, nullptr);//site 1
                Graph.StateNodeMat[0][group]->SetID(0,0);
            }
            else
                Graph.StateNodeMat[k - 1][group]->SetID(prevSiteStateIndex,Graph.StateNodeMat[k - 2][prevSiteStateIndex]->ID);
            //update old state, add in new state
            auto iter = baseWrapper.Graph.StateNodeID2IndexPtr[k - 1].equal_range(Graph.StateNodeMat[k - 1][group]->ID);
            if (iter.first == baseWrapper.Graph.StateNodeID2IndexPtr[k - 1].end())
            {
                Graph.StateNodeMat[k - 1][group]->needMergeUpdate = true;
            }
            else{
                Graph.StateNodeMat[k - 1][group]->needMergeUpdate = true;
                for (auto i = iter.first; i !=iter.second; ++i) {
                    if(Graph.StateNodeMat[k - 1][group]->IsStateIdentical(
                            *(baseWrapper.Graph.StateNodeMat[k - 1][*(i->second)])))
                    {
                        Graph.StateNodeMat[k - 1][group]->needMergeUpdate = false;
                        break;
                    }
                }
            }
            rightCoordinateStat.push_back(rightCoordinate[group]);
            std::sort(dist[group].begin(), dist[group].end());
            clusterMembership.push_back(tmpMem);
        }
    }

    //merge cluster based on KS test
    int merged = 0;
#ifdef DEBUG
    if(k!=0){
        fprintf(stderr, "at site:%d\n", k-1);
        PrintVector(haplotypeCluster[k-1],"haplotype state before merge state");
        PrintVector(clusterAllele[k-1],"state allele before merge allele");
    }
#endif
    if (k >= 1) {
        if (GetNumStates(k-1) != 1) merged = RegressionMergeAtSite(k - 1,false);//TODO:implement this function
        Graph.NormalizeCurrentSiteTransitionProb(k - 2);
//        std::cerr<<"new mat num state:"<<GetNumStates(k-1)<<" baseWrapper:"<<baseWrapper.GetNumStates(k-1)<<std::endl;
#ifdef DEBUG
        if(test) {
            fprintf(stderr, "at site:%d\n", k-1);
            PrintVector(haplotypeCluster[k-1], "haplotype state after merge state");
            PrintVector(clusterAllele[k-1],"state allele after merge allele");
            fprintf(stderr,"\n");
        }
#endif
    }



    //copy haplotypes into forwardCursor->y
    CopyHap(k, forwardCursor);

    int u = 0, v = 0;
    int p = k + 1;
    int q = k + 1;

    for (rank = 0; rank < forwardCursor->M; ++rank) {

        if (forwardCursor->d[rank] > p) p = forwardCursor->d[rank];
        if (forwardCursor->d[rank] > q) q = forwardCursor->d[rank];

        if (forwardCursor->sortedY[rank] == 0) {
            forwardCursor->a[u] = forwardCursor->a[rank];
            forwardCursor->d[u] = p;
            ++u;
            p = 0;

            forwardCursor->c++;
            this->u[k][rank] = u;

        }
        else {
            forwardCursor->b[v] = forwardCursor->a[rank];
            forwardCursor->e[v] = q;
            ++v;
            q = 0;

            this->u[k][rank] = u;
        }
    }

    memcpy(forwardCursor->a + u, forwardCursor->b, v * sizeof(int));
    memcpy(forwardCursor->d + u, forwardCursor->e, v * sizeof(int));
    c[k] = u;

    forwardCursor->d[forwardCursor->M] = k + 2;
    a[k].assign(forwardCursor->a, forwardCursor->a + forwardCursor->M);

    for (int j = 0; j < (int)a[k].size(); ++j) {
        aMap[k][a[k][j]] = j;
    }
    d[k].assign(forwardCursor->d, forwardCursor->d + forwardCursor->M);

    if (k == N - 1)//deal with last columns
    {
        i0 = 0;
        group = 0;

        clusterMembership.clear();

        dist.clear();
        for (rank = 0; rank < forwardCursor->M; ++rank) {

            if (forwardCursor->d[rank] >
                (k - tmpT)) {//new cluster if current sequence and last sequence have common sequence less than T
                if (rank != 0) {
                    tmpNumHap = rank - i0;
                    allele = haplotype[GetHapIDFromFwd(i0)][k];

                    Graph.StateNodeMat[k].push_back(new StateNode(group, tmpNumHap, allele));
                    Graph.StateNodeMat[k][group]->nodeIndex = group;
                    dist.push_back(std::vector<int>(tmpNumHap, 0));//state->rank_occupied
                    std::vector<int> tmpMem;
                    for (ia = i0; ia < rank; ++ia) {
                        hapID = GetHapIDFromFwd(ia);//original ID, by treating rank as backward ranking
                        haplotypeCluster[k][hapID] = group;
                        tmpMem.push_back(ia);
                        dist[group][ia - i0] = GetRankFromBack(k, hapID);
                        if (k >= 1) {
                            prevSiteStateIndex = GetHapStateFromFwd(k - 1, hapID);
                            Graph.StateNodeMat[k][group]->AddParentNode(prevSiteStateIndex,
                                                                        Graph.StateNodeMat[k - 1][prevSiteStateIndex]);
                            Graph.StateNodeMat[k - 1][prevSiteStateIndex]->AddChildNode(allele,
                                                                                        &(Graph.StateNodeMat[k][group]->nodeIndex),
                                                                                        1);
                        }
                    }
                    Graph.StateNodeMat[k][group]->AddChildNode(0, nullptr, 0);//end of the chain

                    Graph.StateNodeMat[k][group]->SetID(prevSiteStateIndex,Graph.StateNodeMat[k - 1][prevSiteStateIndex]->ID);
                    //update old state, add in new state
                    auto iter = baseWrapper.Graph.StateNodeID2IndexPtr[k].equal_range(Graph.StateNodeMat[k][group]->ID);
                    if (iter.first == baseWrapper.Graph.StateNodeID2IndexPtr[k].end())
                    {
                        Graph.StateNodeMat[k][group]->needMergeUpdate = true;
                    }
                    else{
                        Graph.StateNodeMat[k][group]->needMergeUpdate = true;
                        for (auto i = iter.first; i !=iter.second; ++i) {
                            if(Graph.StateNodeMat[k][group]->IsStateIdentical(
                                    *(baseWrapper.Graph.StateNodeMat[k][*(i->second)])))
                            {
                                Graph.StateNodeMat[k][group]->needMergeUpdate = false;
                                break;
                            }
                        }
                    }
                    std::sort(dist[group].begin(), dist[group].end());
                    clusterMembership.push_back(tmpMem);
                    i0 = rank;
                    group++;
                }
            }
        }
        //finish the last segment if i0 didn't reach the end
        if (i0 < forwardCursor->M) {
            tmpNumHap = forwardCursor->M - i0;
            allele = haplotype[GetHapIDFromFwd(i0)][k];

            dist.push_back(std::vector<int>(tmpNumHap, 0));//state->rank_occupied
            Graph.StateNodeMat[k].push_back(new StateNode(group, tmpNumHap, allele));
            Graph.StateNodeMat[k][group]->nodeIndex = group;
            std::vector<int> tmpMem;
            for (ia = i0; ia < forwardCursor->M; ++ia) {
                hapID = GetHapIDFromFwd(ia);//original ID, by treating rank as backward ranking
                haplotypeCluster[k][hapID] = group;
                tmpMem.push_back(ia);
                dist[group][ia - i0] = GetRankFromBack(k, hapID);
                if (k >= 1) {
                    prevSiteStateIndex = GetHapStateFromFwd(k - 1, hapID);
                    Graph.StateNodeMat[k][group]->AddParentNode(prevSiteStateIndex,
                                                                Graph.StateNodeMat[k - 1][prevSiteStateIndex]);
                    Graph.StateNodeMat[k - 1][prevSiteStateIndex]->AddChildNode(allele,
                                                                                &(Graph.StateNodeMat[k][group]->nodeIndex),
                                                                                1);
                }
            }
            Graph.StateNodeMat[k][group]->AddChildNode(0, nullptr, 0);

            Graph.StateNodeMat[k][group]->SetID(prevSiteStateIndex,Graph.StateNodeMat[k - 1][prevSiteStateIndex]->ID);
            //update old state, add in new state
            auto iter = baseWrapper.Graph.StateNodeID2IndexPtr[k].equal_range(Graph.StateNodeMat[k][group]->ID);
            if (iter.first == baseWrapper.Graph.StateNodeID2IndexPtr[k].end())
            {
                Graph.StateNodeMat[k][group]->needMergeUpdate = true;
            }
            else{
                Graph.StateNodeMat[k][group]->needMergeUpdate = true;
                for (auto i = iter.first; i !=iter.second; ++i) {
                    if(Graph.StateNodeMat[k][group]->IsStateIdentical(
                            *(baseWrapper.Graph.StateNodeMat[k][*(i->second)])))
                    {
                        Graph.StateNodeMat[k][group]->needMergeUpdate = false;
                        break;
                    }
                }
            }
            std::sort(dist[group].begin(), dist[group].end());
            clusterMembership.push_back(tmpMem);
        }
        Graph.NormalizeCurrentSiteTransitionProb(k - 1);
    }
    return 0;
}
*/
int PBWTWrapper::CursorBackwards() {

    for (int i = 0; i != N; i++) {
        CursorBackwardsTo(i, 50);
    }
    return 0;
}

int PBWTWrapper::CursorBackwardsTo(int siteBackword, int T) {

    int rank;
    //copy haplotypes into forwardCursor->y
    int siteForward = N - siteBackword - 1;//forward site 0 based
    CopyHap(siteForward, alpha);
    std::vector<int> e(M, 0), b(M, 0);
    int u = 0, v = 0;
    int p = siteBackword + 1;
    int q = siteBackword + 1;
    float cumCoordinate = 0.0f;
    float *tmpD1 = new float[M];
    float *tmpD2 = new float[M];
    for (rank = 0; rank < M; ++rank) {
        if (delta[rank] > p) p = delta[rank];
        if (delta[rank] > q) q = delta[rank];
        if (sortedY[rank] == 0) {
            alpha[u] = alpha[rank];
            delta[u] = p;
            cumCoordinate += (siteBackword - p + 1) > 0 ? 1.0f / (float) (siteBackword - p + 1) : 1.;
            tmpD1[u] = cumCoordinate;
            p = 0;
            ++u;
//            reverseCursor->c++;
//            ultra[siteForward][rank] = u;
        } else {
            b[v] = alpha[rank];
            e[v] = q;
            cumCoordinate += (siteBackword - q + 1) > 0 ? 1.0f / (float) (siteBackword - q + 1) : 1.;
            tmpD2[v] = cumCoordinate;
            q = 0;
            ++v;
//            ultra[siteForward][rank] = u;
        }
    }
//    memcpy(alpha + u, reverseCursor->b, v * sizeof(int));
    std::copy(b.begin(), b.begin() + v, alpha.begin() + u);
//    memcpy(delta + u, reverseCursor->e, v * sizeof(int));
    std::copy(e.begin(), e.begin() + v, delta.begin() + u);
    memcpy(tmpD1 + u, tmpD2, v * sizeof(float));
//    delta[M] = siteBackword + 2; /* sentinels */
//    alpha.assign(alpha.begin(), alpha.end());
    celta[siteForward] = u;
    for (int j = 0; j < (int) alpha.size(); ++j) {
        alphaMap[siteForward][alpha[j]] = j;
    }
    bkDistance[siteForward].assign(tmpD1, tmpD1 + M);
    allDelta[siteForward].assign(delta.begin(), delta.end());
    delete[] tmpD1;
    delete[] tmpD2;
    return 0;
}

int PBWTWrapper::CopyHap(int k, std::vector<int> &tmpA) {//this function has the same effect as forward/backward read
    for (int i = 0; i != M; ++i) {//ensure that alleles are 0/1 values not '0'/'1'
        sortedY[i] = haplotype[tmpA[i]][k];
    }
    //PrintVector(Cursor->sortedY,Cursor->M,"fromCopyHap");
    return 0;
}

void PBWTWrapper::MergeSortedArrayToA(std::vector<int> &a, std::vector<int> &b) {
    int indexA(0), indexB(0), indexTotal(0);
    int sizeA(a.size()), sizeB(b.size());
    std::vector<int> mergedDist(sizeA + sizeB);
    while (indexA < sizeA && indexB < sizeB) {
        if (a[indexA] < b[indexB]) {
            mergedDist[indexTotal++] = a[indexA];
            indexA++;

        } else {
            mergedDist[indexTotal++] = b[indexB];
            indexB++;
        }

    }
    while (indexA < sizeA) {
        mergedDist[indexTotal++] = a[indexA];
        indexA++;
    }
    while (indexB < sizeB) {
        mergedDist[indexTotal++] = b[indexB];
        indexB++;
    }
    a = mergedDist;
}

bool PBWTWrapper::IsRecipricalLengthOK(std::vector<int> &a, std::vector<int> &b) {
    double size = a.size() < b.size() ? a.size() : b.size();
    int cnt = 0;
    if (a.front() < b.front()) {
        if (a.back() < b.back()) {
            for (auto &value: b) {
                if (value < a.back()) cnt++;
                else break;
            }
            return cnt / size > 0.5;
        } else
            return true;
    } else {
        if (a.back() < b.back())
            return true;
        else {
            for (auto &value:a) {
                if (value < b.back()) cnt++;
                else break;
            }
            return cnt / size > 0.5;
        }
    }
}

bool PBWTWrapper::IsEditDistanceOK(int stateA, int stateB, int index, int error_thresh) {
    int thresh_pos = N - 1 - index - error_thresh;
    //ensure stateA has larger size
    if (dist[stateA].size() < dist[stateB].size()) std::swap(stateA, stateB);

    int lower, upper;
    int backRankA, backRankB;
    int numTruth(0);
    const int sizeB = (int) dist[stateB].size();
    for (int i = 0; i < sizeB; ++i) {//for each item in smaller set
        backRankB = dist[stateB][i];//get its rank
        auto lowerRankA = std::lower_bound(dist[stateA].begin(), dist[stateA].end(),
                                           backRankB);//find closest item in larger set
        if (lowerRankA == dist[stateA].end()) lower = backRankB + 1;//if not found, lower bound is rankB + 1
        else if (lowerRankA == dist[stateA].begin()) {//if found as smallest in larger set
            lower = -1;//lower bound is minus 1
        } else {
            lowerRankA--;//otherwise, get the element in larger set immediately smaller than current item in smaller set
            lower = *lowerRankA;
        }
        auto upperRankB = std::upper_bound(dist[stateA].begin(), dist[stateA].end(),
                                           backRankB);//same thing for upper bound
        if (upperRankB == dist[stateA].end()) upper = backRankB - 1;
        else upper = *upperRankB;
        while (lower != -1 and lower < backRankB) {//if lower bound is smaller than current item
            if (allDelta[index][++lower] < thresh_pos)
                continue;//keep going up if shared prefix ended beyound thresh_pos
            else break;
        }
        while (backRankB < upper) {//if upper bound is larger than current item
            if (allDelta[index][--upper] < thresh_pos)
                continue;//keep going down if share prefix ended beyound thresh_pos
            else break;
        }
        if (lower == backRankB || upper == backRankB)
            numTruth++;//if either bound meet current item's rank, that means current item is a member of them
    }
    if (numTruth == sizeB) return true;//if all of the items in smaller set are members of larger set, return true
    return false;
}

/*
bool PBWTWrapper::IsInSameBackCluster(int stateA, int stateB, int site, int error_thresh) {
    //ensure stateA has larger size
    if (dist[stateA].size() < dist[stateB].size()) std::swap(stateA, stateB);

    int lower, upper;
    int hapIDa,hapIDb;
    int backRankA, backRankB;
    int numTruth(0);
    int sizeB=(int)dist[stateB].size();
    for (int i = 0; i < sizeB; ++i) {
        backRankB = dist[stateB][i];
        auto lowerRankA = std::lower_bound(dist[stateA].begin(), dist[stateA].end(), backRankB);
        if (lowerRankA == dist[stateA].end()) lower = backRankB + 1;
        else if (lowerRankA == dist[stateA].begin()) {
            lower = -1;
        }
        else {
            lowerRankA--;
            lower = *lowerRankA;
        }
        auto upperRankB = std::upper_bound(dist[stateA].begin(), dist[stateA].end(), backRankB);
        if (upperRankB == dist[stateA].end()) upper = backRankB - 1;
        else upper = *upperRankB;
        if(lower!=-1&&lower!=backRankB+1)
        {
            hapIDa=GetHapIDFromBack(site,lower);
//            fprintf(stderr,"site:%d\tlower:%d\thapID:%d\tbackRankB:%d\n",site,lower,hapIDa,backRankB);
            if(GetHapStateFromBack(site,hapIDa) == GetHapStateFromBack(site,GetHapIDFromBack(site,backRankB)))
            {
                numTruth++;
                continue;
            }
        }
        if(upper!=backRankB-1)
        {
            hapIDb=GetHapIDFromBack(site,upper);
            if(GetHapStateFromBack(site,hapIDb) == GetHapStateFromBack(site,GetHapIDFromBack(site,backRankB)))
            {
                numTruth++;
                continue;
            }
        }
    }
    if (numTruth == sizeB) return true;
    return false;

}
*/
int PBWTWrapper::CalculateDmax(double &pval, double &Dmax, std::vector<int> &j, std::vector<int> &k) {
    int sizeJ(j.size()), sizeK(k.size());
    double thresh = 1.44 * 1.22 * sqrt(1. / sizeJ + 1. / sizeK);//apha:0.1
    //assume dist has sorted ranks
    double Ddiff = 0;
    pval = 0;
    Dmax = 0;
    int indexA(0), indexB(0);
    while (indexA < sizeJ && indexB < sizeK) {
        if (j[indexA] < k[indexB]) {
            indexA++;
            Ddiff += 1. / sizeJ;
        } else {
            indexB++;
            Ddiff -= 1. / sizeK;
        }
        if (fabs(Ddiff) > Dmax) {
            Dmax = fabs(Ddiff);
            if (sizeJ * sizeK < 10000) {
                if (GetPValue(sizeJ, sizeK, Dmax) < P_thresh) return 1;
            } else {
                if (Dmax > thresh)
                    return 1;
            }
        }
    }
    while (indexA < sizeJ) {
        indexA++;
        Ddiff += 1. / sizeJ;
        if (fabs(Ddiff) > Dmax) {
            Dmax = fabs(Ddiff);
            if (sizeJ * sizeK < 10000) {
                if (GetPValue(sizeJ, sizeK, Dmax) < P_thresh) return 1;
            } else {
                if (Dmax > thresh)
                    return 1;
            }
        }
    }
    while (indexB < sizeK) {
        indexB++;
        Ddiff -= 1. / sizeK;
        if (fabs(Ddiff) > Dmax) {
            Dmax = fabs(Ddiff);
            if (sizeJ * sizeK < 10000) {
                if (GetPValue(sizeJ, sizeK, Dmax) < P_thresh) return 1;
            } else {
                if (Dmax > thresh)
                    return 1;
            }
        }
    }
    return 0;
};

int PBWTWrapper::CalculateDmaxBeta(double &pval, double &Dmax, std::vector<int> &j, std::vector<int> &k) {
    double thresh = 10;
    //assume dist has sorted ranks
    double Ddiff = 0;
    pval = 0;
    Dmax = 0;
    double nx(0), ny(0), alpha(0.5), beta(0.5), pA(0), pB(0);
    int indexA(0), indexB(0);
    int sizeJ(j.size()), sizeK(k.size());
    while (indexA < sizeJ && indexB < sizeK) {
        if (j[indexA] < k[indexB]) {
            indexA++;
            nx++;
        } else {
            indexB++;
            ny++;
        }
        pA = (nx + alpha) / (sizeJ + alpha + beta);
        pB = (ny + alpha) / (sizeK + alpha + beta);
        Ddiff = (pA - pB) * (pA - pB) / (pA * (1 - pA) / sizeJ + pB * (1 - pB) / sizeK);
        if (Ddiff > Dmax) Dmax = Ddiff;
        if (Dmax >= thresh) {
            Dmax = sqrt(Dmax);
            //goto DIST_END;
            return 1;
        }

    }
    while (indexA < sizeJ) {
        indexA++;
        nx++;
        pA = (nx + alpha) / (sizeJ + alpha + beta);
        Ddiff = (pA - pB) * (pA - pB) / (pA * (1 - pA) / sizeJ + pB * (1 - pB) / sizeK);
        if (Ddiff > Dmax) Dmax = Ddiff;
        if (Dmax >= thresh) {
            Dmax = sqrt(Dmax);
            //goto DIST_END;
            return 1;
        }
    }
    while (indexB < sizeK) {
        indexB++;
        ny++;
        pB = (ny + alpha) / (sizeK + alpha + beta);
        Ddiff = (pA - pB) * (pA - pB) / (pA * (1 - pA) / sizeJ + pB * (1 - pB) / sizeK);
        if (Ddiff > Dmax) Dmax = Ddiff;
        if (Dmax >= thresh) {
            Dmax = sqrt(Dmax);
            //goto DIST_END;
            return 1;
        }
    }
    return 0;
};

#include "../SinglePhasing/libStatGen/include/Random.h"

#define MIN_FREQ 10
int PBWTWrapper::RegressionMergeAtSite(int site, bool isBaseWrapper) {

    Random localRandom(31415);
//    if(recomRate[site-1]>1e-4) return 1;
    int ret(0);
    int currentNumCluster = GetNumStates(site);
//    std::cerr << "Enter Site:" << site << " has " << currentNumCluster << " state and P_thresh:" << P_thresh
//              << std::endl;

    std::vector<bool> removeIndicator(currentNumCluster, false);
    std::vector<bool> retainIndicator(currentNumCluster, false);


    tmpOrder = 0;
    pval = 0;
    EXACT = false;
    stateOrder.clear();//mapping oldState to newOrder
    removeMembership.clear();//rankID,state
    std::vector<StateNode *> tmpNodeVec;

//
//    int totalChanged=0;
//    for (int l = 0; l <Graph.StateNodeMat[site].size() ; ++l) {
//        if(Graph.StateNodeMat[site][l]->needMergeUpdate) totalChanged++;
//    }
//    std::cerr<<"site:"<<site<<"\thas "<<totalChanged<<" states out of "<<Graph.StateNodeMat[site].size()<<" states"<<std::endl;


    StateIndex retainState;
    StateIndex removeState;

    double Dmax = 0.;
    double Ddiff = 0.;

    std::vector<StateIndex> stateWithSibs, stateWithoutSibs;
    for (StateIndex i = 0; i < dist.size(); ++i) {
        if (!HasSiblings(site, i)) stateWithoutSibs.push_back(i);
        else stateWithSibs.push_back(i);
    }
    stateWithSibs.insert(stateWithSibs.end(), stateWithoutSibs.begin(), stateWithoutSibs.end());

    StateIndex stateL(0), stateR(0);
    size_t sizeL(0), sizeR(0);
    int totalPair = 0;
//    int skippedPair=0;
    for (int j = 0; j < int(stateWithSibs.size() -
                            stateWithoutSibs.size()); ++j) {//enumerate through all the states, usually retain stateWithSibs
        stateL = stateWithSibs[j];
        for (int k = j + 1; k < (int) stateWithSibs.size(); ++k) {
            stateR = stateWithSibs[k];

            if (GetAllele(site, stateL) != GetAllele(site, stateR))
                continue;

            if ((dist[stateL].size() == 1 or dist[stateR].size() == 1) and
                (dist[stateL].size() != dist[stateR].size())) {
                if (!IsEditDistanceOK(stateL, stateR, site, 5))
                    continue;
                MaxPair mergePair = {stateL, stateR, Dmax, EXACT, 1, localRandom.Next()};
                mergePairList.push(mergePair);
            } else {
                if (!IsRecipricalLengthOK(dist[stateL], dist[stateR]))
                    continue;

//            if (rightCoordinateStat[stateL].Combine(rightCoordinateStat[stateR]).IsSignificant()) {
//                    PrintVector(dist[stateL],dist[stateL].size(),"rank j:");
//                    PrintVector(dist[stateR],dist[stateR].size(),"rank k:");
//                    PrintVector(rightCoordinate[j],rightCoordinate[j].size(),"big beta hat right j:");
//                    PrintVector(rightCoordinate[k],rightCoordinate[k].size(),"big beta hat right k:");
//                continue;
//            }

                if (dist[stateL].size() <= MIN_FREQ) {
                    if (dist[stateR].size() <= MIN_FREQ) {//both rare
                        continue;
                    } else//j rare, k not
                    {
                        if (!IsEditDistanceOK(stateL, stateR, site, 5))
                            continue;
                    }
                } else if (dist[stateR].size() <= MIN_FREQ)//j not, k rare
                {
                    if (!IsEditDistanceOK(stateL, stateR, site, 5))
                        continue;
                }

//            if (dist[stateL].size() <= MIN_FREQ || dist[stateR].size() <= MIN_FREQ) {
//                if (!IsEditDistanceOK(stateL, stateR, site, 5))
//                    continue;
//            }


                sizeL = clusterMembership[stateL].size();
                sizeR = clusterMembership[stateR].size();

                if (sizeL * sizeR < 10000)//exact
                {
                    EXACT = true;
                } else {
                    EXACT = false;
                }

                if (CalculateDmax(pval, Dmax, dist[stateL], dist[stateR]))//return early
                {
                    goto DIST_END;
                }

                if (EXACT) {
                    pval = GetPValue(sizeL, sizeR, Dmax);
                } else {
                    Ddiff = sqrt(double(sizeL * sizeR / (sizeL + sizeR))) * Dmax;
                    pval = 1. - pkstwo_wrapper(1, &Ddiff, 1e-06);
                }

                if (pval > P_thresh)//potentially from same group
                {
                    MaxPair mergePair = {stateL, stateR, Dmax, EXACT, pval, localRandom.Next()};
                    mergePairList.push(mergePair);
                }
            }

            DIST_END:
//            fprintf(stderr, "first:(%d,%d) Dmax:%f and Thresh:%f and Pval:%f, with sample size:%d,%d and %d,%d\n",
//                    j, k, Dmax, P_thresh, pval, dist[j].size(), dist[k].size(), clusterMembership[j].size(),
//                    clusterMembership[k].size());
            continue;

        }

    }

//    std::cerr<<"not in vain:"<<site<<"\t skipped "<<skippedPair<<" out of "<<totalPair<<std::endl;

    StateIndex clusterA(0), clusterB(0);
    double old_p_max = 0;
    while (!mergePairList.empty()) {
        MaxPair iter_pair = mergePairList.top();
        mergePairList.pop();
        old_p_max = iter_pair.pval;
        clusterA = iter_pair.clusterA;
        clusterB = iter_pair.clusterB;
        sizeL = clusterMembership[clusterA].size();
        sizeR = clusterMembership[clusterB].size();


        if (removeIndicator[clusterA] || removeIndicator[clusterB]) continue;//skip
        if (retainIndicator[clusterA] ||
            retainIndicator[clusterB]) {//if either of these two cluster ever merged before, recalculate p value
            if (GetAllele(site, stateL) != GetAllele(site, stateR))
                continue;
            if (dist[clusterA].size() <= MIN_FREQ) {
                if (dist[clusterB].size() <= MIN_FREQ) {//both rare
                    continue;
                } else//j rare, k not
                {
                    if (!IsEditDistanceOK(clusterA, clusterB, site, 5))
                        continue;
                }
            } else if (dist[clusterB].size() <= MIN_FREQ)//j not, k rare
            {
                if (!IsEditDistanceOK(clusterA, clusterB, site, 5))
                    continue;
            }

            {
                if (!IsRecipricalLengthOK(dist[clusterA], dist[clusterB]))//||!IsEditDistanceOK(backBone,j,k,site,0.5))
                    continue;
            }

//            if (dist[stateL].size() <= MIN_FREQ || dist[stateR].size() <= MIN_FREQ) {
//                if (!IsEditDistanceOK(stateL, stateR, site, 5))
//                    continue;
//            }

//            if (rightCoordinateStat[clusterA].Combine(rightCoordinateStat[clusterB]).IsSignificant()) {
//                continue;
//            }

            if (sizeL * sizeR < 10000)//exact
            {
                iter_pair.exact = true;
//                thresh = GetExactThresh(clusterMembership[clusterA].size(), clusterMembership[clusterB].size());
            } else {
                iter_pair.exact = false;
//                thresh = 1.22 * sqrt((clusterMembership[clusterA].size() + clusterMembership[clusterB].size()) /
//                                     (clusterMembership[clusterA].size() * clusterMembership[clusterB].size()));
            }

            if (CalculateDmax(iter_pair.pval, iter_pair.Dmax, dist[clusterA], dist[clusterB]))//return early
            {
                goto END_WHILE;
            }

            iter_pair.Dmax = Dmax;

            if (iter_pair.exact) {
                iter_pair.pval = GetPValue(sizeL, sizeR,
                                           Dmax);//1 - psmirnov2x(iter_pair.Dmax, clusterMembership[clusterA].size(), clusterMembership[clusterB].size());
            } else {
                Ddiff = sqrt(double(sizeL + sizeR)) * iter_pair.Dmax;
                iter_pair.pval = 1 - pkstwo_wrapper(1, &Ddiff, 1e-06);//two sided ks test p value
            }
            if (iter_pair.pval < old_p_max) {//the smaller the p val, the less likely to merge
                MaxPair mergePair = {clusterA, clusterB, iter_pair.Dmax, iter_pair.exact, iter_pair.pval,
                                     localRandom.Next()};
                mergePairList.push(mergePair);
                continue;
            }

        }


#ifdef DEBUG

        {
                    std::cerr << "\nenter rank distribution section, site " << site << ":" << std::endl;
                    for (auto i = 0; i != dist.size(); ++i) {
                        // PrintDistributionAtSite(i,haplotypeCluster[i]);
                        if (dist[i].size() == 0) continue;
                        std::cerr<<"state "<<i<<" :\t";
                        for (int j = 0; j <dist[i].size() ; ++j) {
                            //if(dist[i][j])
                            std::cerr << dist[i][j] << "\t";
                            //std::cerr<<clusterMembership[i][j]<<"\t";
                        }
                        std::cerr<<std::endl;
                    }
                    std::cerr << "exit rank distribution section!\n" << std::endl;
                }
                {
                    std::cerr << "\nenter membership section, site " << site << ":" << std::endl;
                    for (auto i = 0; i != clusterMembership.size(); ++i) {
                        if (clusterMembership[i].size() == 0) continue;
                        std::cerr<<"state "<<i<<" :\t";
                        for (int j = 0; j <clusterMembership[i].size() ; ++j) {
                            std::cerr << GetHapIDFromFwd(site, clusterMembership[i][j]) << "\t";
                            //std::cerr<<clusterMembership[i][j]<<"\t";
                        }
                        std::cerr<<std::endl;
                    }
                    std::cerr << "exit membership section!\n" << std::endl;
                }

#endif

        if (iter_pair.pval > P_thresh) {
//            if(site==10)
//            {
//
//            printf("merge of clusterA %d\tclusterB %d\tDmax:%f\tpval:%f\tret:%d\tsizeL:%d\tsizeR:%d\n",
//                   clusterA, clusterB, iter_pair.Dmax, iter_pair.pval, ret, sizeL, sizeR);
//            }
            retainState = clusterA;
            removeState = clusterB;
//            PrintVector(dist[retainState],"retainState");
//            PrintVector(dist[removeState],"removeState");
            ret = 1;
            if (!HasSiblings(site, clusterA)) {
                retainState = clusterB;
                removeState = clusterA;
            }

            DoMerge(site, retainState, removeState, dist, removeIndicator, retainIndicator, removeMembership);
            rightCoordinateStat[retainState] = rightCoordinateStat[retainState] + rightCoordinateStat[removeState];
            Graph.JoinNodes(site, retainState, removeState);
        }

        END_WHILE:
        continue;
    }

    if (ret) {
        for (StateIndex stateM = 0;
             stateM < (int) dist.size(); ++stateM)//loop through all remained states with the help of mergeIndicator
        {
            if (removeIndicator[stateM]) {
                delete Graph.StateNodeMat[site][stateM];
                Graph.StateNodeMat[site][stateM] = nullptr;
                continue;
            }
//            tmpAllele.push_back(GetAllele(site,stateM));
//            Graph.StateNodeMat[site][stateM]->nodeIndex = tmpOrder;
            for (auto kv:Graph.StateNodeMat[site][stateM]->GetParentIndexSet()) {
                Graph.UpdateChildNodeIndex(site - 1, kv, tmpOrder, GetAllele(site, stateM));
            }

            tmpNodeVec.push_back(Graph.StateNodeMat[site][stateM]);
            stateOrder[stateM] = tmpOrder;
            tmpOrder++;
        }
//        PrintVector(clusterAllele[site],"allele cluster states after");
        //PrintVector(haplotypeCluster[site],"haplotype cluster states after");
        for (auto &value:haplotypeCluster[site]) {
            value = stateOrder[value];//TODO:don't forget index in StateNode
        }
//        clusterAllele[site] = tmpAllele;//update merged cluster allele
        Graph.StateNodeMat[site] = tmpNodeVec;
//        PrintVector(clusterAllele[site],"allele cluster states final");
        //PrintVector(haplotypeCluster[site],"haplotype cluster states final");
        //adjust d array and a array
        MoveSegment(removeMembership, site);
    }
//    std::cerr << "Exit Site:" << site << " has " << GetNumStates(site) << " state" << std::endl;

//    if(site>12) exit(EXIT_FAILURE);
    return ret;
}

int PBWTWrapper::MergeAtSite(int site) {

    Random localRandom(31415);
    int ret(0);
    int currentNumCluster = GetNumStates(site);

    std::cerr << "Enter Site:" << site << " has " << currentNumCluster << " state and P_thresh:" << P_thresh
              << std::endl;

    std::vector<bool> removeIndicator(currentNumCluster, false);
    std::vector<bool> retainIndicator(currentNumCluster, false);


    tmpOrder = 0;
    pval = 0;
    EXACT = false;
    stateOrder.clear();//mapping oldState to newOrder

    removeMembership.clear();//rankID,state
    std::vector<StateNode *> tmpNodeVec;

    StateIndex retainState;
    StateIndex removeState;

    double Dmax = 0.;
    double Ddiff = 0.;

    std::vector<StateIndex> stateWithSibs, stateWithoutSibs;
    for (StateIndex i = 0; i < dist.size(); ++i) {
        if (!HasSiblings(site, i)) stateWithoutSibs.push_back(i);
        else stateWithSibs.push_back(i);
    }
    stateWithSibs.insert(stateWithSibs.end(), stateWithoutSibs.begin(), stateWithoutSibs.end());

    StateIndex stateL(0), stateR(0);
    size_t sizeL(0), sizeR(0);
    int totalPair = 0;
    for (int j = 0; j < int(stateWithSibs.size() -
                            stateWithoutSibs.size()); ++j) {//enumerate through all the states, usually retain stateWithSibs
        stateL = stateWithSibs[j];
        for (int k = j + 1; k < (int) stateWithSibs.size(); ++k) {
            stateR = stateWithSibs[k];

            if (GetAllele(site, stateL) != GetAllele(site, stateR))
                continue;

            if (!IsRecipricalLengthOK(dist[stateL], dist[stateR]))
                continue;

            if (dist[stateL].size() <= 20 || dist[stateR].size() <= 20) {
                if (!IsEditDistanceOK(stateL, stateR, site, 50))
                    continue;
            }

            sizeL = clusterMembership[stateL].size();
            sizeR = clusterMembership[stateR].size();

            if (sizeL * sizeR < 10000)//exact
            {
                EXACT = true;
            } else {
                EXACT = false;
            }

            if (CalculateDmax(pval, Dmax, dist[stateL], dist[stateR]))//return early
            {
                goto DIST_END;
            }

            if (EXACT) {
                pval = GetPValue(sizeL, sizeR, Dmax);
            } else {
                Ddiff = sqrt(double(sizeL + sizeR)) * Dmax;
                pval = 1. - pkstwo_wrapper(1, &Ddiff, 1e-06);
            }

            if (pval > P_thresh)//potentially from same group
            {
                MaxPair mergePair = {stateL, stateR, Dmax, EXACT, pval, localRandom.Next()};
                mergePairList.push(mergePair);
            }

            DIST_END:
//            fprintf(stderr, "first:(%d,%d) Dmax:%f and Thresh:%f and Pval:%f, with sample size:%d,%d and %d,%d\n",
//                    j, k, Dmax, P_thresh, pval, dist[j].size(), dist[k].size(), clusterMembership[j].size(),
//                    clusterMembership[k].size());
            continue;

        }

    }

    //process merging
    StateIndex clusterA(0), clusterB(0);
    double old_p_max = 0;
    while (!mergePairList.empty()) {
        MaxPair iter_pair = mergePairList.top();
        mergePairList.pop();
        old_p_max = iter_pair.pval;
        clusterA = iter_pair.clusterA;
        clusterB = iter_pair.clusterB;
        sizeL = clusterMembership[clusterA].size();
        sizeR = clusterMembership[clusterB].size();

        if (removeIndicator[clusterA] || removeIndicator[clusterB]) continue;//skip
        if (retainIndicator[clusterA] ||
            retainIndicator[clusterB]) {//if either of these two cluster ever merged before, recalculate p value

            if (GetAllele(site, stateL) != GetAllele(site, stateR))
                continue;
            if (dist[clusterA].size() <= 4) {
                if (dist[clusterB].size() <= 4) {//both rare
                    continue;
                } else//j rare, k not
                {
                    if (!IsEditDistanceOK(clusterA, clusterB, site, 50))
                        continue;
                }
            } else if (dist[clusterB].size() <= 4)//j not, k rare
            {
                if (!IsEditDistanceOK(clusterA, clusterB, site, 50))
                    continue;
            }

            {
                if (!IsRecipricalLengthOK(dist[clusterA], dist[clusterB]))
                    continue;
            }

            if (sizeL * sizeR < 10000)//exact
            {
                iter_pair.exact = true;
            } else {
                iter_pair.exact = false;
            }

            if (CalculateDmax(iter_pair.pval, iter_pair.Dmax, dist[clusterA], dist[clusterB]))//return early
            {
                goto END_WHILE;
            }

            iter_pair.Dmax = Dmax;

            if (iter_pair.exact) {
                iter_pair.pval = GetPValue(sizeL, sizeR,
                                           Dmax);//1 - psmirnov2x(iter_pair.Dmax, clusterMembership[clusterA].size(), clusterMembership[clusterB].size());
            } else {
                Ddiff = sqrt(double(sizeL + sizeR)) * iter_pair.Dmax;
                iter_pair.pval = 1 - pkstwo_wrapper(1, &Ddiff, 1e-06);//two sided ks test p value
            }
            if (iter_pair.pval < old_p_max) {//the smaller the p val, the less likely to merge
                MaxPair mergePair = {clusterA, clusterB, iter_pair.Dmax, iter_pair.exact, iter_pair.pval,
                                     localRandom.Next()};
                mergePairList.push(mergePair);
                continue;
            }

        }


#ifdef DEBUG

        {
                    std::cerr << "\nenter rank distribution section, site " << site << ":" << std::endl;
                    for (auto i = 0; i != dist.size(); ++i) {
                        // PrintDistributionAtSite(i,haplotypeCluster[i]);
                        if (dist[i].size() == 0) continue;
                        std::cerr<<"state "<<i<<" :\t";
                        for (int j = 0; j <dist[i].size() ; ++j) {
                            //if(dist[i][j])
                            std::cerr << dist[i][j] << "\t";
                            //std::cerr<<clusterMembership[i][j]<<"\t";
                        }
                        std::cerr<<std::endl;
                    }
                    std::cerr << "exit rank distribution section!\n" << std::endl;
                }
                {
                    std::cerr << "\nenter membership section, site " << site << ":" << std::endl;
                    for (auto i = 0; i != clusterMembership.size(); ++i) {
                        if (clusterMembership[i].size() == 0) continue;
                        std::cerr<<"state "<<i<<" :\t";
                        for (int j = 0; j <clusterMembership[i].size() ; ++j) {
                            std::cerr << GetHapIDFromFwd(site, clusterMembership[i][j]) << "\t";
                            //std::cerr<<clusterMembership[i][j]<<"\t";
                        }
                        std::cerr<<std::endl;
                    }
                    std::cerr << "exit membership section!\n" << std::endl;
                }

#endif

        if (iter_pair.pval > P_thresh) {
//            if(site==10)
//            {
//
//            printf("merge of clusterA %d\tclusterB %d:%f\tpval:%f\tret:%d\tsizeL:%d\tsizeR:%d\n",
//                   clusterA, clusterB, iter_pair.Dmax, iter_pair.pval, ret, sizeL, sizeR);
//            }
            retainState = clusterA;
            removeState = clusterB;
//            PrintVector(dist[retainState],"retainState");
//            PrintVector(dist[removeState],"removeState");
            ret = 1;
            if (!HasSiblings(site, clusterA)) {
                retainState = clusterB;
                removeState = clusterA;
            }

            DoMerge(site, retainState, removeState, dist, removeIndicator, retainIndicator, removeMembership);
            rightCoordinateStat[retainState] = rightCoordinateStat[retainState] + rightCoordinateStat[removeState];
            Graph.JoinNodes(site, retainState, removeState);
        }

        END_WHILE:
        continue;
    }

    if (ret) {
        for (StateIndex stateM = 0;
             stateM < (int) dist.size(); ++stateM)//loop through all remained states with the help of mergeIndicator
        {
            if (removeIndicator[stateM]) {
                delete Graph.StateNodeMat[site][stateM];
                Graph.StateNodeMat[site][stateM] = nullptr;
                continue;
            }

            for (auto kv:Graph.StateNodeMat[site][stateM]->GetParentIndexSet()) {
                Graph.UpdateChildNodeIndex(site - 1, kv, tmpOrder, GetAllele(site, stateM));
            }

            tmpNodeVec.push_back(Graph.StateNodeMat[site][stateM]);
            stateOrder[stateM] = tmpOrder;
            tmpOrder++;
        }

        for (auto &value:haplotypeCluster[site]) {
            value = stateOrder[value];//TODO:don't forget index in StateNode
        }
        Graph.StateNodeMat[site] = tmpNodeVec;

        //adjust d array and a array
        MoveSegment(removeMembership, site);
    }
    std::cerr << "Exit Site:" << site << " has " << GetNumStates(site) << " state" << std::endl;

    return ret;
}

void PBWTWrapper::DoMerge(int site, StateIndex retainState, StateIndex removeState, std::vector<std::vector<int>> &dist,
                          std::vector<bool, std::allocator<bool>> &removeIndicator,
                          std::vector<bool, std::allocator<bool>> &retainIndicator,
                          std::unordered_map<int, int> &removeMembership) {//move dist occupation from stateB to stateA
//    fprintf(stderr,"now merge state %d,%d at site %d\n",retainState,removeState,site);
//    for (int t = 0; t != dist[removeState].size(); ++t) {
//                dist[retainState][t] += dist[removeState][t];
//            }
    MergeSortedArrayToA(dist[retainState], dist[removeState]);

    dist[removeState].clear();
    //Merge Action, change removeIndicator
    retainIndicator[retainState] = true;
    removeIndicator[removeState] = true;

    int removeRankID = 0;
    for (int t = 0; t < (int) clusterMembership[removeState].size(); ++t) {
        removeRankID = clusterMembership[removeState][t];
        haplotypeCluster[site][GetHapIDFromFwd(removeRankID)] = retainState;
        clusterMembership[retainState].push_back(removeRankID);
        removeMembership[removeRankID] = removeState;
    }
    clusterMembership[removeState].clear();
//    for(auto kv:inEdges[site][removeState])
//    {
//        inEdges[site][retainState][kv.first]=true;
//    }
}


int PBWTWrapper::SetHaps(char **haps, int CopyStart, int CopyEnd, char **sampledHaps, int CopyStart2, int CopyEnd2,
                         float *rate, int phased) {

    if (rate)
        recomRate.assign(rate, rate + nMarkers - 1);
//    for (int k = 0; k <nMarkers-1 ; ++k) {
//        std::cerr<<"marker "<<k+1<<":"<<recomRate[k]<<std::endl;
//    }
    haplotype = new char *[M];
    int i = 0;
    for (int j = CopyStart; i < M && j < CopyEnd; ++i, ++j) {
        haplotype[i] = haps[j];
    }
    //only shuffle phased ref panel
    CopyStart = CopyStart > phased * 2 ? CopyStart : phased * 2;
    std::random_shuffle(haplotype + phased * 2, haplotype + (CopyEnd - CopyStart));
    if (sampledHaps)
        for (int j = CopyStart2; i < M && j < CopyEnd2; ++i, ++j) {
            haplotype[i] = sampledHaps[j];
        }

    return 0;
}

int PBWTWrapper::PrintDistributionAtSite(int state, std::vector<int> &dist) {
    //std::ofstream fout("/Users/fanzhang/Downloads/PlutoTest/rank.txt",std::ofstream::app);
    std::cerr << "state:" << state << ":\t";
    // fout<<"state:"<<state<<":\t";
    for (auto k:dist) {
        std::cerr << k << "\t";
        //    fout << k << "\t";

    }
    std::cerr << std::endl;
    //fout<<std::endl;

    return 0;
}

int PBWTWrapper::PrintSummary() {
//    mean nodes/level =  70.38  max nodes/level = 111  nodes = 10135
//    mean edges/level =  94.03  max edges/level = 177  edges = 13541
//    mean edges/node  =   1.34  mean count/node =  53.45
    int totalNodes(0), maxNodes(0);
    int totalEdges(0), maxEdges(0);
    float meanEdges(0.0), meanNodes(0.0);
    for (int i = 0; i < N; ++i) {
        totalNodes += GetNumStates(i);//clusterAllele[i].size();
        if (GetNumStates(i) > maxNodes) maxNodes = GetNumStates(i);//clusterAllele[i].size();
    }
    meanNodes = (float) totalNodes / N;
//
//    for (int j = 0; j < transVector.size(); ++j) {
//        int tmpEdges(0);
//        for (int i = 0; i < transVector[j].size(); ++i) {
//            for (int k = 0; k < transVector[j][i].size(); ++k) {
//                if (transVector[j][i][k] != 0) tmpEdges++;
//            }
//        }
//
//        if (tmpEdges > maxEdges) maxEdges = tmpEdges;
//        totalEdges += tmpEdges;
//    }
    meanEdges = (float) totalEdges / (N - 1);

    fprintf(stderr, "mean nodes/level = %f\tmax nodes/level = %d\tnodes = %d\n", meanNodes, maxNodes, totalNodes);
    fprintf(stderr, "mean edges/level = %f\tmax edges/level = %d\tedges = %d\n", meanEdges, maxEdges, totalEdges);
    fprintf(stderr, "mean edges/node = %f\tmean count/node = %f\n", (float) totalEdges / totalNodes,
            (float) (M) * (N) / totalNodes);

    return 0;
}

int PBWTWrapper::MoveSegment(const std::unordered_map<int, int> &mergedMembership, int site) {
//    fprintf(stderr,"before d:\n");
//    for (int k = 0; k <M ; ++k) {
//        fprintf(stderr,"%d\t",forwardCursor->d[k]);
//    }
//    fprintf(stderr,"\n");
//    fprintf(stderr,"before a:\n");
//    for (int k = 0; k <M ; ++k) {
//        fprintf(stderr,"%d\t",forwardCursor->a[k]);
//    }
//    fprintf(stderr,"\n");

    std::vector<int> tmpD, tmpA;
    int prevState(0), newD(-1);
    int lastHapID;
    int firstHapID;
    int index;
    for (int i = 0; i < (int) clusterMembership.size(); ++i) {
        if (clusterMembership[i].size() != 0) {
            if (prevState == 0 && clusterMembership[prevState].size() == 0) {
                newD = site + 1;
            } else if (i != 0 && prevState != int(i -
                                                  1)) {//if current state is not the first state, and prevState is not the immediate preious state
                //fprintf(stderr,"site:%d\tclusterSize:%d\ti:%d\tprevState:%d\tprev Size:%d\n",site,clusterMembership[i].size(),i,prevState,clusterMembership[prevState].size());
                lastHapID = GetHapIDFromFwd(clusterMembership[prevState].back());//lastHap of previous state
                firstHapID = GetHapIDFromFwd(clusterMembership[i].front());//firstHap of current state
                index = site;
                while (index--) {
                    if (haplotype[lastHapID][index] != haplotype[firstHapID][index]) break;
                }
                newD = index;
            }

            for (int j = 0; j < (int) clusterMembership[i].size(); ++j) {
                if (mergedMembership.find(clusterMembership[i][j]) != mergedMembership.end()) {
                    tmpD.push_back(0);
                } else if (j == 0 && newD != -1) {
                    tmpD.push_back(newD);
                    newD = -1;
                } else {// newD==-1 means no new D value, use the old d value
                    tmpD.push_back(d[clusterMembership[i][j]]);
                }
                tmpA.push_back(a[clusterMembership[i][j]]);
            }
            prevState = i;
        }
    }

//    std::copy(tmpD.begin(), tmpD.end(), d);
//    std::copy(tmpA.begin(), tmpA.end(), a);
    d.assign(tmpD.begin(), tmpD.end());
    a.assign(tmpA.begin(), tmpA.end());

//    fprintf(stderr,"after d:\n");
//    for (int k = 0; k <M ; ++k) {
//        fprintf(stderr,"%d\t",forwardCursor->d[k]);
//    }
//    fprintf(stderr,"\n");
//    fprintf(stderr,"after a:\n");
//    for (int k = 0; k <M ; ++k) {
//        fprintf(stderr,"%d\t",forwardCursor->a[k]);
//    }
//    fprintf(stderr,"\n");
//    for (int j = 0; j < (int)tmpA.size(); ++j) {
//        aMap[site][forwardCursor->a[j]] = j;//at site k, hapID a[k] is stored at jth row/rank
//    }
    return 0;

}


