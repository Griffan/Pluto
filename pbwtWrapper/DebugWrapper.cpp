//
// Created by Fan Zhang on 1/1/16.
//

#include "DebugWrapper.h"
#include <fstream>
#include "ks.h"
#define prefixLength 1200
#define ROUND 10
DebugWrapper::DebugWrapper(int a, int b):PBWTWrapper(a,b) {
    numRight=0;
    numAltRight=0;
}

int DebugWrapper::MergeCluster(int site, RESULT* result) {
    int ret(0);
    int oldNumCluster = GetNumStates(site);
    int numHaps = haplotypeCluster[site].size();

    //std::vector<std::vector<int> > clusterMemberShip(oldNumCluster,std::vector<int>());// state->position in current M array
    std::vector<std::vector<int> > dist(oldNumCluster,std::vector<int>(numHaps,0));//state->rank_occupied
    std::vector<bool> mergeIndicator(oldNumCluster,false);

    std::vector<std::vector<double> > Dmax(oldNumCluster,std::vector<double>(oldNumCluster,0));
    std::vector<std::vector<float> > rankSum(oldNumCluster,std::vector<float>(oldNumCluster,0));
    std::vector<std::vector<std::pair<double,double> > > hapsCounted(oldNumCluster,std::vector<std::pair<double,double> >(oldNumCluster,std::make_pair<int,int>(0,0)));

    if(site==11) {
//        std::ofstream fout("scatter.txt",std::fstream::app);
//        if(!fout.is_open()){fprintf(stderr,"open scatter.txt failed\n");exit(EXIT_FAILURE);}

        int ID;
        int TotalNum=0;
        std::vector<std::vector<int> > tmpRank(clusterMembership.size(), std::vector<int>());
        fprintf(stderr, "\n");
        fprintf(stderr, "membership for site %d:\n", site);
        for (int k = 0; k < clusterMembership.size(); ++k) {
            //   TotalNum+=clusterMembership[k].size();
            fprintf(stderr, "state:%d\n", k);
            for (int i = 0; i < clusterMembership[k].size(); ++i) {
                fprintf(stderr, "%d\t", alphaMap[site + 1][a[site][clusterMembership[k][i]]]);
                tmpRank[k].push_back(alphaMap[site + 1][a[site][clusterMembership[k][i]]]);
            }
            fprintf(stderr, "\n");
            fprintf(stderr, "seqID:%d\n", k);
            for (int i = 0; i < clusterMembership[k].size(); ++i) {
                fprintf(stderr, "%d\t", a[site][clusterMembership[k][i]]);
            }
            fprintf(stderr, "\n");

//            for (int j = 0; j < clusterMembership[k].size(); ++j) {
//                ID=a[site][clusterMembership[k][j]];
//               // if(site==result->indexToPut) {
//
//                    if (result->branchA.find(ID) != result->branchA.end())
//                        fout <<site<<"\t"<<ID<<"\t"<< clusterMembership[k][j] << "\t" <<
//                        alphaMap[site + 1][ID] << "\t" << 1 << "\t"<<site-result->indexToPut<<std::endl;
//                    else if (result->branchB.find(ID) != result->branchB.end())
//                        fout << site<<"\t"<<ID<<"\t"<<clusterMembership[k][j] << "\t" <<
//                        alphaMap[site + 1][ID] << "\t" << 2 <<"\t"<<site-result->indexToPut<< std::endl;
//                    else
//                        fout << site<<"\t"<<ID<<"\t"<<clusterMembership[k][j] << "\t" << alphaMap[site + 1][ID] <<"\t"<<3<<"\t"<<site-result->indexToPut<< std::endl;
//               // }
//               // else
//               //     fout <<site<<"\t"<<ID<<"\t"<< clusterMembership[k][j] << "\t" << alphaMap[site + 1][ID] <<"\t"<<3<< "\t"<<site-result->indexToPut<<std::endl;
//
//            }
        }
        //fprintf(stderr,"TotalNum:%d\n",TotalNum);
        //       fprintf(stderr, "\n");
//        fprintf(stderr, "calculating Dmax:\n");
//        for (int l = 0; l < clusterMembership.size(); ++l) {
//            for (int i = l + 1; i < clusterMembership.size(); ++i) {
//                fprintf(stderr, "state %d and %d: %f and thresh:%f and sample size:%d and %d\n", l, i,
//                        CalDmax(tmpRank[l], tmpRank[i]), 1.36 * std::sqrt(
//                        double(clusterMembership[l].size() + clusterMembership[i].size()) /
//                        (clusterMembership[l].size() * clusterMembership[i].size())), clusterMembership[l].size(),
//                        clusterMembership[i].size());
//            }
//        }
    }

    double tmpABS=0;
    int currentNumCluster=oldNumCluster-1;
    std::unordered_map<int,int> stateOrder;//mapping oldState to newOrder
    int tmpOrder(0);
    std::vector<uchar> tmpAllele;
    double pval(0);

    for (int i = 0; i != numHaps; ++i)//from rank 0 to rank numHaps-1
    {
        //alpha[site][i]:original index of haps at i th place, e.g. david is the ith haplotype
        //haplotypeCluster[site][alpha[site][i]]: david's state status, saying state is michigan
        //dist[haplotypeCluster[site][alpha[site][i]]][i]: state michigan's ith slot is occupied, or where is the ith haplotype located, michigan state's ith slot
        //the structure looks like:
        //state 1:0000001000100
        //state 2:0100100000010
        //state 3:1011010111001
        int hapID = alpha[site + 1][i];//original ID
        int hapState = haplotypeCluster[site][aMap[site][hapID]];
        //std::cerr<<"rank:hapID:hapState\t"<<i<<"\t"<<hapID<<"\t"<<hapState<<std::endl;

        dist[hapState][i] = 1;//record rank occupation indicator for each cluster, haplotypeCluster record state status for haplotype alpha[site][i]]
//    }
//    for (int i = 0; i != numHaps; ++i)//from rank 0 to rank numHaps-1
//    {
        for (int j = 0; j < dist.size(); ++j) {//enumerate through all the states
            if(i==numHaps-1)
            {
                if(mergeIndicator[j])//ToDO: need adjustment for trangle situation
                {
                    continue;
                }
                else
                {
                    tmpAllele.push_back(clusterAllele[site][j]);
                    stateOrder[j]=tmpOrder;
                    tmpOrder++;
                }
            }
            std::vector<int> tmpDist(dist[j]);//record the pivot state
            int tmpMembershipSize=clusterMembership[j].size();

            //if(tmpMembershipSize<5) continue;

            for (int k = j+1; k < dist.size(); ++k)//comparing state j and state k
            {
                if (/*clusterAllele[site][j]!=clusterAllele[site][k] || */mergeIndicator[k]) continue;//if alleles are different or merged once

                double total = tmpMembershipSize + clusterMembership[k].size();
                //double sizeRatio = tmpMembershipSize / total;
                //if (sizeRatio > 0.99 || sizeRatio < 0.01 || clusterMembership[k].size() < 5) continue;


                if (tmpDist[i] == 1)//if j's ith slot is occupied
                {
                    rankSum[j][k] += 1;//rankSum between j and k increase by 1
                    hapsCounted[j][k].first = rankSum[j][k];//local rank between j and k updated and normalized by total rankSum
                    //hapsCounted[j][k].first++;
                }
                else if (dist[k][i] == 1)//if k's ith slot is occupied
                {
                    rankSum[j][k] += 1;
                    hapsCounted[j][k].second = rankSum[j][k];
                    //hapsCounted[j][k].second++;
                }
                else if (i != numHaps - 1)
                    continue;

                tmpABS = fabs(hapsCounted[j][k].first - hapsCounted[j][k].second);//tmp Dmax
                //if(j==0 and k==1)std::cerr<<"hapsCounted["<<j<<"]["<<k<<"]:"<<hapsCounted[j][k].first<<"\t"<<hapsCounted[j][k].second<<"\t"<<tmpABS<<std::endl;

                //_LIBCPP_ASSERT(total,rankSum[j][k]);
                if (tmpABS > Dmax[j][k]) Dmax[j][k] = tmpABS;

#ifdef DEBUG
                if (0  && i == numHaps - 1) {
                    std::cerr << "\nenter rank distribution section, site " << site << ":" << std::endl;
                    for (auto i = 0; i != dist.size(); ++i) {
                        // PrintDistributionAtSite(i,haplotypeCluster[i]);
                        PrintDistributionAtSite(i, dist[i]);
                    }
                    std::cerr << "exit rank distribution section!\n" << std::endl;
                    //exit(0);
                }
#endif
//                if(i==numHaps-1)
//                {
//                    std::cerr<<"Testing state "<<j<<" and "<<k<<std::endl;
//                    //std::cerr<<"hapsCounted["<<j<<"]["<<k<<"]:"<<hapsCounted[j][k].first<<"\t"<<hapsCounted[j][k].second<<std::endl;
//                }

                if (i == numHaps - 1)
                {
                    Dmax[j][k] /= total;

                    pval = 1-psmirnov2x(&Dmax[j][k], tmpMembershipSize, clusterMembership[k].size());

                    fprintf(stderr,"state %d and %d: Pvalue:%f,Dmax:%f,sizeA:%d, sizeB:%d\n",j,k,pval,Dmax[j][k],clusterMembership[j].size(),clusterMembership[k].size());
                }
                if(i==numHaps-1 && 0)//KStest(Dmax[j][k]/total,tmpMembershipSize,clusterMembership[k].size()))//last haplotypes, deal with merging test
                {
                    ret = 1;
                    //PrintVector(dist[i],"state i");
                    //PrintVector(dist[j],"state j");
                    currentNumCluster--;
                    //Merge Action, change mergeIndicator
                    for(int t=0;t!=dist[j].size();++t)
                    {
                        dist[j][t]+=dist[k][t];
                        dist[k][t]=-65534;
                    }
                    mergeIndicator[k]=true;// j th cluster has been merged into i th cluster

                    for(int t=0;t!=clusterMembership[k].size();++t)
                    {
                        haplotypeCluster[site][clusterMembership[k][t]]=j;
                        clusterMembership[j].push_back(clusterMembership[k][t]);
                    }
                    clusterMembership[k].clear();

                }
            }
        }

    }
    if(ret) {
        //adjust d array and a array
       //MoveSegment(clusterMembership);
        //fprintf(stderr,"site:%d merged...\n",site);

        //PrintVector(clusterAllele[site],"allele cluster states after");
        //PrintVector(haplotypeCluster[site],"haplotype cluster states after");
        for (int k = 0; k < haplotypeCluster[site].size(); ++k) {
            haplotypeCluster[site][k] = stateOrder[haplotypeCluster[site][k]];
        }
        clusterAllele[site] = tmpAllele;//update merged cluster allele
        //PrintVector(clusterAllele[site],"allele cluster states final");
        //PrintVector(haplotypeCluster[site],"haplotype cluster states final");
    }
    return ret;
}

int DebugWrapper::MergeCluster(int site, Index2ID& MAP1, ID2POP& MAP2) {
    int ret(0);
    int oldNumCluster = GetNumStates(site);
    int numHaps = haplotypeCluster[site].size();
    int indexToPut=4000;

    if(site%1000==0)fprintf(stderr,"processing site:%d\n",site);
    if(site%5==0&&site>=indexToPut-50&&site<indexToPut+50) {

        std::ofstream fout("/Users/fanzhang/Downloads/PlutoTest/scatter.txt",std::iostream::app);
        if(!fout.is_open()){fprintf(stderr,"open scatter.txt failed\n");exit(EXIT_FAILURE);}

        int ID;
        int TotalNum=0;
        int leftD=0;
        int rightD=0;
        std::vector<std::vector<int> > tmpRank(clusterMembership.size(), std::vector<int>());
//        fprintf(stderr, "\n");
//        fprintf(stderr, "membership for site %d:\n", site);
        for (int k = 0; k < clusterMembership.size(); ++k) {
            TotalNum+=clusterMembership[k].size();
            //fprintf(stderr, "state:%d\n", k);
//            for (int i = 0; i < clusterMembership[k].size(); ++i) {
//                fprintf(stderr, "%d\t", alphaMap[site + 1][a[site][clusterMembership[k][i]]]);
//                tmpRank[k].push_back(alphaMap[site + 1][a[site][clusterMembership[k][i]]]);
//            }
//            fprintf(stderr, "\n");
//            fprintf(stderr, "seqID:%d\n", k);
//            for (int i = 0; i < clusterMembership[k].size(); ++i) {
//                fprintf(stderr, "%d\t", a[site][clusterMembership[k][i]]);
//            }
//            fprintf(stderr, "\n");

            for (int j = 0; j < clusterMembership[k].size(); ++j) {
                ID=a[site][clusterMembership[k][j]];
                leftD=d[site][clusterMembership[k][j]];
                rightD=delta[site][alphaMap[site ][ID]];
                fout << site<<"\t"<<ID<<"\t"<<clusterMembership[k][j] << "\t" << alphaMap[site][ID] <<"\t"<<leftD<<"\t"<<rightD<<"\t"<<3<<"\t"<<site-indexToPut<<"\t"<<MAP1[ID/2]<<"\t";
                if(MAP2.find(MAP1[ID/2])!=MAP2.end()) {
                    fout << MAP2[MAP1[ID / 2]] << std::endl;
                }
                else
                    fout<< "NA"<<std::endl;

            }
        }
        fout.close();
        // fprintf(stderr,"TotalNum:%d\n",TotalNum);
//        fprintf(stderr, "\n");
//        fprintf(stderr, "calculating Dmax:\n");
//        for (int l = 0; l < clusterMembership.size(); ++l) {
//            for (int i = l + 1; i < clusterMembership.size(); ++i) {
//                fprintf(stderr, "state %d and %d: %f and thresh:%f and sample size:%d and %d\n", l, i,
//                        CalDmax(tmpRank[l], tmpRank[i]), 1.36 * std::sqrt(
//                                double(clusterMembership[l].size() + clusterMembership[i].size()) /
//                                (clusterMembership[l].size() * clusterMembership[i].size())), clusterMembership[l].size(),
//                        clusterMembership[i].size());
//            }
//        }
    }


    return 0;//debug
}

int DebugWrapper::CursorForwards(RESULT* result) {//so far only implemented for test purpose


    //PrintVector(forwardCursor->a,M,"end arrary aFend check 0");

    for (int k = 0; k != pbwtCore->N; ++k) {
        //fprintf(stderr,"at site %d\n",k);
        CursorForwardsTo(k, prefixLength, result);
    }
    //copy end of a to PBWT
    //PrintVector(forwardCursor->a,M,"end arrary aFend check 1");

    pbwtCursorToAFend(forwardCursor, pbwtCore);

    for (int i=0;i != pbwtCore->N; i++) {

        UpdateTransVector(i);
    }
    PrintSummary();
    //update crossover rate?
    return 0;
}


int DebugWrapper::CursorForwards(Index2ID& MAP1,ID2POP& MAP2) {//so far only implemented for test purpose


    //PrintVector(forwardCursor->a,M,"end arrary aFend check 0");

    for (int k = 0; k != pbwtCore->N; ++k) {
        //fprintf(stderr,"at site %d\n",k);
        CursorForwardsTo(k,  MAP1, MAP2,prefixLength);
    }
    //copy end of a to PBWT
    //PrintVector(forwardCursor->a,M,"end arrary aFend check 1");

    pbwtCursorToAFend(forwardCursor, pbwtCore);

//    for (int i=0;i != pbwtCore->N; i++) {
//
//        UpdateTransVector(i);
//    }
//    PrintSummary();
    //update crossover rate?
    return 0;
}

int DebugWrapper::CursorForwardsTo(int k, int T, RESULT *result) {
/*T is the length that how far you look back
 *This function must be called along the sites, no skip permitted;
 *Mask the site you want to skip at the begining if you have to.
 */
    int i, i0 = 0, ia;
    int group = 0;
    int *tmpA,*tmpD;
    tmpA=new int [M];
    tmpD=new int [M];
    clusterMembership.clear();
    /*coppy array d*/
    // int *lastD = new int[forwardCursor->M + 1];
    // uchar *lastY = new uchar[forwardCursor->M +1];
    // memcpy(lastD, forwardCursor->d, (forwardCursor->M + 1) * sizeof(int));
    //memcpy(lastY, forwardCursor->sortedY, (forwardCursor->M + 1) * sizeof(uchar));
#ifdef DEBUG
    if((k>=64||k==63||k==62||k==61)) {
        fprintf(stderr, "last a array %d:\n",k);
        for (int kt = 0; kt < M; ++kt) {
            fprintf(stderr, "%d\t", forwardCursor->a[kt]);
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "last Y array:%d\n",k);
        for (int kt = 0; kt < M; ++kt) {
            fprintf(stderr, "%d\t", forwardCursor->sortedY[kt]);
        }
        fprintf(stderr, "\n");
    }
#endif
    //copy haplotypes into forwardCursor->y
    CopyHap(k, forwardCursor);

    if (k==pbwtCore->N-1)//deal with last columns
    {
        for (i = 0; i < forwardCursor->M; ++i) {
            haplotypeCluster[k][i] = 0;
        }
        clusterAllele[k].push_back(0);

    }

    int tmpT= k > T ? T : k;

    /*reprot haolotype cluster based on prefix, so current site not included*/
    //fprintf(stderr,"k:%d,T:%d\t",k,T);PrintVector(forwardCursor->d,forwardCursor->M,"before tmpD");
    //fprintf(stderr,"k:%d,T:%d\t",k,T);PrintVector(forwardCursor->a,forwardCursor->M,"olda");
    //fprintf(stderr,"k:%d,T:%d\t",k,T);PrintVector(forwardCursor->sortedY,forwardCursor->M,"sortedY");
    int u = 0, v = 0;
    int p = k + 1;
    int q = k + 1;

    for (i = 0; i < forwardCursor->M; ++i) {

        if (forwardCursor->d[i] > p) p = forwardCursor->d[i];
        if (forwardCursor->d[i] > q) q = forwardCursor->d[i];

        /*assign states of last column based on previous d and sortedY*/
        if (forwardCursor->d[i] > (k - tmpT)&& k !=0) {//if current sequence and last sequence have common sequence longer than T
            //if (na && nb)        /* then there is something to report */
            if(i!=0) {
                //fprintf(stderr,"d:%d\tk-tmpT:%d\n",forwardCursor->d[i],k-tmpT);
                std::vector<int> tmpMem;
                for (ia = i0; ia < i; ++ia) {
                    haplotypeCluster[k - 1][ia] = group;
                    tmpMem.push_back(ia);
                }

                clusterAllele[k - 1].push_back(haplotype[forwardCursor->a[i0]][k-1]);
                clusterMembership.push_back(tmpMem);
                // na = 0;
                // nb = 0;
                i0 = i;
                group++;

            }

        }

        if (forwardCursor->sortedY[i] == 0) {
            //forwardCursor->a[u] = forwardCursor->a[i];
            //forwardCursor->d[u] = p;
            tmpA[u] = forwardCursor->a[i];
            tmpD[u] = p;
            ++u;
            p = 0;
            // na++;
            forwardCursor->c++;
        }
        else {
            forwardCursor->b[v] = forwardCursor->a[i];
            forwardCursor->e[v] = q;
            ++v;
            q = 0;
            //nb++;
        }
    }

    if( i0 < forwardCursor->M)
    {
        if (k != 0)
        {
            std::vector<int> tmpMem;
            for (ia = i0; ia < forwardCursor->M; ++ia)
            {
                haplotypeCluster[k-1][ia] = group;
                tmpMem.push_back(ia);
            }
            clusterAllele[k-1].push_back(haplotype[forwardCursor->a[i0]][k-1]);
            clusterMembership.push_back(tmpMem);
            //fprintf(stderr,"site:%d\tnumStates:%d\n",k-1,clusterAllele[k-1].size());

        }

    }
    int test = 0;
#ifdef DEBUG
    if(k!=0){
        fprintf(stderr, "at site:%d\n", k-1);
        PrintVector(haplotypeCluster[k-1],"haplotype state before merge state");
        PrintVector(clusterAllele[k-1],"state allele before merge allele");
    }
#endif
    if(k!=0&&clusterAllele[k-1].size()==0) {fprintf(stderr,"0 states, abort!");abort();}
    if(k!=0&&clusterAllele[k-1].size()!=1) test=MergeCluster(k-1, result);//TODO:implement this function
#ifdef DEBUG
    if(k!=0&&test) {
        fprintf(stderr, "at site:%d\n", k-1);

        PrintVector(haplotypeCluster[k-1], "haplotype state after merge state");
        PrintVector(clusterAllele[k-1],"state allele after merge allele");
        fprintf(stderr,"\n");
    }
#endif
    //numCluster[k] = group;
    //forwardCursor->c = na;
    //numZero[k]=na;
    memcpy(forwardCursor->a , tmpA, u * sizeof(int));
    memcpy(forwardCursor->d , tmpD, u * sizeof(int));
    memcpy(forwardCursor->a + u, forwardCursor->b, v * sizeof(int));
    memcpy(forwardCursor->d + u, forwardCursor->e, v * sizeof(int));
    //forwardCursor->d[0] = k + 2;
    forwardCursor->d[forwardCursor->M] = k + 2; /* sentinels */
    a[k].assign(forwardCursor->a,forwardCursor->a+forwardCursor->M);

    for (int j = 0; j <a[k].size() ; ++j) {
        //std::cerr<<"a size:"<<a[k].size()<<" and "<<a[k][j]<<std::endl;
        //aMap[k].insert(std::make_pair(a[k][j],j));
        aMap[k][a[k][j]]=j;
    }
    d[k].assign(forwardCursor->d,forwardCursor->d+forwardCursor->M);
    //sortedY[k].assign(forwardCursor->sortedY,forwardCursor->sortedY+forwardCursor->M);
    //delete [] lastD;
    //delete [] lastY;

    //pbwtCursorForwardsReadAD(forwardCursor, k);
    // updateCursorForwards();//
    //fprintf(stderr,"k:%d,T:%d\t",k,T);PrintVector(forwardCursor->d,forwardCursor->M,"after tmpD");

    return 0;
}


int DebugWrapper::CursorForwardsTo(int k, Index2ID &MAP1, ID2POP& MAP2, int T) {
/*T is the length that how far you look back
 *This function must be called along the sites, no skip permitted;
 *Mask the site you want to skip at the begining if you have to.
 */
    int i, i0 = 0, ia;
    int group = 0;
    int *tmpA,*tmpD;
    tmpA=new int [M];
    tmpD=new int [M];
    clusterMembership.clear();

    //copy haplotypes into forwardCursor->y
    CopyHap(k, forwardCursor);

    if (k==pbwtCore->N-1)//deal with last columns
    {
        for (i = 0; i < forwardCursor->M; ++i) {
            haplotypeCluster[k][i] = 0;
        }
        clusterAllele[k].push_back(0);

    }

    int tmpT= k > T ? T : k;

    int u = 0, v = 0;
    int p = k + 1;
    int q = k + 1;

    for (i = 0; i < forwardCursor->M; ++i) {

        if (forwardCursor->d[i] > p) p = forwardCursor->d[i];
        if (forwardCursor->d[i] > q) q = forwardCursor->d[i];

        /*assign states of last column based on previous d and sortedY*/
        if (forwardCursor->d[i] > (k - tmpT)&& k !=0) {//if current sequence and last sequence have common sequence longer than T

            if(i!=0) {
                std::vector<int> tmpMem;
                for (ia = i0; ia < i; ++ia) {
                    haplotypeCluster[k - 1][ia] = group;
                    tmpMem.push_back(ia);
                }

                clusterAllele[k - 1].push_back(haplotype[forwardCursor->a[i0]][k-1]);
                clusterMembership.push_back(tmpMem);
                i0 = i;
                group++;

            }

        }

        if (forwardCursor->sortedY[i] == 0) {
            //forwardCursor->a[u] = forwardCursor->a[i];
            //forwardCursor->d[u] = p;
            tmpA[u] = forwardCursor->a[i];
            tmpD[u] = p;
            ++u;
            p = 0;
            // na++;
            forwardCursor->c++;
        }
        else {
            forwardCursor->b[v] = forwardCursor->a[i];
            forwardCursor->e[v] = q;
            ++v;
            q = 0;
            //nb++;
        }
    }

    if( i0 < forwardCursor->M)
    {
        if (k != 0)
        {
            std::vector<int> tmpMem;
            for (ia = i0; ia < forwardCursor->M; ++ia)
            {
                haplotypeCluster[k-1][ia] = group;
                tmpMem.push_back(ia);
            }
            clusterAllele[k-1].push_back(haplotype[forwardCursor->a[i0]][k-1]);
            clusterMembership.push_back(tmpMem);
        }

    }
    int test = 0;
#ifdef DEBUG
    if(k!=0){
        fprintf(stderr, "at site:%d\n", k-1);
        PrintVector(haplotypeCluster[k-1],"haplotype state before merge state");
        PrintVector(clusterAllele[k-1],"state allele before merge allele");
    }
#endif
    if(k!=0&&clusterAllele[k-1].size()==0) {fprintf(stderr,"0 states, abort!");abort();}
    if(k!=0&&clusterAllele[k-1].size()!=1) test=MergeCluster(k-1, MAP1, MAP2);//TODO:implement this function
#ifdef DEBUG
    if(k!=0&&test) {
        fprintf(stderr, "at site:%d\n", k-1);
        PrintVector(haplotypeCluster[k-1], "haplotype state after merge state");
        PrintVector(clusterAllele[k-1],"state allele after merge allele");
        fprintf(stderr,"\n");
    }
#endif
    memcpy(forwardCursor->a , tmpA, u * sizeof(int));
    memcpy(forwardCursor->d , tmpD, u * sizeof(int));
    memcpy(forwardCursor->a + u, forwardCursor->b, v * sizeof(int));
    memcpy(forwardCursor->d + u, forwardCursor->e, v * sizeof(int));
    //forwardCursor->d[0] = k + 2;
    forwardCursor->d[forwardCursor->M] = k + 2; /* sentinels */
    a[k].assign(forwardCursor->a,forwardCursor->a+forwardCursor->M);

    for (int j = 0; j <a[k].size() ; ++j) {
        //std::cerr<<"a size:"<<a[k].size()<<" and "<<a[k][j]<<std::endl;
        //aMap[k].insert(std::make_pair(a[k][j],j));
        aMap[k][a[k][j]]=j;
    }
    d[k].assign(forwardCursor->d,forwardCursor->d+forwardCursor->M);

    return 0;
}


int DebugWrapper::ConfidentOrNot(char** individual, int siteA, int siteB) {
    if(siteA>siteB) swap(siteA,siteB);
    //find maximal length of shared prefix
    int length=siteB+1;
    char* haps[2];
    haps[0] = new char [length];//start from 0, hence siteB is length of prefix

    memcpy(haps[0],individual[0],length);

    char* altHaps[2];
    altHaps[0] = new char [length];
    //altHaps[1] = new char [length];
    memcpy(altHaps[0],individual[0],length);

    altHaps[0][siteA]=individual[1][siteA];

    int leftLengthA=FindLengthOfPrefix(haps[0],siteB);

    int altLeftLengthA=FindLengthOfPrefix(altHaps[0],siteB);

    delete haps[0];

    delete altHaps[0];

    //find maximal length of shared suffix
    length=N-siteA;
    haps[0] = new char [length];//start from 0, hence siteB is length of prefix

    memcpy(haps[0],individual[0]+siteA,length);

    altHaps[0] = new char [length];

    memcpy(altHaps[0],&individual[0][siteA],length);

    altHaps[0][siteB-siteA] = individual[1][siteB];
    int rightLengthA=FindLengthOfSuffix(haps[0],siteA);

    int altRightLengthA=FindLengthOfSuffix(altHaps[0],siteA);


    int ret(0);
//    if(leftLengthA >= altLeftLengthA && rightLengthA >= altRightLengthA) {/*numRight++; cerr<<"so far numRight:"<<numRight<<endl;*/ ret|=0x11;}//confident and right
//    if(altLeftLengthA > leftLengthA && altRightLengthA > rightLengthA) {/*numAltRight++;cerr<<"so far numAltRight:"<<numAltRight<<endl;*/ ret|=0x1;}//confident  and wrong
//    if(leftLengthA * rightLengthA >= altLeftLengthA*altRightLengthA) ret|=0x1100;//ambiguous and right
//    if(altLeftLengthA * altRightLengthA > leftLengthA * rightLengthA) ret|=0x100;//ambiguous and wrong
    //cerr<<siteA<<"\t"<<siteB<<"\t"<<leftLengthA<<"\t"<<altLeftLengthA<<"\t"<<rightLengthA<<"\t"<<altRightLengthA<<endl;
    if (leftLengthA == altLeftLengthA && rightLengthA == altRightLengthA)
    {
        ret=rand()%2;
        if(ret)
        {
            ret=0;
            ret|=0x11;
        }
        else
        {
            ret=0;
            ret|=0x1;
        }
    } //then it is ambiguous and use random guess or fractional counts
    else if (leftLengthA >= altLeftLengthA && rightLengthA >= altRightLengthA)
    {
        if((haps[0][siteA-siteA]==this->haplotype[M][siteA]&&haps[0][siteB-siteA]==this->haplotype[M][siteB])||(haps[0][siteA-siteA]==this->haplotype[M+1][siteA]&&haps[0][siteB-siteA]==this->haplotype[M+1][siteB]))//haps is the same with ref
        {
            ret|=0x11;
            //cerr<<siteA<<"\t"<<siteB<<"\tphasing via hap"<<endl;
        }//confident right
        else
            ret|=0x1;//confident wrong
    }
    else if (leftLengthA <= altLeftLengthA && rightLengthA <= altRightLengthA)
    {
        if((altHaps[0][siteA-siteA]==this->haplotype[M][siteA]&&altHaps[0][siteB-siteA]==this->haplotype[M][siteB])||(altHaps[0][siteA-siteA]==this->haplotype[M+1][siteA]&&altHaps[0][siteB-siteA]==this->haplotype[M+1][siteB]))
        {
            ret|=0x11;
            //cerr<<siteA<<"\t"<<siteB<<"\tphasing via althap"<<endl;

        }//confident right
        else
            ret|=0x1;//confident wrong
    }
    else if (leftLengthA * rightLengthA == altLeftLengthA * altRightLengthA )
    {
        ret=rand()%2;
        if(ret)
        {
            ret=0;
            ret|=0x1100;
        }
        else
        {
            ret=0;
            ret|=0x100;
        }
    }//then it is ambiguous and use random guess or fractional counts
    else if (leftLengthA * rightLengthA > altLeftLengthA * altRightLengthA)
    {
        if((haps[0][siteA-siteA]==this->haplotype[M][siteA]&&haps[0][siteB-siteA]==this->haplotype[M][siteB])||(haps[0][siteA-siteA]==this->haplotype[M+1][siteA]&&haps[0][siteB-siteA]==this->haplotype[M+1][siteB]))//haps is the same with ref
        {
            ret|=0x1100;
            //cerr<<siteA<<"\t"<<siteB<<"\tphasing via ambiguous hap"<<endl;
        }
        else
            ret|=0x100;
    }// then it is ambiguous right
    else if (leftLengthA * rightLengthA < altLeftLengthA * altRightLengthA)
    {
        if((altHaps[0][siteA-siteA]==this->haplotype[M][siteA]&&altHaps[0][siteB-siteA]==this->haplotype[M][siteB])||(altHaps[0][siteA-siteA]==this->haplotype[M+1][siteA]&&altHaps[0][siteB-siteA]==this->haplotype[M+1][siteB]))
        {
            ret |= 0x1100;
            //cerr<<siteA<<"\t"<<siteB<<"\tphasing via ambiguous althap"<<endl;
        }
        else
            ret|=0x100;
    }// then it is ambiguous wrong
    else { std::cerr<<"matching length error"<<std::endl;exit(EXIT_FAILURE);}//spit error message.


    delete haps[0];

    delete altHaps[0];

    return ret;
}


int DebugWrapper::Phase(char** individual, int siteA, int siteB) {
    if(siteA>siteB) swap(siteA,siteB);
    //find maximal length of shared prefix
    int length=siteB+1;
    char* haps[2];
    haps[0] = new char [length];//start from 0, hence siteB is length of prefix

    memcpy(haps[0],individual[0],length);

    char* altHaps[2];
    altHaps[0] = new char [length];
    //altHaps[1] = new char [length];
    memcpy(altHaps[0],individual[0],length);

    altHaps[0][siteA]=individual[1][siteA];

    int leftLengthA=FindLengthOfPrefix(haps[0],siteB);

    int altLeftLengthA=FindLengthOfPrefix(altHaps[0],siteB);

    delete haps[0];

    delete altHaps[0];

    //find maximal length of shared suffix
    length=N-siteA;
    haps[0] = new char [length];//start from 0, hence siteB is length of prefix

    memcpy(haps[0],individual[0]+siteA,length);

    altHaps[0] = new char [length];

    memcpy(altHaps[0],&individual[0][siteA],length);

    altHaps[0][siteB-siteA] = individual[1][siteB];
    int rightLengthA=FindLengthOfSuffix(haps[0],siteA);

    int altRightLengthA=FindLengthOfSuffix(altHaps[0],siteA);


    int ret(0);
//    if(leftLengthA >= altLeftLengthA && rightLengthA >= altRightLengthA) {/*numRight++; cerr<<"so far numRight:"<<numRight<<endl;*/ ret|=0x11;}//confident and right
//    if(altLeftLengthA > leftLengthA && altRightLengthA > rightLengthA) {/*numAltRight++;cerr<<"so far numAltRight:"<<numAltRight<<endl;*/ ret|=0x1;}//confident  and wrong
//    if(leftLengthA * rightLengthA >= altLeftLengthA*altRightLengthA) ret|=0x1100;//ambiguous and right
//    if(altLeftLengthA * altRightLengthA > leftLengthA * rightLengthA) ret|=0x100;//ambiguous and wrong

    if (leftLengthA == altLeftLengthA && rightLengthA == altRightLengthA)
    {

    } //then it is ambiguous and use random guess or fractional counts
    else if (leftLengthA >= altLeftLengthA && rightLengthA >= altRightLengthA)
    {

            ret|=0x11;

    }
    else if (leftLengthA <= altLeftLengthA && rightLengthA <= altRightLengthA)
    {

            swap(individual[1][siteB],individual[0][siteB]);

            ret|=0x1;//confident wrong
    }
    else if (leftLengthA * rightLengthA == altLeftLengthA * altRightLengthA )
    {

    }//then it is ambiguous and use random guess or fractional counts
    else if (leftLengthA * rightLengthA > altLeftLengthA * altRightLengthA)
    {

            ret|=0x1100;

    }// then it is ambiguous right
    else if (leftLengthA * rightLengthA < altLeftLengthA * altRightLengthA)
    {

            swap(individual[1][siteB],individual[0][siteB]);
            ret|=0x100;
    }// then it is ambiguous wrong
    else { std::cerr<<"matching length error"<<std::endl;exit(EXIT_FAILURE);}//spit error message.


    delete haps[0];

    delete altHaps[0];

    return ret;
}

int DebugWrapper::FindLengthOfSuffix(char *haplotypeX,int siteA) {

    int maxLen(0);

    /* match query in turn */


        int f = 0;
        int g = M-1;


        int f1 = 0;
        int g1 = M-1;

        for (int k = siteA; k < N; ++k) {               /* use classic FM updates to extend [f,g) interval to next position */
            f1 = haplotypeX[k-siteA] ? c[k] + (f - u[k][f]) : u[k][f];
            g1 = haplotypeX[k-siteA] ? c[k] + (g - u[k][g]) : u[k][g];
            //cerr<<"suffix->k:"<<k<<"\tg:"<<g<<endl;
            //if(k>=c.size()||k<0) cerr<<"c size:"<<c.size()<<"\tk:"<<k<<endl;
            /*if(k>u.size()||f>u[k].size()||g1>u[k].size())*/
            //cerr<<"M:"<<M<<"\tSuffix->k:"<<k<<"\t"<<u[k].size()<<"\tf:"<<f<<"\tg1:"<<g1<<"\tg:"<<g<<"\t1st:"<<c[k]<<"\t2st:"<<u[k][g]<<endl;
            if(f1<0) f1=0;
            if(g1>M-1) g1=M-1;
            /* if the interval is non-zero we can just proceed */
            if (g1 > f1) {
                f = f1;
                g = g1;
            } /* no change to e */
            else        /* we have reached a maximum - need to report and update e, f*,g* */
            {
                break;
            }
            maxLen++;
        }

    return maxLen;
}

int DebugWrapper::FindLengthOfPrefix(char *haplotypeX, int siteB) {
    int maxLen(0);

    /* match query in turn */

    int f = 0;
    int g = M-1;

    //cerr<<"in DebugWraooer M:"<<M<<endl;
    int f1 = 0;
    int g1 = M-1;


    for (int k = siteB; k > -1; --k) {               /* use classic FM updates to extend [f,g) interval to next position */
        f1 = haplotypeX[k] ? celta[k] + (f - ultra[k][f]) : ultra[k][f];
        g1 = haplotypeX[k] ? celta[k] + (g - ultra[k][g]) : ultra[k][g];
        //cerr<<"prefix->k:"<<k<<"\tg:"<<g<<endl;
        //cerr<<celta[k]<<"\t"<<ultra[k][f]<<"\t"<<ultra[k][g]<<endl;
        //if(k>=celta.size()||k<0) cerr<<"celta size:"<<celta.size()<<"\tk:"<<k<<endl;
        /*if(k>ultra.size()||f>ultra[k].size()||g1>ultra[k].size())*/
        //cerr<<"M"<<M<<"Prefix->k:"<<k<<"\t"<<ultra[k].size()<<"\tf:"<<f<<"\tg:"<<g<<"\tg1:"<<g1<<"\t1st:"<<celta[k]<<"\t2nd:"<<ultra[k][g]<<endl;
        //if(f1<0||g1>M-1) continue;
        if(f1<0) f1=0;
        if(g1>M-1) g1=M-1;
        /* if the interval is non-zero we can just proceed */
        if (g1 > f1) {
            f = f1;
            g = g1;
        } /* no change to e */
        else        /* we have reached a maximum - need to report and update e, f*,g* */
        {
            break;
        }

        maxLen++;

    }

    return maxLen;
}

char **DebugWrapper::ExtractSubset(int individual) {
    char* tmpHapA,*tmpHapB;
    tmpHapA = this->haplotype[individual*2];
    tmpHapB = this->haplotype[individual*2+1];
    this->haplotype[individual*2] = this->haplotype[M-2];
    this->haplotype[individual*2+1] = this->haplotype[M-1];
    this->haplotype[M-2] = tmpHapA;
    this->haplotype[M-1] = tmpHapB;
    char** haplotypeX = new char*[2];
    //haplotypeX[0] = tmpHapA;
    //haplotypeX[1] = tmpHapB;
    haplotypeX[0]=new char [N];
    haplotypeX[1]=new char [N];
    memcpy(haplotypeX[0],this->haplotype[M-2],N);
    memcpy(haplotypeX[1],this->haplotype[M-1],N);

//    srand(time(NULL));
//    for (int i = 0; i <N ; ++i) {
//        if(rand()%2==1)
//        swap(haplotypeX[0][i],haplotypeX[1][i]);
//    }
    M-=2;
    resetWrapper();
    return haplotypeX;
}

int DebugWrapper::SubsetResume(char **haplotypeX, int individual) {

    M+=2;
    char* tmpHapA,*tmpHapB;
    tmpHapA = this->haplotype[individual*2];
    tmpHapB = this->haplotype[individual*2+1];
    this->haplotype[individual*2] = this->haplotype[M-2];
    this->haplotype[individual*2+1] = this->haplotype[M-1];
    this->haplotype[M-2] = tmpHapA;
    this->haplotype[M-1] = tmpHapB;
    delete haplotypeX[0];
    delete haplotypeX[1];
    delete haplotypeX;

    return 0;
}

struct CONRE
{
    int confidentRight;
    int ambiguousRight;
    int confidentTotal;
    int ambiguousTotal;
    CONRE()
    {
        confidentRight=0;
        ambiguousRight=0;
        confidentTotal=0;
        ambiguousTotal=0;
    }
};


int DebugWrapper::Process(int nMarkers, int nSamples, char** haps) {

    std::vector<CONRE> conreIndividual(nSamples,CONRE());
    std::unordered_map<int,CONRE> conreSNP;
    char*hapInUse[2];
    hapInUse[0]= new char [N];
    hapInUse[1]= new char [N];

    for (int individual = 0; individual < 5/*nSamples*/; ++individual) {

            char **haplotypeX = ExtractSubset(individual);

            fprintf(stderr, "finished initializing graph\n");

        CursorBackwards();
            fprintf(stderr, "finished backward procedure\n");
        PBWTWrapper::CursorForwards();
            fprintf(stderr, "finished forward procedure\n");

            std::vector<int> heterIndex;
            for (int i = 0; i < N; ++i) {
                if (haplotypeX[0][i] != haplotypeX[1][i]) heterIndex.push_back(i);
            }
            fprintf(stderr, "number of heterozygous site for individual %d :%d\n", individual, heterIndex.size());

        //phasing multiple rounds
        char*** voteHapVec = new char** [ROUND];
        int seed=time(NULL);
        srand(seed);
        //std::cerr<<"Seed:"<<seed<<std::endl;
        for (int k = 0; k <ROUND; ++k) {//phasing round
            voteHapVec[k]=new char* [2];
            voteHapVec[k][0]=new char [N];
            voteHapVec[k][1]=new char [N];

            memcpy(hapInUse[0], haplotypeX[0], N);
            memcpy(hapInUse[1], haplotypeX[1], N);

            for (int i = 0; i <N ; ++i) {
                if(rand()%2==1)
                    swap(hapInUse[0][i],hapInUse[1][i]);
            }
            for (int j = 0; j < heterIndex.size() - 1; ++j) {
                Phase(hapInUse, heterIndex[j], heterIndex[j + 1]);
            }
            memcpy(voteHapVec[k][0],hapInUse[0],N);
            memcpy(voteHapVec[k][1],hapInUse[1],N);
        }
////iterative updating
//        for (int i = 0; i <N ; ++i) {
//                if(rand()%2==1)
//                    swap(haplotypeX[0][i],haplotypeX[1][i]);
//        }
//        for (int k = 0; k <ROUND; ++k) {//phasing round
//            voteHapVec[k]=new char* [2];
//            voteHapVec[k][0]=new char [N];
//            voteHapVec[k][1]=new char [N];
//
//            for (int j = 0; j < heterIndex.size() - 1; ++j) {
//                Phase(haplotypeX, heterIndex[j], heterIndex[j + 1]);
//            }
//            memcpy(voteHapVec[k][0],haplotypeX[0],N);
//            memcpy(voteHapVec[k][1],haplotypeX[1],N);
//        }
        //major voting
        MajorVoting(haplotypeX,voteHapVec,ROUND,heterIndex);

        for (int l = 0; l < ROUND; ++l) {
            delete [] voteHapVec[l][0];
            delete [] voteHapVec[l][1];
        }

        for (int j = 0; j < heterIndex.size()-1; ++j) {
            int ret = ConfidentOrNot(haplotypeX,heterIndex[j],heterIndex[j+1]);
            if((ret&0x10)&&(ret&0x1))
            {
                conreIndividual[individual].confidentTotal+=1;
                conreIndividual[individual].confidentRight+=1;
                if(conreSNP.find(heterIndex[j])!=conreSNP.end())
                {
                    conreSNP[heterIndex[j]].confidentTotal+=1;
                    conreSNP[heterIndex[j]].confidentRight+=1;
                }
                else
                {
                    CONRE conre;
                    conre.confidentTotal=1;
                    conre.confidentRight=1;
                    conreSNP.insert(make_pair(heterIndex[j],conre));
                }
            }
            else if(ret&0x1)
            {
                conreIndividual[individual].confidentTotal+=1;
                if(conreSNP.find(heterIndex[j])!=conreSNP.end())
                {
                    conreSNP[heterIndex[j]].confidentTotal+=1;
                }
                else
                {
                    CONRE conre;
                    conre.confidentTotal=1;
                    conreSNP.insert(make_pair(heterIndex[j],conre));
                }
            }
            if((ret&0x1000)&&(ret&0x100))
            {
                conreIndividual[individual].ambiguousRight+=1;
                conreIndividual[individual].ambiguousTotal+=1;
                if(conreSNP.find(heterIndex[j])!=conreSNP.end())
                {
                    conreSNP[heterIndex[j]].ambiguousRight+=1;
                    conreSNP[heterIndex[j]].ambiguousTotal+=1;
                }
                else
                {
                    CONRE conre;
                    conre.ambiguousRight=1;
                    conre.ambiguousTotal=1;
                    conreSNP.insert(make_pair(heterIndex[j],conre));
                }
            }
            else if(ret&0x100)
            {
                conreIndividual[individual].ambiguousTotal+=1;
                if(conreSNP.find(heterIndex[j])!=conreSNP.end())
                {
                    conreSNP[heterIndex[j]].ambiguousTotal+=1;
                }
                else
                {
                    CONRE conre;
                    conre.ambiguousTotal=1;
                    conreSNP.insert(make_pair(heterIndex[j],conre));
                }
            }

        }

        SubsetResume(haplotypeX,individual);

    }

    ofstream indvOUT("/Users/fanzhang/Downloads/PlutoTest/individual.out");
    ofstream SNPOUT("/Users/fanzhang/Downloads/PlutoTest/snp.out");
    for (int k = 0; k <nSamples ; ++k) {
        indvOUT<<k<<"\t"<<conreIndividual[k].confidentTotal<<"\t"<<conreIndividual[k].confidentRight<<"\t"<<conreIndividual[k].ambiguousTotal<<"\t"<<conreIndividual[k].ambiguousRight<<std::endl;
    }
    for (std::unordered_map<int,CONRE>::iterator iter = conreSNP.begin(); iter!=conreSNP.end(); ++iter) {
        SNPOUT<<iter->first<<"\t"<<iter->second.confidentTotal<<"\t"<<iter->second.confidentRight<<"\t"<<iter->second.ambiguousTotal<<"\t"<<iter->second.ambiguousRight<<endl;
    }

    conreSNP.clear();
    delete [] hapInUse[0];
    delete [] hapInUse[1];
    return 0;
}

int DebugWrapper::MajorVoting(char **individual, char ***voteHapVec, int round,std::vector<int>& heterIndex) {
    std::cerr<<"Enter Major Voting..."<<std::endl;
    int vote(0);
    int tmp1(0),tmp2(0);
    for (int j = 0; j <heterIndex.size()-1; ++j) {
        vote=0;
        for (int i = 0; i < round; ++i) {
           // std::cerr<<"Site:"<<heterIndex[j]<<"\ttmp1-tmp2:"<<tmp1<<tmp2<<"\tRound:"<<i<<"\t hap1:"<<(int)voteHapVec[i][0][heterIndex[j]]<<(int)voteHapVec[i][0][heterIndex[j+1]]<<"\thap2:"<<(int)voteHapVec[i][1][heterIndex[j]]<<(int)voteHapVec[i][1][heterIndex[j+1]]<<"\tvote:"<<vote<<std::endl;
            if(i==0)
            {
                tmp1=voteHapVec[i][0][heterIndex[j]];
                tmp2=voteHapVec[i][0][heterIndex[j+1]];
                vote=1;
            }
            else if(tmp1==voteHapVec[i][0][heterIndex[j]] && tmp2==voteHapVec[i][0][heterIndex[j+1]])
            {
                vote++;
            }
            else if(tmp1==voteHapVec[i][1][heterIndex[j]] && tmp2==voteHapVec[i][1][heterIndex[j+1]])
            {
                vote++;
            }
        }
        //std::cerr<<"vote:"<<vote<<"\tM/2:"<<round/2<<std::endl;
        if(vote > round/2)
        {
            //std::cerr<<"voting success at site "<<heterIndex[j]<<std::endl;
            individual[0][heterIndex[j]]=voteHapVec[0][0][heterIndex[j]];
            individual[0][heterIndex[j+1]]=voteHapVec[0][0][heterIndex[j+1]];
            individual[1][heterIndex[j]]=voteHapVec[0][1][heterIndex[j]];
            individual[1][heterIndex[j+1]]=voteHapVec[0][1][heterIndex[j+1]];
        }
        else
        {
            individual[0][heterIndex[j]]=voteHapVec[0][0][heterIndex[j]];
            individual[0][heterIndex[j+1]]=voteHapVec[0][1][heterIndex[j+1]];
            individual[1][heterIndex[j]]=voteHapVec[0][0][heterIndex[j]];
            individual[1][heterIndex[j+1]]=voteHapVec[0][1][heterIndex[j+1]];
        }
    }
    return 0;
}
