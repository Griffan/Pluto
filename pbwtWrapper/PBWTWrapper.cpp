//
// Created by Fan Zhang on 7/20/15.
//

#include <unordered_map>
#include <map>
#include <cmath>
#include "PBWTWrapper.h"
#include "iostream"
#include "pbwt/pbwt.h"
#include <fstream>
#include "math.h"
#include "ks.h"
//#include <numeric>

//PBWTWrapper::PBWTWrapper(const char **haplotype, int nhaps, int nsnps) {
//
//    int version = 2;
//    printf("Read %i SNPs %i haplotypes and %i individuals from PHASE format version %i\n", nsnps, nhaps, nhaps / 2,
//           version);
//    PBWT *p = pbwtCreate(nhaps, nsnps);
//    //p->chrom = strdup (fgetword(fp)) ; /* example 4th line is P followed by site positions */
//    //p->sites = arrayCreate (4096, Site) ;
//    //int i ; for (i = 0 ; i < p->N ; ++i) arrayp(p->sites,i,Site)->x = atoi (fgetword(fp)) ;
//
//    int i, j;
//    p->CompressedAllele = arrayCreate(4096 * 32, uchar);
//    PbwtCursor *u = pbwtCursorCreate(p, TRUE, TRUE);
//    for (i = 0; i < p->N; ++i) {
//        for (j = 0; j < p->M; ++j)
//            u->sortedY[j] = haplotype[i][u->a[j]];
//        pbwtCursorWriteForwards(u);
//        if (nCheckPoint && !((i + 1) % nCheckPoint)) pbwtCheckPoint(u, p);
//    }
//    pbwtCursorToAFend(u, p);
//
//    fprintf(stderr, "transform phase file");
//    //if (p->chrom) fprintf(stderr, " for chromosome %s", p->chrom);
//    fprintf(stderr, ": M, N are\t%d\t%d; yz length is %ld\n", p->M, p->N, arrayMax(p->CompressedAllele));
//
//    pbwtCursorDestroy(u);
//    pbwtCore=p;
//
//}


PBWTWrapper::PBWTWrapper(int nhaps, int nsnps):a(nsnps,std::vector<int>(nhaps,0)), alpha(a),
                                               alphaMap(nsnps,std::unordered_map<int,int>()),
                                               aMap(nsnps,std::unordered_map<int,int>()),
                                               d(a), /*sortedY(nsnps,std::vector<uchar>(nhaps,0)),*/
                                               c(nsnps,0),celta(c),u(a),ultra(a),
                                               haplotypeCluster(a),
                                               clusterAllele(nsnps,std::vector<uchar>())
{
    nSamples=nhaps/2;
    nMarkers=nsnps;
    N=nsnps;
    M=nhaps;//last two haps are slots for current individual need to be phased
    //cerr<<"Inside PBWTWrapper M:"<<M<<endl;
    pbwtCore = pbwtCreate(nhaps, nsnps);
    //pbwtCore->CompressedAllele = arrayCreate(4096 * 32, uchar);

    forwardCursor = pbwtCursorCreate(pbwtCore, TRUE, TRUE);

    reverseCursor = pbwtCursorCreate(pbwtCore, FALSE, TRUE);

    haplotype= nullptr;
    //cerr<<"Inside PBWTWrapper M:"<<M<<endl;
}

#define PREFIX_LENGTH 1200

int PBWTWrapper::CursorForwards() {//so far only implemented for test purpose


    //PrintVector(forwardCursor->a,M,"end arrary aFend check 0");

    for (int k = 0; k != pbwtCore->N; ++k) {
        //fprintf(stderr,"at site %d\n",k);
        CursorForwardsTo(k, PREFIX_LENGTH);
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

int PBWTWrapper::CursorForwardsTo(int k, int T) {
/*T is the length that how far you look back
 *This function must be called along the sites, no skip permitted;
 *Mask the site you want to skip at the begining if you have to.
 */
    int i, i0 = 0, ia;
    int group = 0;

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

    /*reprot haolotype cluster based on prefix, so current site not included*/
    //cluster of the previous site k-1
    for (i = 0; i < forwardCursor->M; ++i) {
        /*assign states of last column based on previous d and sortedY*/
        if (forwardCursor->d[i] > (k - tmpT)&& k !=0) {//new cluster if current sequence and last sequence have common sequence less than T
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
                if(DEBUG && k==64) fprintf(stderr,"i0:%d\tallele:%d\n",i0,haplotype[forwardCursor->a[i0]][k-1]);
                // na = 0;
                // nb = 0;
                i0 = i;
                group++;
            }
        }

    }
    //finish the last segment if i0 didn't reach the end
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
            if(DEBUG && k==64) fprintf(stderr,"i0:%d\tallele:%d\n",i0,haplotype[forwardCursor->a[i0]][k-1]);
            clusterAllele[k-1].push_back(haplotype[forwardCursor->a[i0]][k-1]);
            clusterMembership.push_back(tmpMem);
            //fprintf(stderr,"site:%d\tnumStates:%d\n",k-1,clusterAllele[k-1].size());

        }
    }
    //merge cluster based on KS test
    int test = 0;
    if(DEBUG&&k!=0){
        fprintf(stderr, "at site:%d\n", k-1);
        PrintVector(haplotypeCluster[k-1],"haplotype state before merge state");
        PrintVector(clusterAllele[k-1],"state allele before merge allele");
    }
    if(k!=0&&clusterAllele[k-1].size()==0) {fprintf(stderr,"0 states, abort!");abort();}
    //if(k!=0&&clusterAllele[k-1].size()!=1) test= MergeAtSite(k - 1);//TODO:implement this function
    if(k!=0&&clusterAllele[k-1].size()!=1) test = MergeAtSiteExperiment(k-1);
    if(k!=0&&DEBUG&&test) {
        fprintf(stderr, "at site:%d\n", k-1);
        PrintVector(haplotypeCluster[k-1], "haplotype state after merge state");
        PrintVector(clusterAllele[k-1],"state allele after merge allele");
        fprintf(stderr,"\n");
    }

    //now use haplotype alleles on current site k, to update array a and array d
    //fprintf(stderr,"k:%d,T:%d\t",k,T);PrintVector(forwardCursor->d,forwardCursor->M,"before tmpD");
    //fprintf(stderr,"k:%d,T:%d\t",k,T);PrintVector(forwardCursor->a,forwardCursor->M,"olda");
    //fprintf(stderr,"k:%d,T:%d\t",k,T);PrintVector(forwardCursor->sortedY,forwardCursor->M,"sortedY");
    int u = 0, v = 0;
    int p = k + 1;
    int q = k + 1;
    int *tmpA,*tmpD;
    tmpA=new int [M];
    tmpD=new int [M];
    for (i = 0; i < forwardCursor->M; ++i) {

        if (forwardCursor->d[i] > p) p = forwardCursor->d[i];
        if (forwardCursor->d[i] > q) q = forwardCursor->d[i];

        if (forwardCursor->sortedY[i] == 0) {
            //forwardCursor->a[u] = forwardCursor->a[i];
            //forwardCursor->d[u] = p;
            tmpA[u] = forwardCursor->a[i];
            tmpD[u] = p;
            ++u;
            p = 0;
            // na++;
            forwardCursor->c++;
            this->u[k][i]=u;

        }
        else {
            forwardCursor->b[v] = forwardCursor->a[i];
            forwardCursor->e[v] = q;
            ++v;
            q = 0;
            //nb++;
            this->u[k][i]=u;
        }
    }
    memcpy(forwardCursor->a , tmpA, u * sizeof(int));
    memcpy(forwardCursor->d , tmpD, u * sizeof(int));
    delete [] tmpA;
    delete [] tmpD;
    memcpy(forwardCursor->a + u, forwardCursor->b, v * sizeof(int));
    memcpy(forwardCursor->d + u, forwardCursor->e, v * sizeof(int));
    c[k] = u;

    forwardCursor->d[forwardCursor->M] = k + 2; /* sentinels */
    a[k].assign(forwardCursor->a,forwardCursor->a+forwardCursor->M);

    for (int j = 0; j <a[k].size() ; ++j) {
        //std::cerr<<"a size:"<<a[k].size()<<" and "<<a[k][j]<<std::endl;
        aMap[k].insert(std::make_pair(a[k][j],j));
    }
    d[k].assign(forwardCursor->d,forwardCursor->d+forwardCursor->M);


    return 0;
}

int PBWTWrapper::CursorBackwards() {
    for (int i = pbwtCore->N-1; i!=-1; i--) {

        CursorBackwardsTo(i, 5);
    }
    return 0;
}

int PBWTWrapper::CursorBackwardsTo(int k, int T) {

    int i;
    //copy haplotypes into forwardCursor->y
    CopyHap(k, reverseCursor);


    /*reprot haolotype cluster based on prefix, so current site not included*/
    //fprintf(stderr,"k:%d,T:%d\t",k,T);PrintVector(forwardCursor->d,forwardCursor->M,"before tmpD");
    //fprintf(stderr,"k:%d,T:%d\t",k,T);PrintVector(forwardCursor->a,forwardCursor->M,"olda");
    //fprintf(stderr,"k:%d,T:%d\t",k,T);PrintVector(forwardCursor->sortedY,forwardCursor->M,"sortedY");
    int u = 0, v = 0;
    for (i = 0; i < reverseCursor->M; ++i) {


        if (reverseCursor->sortedY[i] == 0) {
            reverseCursor->a[u] = reverseCursor->a[i];
            ++u;
            reverseCursor->c++;
            ultra[k][i]=u;
        }
        else {
            reverseCursor->b[v] = reverseCursor->a[i];
            ++v;
            ultra[k][i]=u;
        }
    }


    memcpy(reverseCursor->a + u, reverseCursor->b, v * sizeof(int));
    alpha[k].assign(reverseCursor->a,reverseCursor->a+reverseCursor->M);
    celta[k] = u;
    for (int j = 0; j <alpha[k].size() ; ++j) {
        alphaMap[k].insert(std::make_pair(alpha[k][j],j));
    }
    return 0;
}
int PBWTWrapper::ObtainHapFromSinglePhasing(char **haps) {
    haplotype = haps;
    pbwtCore->CompressedAllele = arrayCreate(4096 * 32, uchar);
    forwardCursor = pbwtCursorCreate(pbwtCore, TRUE, TRUE);
    for (int i = 0; i < pbwtCore->N; ++i) {
        for (int j = 0; j < pbwtCore->M; ++j) forwardCursor->sortedY[j] = haplotype[forwardCursor->a[j]][i];
        pbwtCursorWriteForwards(forwardCursor);
        if (nCheckPoint && !((i + 1) % nCheckPoint)) pbwtCheckPoint(forwardCursor, pbwtCore);
    }
    pbwtCursorToAFend(forwardCursor, pbwtCore);
    return 0;
}

int PBWTWrapper::CopyHap(int k, PbwtCursor *Cursor) {//this function has the same effect as forward/backward read
    for (int i = 0; i != Cursor->M; ++i) {
        if (haplotype[Cursor->a[i]][k] >= '0')
            Cursor->sortedY[i] = haplotype[Cursor->a[i]][k] -'0';
        else //fprintf(stderr,"alert!!!! %d,%d,%d,%d\n",haplotype[Cursor->a[i]][k],k,i,Cursor->a[i]);
            Cursor->sortedY[i] = haplotype[Cursor->a[i]][k];
    }

    //PrintVector(Cursor->sortedY,Cursor->M,"fromCopyHap");
    return 0;
}

int PBWTWrapper::UpdateTransVector(int site)//calculate trans probability of site to-1 after site to
{
	if (site > pbwtCore->N||site<0) die((char*)"Site is out of range!");

	if (site == 0)
    {
        return 0;
    }

    int prevSite = site - 1;

    //fprintf(stderr,"site:%d\tprevSite:%d\tsite:%d\n",site,GetNumStates(prevSite),GetNumStates(site));

    //transVector.insert(std::make_pair(prevSite,std::vector<std::vector<float> >(GetNumStates(prevSite),std::vector<float>(GetNumStates(site), 0.))));
    transVector.push_back(std::vector<std::vector<float> >(GetNumStates(prevSite),std::vector<float>(GetNumStates(site), 0.)));
    if(prevSite!=transVector.size()-1) {
        fprintf(stderr,"transvector index error\n");
        exit(1);
    }
	std::vector<float> marginal(GetNumStates(prevSite), 0.000001);
    //int maxi=0;
    //int maxj=0;
	for (int i = 0; i != haplotypeCluster[site].size(); ++i)
	{
        //fprintf(stderr,"prevSite:%d,\ttransVector[from].size():%d\t[to].size():%d\tprevstates:%d\tstates:%d\n",prevSite,transVector[prevSite].size(),transVector[prevSite][haplotypeCluster[prevSite][i]].size(),haplotypeCluster[prevSite][i],haplotypeCluster[site][i]);
		transVector[prevSite][haplotypeCluster[prevSite][i]][haplotypeCluster[site][i]]++;
		//fprintf(stderr,"sum size:%d\ta:%d\n",sum.size(),haplotypeCluster[from][i]);
        marginal[haplotypeCluster[prevSite][i]]++;
        //if(haplotypeCluster[prevSite][i]>maxi) maxi=haplotypeCluster[prevSite][i];
        //if(haplotypeCluster[site][i]>maxj) maxj=haplotypeCluster[site][i];
	}
	for (int i = 0; i != GetNumStates(prevSite); ++i)
	{
		for (int j = 0; j != GetNumStates(site); ++j) {
            transVector[prevSite][i][j] /= marginal[i];
            //fprintf(stderr,"i:%d to j:%d is %f\t",i,j,transVector[prevSite][i][j]);
        }
       //fprintf(stderr,"\n");
	}

    //fprintf(stderr,"finish %d and:prevStates:%d,States:%d\n",prevSite,getNumStates(prevSite),GetNumStates(site));

	return 0;
}

int PBWTWrapper::MergeAtSite(int site) {
    int ret(0);
    int oldNumCluster = GetNumStates(site);
    int numHaps = haplotypeCluster[site].size();
    //std::vector<std::vector<int> > clusterMemberShip(oldNumCluster,std::vector<int>());// state->position in current M array
    std::vector<std::vector<int> > dist(oldNumCluster,std::vector<int>(numHaps,0));//state->rank_occupied
    std::vector<bool> mergeIndicator(oldNumCluster,false);

    std::vector<std::vector<double> > Dmax(oldNumCluster,std::vector<double>(oldNumCluster,0));
    std::vector<std::vector<float> > rankSum(oldNumCluster,std::vector<float>(oldNumCluster,0));
    std::vector<std::vector<std::pair<double,double> > > hapsCounted(oldNumCluster,std::vector<std::pair<double,double> >(oldNumCluster,std::make_pair<int,int>(0,0)));



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
                if (clusterAllele[site][j]!=clusterAllele[site][k] || mergeIndicator[k]) continue;//if alleles are different or merged once

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


                if (0 && DEBUG && i == numHaps - 1) {
                    std::cerr << "\nenter rank distribution section, site " << site << ":" << std::endl;
                    for (auto i = 0; i != dist.size(); ++i) {
                        // PrintDistributionAtSite(i,haplotypeCluster[i]);
                        PrintDistributionAtSite(i, dist[i]);
                    }
                    std::cerr << "exit rank distribution section!\n" << std::endl;
                    //exit(0);
                }
//                if(i==numHaps-1)
//                {
//                    std::cerr<<"Testing state "<<j<<" and "<<k<<std::endl;
//                    //std::cerr<<"hapsCounted["<<j<<"]["<<k<<"]:"<<hapsCounted[j][k].first<<"\t"<<hapsCounted[j][k].second<<std::endl;
//                }

                if (i == numHaps - 1)
                {
                    Dmax[j][k] /= total;

                    pval = 1-psmirnov2x(&Dmax[j][k], tmpMembershipSize, clusterMembership[k].size());

                    //fprintf(stderr,"state %d and %d: %f\n",j,k,pval);
                }
                if(i==numHaps-1 && pval< 0.1)//KStest(Dmax[j][k]/total,tmpMembershipSize,clusterMembership[k].size()))//last haplotypes, deal with merging test
                {
                    ret = 1;
                    //PrintVector(dist[i],"state i");
                    //PrintVector(dist[j],"state j");

                    if(DEBUG)fprintf(stderr,"merge state:%d and state:%d\n",j,k);
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
        MoveSegment(clusterMembership);
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

int PBWTWrapper::MergeCluster(int indexA, int indexB) {
    return 0;
}
//int PBWTWrapper::MergeAtSiteExperiment(int site) {
//    int ret(0);
//    int oldNumCluster = GetNumStates(site);
//    int numHaps = GetNumHaps(site);
//
//    std::vector<std::vector<int> > dist(oldNumCluster,std::vector<int>(numHaps,0));//state->rank_occupied
//    std::vector<bool> mergeIndicator(oldNumCluster,false);
//
//    std::vector<std::vector<double> > Dmax(oldNumCluster,std::vector<double>(oldNumCluster,0));
//    std::vector<std::vector<float> > rankSum(oldNumCluster,std::vector<float>(oldNumCluster,0));
//    std::vector<std::vector<std::pair<double,double> > >
//            hapsCounted(oldNumCluster,std::vector<std::pair<double,double> >(oldNumCluster,std::make_pair<double,double>(0,0)));//cluster-cluster-<double,double>
//
//
//
//    double tmpABS=0;
//    int currentNumCluster=oldNumCluster-1;
//    std::unordered_map<int,int> stateOrder;//mapping oldState to newOrder
//    int tmpOrder(0);
//    std::vector<uchar> tmpAllele;
//    double pVal(0),maxPVal(-1.0),prevPVal(maxPVal);
//
//    std::pair<int,int> pairToBeMerged(0,0);
//
//    for (int ranki = 0; ranki != numHaps; ++ranki)//from rank 0 to rank numHaps-1
//    {
//        //alpha[site][i]:original index of haps at i th place, e.g. david is the ith haplotype
//        //haplotypeCluster[site][alpha[site][i]]: david's state status, saying state is michigan
//        //dist[haplotypeCluster[site][alpha[site][i]]][i]: state michigan's ith slot is occupied, or where is the ith haplotype located, michigan state's ith slot
//        //the structure looks like:
//        //state 1:0000001000100
//        //state 2:0100100000010
//        //state 3:1011010111001
//        int hapID = GetOriginalHapIDFromBack(site, ranki);//original ID
//        int hapState = GetHapState(site, hapID);
//        //std::cerr<<"rank:hapID:hapState\t"<<i<<"\t"<<hapID<<"\t"<<hapState<<std::endl;
//
//        //record rank occupation indicator for each cluster, haplotypeCluster record state status for haplotype alpha[site][i]]
//        //which means hapID occupied ranki slot, and so that hapState occupied ranki slot
//        dist[hapState][ranki] = 1;
//
//        for (int statej = 0; statej < dist.size(); ++statej)
//        {//enumerate through all the states
//            if(mergeIndicator[statej])//ToDO: need adjustment for trangle situation
//            {
//                continue;
//            }
//            else if(currentNumCluster==1) break;
//
//            //if(tmpMembershipSize<5) continue;
//
//            for (int statek = statej + 1; statek < dist.size(); ++statek)//comparing state j and state k
//            {
//                if (clusterAllele[site][statej] != clusterAllele[site][statek] || mergeIndicator[statek])
//                    continue;//if alleles are different or merged once
//
//                double localTotal = clusterMembership[statej].size() + clusterMembership[statek].size();
//                //double sizeRatio = tmpMembershipSize / localTotal;
//                //if (sizeRatio > 0.99 || sizeRatio < 0.01 || clusterMembership[k].size() < 5) continue;
//
//
//                if (dist[statej][ranki] == 1)//if j's ith slot is occupied
//                {
//                    rankSum[statej][statek] += 1;//rankSum between j and k increase by 1
//                    hapsCounted[statej][statek].first = rankSum[statej][statek];//local rank between j and k updated and normalized by localTotal rankSum
//                    //hapsCounted[j][k].first++;
//                }
//                else if (dist[statek][ranki] == 1)//if k's ith slot is occupied
//                {
//                    rankSum[statej][statek] += 1;
//                    hapsCounted[statej][statek].second = rankSum[statej][statek];
//                    //hapsCounted[j][k].second++;
//                }
//                else if (ranki != numHaps - 1)
//                    continue;
//
//                tmpABS = fabs(hapsCounted[statej][statek].first - hapsCounted[statej][statek].second);//tmp Dmax
//
//                if (tmpABS > Dmax[statej][statek]) Dmax[statej][statek] = tmpABS;
//
//
//                Dmax[statej][statek] /= localTotal;
//
//                pVal = 1 - psmirnov2x(&Dmax[statej][statek], clusterMembership[statej].size(), clusterMembership[statek].size());
//
//                if(pVal > maxPVal)
//                {
//                    maxPVal = pVal;
//                    pairToBeMerged=std::make_pair(statej,statek);
//                }
//
//
//            }//inner loop end
//
//
//            if(statej!= dist.size()-1 || ranki != numHaps - 1) continue;//if not the end where both loop finished
//            if(DEBUG) std::cerr<<"state:"<<statej<<" min Pvalue:"<<maxPVal<<" Dmax:"<<Dmax[pairToBeMerged.first][pairToBeMerged.second]<<std::endl;
//            std::cerr << "\nenter rank distribution section, site " << site << ":" << std::endl;
//            for (auto i = 0; i != dist.size(); ++i) {
//                // PrintDistributionAtSite(i,haplotypeCluster[i]);
//                PrintDistributionAtSite(i, dist[i]);
//            }
//            std::cerr << "exit rank distribution section!\n" << std::endl;
//
//            if ( maxPVal!=prevPVal && maxPVal > 0.01)//KStest(Dmax[j][k]/total,tmpMembershipSize,clusterMembership[k].size()))//last haplotypes, deal with merging test
//            {
//                prevPVal = maxPVal;
//                ret = 1;
//                int stateA = pairToBeMerged.first;
//                int stateB = pairToBeMerged.second;
//                //PrintVector(dist[stateA],"stateA");
//                //PrintVector(dist[stateB],"stateB");
//
//                if (DEBUG)fprintf(stderr, "merge state:%d and state:%d, num of states remained %d at site:%d with P value:%f and Dmax:%f\n", stateA, stateB, currentNumCluster,site,maxPVal,Dmax[stateA][stateB]);
//                currentNumCluster--;
//
//                //move dist occupation from stateB to stateA
//                for (int t = 0; t != dist[stateA].size(); ++t) {
//                    dist[stateA][t] += dist[stateB][t];
//                    dist[stateB][t] = -65534;
//                }
//                //Merge Action, change mergeIndicator
//                mergeIndicator[stateB] = true;// j th cluster has been merged into i th cluster
//
//
//                for (int stateM = 0; stateM < dist.size(); ++stateM)//loop through all remained states with the help of mergeIndicator
//                {
//                    if(mergeIndicator[stateM]) continue;
//                    tmpAllele.push_back(clusterAllele[site][stateM]);
//                    stateOrder[stateM] = tmpOrder;
//                    tmpOrder++;
//                    for (int stateN = stateM+1; stateN < dist.size(); ++stateN) {
//                        Dmax[stateM][stateN]=0;
//                        hapsCounted[stateM][stateN]={0,0};
//                        rankSum[stateM][stateN]=0;
//                    }
//                }
//
//                int HapID=0;
//                for (int t = 0; t != clusterMembership[stateB].size(); ++t) {
//                    HapID = clusterMembership[stateB][t];
//                    haplotypeCluster[site][HapID] = stateA;
//                    clusterMembership[stateA].push_back(HapID);
//                }
//                clusterMembership[stateB].clear();
//
//                //adjust d array and a array
//                MoveSegment(clusterMembership);
//                //fprintf(stderr,"site:%d merged...\n",site);
//
//                //PrintVector(clusterAllele[site],"allele cluster states after");
//                //PrintVector(haplotypeCluster[site],"haplotype cluster states after");
//                for (int k = 0; k < haplotypeCluster[site].size(); ++k) {
//                    haplotypeCluster[site][k] = stateOrder[haplotypeCluster[site][k]];
//                }
//                clusterAllele[site] = tmpAllele;//update merged cluster allele
//                //PrintVector(clusterAllele[site],"allele cluster states final");
//                //PrintVector(haplotypeCluster[site],"haplotype cluster states final");
//                tmpAllele.clear();
//                stateOrder.clear();
//                tmpOrder = 0;
//                maxPVal=-1;
//                statej=0;//restart outer loop
//            }
//
//        }
//
//    }
//
//    return ret;
//}

int PBWTWrapper::MergeAtSiteExperiment(int site) {
    int ret(0);
    int oldNumCluster = GetNumStates(site);
    std::cerr<<"Enter Site:"<<site<<" has "<< oldNumCluster<<" state"<<std::endl;
    int numHaps = GetNumHaps(site);
    //if(DEBUG&&site<20) return 0;
    //std::vector<std::vector<int> > dist(oldNumCluster,std::vector<int>(numHaps,0));//state->rank_occupied
    std::vector<std::vector<int> > dist(oldNumCluster,std::vector<int>(0,0));//state->rank_occupied
    std::vector<bool> mergeIndicator(oldNumCluster,false);

//    std::vector<std::vector<double> > Dmax(oldNumCluster,std::vector<double>(oldNumCluster,0));
//    std::vector<std::vector<float> > rankSum(oldNumCluster,std::vector<float>(oldNumCluster,0));
//    std::vector<std::vector<std::pair<double,double> > >
//            hapsCounted(oldNumCluster,std::vector<std::pair<double,double> >(oldNumCluster,std::make_pair<double,double>(0,0)));//cluster-cluster-<double,double>



    double tmpABS=0;
    int currentNumCluster=oldNumCluster+1;
    std::unordered_map<int,int> stateOrder;//mapping oldState to newOrder
    int tmpOrder(0);
    std::vector<uchar> tmpAllele;
    double pVal(0),maxPVal(-1.0),prevPVal(maxPVal);

    std::pair<int,int> pairToBeMerged(0,0);

    for (int ranki = 0; ranki != numHaps; ++ranki)//from rank 0 to rank numHaps-1
    {
        //alpha[site][i]:original index of haps at i th place, e.g. david is the ith haplotype
        //haplotypeCluster[site][alpha[site][i]]: david's state status, saying state is michigan


        int hapID = GetOriginalHapIDFromBack(site, ranki);//backward ID, because we use backward rank to test
        //int hapIDFwd = GetOriginalHapIDFromFwd(site-1, ranki);//backward and forward hapID are the same
        int hapState = GetHapState(site, hapID);
        //std::cerr<<"rank:hapID:hapIDFwd:hapState\t"<<ranki<<"\t"<<hapID<<"\t"<<hapIDFwd<<"\t"<<hapState<<std::endl;

        //record rank occupation indicator for each cluster, haplotypeCluster record state status for haplotype alpha[site][i]]
        //which means hapID occupied ranki slot, and so that hapState occupied ranki slot
        dist[hapState].push_back(ranki);
    }

    {
        int mergeSingleton=false;
        FOR_SINGLETON:
            while(currentNumCluster!=oldNumCluster) {
                oldNumCluster = currentNumCluster;
                for (int statej = 0; statej < dist.size(); ++statej) {//enumerate through all the states
                    if (mergeIndicator[statej] &&
                        clusterMembership[statej].size() == 0)//ToDO: need adjustment for trangle situation
                    {
                        continue;
                    }
                    else if (currentNumCluster == 1) break;

                    //if(tmpMembershipSize<5) continue;

                    for (int statek = statej + 1; statek < dist.size(); ++statek)//comparing state j and state k
                    {
                        if (clusterAllele[site][statej] != clusterAllele[site][statek] || (mergeIndicator[statek] &&
                                                                                               clusterMembership[statek].size() ==
                                                                                               0))
                            continue;//if alleles are different or merged once
                        //std::cerr<<"Now checkout state:"<<statej<<"\tstate:"<<statek<<std::endl;

                        if (!mergeSingleton) {
                            if (dist[statej].size() > 1 && dist[statek].size() > 1) {//both states are not singleton

                                pVal = ks_test(dist[statej], dist[statek]);
                                if (pVal > maxPVal) {
                                    maxPVal = pVal;
                                    pairToBeMerged = std::make_pair(statej, statek);
                                }
                            }
                            else continue;
                        }
                        else {
                            if (dist[statej].size() == 1 || dist[statek].size() == 1) {

                                pVal = ks_test(dist[statej], dist[statek]);

                                if (pVal > maxPVal) {
                                    maxPVal = pVal;
                                    pairToBeMerged = std::make_pair(statej, statek);
                                }
                            }
                            else continue;
                        }
                        if (DEBUG)
                            std::cerr << "out of " << GetNumStates(site) << " states," << statej << "\t" << statek <<
                            " p value:" << pVal << "\twhile max P value:" << maxPVal << "\twith signal:" <<
                            mergeSingleton << std::endl;
                    }//inner loop end
                }

                if (DEBUG) {
                    std::cerr << "\nenter rank distribution section, site " << site << ":" << std::endl;
                    for (auto i = 0; i != dist.size(); ++i) {
                        // PrintDistributionAtSite(i,haplotypeCluster[i]);
                        if (dist[i].size() == 0) continue;
                        PrintDistributionAtSite(i, dist[i]);
                    }
                    std::cerr << "exit rank distribution section!\n" << std::endl;
                }
                if (DEBUG) {
                    std::cerr << "\nenter membership section, site " << site << ":" << std::endl;
                    for (auto i = 0; i != clusterMembership.size(); ++i) {
                        if (clusterMembership[i].size() == 0) continue;
                        std::cerr<<"state "<<i<<" :\t";
                        for (int j = 0; j <clusterMembership[i].size() ; ++j) {
                            std::cerr << GetOriginalHapIDFromFwd(site, clusterMembership[i][j]) << "\t";
                            //std::cerr<<clusterMembership[i][j]<<"\t";
                        }
                        std::cerr<<std::endl;
                    }
                    std::cerr << "exit membership section!\n" << std::endl;
                }

                if ( maxPVal > 0.1)//KStest(Dmax[j][k]/total,tmpMembershipSize,clusterMembership[k].size()))//last haplotypes, deal with merging test
                {
                    prevPVal = maxPVal;
                    ret = 1;
                    int stateA = pairToBeMerged.first;
                    int stateB = pairToBeMerged.second;
                    //PrintVector(dist[stateA],"stateA");
                    //PrintVector(dist[stateB],"stateB");

                    if (DEBUG)fprintf(stderr, "merge state:%d and state:%d, num of states remained %d at site:%d with P value:%f\n", stateA, stateB, currentNumCluster-1,site,maxPVal);
                    currentNumCluster--;

                    //move dist occupation from stateB to stateA
                    for (int t = 0; t != dist[stateB].size(); ++t) {
                        dist[stateA].push_back(dist[stateB][t]);
                    }
                    dist[stateB].clear();
                    //Merge Action, change mergeIndicator
                    mergeIndicator[stateB] = true;// j th cluster has been merged into i th cluster


                    for (int stateM = 0; stateM < dist.size(); ++stateM)//loop through all remained states with the help of mergeIndicator
                    {
                        if(mergeIndicator[stateM]) continue;
                        tmpAllele.push_back(clusterAllele[site][stateM]);
                        stateOrder[stateM] = tmpOrder;
                        tmpOrder++;
                    }

                    int HapID=0;
                    for (int t = 0; t != clusterMembership[stateB].size(); ++t) {
                        HapID = clusterMembership[stateB][t];
                        haplotypeCluster[site][HapID] = stateA;
                        clusterMembership[stateA].push_back(HapID);
                    }
                    clusterMembership[stateB].clear();

                    //adjust d array and a array
                    MoveSegment(clusterMembership);
                    //fprintf(stderr,"site:%d merged...\n",site);

                    //PrintVector(clusterAllele[site],"allele cluster states after");
                    //PrintVector(haplotypeCluster[site],"haplotype cluster states after");
                    for (int k = 0; k < haplotypeCluster[site].size(); ++k) {
                        haplotypeCluster[site][k] = stateOrder[haplotypeCluster[site][k]];
                    }
                    clusterAllele[site] = tmpAllele;//update merged cluster allele
                    //PrintVector(clusterAllele[site],"allele cluster states final");
                    //PrintVector(haplotypeCluster[site],"haplotype cluster states final");
                    tmpAllele.clear();
                    stateOrder.clear();
                    tmpOrder = 0;
                    maxPVal=-1;
                    //statej=-1;//restart outer loop, set as -1 because for loop automatically ++ at the last round of previous loop
                    if(DEBUG)std::cerr<<"Keep shrinking at site:"<<site<<std::endl;
                }
            }
        if(DEBUG)std::cerr<<"End of shrinking at site:"<<site<<std::endl;
        if(!mergeSingleton)//if singleton not processed
        {
            mergeSingleton = true;
            currentNumCluster+=1;
            tmpAllele.clear();
            stateOrder.clear();
            tmpOrder = 0;
            maxPVal=-1;
            goto FOR_SINGLETON;
        }

    }
    std::cerr<<"Exit Site:"<<site<<" has "<< GetNumStates(site)<<" state"<<std::endl;
    return ret;
}


//bool PBWTWrapper::KStest(const double& Dmax, const int & sizeA, const int & sizeB) {
//
//    double thresh=1.36*std::sqrt(double(sizeA+sizeB)/(sizeA*sizeB));
//    if(DEBUG)fprintf(stderr,"Dmax:%lf and threshold:%f\n",Dmax,thresh);
//    if(isinf(thresh)||isnan(thresh))
//        return false;
//    else if(Dmax > thresh)//1.36 is 0.05 significance parameter
//        return false;//reject null hypo, they are different
//    else
//        return true;//accept null hypo, they are the same
//}

int PBWTWrapper::SetHaps(char **haps) {
    haplotype=haps;
    return 0;
}

int PBWTWrapper::PrintDistributionAtSite(int state,std::vector<int> &dist) {
    //std::ofstream fout("/Users/fanzhang/Downloads/PlutoTest/rank.txt",std::ofstream::app);
    std::cerr<<"state:"<<state<<":\t";
   // fout<<"state:"<<state<<":\t";
    for(auto k:dist) {
        std::cerr << k << "\t";
    //    fout << k << "\t";

    }
    std::cerr<<std::endl;
    //fout<<std::endl;

    return 0;
}

int PBWTWrapper::PrintSummary() {
//    mean nodes/level =  70.38  max nodes/level = 111  nodes = 10135
//    mean edges/level =  94.03  max edges/level = 177  edges = 13541
//    mean edges/node  =   1.34  mean count/node =  53.45
    int totalNodes(0),maxNodes(0);
    int totalEdges(0),maxEdges(0);
    float meanEdges(0.0),meanNodes(0.0);
    for (int i = 0; i <clusterAllele.size(); ++i) {
        totalNodes+=clusterAllele[i].size();
        if(clusterAllele[i].size()>maxNodes) maxNodes=clusterAllele[i].size();
    }
    meanNodes=(float) totalNodes/clusterAllele.size();

    for (int j = 0; j < transVector.size(); ++j) {
        int tmpEdges(0);
        for (int i = 0; i < transVector[j].size(); ++i) {
            for (int k = 0; k <transVector[j][i].size() ; ++k) {
                if(transVector[j][i][k]!=0) tmpEdges++;
            }
        }

        if(tmpEdges>maxEdges) maxEdges=tmpEdges;
        totalEdges+=tmpEdges;
    }
    meanEdges=(float) totalEdges/(clusterAllele.size()-1);

    printf("mean nodes/level = %f\tmax nodes/level = %d\tnodes = %d\n",meanNodes,maxNodes,totalNodes);
    printf("mean edges/level = %f\tmax edges/level = %d\tedges = %d\n",meanEdges,maxEdges,totalEdges);
    printf("mean edges/node = %f\tmean count/node = %f\n",(float)totalEdges/totalNodes,(float)(pbwtCore->M)*(pbwtCore->N)/totalNodes);

    return 0;
}

int PBWTWrapper:: MoveSegment(std::vector<std::vector<int> >&Membership) {//fromEnd don't include
//
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
    std::vector<int> tmpD,tmpA;
    for (int i = 0; i < Membership.size() ; ++i) {
        if(Membership[i].size() > 0) {
            tmpD.push_back(forwardCursor->d[Membership[i][0]]);//keep original d
            tmpA.push_back(forwardCursor->a[Membership[i][0]]);
        }
        else
            continue;
        //int d=forwardCursor->d[MemberShip[i][0]];
        for (int j = 1; j < Membership[i].size() ; ++j) {
            tmpD.push_back(0);
            tmpA.push_back(forwardCursor->a[Membership[i][j]]);//MemberShip
        }
    }
    std::copy(tmpD.begin(),tmpD.end(),forwardCursor->d);
    std::copy(tmpA.begin(),tmpA.end(),forwardCursor->a);

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
    return 0;
}

int PBWTWrapper::UpdateRankWithinState(std::vector<std::vector<int> > &dist,int stateA,int stateB) {

    for(int i=0;i!= dist.size();++i)
    {
        for(int j(0),v(0);j!=dist[i].size();++j)
        {
            if(dist[i][j]==1) v++;
            dist[i][j]=v;
        }
    }
    return 0;
}

int PBWTWrapper::RemoveIndividualFromPBWT(int individualIndex) {
    char * haps[2];
    haps[0]=haplotype[2*individualIndex];
    haps[1]=haplotype[2*individualIndex+1];

    char* tmpHapA,*tmpHapB;
    tmpHapA = this->haplotype[individualIndex*2];
    tmpHapB = this->haplotype[individualIndex*2+1];
    this->haplotype[individualIndex*2] = this->haplotype[M-2];
    this->haplotype[individualIndex*2+1] = this->haplotype[M-1];
    this->haplotype[M-2] = tmpHapA;
    this->haplotype[M-1] = tmpHapB;

    for (int i = 0; i <N; ++i) {



    }



    M-=2;
    return 0;
}

int PBWTWrapper::InsertIndividualBackToPBWT(int individualIndex, char **haps) {
    return 0;
}


