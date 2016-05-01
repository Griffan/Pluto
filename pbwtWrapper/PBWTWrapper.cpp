//
// Created by Fan Zhang on 7/20/15.
//


#include <cmath>
#include "PBWTWrapper.h"
#include "iostream"
#include "pbwt/pbwt.h"
#include <fstream>
#include "math.h"
#include <functional>

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
int lmax;
int lmin;
int rmax;
int rmin;
bool comparator(const max_pair_t& lhs, const max_pair_t& rhs) {

//    if( rhs.sizeA <= 5 && rhs.sizeB <=5)//pair with both clusters have size 1 will has lowest priority
//    {
//        return false;//not allow two near singletons to merge
//    }
//    if( lhs.sizeA <=5 && lhs.sizeB <=5)
//    {
//        return true;
//    }
//    else
//    {
//        if(lhs.sizeA>lhs.sizeB)
//        {
//            lmax=lhs.sizeA;
//            lmin=lhs.sizeB;
//        }
//        else
//        {
//            lmax=lhs.sizeB;
//            lmin=lhs.sizeA;
//        }
//        if(rhs.sizeA>rhs.sizeB)
//        {
//            rmax=rhs.sizeA;
//            rmin=rhs.sizeB;
//        }
//        else
//        {
//            rmax=rhs.sizeB;
//            rmin=rhs.sizeA;
//        }
//
//        if (lmax < rmax )//pair with both clusters larger has higher priority
//        {
//            return true;
//        }
//        else {
//            if (lmin < rmin) {
//
//                return lhs.pval < rhs.pval;
//            }
//            else// each has a larger cluster
//            {
//                    return false;
//            }
//        }


//    }

    return lhs.pval<rhs.pval;
}

PBWTWrapper::PBWTWrapper(int nhaps, int nsnps):a(nsnps,std::vector<int>(nhaps,0)), alpha(a),alphaMap(a),aMap(a),
                                               /*alphaMap(nsnps,std::unordered_map<int,int>()),
                                               aMap(nsnps,std::unordered_map<int,int>()),*/
                                               d(a),delta(a), /*sortedY(nsnps,std::vector<uchar>(nhaps,0)),*/
                                               c(nsnps,0),celta(c),u(a),ultra(a),
                                               haplotypeCluster(a),
                                               clusterAllele(nsnps,std::vector<uchar>()),
                                               mergePairList(std::function<bool(const max_pair_t& , const max_pair_t&)>(comparator))
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

//    exact_ks_test_p_val=0.05;
    //CalculateDvalueMatrix();
    //PrintMatrix(DvalueMatrix,"Dvalue:");
    //cerr<<"Inside PBWTWrapper M:"<<M<<endl;
}

#define PREFIX_LENGTH 1200

int PBWTWrapper::CursorForwards()
{//so far only implemented for test purpose


    //PrintVector(forwardCursor->a,M,"end arrary aFend check 0");

    for (int k = 0; k != pbwtCore->N; ++k) {
        //fprintf(stderr,"at site %d\n",k);
        CursorForwardsTo(k, PREFIX_LENGTH);
    }
    //copy end of a to PBWT
    //PrintVector(forwardCursor->a,M,"end arrary aFend check 1");

    pbwtCursorToAFend(forwardCursor, pbwtCore);

//    for (int i=0;i != pbwtCore->N; i++) {
//
//        UpdateTransVector(i);
//    }
    PrintSummary();
    //update crossover rate?
    return 0;
}

int PBWTWrapper::CursorForwardsTo(int k, int T)
{
/*T is the length that how far you look back
 *This function must be called along the sites, no skip permitted;
 *Mask the site you want to skip at the begining if you have to.
 */
    int rank, i0 = 0, ia;
    int group = 0;

    int hapID(0);

    clusterMembership.clear();
    hasSiblings.clear();
    dist.clear();


    int tmpT= k > T ? T : k;

    /*reprot haolotype cluster based on prefix, so current site not included*/
    //cluster of the previous site k-1
    for (rank = 0; rank < forwardCursor->M; ++rank) {
        /*assign states of last column based on previous d and sortedY*/
        if (forwardCursor->d[rank] > (k - tmpT) && k != 0) {//new cluster if current sequence and last sequence have common sequence less than T
            //if (na && nb)        /* then there is something to report */
            if(rank != 0) {
                //fprintf(stderr,"d:%d\tk-tmpT:%d\n",forwardCursor->d[i],k-tmpT);
                dist.push_back(std::vector<int>(rank-i0, 0));//state->rank_occupied
                std::vector<int> tmpMem;
                for (ia = i0; ia < rank; ++ia) {
                    hapID= GetHapIDFromFwd(k-1, ia);//original ID, by treating rank as backward ranking
                    haplotypeCluster[k - 1][hapID] = group;
                    tmpMem.push_back(ia);
                    dist[group][ia-i0] = GetRankFromBack(k-1,hapID);
                }
                std::sort(dist[group].begin(),dist[group].end());
                clusterAllele[k - 1].push_back(haplotype[GetHapIDFromFwd(k-1,i0)][k-1]);
//                if(k==6)fprintf(stderr,"site 5, GetHapID %d or %d and allele %d\n",GetHapIDFromFwd(k-2,i0),GetHapIDFromFwd(k-1,i0),haplotype[GetHapIDFromFwd(k-1,i0)][k-1]);
                clusterMembership.push_back(tmpMem);
                hasSiblings.push_back(false);
                // na = 0;
                // nb = 0;
                i0 = rank;
                group++;
            }
        }

    }
    //finish the last segment if i0 didn't reach the end
    if( i0 < forwardCursor->M)
    {
        if (k != 0)
        {
            dist.push_back(std::vector<int>(forwardCursor->M-i0, 0));//state->rank_occupied
            std::vector<int> tmpMem;
            for (ia = i0; ia < forwardCursor->M; ++ia)
            {
                hapID= GetHapIDFromFwd(k-1, ia);//original ID, by treating rank as backward ranking
                haplotypeCluster[k-1][hapID] = group;
                tmpMem.push_back(ia);
                dist[group][ia-i0] = GetRankFromBack(k-1,hapID);
            }
            std::sort(dist[group].begin(),dist[group].end());
            clusterAllele[k-1].push_back(haplotype[GetHapIDFromFwd(k-1,i0)][k-1]);
//            if(k==6)fprintf(stderr,"site 5, GetHapID %d or %d and allele %d\n",GetHapIDFromFwd(k-2,i0),GetHapIDFromFwd(k-1,i0),haplotype[GetHapIDFromFwd(k-1,i0)][k-1]);
            clusterMembership.push_back(tmpMem);
            hasSiblings.push_back(false);
            //fprintf(stderr,"site:%d\tnumStates:%d\n",k-1,clusterAllele[k-1].size());

        }
    }
    //merge cluster based on KS test
    int test = 0;
#ifdef DEBUG
    if(k!=0){
        fprintf(stderr, "at site:%d\n", k-1);
        PrintVector(haplotypeCluster[k-1],"haplotype state before merge state");
        PrintVector(clusterAllele[k-1],"state allele before merge allele");
    }
#endif
    if(k!=0) {
        LabelNoSiblingCluster(k - 1);
        if (clusterAllele[k - 1].size() != 1) test = MergeAtSite(k - 1);//TODO:implement this function
        //if(k!=0&&clusterAllele[k-1].size()!=1) test = MergeAtSiteExperiment(k-1);
#ifdef DEBUG
        if(test) {
            fprintf(stderr, "at site:%d\n", k-1);
            PrintVector(haplotypeCluster[k-1], "haplotype state after merge state");
            PrintVector(clusterAllele[k-1],"state allele after merge allele");
            fprintf(stderr,"\n");
        }
#endif
        UpdateTransVector(k - 1);
    }


    //copy haplotypes into forwardCursor->y
    CopyHap(k, forwardCursor);
    //now use haplotype alleles on current site k, to update array a and array d
    //fprintf(stderr,"k:%d,T:%d\t",k,T);PrintVector(forwardCursor->d,forwardCursor->M,"before tmpD");
    //fprintf(stderr,"k:%d,T:%d\t",k,T);PrintVector(forwardCursor->a,forwardCursor->M,"olda");
    //fprintf(stderr,"k:%d,T:%d\t",k,T);PrintVector(forwardCursor->sortedY,forwardCursor->M,"sortedY");
    int u = 0, v = 0;
    int p = k + 1;
    int q = k + 1;
//    int *tmpA,*tmpD;
//    tmpA=new int [M];
//    tmpD=new int [M];
    for (rank = 0; rank < forwardCursor->M; ++rank) {

        if (forwardCursor->d[rank] > p) p = forwardCursor->d[rank];
        if (forwardCursor->d[rank] > q) q = forwardCursor->d[rank];

        if (forwardCursor->sortedY[rank] == 0) {
            forwardCursor->a[u] = forwardCursor->a[rank];
            forwardCursor->d[u] = p;
//            tmpA[u] = forwardCursor->a[rank];
//            tmpD[u] = p;
            ++u;
            p = 0;
            // na++;
            forwardCursor->c++;
            this->u[k][rank]=u;

        }
        else {
            forwardCursor->b[v] = forwardCursor->a[rank];
            forwardCursor->e[v] = q;
            ++v;
            q = 0;
            //nb++;
            this->u[k][rank]=u;
        }
    }
//    memcpy(forwardCursor->a , tmpA, u * sizeof(int));
//    memcpy(forwardCursor->d , tmpD, u * sizeof(int));
//    delete [] tmpA;
//    delete [] tmpD;
    memcpy(forwardCursor->a + u, forwardCursor->b, v * sizeof(int));
    memcpy(forwardCursor->d + u, forwardCursor->e, v * sizeof(int));
    c[k] = u;

    forwardCursor->d[forwardCursor->M] = k + 2; /* sentinels */
    a[k].assign(forwardCursor->a,forwardCursor->a+forwardCursor->M);

    for (int j = 0; j <a[k].size() ; ++j) {
        //std::cerr<<"a size:"<<a[k].size()<<" and "<<a[k][j]<<std::endl;
        //aMap[k].insert(std::make_pair(a[k][j],j));
        aMap[k][a[k][j]]=j;
    }
    d[k].assign(forwardCursor->d,forwardCursor->d+forwardCursor->M);

    if (k==pbwtCore->N-1)//deal with last columns
    {
        i0 = 0;
        group = 0;

        clusterMembership.clear();
        hasSiblings.clear();
        dist.clear();
        for (rank = 0; rank < forwardCursor->M; ++rank) {
            /*assign states of last column based on previous d and sortedY*/
            if (forwardCursor->d[rank] > (k - tmpT)) {//new cluster if current sequence and last sequence have common sequence less than T
                //if (na && nb)        /* then there is something to report */
                if(rank != 0) {
                    //fprintf(stderr,"d:%d\tk-tmpT:%d\n",forwardCursor->d[i],k-tmpT);
                    dist.push_back(std::vector<int>(rank-i0, 0));//state->rank_occupied
                    std::vector<int> tmpMem;
                    for (ia = i0; ia < rank; ++ia) {
                        hapID= GetHapIDFromFwd(k, ia);//original ID, by treating rank as backward ranking
                        haplotypeCluster[k][hapID] = group;
                        tmpMem.push_back(ia);
                        dist[group][ia-i0] = GetRankFromBack(k,hapID);
                    }
                    std::sort(dist[group].begin(),dist[group].end());
                    clusterAllele[k].push_back(haplotype[GetHapIDFromFwd(k,i0)][k]);
                    clusterMembership.push_back(tmpMem);
                    hasSiblings.push_back(false);
                    // na = 0;
                    // nb = 0;
                    i0 = rank;
                    group++;
                }
            }

        }
        //finish the last segment if i0 didn't reach the end
        if( i0 < forwardCursor->M)
        {
                dist.push_back(std::vector<int>(forwardCursor->M-i0, 0));//state->rank_occupied
                std::vector<int> tmpMem;
                for (ia = i0; ia < forwardCursor->M; ++ia)
                {
                    hapID= GetHapIDFromFwd(k, ia);//original ID, by treating rank as backward ranking
                    haplotypeCluster[k][hapID] = group;
                    tmpMem.push_back(ia);
                    dist[group][ia-i0] = GetRankFromBack(k,hapID);
                }
                std::sort(dist[group].begin(),dist[group].end());
                clusterAllele[k].push_back(haplotype[GetHapIDFromFwd(k,i0)][k]);
                clusterMembership.push_back(tmpMem);
                hasSiblings.push_back(false);
                //fprintf(stderr,"site:%d\tnumStates:%d\n",k-1,clusterAllele[k-1].size());
        }

            LabelNoSiblingCluster(k);
            UpdateTransVector(k);//for the last second site

    }

    return 0;
}

int PBWTWrapper::CursorBackwards()
{
    //for (int i = pbwtCore->N-1; i!=-1; i--) {
    for (int i = 0; i!=pbwtCore->N; i++) {

        CursorBackwardsTo(i, 5);
    }
    return 0;
}

int PBWTWrapper::CursorBackwardsTo(int k, int T)
{

    int rank;
    //copy haplotypes into forwardCursor->y
    int site=N-k-1;
    CopyHap(site, reverseCursor);


    /*reprot haolotype cluster based on prefix, so current site not included*/
    //fprintf(stderr,"k:%d,T:%d\t",k,T);PrintVector(forwardCursor->d,forwardCursor->M,"before tmpD");
    //fprintf(stderr,"k:%d,T:%d\t",k,T);PrintVector(forwardCursor->a,forwardCursor->M,"olda");
    //fprintf(stderr,"k:%d,T:%d\t",k,T);PrintVector(forwardCursor->sortedY,forwardCursor->M,"sortedY");
    int u = 0, v = 0;
    int p = k + 1;
    int q = k + 1;
    for (rank = 0; rank < reverseCursor->M; ++rank) {
        if (reverseCursor->d[rank] > p) p = reverseCursor->d[rank];
        if (reverseCursor->d[rank] > q) q = reverseCursor->d[rank];

        if (reverseCursor->sortedY[rank] == 0) {
            reverseCursor->a[u] = reverseCursor->a[rank];
            reverseCursor->d[u] = p;
            p=0;
            ++u;
            reverseCursor->c++;
            ultra[site][rank]=u;
        }
        else {
            reverseCursor->b[v] = reverseCursor->a[rank];
            reverseCursor->e[v] = q;
            q=0;
            ++v;
            ultra[site][rank]=u;
        }
    }


    memcpy(reverseCursor->a + u, reverseCursor->b, v * sizeof(int));
    memcpy(reverseCursor->d + u, reverseCursor->e, v * sizeof(int));
    reverseCursor->d[reverseCursor->M] = k + 2; /* sentinels */
    alpha[site].assign(reverseCursor->a,reverseCursor->a+reverseCursor->M);
    celta[site] = u;
    for (int j = 0; j <alpha[site].size() ; ++j) {
        //alphaMap[k].insert(std::make_pair(alpha[k][j],j));//hapID,rank
        alphaMap[site][alpha[site][j]]=j;
    }
    delta[site].assign(reverseCursor->d,reverseCursor->d+reverseCursor->M);
    return 0;
}

int PBWTWrapper::CopyHap(int k, PbwtCursor *Cursor)
{//this function has the same effect as forward/backward read
    for (int i = 0; i != Cursor->M; ++i) {
        if (haplotype[Cursor->a[i]][k] >= '0')
            Cursor->sortedY[i] = haplotype[Cursor->a[i]][k] -'0';
        else //fprintf(stderr,"alert!!!! %d,%d,%d,%d\n",haplotype[Cursor->a[i]][k],k,i,Cursor->a[i]);
            Cursor->sortedY[i] = haplotype[Cursor->a[i]][k];
    }

    //PrintVector(Cursor->sortedY,Cursor->M,"fromCopyHap");
    return 0;
}

int PBWTWrapper::LabelNoSiblingCluster(int site)
{
    if (site == 0)
    {
        inEdges.push_back(EDGE());
        return 0;
    }
    int prevSite = site - 1;
    inEdges.push_back(EDGE());
    outEdges.push_back(EDGE());

    for (int hapID = 0; hapID != haplotypeCluster[site].size(); ++hapID)//loop through each hapID
    {
        inEdges[site][GetHapState(site,hapID)][GetHapState(prevSite,hapID)]=true;
        outEdges[prevSite][GetHapState(prevSite,hapID)][GetHapState(site,hapID)]=true;
    }
    int cnt=0;
    for (int state = 0; state <hasSiblings.size(); ++state) {
        for (auto kv:inEdges[site][state]) {
            int parentState=kv.first;
            if(outEdges[prevSite][parentState].size()>1) {
                hasSiblings[state]=true;
                cnt++;
                break;
            }
        }
    }
    std::cerr<<"site:"<<site<<" has "<<cnt<<" states has Sibs"<<std::endl;
    return 0;
}

int PBWTWrapper::UpdateTransVector(int site)//calculate trans probability of site to-1 after site to
{
	//if (site > pbwtCore->N||site<0) die((char*)"Site is out of range!");

	if (site == 0)
    {
        //inEdges.push_back(EDGE());
        return 0;
    }

    int prevSite = site - 1;
    transVector.push_back(std::vector<std::vector<float> >(GetNumStates(prevSite),std::vector<float>(GetNumStates(site),0)));
	std::vector<float> marginal(GetNumStates(prevSite), 0.0);
    //inEdges.push_back(EDGE());
    inEdges[site].clear();
    //outEdges[prevSite].clear();
	for (int hapID = 0; hapID != haplotypeCluster[site].size(); ++hapID)//loop through each hapID
	{
//        fprintf(stderr,"prevSite:%d,\ttransVector[from].size():%d\t[to].size():%d\tprevstates:%d\tstates:%d\n",
//                prevSite,transVector[prevSite].size(),transVector[prevSite][GetHapState(prevSite, hapID)].size(),
//                GetHapState(prevSite, hapID),GetHapState(site, hapID));
        transVector[prevSite][GetHapState(prevSite, hapID)][GetHapState(site, hapID)]+=1.;
		//fprintf(stderr,"sum size:%d\ta:%d\n",sum.size(),haplotypeCluster[from][i]);
        marginal[GetHapState(prevSite, hapID)]+=1.;
        //if(haplotypeCluster[prevSite][i]>maxi) maxi=haplotypeCluster[prevSite][i];
        //if(haplotypeCluster[site][i]>maxj) maxj=haplotypeCluster[site][i];

	}
	for (int i = 0; i != GetNumStates(prevSite); ++i)
	{
		for (int j = 0; j != GetNumStates(site); ++j) {
            if(transVector[prevSite][i][j]>0)
            {
                transVector[prevSite][i][j] /= marginal[i];
                inEdges[site][j][i] = true;
                //outEdges[prevSite][i][j] = true;
            }
//            fprintf(stderr,"i:%d to j:%d is %f\t",i,j,transVector[prevSite][i][j]);
        }
//       fprintf(stderr,"\n");
	}
//    fprintf(stderr,"finish %d and:prevStates:%d,States:%d\n",prevSite,GetNumStates(prevSite),GetNumStates(site));

	return 0;
}

void PBWTWrapper::MergeSortedArrayToA(std::vector<int> &a, std::vector<int> &b)
{
    int indexA(0),indexB(0),indexTotal(0);
    std::vector<int> mergedDist(a.size()+b.size());
    while(indexA < a.size() && indexB <b.size())
    {
        if (a[indexA] < b[indexB]) {
            mergedDist[indexTotal++]=a[indexA];
            indexA++;

        }
        else {
            mergedDist[indexTotal++] = b[indexB];
            indexB++;
        }

    }
    while(indexA < a.size()) {
        mergedDist[indexTotal++] = a[indexA];
        indexA++;
    }
    while(indexB < b.size())
    {
        mergedDist[indexTotal++]=b[indexB];
        indexB++;
    }
        a=mergedDist;
}

bool PBWTWrapper::IsRecipricalLengthOK(std::vector<int> &a, std::vector<int> &b)
{
    double size=a.size()<b.size()?a.size():b.size();
    int cnt=0;
       if(a.front()<b.front())
       {
           if(a.back()<b.back()) {
               for (int i = 0; i <b.size() ; ++i) {
                   if(b[i]<a.back()) cnt++;
                   else break;
               }
               return cnt / size > 0.8;
           }
           else
               return true;
       }
        else
       {
            if(a.back()<b.back())
                return true;
           else
            {
                for (int i = 0; i <a.size() ; ++i) {
                    if(a[i]<b.back()) cnt++;
                    else break;
                }
                return cnt / size > 0.8;
            }
       }
}

bool PBWTWrapper::IsEditDistanceOK(const std::vector< std::vector<char*> >& backBone,int stateA, int stateB, int index, int eror_thresh)
{
    int size=N-index>500?500:N-index;
    int end=index+size;
    int start=index;

    for (int i = 0; i <backBone[stateA].size(); ++i) {
        for (int j = 0; j <backBone[stateB].size(); ++j) {
            start=index;
            while (start < end) {
                if (backBone[stateA][i][start] != backBone[stateB][j][start]) break;
                else if (start > eror_thresh) return true;
                start++;
            }

        }
    }
    return false;
}

int PBWTWrapper::CalculateDmax(double & pval, double & Dmax, std::vector<int> & j, std::vector<int>& k) {
    double thresh = 1.22 * sqrt(1. / j.size() + 1. / k.size());
    //assume dist has sorted ranks
    double Ddiff = 0;
    pval = 0;
    Dmax = 0;
    int indexA(0), indexB(0);
    while (indexA < j.size() && indexB < k.size()) {
        if (j[indexA] < k[indexB]) {
            indexA++;
            Ddiff += 1. / j.size();

        }
        else {
            indexB++;
            Ddiff -= 1. / k.size();
        }
        if (Ddiff > Dmax) Dmax = Ddiff;
        if (fabs(Dmax) > thresh) {
            //goto DIST_END;
            return 1;
        }

    }
    while (indexA < j.size()) {
        indexA++;
        Ddiff += 1. / j.size();
        if (Ddiff > Dmax) Dmax = Ddiff;
        if (fabs(Dmax) > thresh) {
            //goto DIST_END;
            return 1;
        }
    }
    while (indexB < k.size()) {
        indexB++;
        Ddiff -= 1. / k.size();
        if (Ddiff > Dmax) Dmax = Ddiff;
        if (fabs(Dmax) > thresh) {
            //goto DIST_END;
            return 1;
        }
    }
    return 0;
};
int PBWTWrapper::MergeAtSite(int site)
{
//    if(dist.size()<400) return 1;
//    if(freq1s[site]<0.1 or freq1s[site]>0.9) return 1;
    int ret(0);
    int currentNumCluster = GetNumStates(site);
    std::cerr<<"Enter Site:"<<site<<" has "<< currentNumCluster<<" state and List size:"<<mergePairList.size()<<std::endl;

//    unsigned long numHaps = haplotypeCluster[site].size();


    std::vector<bool> removeIndicator(currentNumCluster, false);
    std::vector<bool> retainIndicator(currentNumCluster, false);

    std::vector< std::vector<char*> > backBone(currentNumCluster,std::vector<char*>(0,nullptr));

    for (int i = 0; i <currentNumCluster ; ++i) {
        backBone[i].push_back(haplotype[GetHapIDFromFwd(site,clusterMembership[i].front())]);
    }

    tmpABS = 0;
    tmpOrder=0;
    pval=0;
    EXACT=false;
    stateOrder.clear();//mapping oldState to newOrder
    tmpAllele.clear();
    removeMembership.clear();//rankID,state


    int retainState;
    int removeState;

    double Dmax=0.;
    double Ddiff=0.;

    for (int j = 0; j < dist.size(); ++j) {//enumerate through all the states
        if (!hasSiblings[j]/* || clusterMembership[j].size() == 0*/) continue;
        for (int k = 0; k < j; ++k) {

            if (hasSiblings[k]||clusterAllele[site][j] != clusterAllele[site][k] /*|| clusterMembership[k].size() == 0*/) continue;
            if(dist[j].size()<= 2 )
            {
                if(dist[k].size() <= 2) {//both rare
                    if (!IsEditDistanceOK(backBone, j, k, site, 40))
                        continue;
                }
                else//j rare, k not
                {
                    if (!IsEditDistanceOK(backBone, j, k, site, 10))
                        continue;
                }
            }
            else if(dist[k].size() <= 2)//j not, k rare
            {
                if (!IsEditDistanceOK(backBone, j, k, site, 10))
                    continue;
            }

            {
                if(!IsRecipricalLengthOK(dist[j], dist[k]))//||!IsEditDistanceOK(backBone,j,k,site,0.5))
                    continue;
            }
//            if (clusterAllele[site][j] != clusterAllele[site][k] || clusterMembership[k].size() == 0) { std::cerr<<"skip, diff allele"<<std::endl;continue;}
//               if( !IsRecipricalLengthOK(dist[j],dist[k])){ std::cerr<<"skip, not reciprical 80%"<<std::endl; continue;}
//                if(!IsEditDistanceOK(haplotype[clusterMembership[j].front()],haplotype[clusterMembership[k].front()],site,N)){ std::cerr<<"skip, edit distance"<<std::endl;
//                    continue;}

            if (clusterMembership[j].size() * clusterMembership[k].size() < 10000)//exact
            {
                EXACT = true;
                //thresh = GetExactThresh(clusterMembership[j].size(), clusterMembership[k].size());
            }
            else {
                EXACT = false;
//                thresh = 1.22 * sqrt((clusterMembership[j].size() + clusterMembership[k].size()) /
//                                     (clusterMembership[j].size() * clusterMembership[k].size()));
            }

            if(CalculateDmax(pval,Dmax,dist[j],dist[k]))//return early
            {
                goto DIST_END;
            }
            if(EXACT)
                pval = 1 - psmirnov2x(&Dmax, clusterMembership[j].size(), clusterMembership[k].size());
            else {
                Ddiff = sqrt(double(clusterMembership[j].size() + clusterMembership[k].size())) * Dmax;
                pval = 1 - pkstwo_wrapper(1, &Ddiff, 1e-06);
            }

            {
                max_pair_t mergePair={j,k,Dmax,EXACT,pval,dist[j].size(),dist[k].size()};
                mergePairList.push(mergePair);
            }

            DIST_END:
//            fprintf(stderr,"first:(%d,%d) Dmax:%f and Thresh:%f and Pval:%f, with sample size:%d,%d and %d,%d\n",
//                        j,k,Dmax,thresh,pval,dist[j].size(),dist[k].size(),clusterMembership[j].size(),clusterMembership[k].size());
             continue;
        }
    }

    int clusterA(0), clusterB(0);
    double old_p_max =0;
    while(!mergePairList.empty()) {

        max_pair_t iter_pair=mergePairList.top();
        mergePairList.pop();
        old_p_max =iter_pair.pval;
        clusterA = iter_pair.clusterA;
        clusterB = iter_pair.clusterB;
//        fprintf(stderr,"second:(%d,%d) Dmax:%f and Thresh:%f and Pval:%f, with sample size:%d,%d \n", clusterA,clusterB,iter_pair.Dmax,thresh,iter_pair.pval,dist[clusterA].size(),dist[clusterB].size());
//        continue;
        if (removeIndicator[clusterA] || removeIndicator[clusterB]) continue;
        if (retainIndicator[clusterA] || retainIndicator[clusterB]) {
            if (clusterMembership[clusterA].size() * clusterMembership[clusterB].size() < 10000)//exact
            {
                iter_pair.exact = true;
//                thresh = GetExactThresh(clusterMembership[clusterA].size(), clusterMembership[clusterB].size());
            }
            else {
                iter_pair.exact = false;
//                thresh = 1.22 * sqrt((clusterMembership[clusterA].size() + clusterMembership[clusterB].size()) /
//                                     (clusterMembership[clusterA].size() * clusterMembership[clusterB].size()));
            }
            if(CalculateDmax(iter_pair.pval,iter_pair.Dmax,dist[clusterA],dist[clusterB]))//return early
            {
                goto END_WHILE;
            }
            iter_pair.Dmax=Dmax;

            if (iter_pair.exact) {
                iter_pair.pval = 1 - psmirnov2x(&iter_pair.Dmax, clusterMembership[clusterA].size(), clusterMembership[clusterB].size());
            }
            else {
                Ddiff = sqrt(double(clusterMembership[clusterA].size() + clusterMembership[clusterB].size())) * iter_pair.Dmax;
                iter_pair.pval = 1 - pkstwo_wrapper(1, &Ddiff, 1e-06);
            }
            if(iter_pair.pval < old_p_max)
            {
                max_pair_t mergePair={clusterA,clusterB,iter_pair.Dmax,iter_pair.exact,iter_pair.pval,dist[clusterA].size(),dist[clusterB].size()};
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
//        fprintf(stderr, "fail to merge state:%d and state:%d, num of states remained %d at site:%d with P value:%f\t and Dmax:%f\n", clusterA, clusterB, currentNumCluster-1,site,pval,iter_pair.Dmax);

//        fprintf(stderr,"second:(%d,%d) Dmax:%f and Pval:%f, with sample size:%d,%d \n", clusterA,clusterB,iter_pair.Dmax,iter_pair.pval,dist[clusterA].size(),dist[clusterB].size());


        if (iter_pair.pval > 0.01)
        {

            retainState=clusterA;
            removeState=clusterB;
//            PrintVector(dist[retainState],"retainState");
//            PrintVector(dist[removeState],"removeState");
            ret = 1;
//            if(!hasSiblings[clusterA])
//            {
//                retainState=clusterB;
//                removeState=clusterA;
//            }
//            fprintf(stderr,"third:(%d,%d) Dmax:%f and Thresh:%f and Pval:%f, with sample size:%d,%d and %d,%d\n", retainState,removeState,iter_pair.Dmax,thresh,iter_pair.pval,dist[retainState].size(),dist[removeState].size(),clusterMembership[retainState].size(),clusterMembership[removeState].size());

            DoMerge(site, retainState, removeState, dist, removeIndicator, retainIndicator, removeMembership);

            for (int i = 0; i <backBone[removeState].size() ; ++i) {
                backBone[retainState].push_back(backBone[removeState][i]);
            }
            backBone[removeState].clear();

            //finish merge, look for next candidate pair
        }
        else// if(clusterMembership[clusterA].size()<=2 && clusterMembership[clusterA].size()<=2)
        {
            mergePairList= std::priority_queue<max_pair_t,std::vector<max_pair_t>, std::function<bool(const max_pair_t&,const max_pair_t&)> >(comparator);
            break;
        }
        END_WHILE:
        continue;
    }

    if (ret) {
        for (int stateM = 0;
             stateM < dist.size(); ++stateM)//loop through all remained states with the help of mergeIndicator
        {
            if (removeIndicator[stateM]) continue;
            tmpAllele.push_back(clusterAllele[site][stateM]);
            stateOrder[stateM] = tmpOrder;
            tmpOrder++;
        }
//        PrintVector(clusterAllele[site],"allele cluster states after");
        //PrintVector(haplotypeCluster[site],"haplotype cluster states after");
        for (int k = 0; k < haplotypeCluster[site].size(); ++k) {
            haplotypeCluster[site][k] = stateOrder[haplotypeCluster[site][k]];
        }
        clusterAllele[site] = tmpAllele;//update merged cluster allele
//        PrintVector(clusterAllele[site],"allele cluster states final");
        //PrintVector(haplotypeCluster[site],"haplotype cluster states final");
        //adjust d array and a array
        MoveSegment(removeMembership, site);
    }
   std::cerr<<"Exit Site:"<<site<<" has "<< GetNumStates(site)<<" state"<<std::endl;
    return ret;
}

void PBWTWrapper::DoMerge(int site, int retainState, int removeState, std::vector<std::vector<int>> &dist,
                          std::vector<bool, std::allocator<bool>> &removeIndicator,
                          std::vector<bool, std::allocator<bool>> &retainIndicator,
                          std::unordered_map<int, int> &removeMembership)
{//move dist occupation from stateB to stateA
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
    for (int t = 0; t != clusterMembership[removeState].size(); ++t) {
                removeRankID = clusterMembership[removeState][t];
                haplotypeCluster[site][GetHapIDFromFwd(site, removeRankID)] = retainState;
                clusterMembership[retainState].push_back(removeRankID);
                removeMembership[removeRankID] = removeState;
            }
    clusterMembership[removeState].clear();
//    for(auto kv:inEdges[site][removeState])
//    {
//        inEdges[site][retainState][kv.first]=true;
//    }
}

int PBWTWrapper::SetHaps(char **haps,double * freq)
{
    haplotype=haps;
    freq1s=freq;
    return 0;
}

int PBWTWrapper::PrintDistributionAtSite(int state,std::vector<int> &dist)
{
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

int PBWTWrapper::PrintSummary()
{
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

/*int PBWTWrapper:: MoveSegment(const std::vector<std::vector<int> >&Membership, int stateA, int stateB, int site) {//fromEnd don't include

    fprintf(stderr,"before d:\n");
    for (int k = 0; k <M ; ++k) {
        fprintf(stderr,"%d\t",forwardCursor->d[k]);
    }
    fprintf(stderr,"\n");
    fprintf(stderr,"before a:\n");
    for (int k = 0; k <M ; ++k) {
        fprintf(stderr,"%d\t",forwardCursor->a[k]);
    }
    fprintf(stderr,"\n");

    //find the first preceding state M of stateB and the first following state N of stateB
    //find the d value between the first element of N and the last element of M
    int stateN(stateB),stateM(stateB),stateX(stateA);
    int newD1(-1),newD2(-1);

    while(stateN--)
    {
        if(stateN<0) {stateN=-1;break;}//stateB is the last state
        if(Membership[stateN].size()!=0) break;
    }
    std::vector<int> tmpD, tmpA;
    if(stateN>stateA) {//stateA and stateB are not adjacent
        std::cerr<<"Not adjacent"<<std::endl;
        stateN=stateB;
        while(stateN++)
        {
            if(stateN>=Membership.size()) {stateN=-1;break;}//stateB is the last state
            if(Membership[stateN].size()!=0) break;
        }
        while(stateM--)
        {
            if(stateM<0) {break;}//stateB is the first state, which is unlikely
            if(Membership[stateM].size()!=0) break;
        }

        //find the new d value between last element of stateB and the first element of the following state of stateA
        while(stateX++)
        {
            if(stateX>=Membership.size()) {stateX=-1;break;}//stateA is the last state
            if(Membership[stateX].size()!=0) break;
        }
        if(stateM>0&&stateN>0)//find the d value for these two adjacent states
        {
            int p=GetHapIDFromFwd(site,Membership[stateM].back());
            int q=GetHapIDFromFwd(site,Membership[stateN].front());

            int index=site;
            while(index--)
            {
                if(haplotype[p][index]!=haplotype[q][index]) break;
            }
            newD2=index;
        }
        if(stateX>0)
        {
            int p = GetHapIDFromFwd(site, Membership[stateB].back());
            int q = GetHapIDFromFwd(site,Membership[stateX].front());

            int index=site;
            while(index--)
            {
                if(haplotype[p][index]!=haplotype[q][index]) break;
            }
            newD1=index;

        }

        for (int i = 0; i < Membership.size(); ++i) {
            if (i == stateA)//put d in stateB behind of stateA
            {
                for (int j = 0; j < Membership[i].size(); ++j) {
                    tmpD.push_back(forwardCursor->d[Membership[i][j]]);
                    tmpA.push_back(forwardCursor->a[Membership[i][j]]);//MemberShip
                }
                tmpD.push_back(0);//the first element of original stateB
                tmpA.push_back(forwardCursor->a[Membership[stateB][0]]);
                for (int j = 1; j < Membership[stateB].size(); ++j) {
                    tmpD.push_back(forwardCursor->d[Membership[stateB][j]]);
                    tmpA.push_back(forwardCursor->a[Membership[stateB][j]]);
                }
            }
            else if (i == stateX) {
                tmpD.push_back(newD1);//new d value
                tmpA.push_back(forwardCursor->a[Membership[i][0]]);
                for (int j = 1; j < Membership[i].size(); ++j) {
                    tmpD.push_back(forwardCursor->d[Membership[i][j]]);
                    tmpA.push_back(forwardCursor->a[Membership[i][j]]);
                }
            }
            else if (i == stateB) continue;
            else if (i == stateN)//put the newly found d value to the first element of stateN
            {
                tmpD.push_back(newD2);//new d value
                tmpA.push_back(forwardCursor->a[Membership[i][0]]);
                for (int j = 1; j < Membership[i].size(); ++j) {
                    tmpD.push_back(forwardCursor->d[Membership[i][j]]);
                    tmpA.push_back(forwardCursor->a[Membership[i][j]]);
                }
            }
            else//remained unchanged
                for (int j = 0; j < Membership[i].size(); ++j) {
                    tmpD.push_back(forwardCursor->d[Membership[i][j]]);
                    tmpA.push_back(forwardCursor->a[Membership[i][j]]);//MemberShip
                }
        }
    }
    else//stateA and stateB are directly adjacent
    {
        for (int i = 0; i < Membership.size(); ++i) {
            if (i == stateA)//put d in stateB behind of stateA
            {
                for (int j = 0; j < Membership[i].size(); ++j) {
                    tmpD.push_back(forwardCursor->d[Membership[i][j]]);
                    tmpA.push_back(forwardCursor->a[Membership[i][j]]);//MemberShip
                }
                tmpD.push_back(0);//the first element of original stateB
                tmpA.push_back(forwardCursor->a[Membership[stateB][0]]);
                for (int j = 1; j < Membership[stateB].size(); ++j) {
                    tmpD.push_back(forwardCursor->d[Membership[stateB][j]]);
                    tmpA.push_back(forwardCursor->a[Membership[stateB][j]]);
                }
            } else if(i==stateB) continue;
            else
                for (int j = 0; j < Membership[i].size(); ++j) {
                    tmpD.push_back(forwardCursor->d[Membership[i][j]]);
                    tmpA.push_back(forwardCursor->a[Membership[i][j]]);//MemberShip
                }
        }
    }


    std::copy(tmpD.begin(),tmpD.end(),forwardCursor->d);
    std::copy(tmpA.begin(),tmpA.end(),forwardCursor->a);
    fprintf(stderr,"newD1:%d\tstateM:%d,\tstateN:%d,\tstateX:%d,\tnewD2:%d,\tGetNumStates:%d\n",newD1,stateM,stateN,stateX,newD2,GetNumStates(site));
    fprintf(stderr,"after d:\n");
    for (int k = 0; k <M ; ++k) {
        fprintf(stderr,"%d\t",forwardCursor->d[k]);
    }
    fprintf(stderr,"\n");
    fprintf(stderr,"after a:\n");
    for (int k = 0; k <M ; ++k) {
        fprintf(stderr,"%d\t",forwardCursor->a[k]);
    }
    fprintf(stderr,"\n");
    return 0;
}*/
int PBWTWrapper::MoveSegment(const std::unordered_map<int, int> &mergedMembership,int site)
{
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
    int prevState(0),newD(-1);
    int lastHapID;
    int firstHapID;
    int index;
    for (int i = 0; i < clusterMembership.size(); ++i) {
        if(clusterMembership[i].size()!= 0) {
            if(prevState==0 && clusterMembership[prevState].size()==0)
            {
                newD=site+1;
            }
            else if (i != 0 && prevState != (i - 1)) {//if current state is not the first state, and prevState is not the immediate preious state
                //fprintf(stderr,"site:%d\tclusterSize:%d\ti:%d\tprevState:%d\tprev Size:%d\n",site,clusterMembership[i].size(),i,prevState,clusterMembership[prevState].size());
                lastHapID = GetHapIDFromFwd(site, clusterMembership[prevState].back());//lastHap of previous state
                firstHapID = GetHapIDFromFwd(site, clusterMembership[i].front());//firstHap of current state
                index = site;
                while (index--) {
                    if (haplotype[lastHapID][index] != haplotype[firstHapID][index]) break;
                }
                newD = index;
            }

            for (int j = 0; j < clusterMembership[i].size(); ++j) {
                if (mergedMembership.find(clusterMembership[i][j]) != mergedMembership.end()) {
                    tmpD.push_back(0);
                }
                else if (j == 0 && newD != -1) {
                    tmpD.push_back(newD);
                    newD = -1;
                }
                else {// newD==-1 means no new D value, use the old d value
                    tmpD.push_back(forwardCursor->d[clusterMembership[i][j]]);
                }
                tmpA.push_back(forwardCursor->a[clusterMembership[i][j]]);
            }
            prevState=i;
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

int PBWTWrapper::RemoveIndividualFromPBWT(int individualToProcess)
{
    char * haps[2];
    haps[0]=haplotype[2 * individualToProcess];
    haps[1]=haplotype[2 * individualToProcess + 1];

    char* tmpHapA,*tmpHapB;
    tmpHapA = this->haplotype[individualToProcess * 2];
    tmpHapB = this->haplotype[individualToProcess * 2 + 1];
    this->haplotype[individualToProcess * 2] = this->haplotype[M - 2];
    this->haplotype[individualToProcess * 2 + 1] = this->haplotype[M - 1];
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


