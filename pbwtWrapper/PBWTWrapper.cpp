//
// Created by Fan Zhang on 7/20/15.
//

#include <unordered_map>
#include <map>
#include "PBWTWrapper.h"
#include "iostream"
#include "pbwt/pbwt.h"

PBWTWrapper::PBWTWrapper(const char **haplotype, int nhaps, int nsnps) {

    int version = 2;
    printf("Read %i SNPs %i haplotypes and %i individuals from PHASE format version %i\n", nsnps, nhaps, nhaps / 2,
           version);
    PBWT *p = pbwtCreate(nhaps, nsnps);
    //p->chrom = strdup (fgetword(fp)) ; /* example 4th line is P followed by site positions */
    //p->sites = arrayCreate (4096, Site) ;
    //int i ; for (i = 0 ; i < p->N ; ++i) arrayp(p->sites,i,Site)->x = atoi (fgetword(fp)) ;

    int i, j;
    p->CompressedAllele = arrayCreate(4096 * 32, uchar);
    PbwtCursor *u = pbwtCursorCreate(p, TRUE, TRUE);
    for (i = 0; i < p->N; ++i) {
        for (j = 0; j < p->M; ++j)
            u->sortedY[j] = haplotype[i][u->a[j]];
        pbwtCursorWriteForwards(u);
        if (nCheckPoint && !((i + 1) % nCheckPoint)) pbwtCheckPoint(u, p);
    }
    pbwtCursorToAFend(u, p);

    fprintf(stderr, "transform phase file");
    //if (p->chrom) fprintf(stderr, " for chromosome %s", p->chrom);
    fprintf(stderr, ": M, N are\t%d\t%d; yz length is %ld\n", p->M, p->N, arrayMax(p->CompressedAllele));

    pbwtCursorDestroy(u);
    pbwtCore=p;

}

PBWTWrapper::PBWTWrapper(int nhaps, int nsnps):a(nsnps,std::vector<int>(nhaps,0)),alpha(a),d(a), haplotypeCluster(a),numZero(nsnps,0),numCluster(numZero){
    pbwtCore = pbwtCreate(nhaps, nsnps);
    //pbwtCore->CompressedAllele = arrayCreate(4096 * 32, uchar);
    InitializeCursor(TRUE, TRUE);
    reverseCursor=0;

}


int PBWTWrapper::InitializeCursor(BOOL isForwards, BOOL isStart) {

    forwardCursor = pbwtCursorCreate(pbwtCore, TRUE, TRUE);
    return 0;
}


int PBWTWrapper::InitializeReverseCursor(BOOL isForwards, BOOL isStart) {
    // TODO:initialize it with information already get from forward loop
    reverseCursor = pbwtCursorCreate(pbwtCore, FALSE, TRUE);
    return 0;
}


int PBWTWrapper::CursorForwards() {//so far only implemented for test purpose



    for (int k = 0; k != pbwtCore->N; ++k) {
        //copy haplotypes into forwardCursor->y
        CopyHap(k, forwardCursor);
        CursorForwardsTo(k, 2);
    }
    //copy end of a to PBWT
    pbwtCursorToAFend(forwardCursor, pbwtCore);
    //update crossover rate?
    return 0;
}


int PBWTWrapper::CursorForwardsTo(int k, int T) {
/*T is the length that how far you look back
 *This function must be called along the sites, no skip permitted;
 *Mask the site you want to skip at the begining if you have to.
 */
    int i, i0 = 0, ia, ib, na = 0, nb = 0, dmin;
    int group = 0;

    /*coppy array d*/
    int *tmpD = new int[forwardCursor->M + 1];
    memcpy(tmpD, forwardCursor->d, (forwardCursor->M + 1) * sizeof(int));
    /*reprot haolotype cluster based on prefix, so current site not included*/
    fprintf(stderr,"k:%d,T:%d\t",k,T);PrintVector(tmpD,forwardCursor->M,"tmpD");
    fprintf(stderr,"k:%d,T:%d\t",k,T);PrintVector(forwardCursor->a,forwardCursor->M,"olda");
    //fprintf(stderr,"k:%d,T:%d\t",k,T);PrintVector(forwardCursor->sortedY,forwardCursor->M,"sortedY");
    int u = 0, v = 0;
    int p = k + 1;
    int q = k + 1;

    for (i = 0; i < forwardCursor->M; ++i) {

        if (forwardCursor->d[i] > p) p = forwardCursor->d[i];
        if (forwardCursor->d[i] > q) q = forwardCursor->d[i];

        if (tmpD[i] > k - T) {//if current sequence and last sequence have common sequence longer than T
            //if (na && nb)        /* then there is something to report */
            {
                for (ia = i0; ia < i; ++ia)
                {
                    haplotypeCluster[k][ia] = group;
                }
            }
            na = 0;
            nb = 0;
            i0 = i;
            group++;
        }
        if (forwardCursor->sortedY[i] == 0) {
            forwardCursor->a[u] = forwardCursor->a[i];
            forwardCursor->d[u] = p;
            ++u;
            p = 0;
            na++;
        }
        else {
            forwardCursor->b[v] = forwardCursor->a[i];
            forwardCursor->e[v] = q;
            ++v;
            q = 0;
            nb++;
        }
    }
    numCluster[k] = group;
    forwardCursor->c = na;
    numZero[k]=na;
    memcpy(forwardCursor->a + u, forwardCursor->b, v * sizeof(int));
    memcpy(forwardCursor->d + u, forwardCursor->e, v * sizeof(int));
    forwardCursor->d[0] = k + 2;
    forwardCursor->d[forwardCursor->M] = k + 2; /* sentinels */
    a[k].assign(forwardCursor->a,forwardCursor->a+forwardCursor->M+1);
    d[k].assign(forwardCursor->d,forwardCursor->d+forwardCursor->M+1);
    //pbwtCursorForwardsReadAD(forwardCursor, k);
    // updateCursorForwards();//

    return 0;
}


int PBWTWrapper::CursorBackwards() {
    //in the following code, we assume forwardCursor is ready and just finished forward loop
    int i, j, M = pbwtCore->M;

    if (pbwtCore->aFend)
        forwardCursor = pbwtCursorCreate(pbwtCore, TRUE, FALSE);
    else {
        forwardCursor = pbwtCursorCreate(pbwtCore, TRUE, TRUE);
        for (i = 0; i < pbwtCore->N; ++i)    /* first run forwards to the end */
            pbwtCursorForwardsRead(forwardCursor);
        pbwtCursorToAFend(forwardCursor, pbwtCore);
        //error("Please double check the completeness of PBWT structure, I can't find aFend!");
    }

    /* use p->aFend also to start the reverse cursor - this gives better performance */
    if (!pbwtCore->aRstart) pbwtCore->aRstart = new int[M];
    memcpy(pbwtCore->aRstart, forwardCursor->a, M * sizeof(int));// is Rstart the same as Fend and uF->a?
    //pbwtCore->ReverseCompressedAllele = arrayReCreate (pbwtCore->ReverseCompressedAllele, arrayMax(pbwtCore->CompressedAllele),uchar);// I didn't actually use this array
    reverseCursor = pbwtCursorCreate(pbwtCore, FALSE, TRUE); /* will pick up aRstart */

    //isolated from context
    for (i = pbwtCore->N; i--;) {

        CursorBackwardsTo(i, 5);
    }
    //isolated from context


    /* save uR->a, which is the lexicographic order of the sequences */
    if (!pbwtCore->aRend) pbwtCore->aRend = myalloc (M, int);
    memcpy(pbwtCore->aRend, reverseCursor->a, M * sizeof(int));//the end when loop from back to the original first

    // fprintf(stderr, "built reverse PBWT - size %ld\n", arrayMax(p->ReverseCompressedAllele));

    if (isCheck)            /* print out the reversed haplotypes */
    {
        FILE *fp = fopen("rev.haps", "w");
        Array tz = pbwtCore->CompressedAllele;
        pbwtCore->CompressedAllele = pbwtCore->ReverseCompressedAllele;
        int *ta = pbwtCore->aFstart;
        pbwtCore->aFstart = pbwtCore->aRstart;
        pbwtWriteHaplotypes(fp, pbwtCore);
        pbwtCore->CompressedAllele = tz;
        pbwtCore->aFstart = ta;
    }

    return 0;
}

int PBWTWrapper::CursorBackwardsTo(int k, int T) {//this function must be call in order, no skip allowed
   // int j;
    //current status: forwardCursor's sortedY and a both stopped at the final site
    //
    //uchar *x = new uchar[forwardCursor->M];
    //pbwtCursorReadBackwards(forwardCursor);
//    pbwtCursorBackwardsA(forwardCursor);
//    for (j = 0; j < forwardCursor->M; ++j)
//        x[forwardCursor->a[j]] = forwardCursor->sortedY[j];//x has the order same as haplotype
//    for (j = 0; j < forwardCursor->M; ++j)
//        reverseCursor->sortedY[j] = x[reverseCursor->a[j]];// I think uR->a is the same as uF->a
    //pbwtCursorWriteForwards(reverseCursor);
    // delete[] x;
    CopyHap(k,reverseCursor);
    pbwtCursorForwardsA(reverseCursor);
    MergeCluster(k);
    alpha[k].assign(reverseCursor->a,reverseCursor->a+reverseCursor->M);
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
        if (haplotype[Cursor->a[i]][k] != 0 && haplotype[Cursor->a[i]][k] != 1) Cursor->sortedY[i] =
                                                                                        haplotype[Cursor->a[i]][k] -
                                                                                        '0';
        else fprintf(stderr,"alert!!!! %d,%d,%d,%d\n",haplotype[k][Cursor->a[i]],k,i,Cursor->a[i]);
    }

    //PrintVector(Cursor->sortedY,Cursor->M,"fromCopyHap");
    return 0;
}

int PBWTWrapper::MergeCluster(int site) {
    int lastNumCluster=numCluster[site];
    std::vector<std::vector<int> > clusterMemberShip(lastNumCluster,std::vector<int>());
    std::vector<std::vector<int> > dist(lastNumCluster,std::vector<int>(haplotypeCluster[site].size(),0));
    std::unordered_map<int, bool> mergeIndicator;
    std::vector<unsigned long> order(forwardCursor->M,0);
    for(int i=0;i!=reverseCursor->M;++i)
    {
        order[reverseCursor->a[i]]=i;//record where the ith sequence now is
    }
    for (int i = 0; i != haplotypeCluster[site].size(); ++i) {
        mergeIndicator[haplotypeCluster[site][i]]=false;
        clusterMemberShip[haplotypeCluster[site][i]].push_back(i);
        dist[haplotypeCluster[site][i]][order[i]]=1;//record rank distribution for each cluster
    }
    for(auto kv:dist)
    {
        for(int i(0),v(0);i!=kv.size();++i)
        {
            if(kv[i]==1) v++;
            kv[i]=v;
        }
    }
    if(DEBUG)
    {
        std::cerr<<"enter debug section:"<<std::endl;
        for (auto i = 0; i != haplotypeCluster.size(); ++i) {
            PrintDistributionAtSite(i,haplotypeCluster[i]);
        }
        std::cerr<<"exit debug section!"<<std::endl;
        exit(0);
    }
    int currentNumCluster=lastNumCluster-1;
    while(currentNumCluster!=lastNumCluster) {
        lastNumCluster=currentNumCluster;

        for (auto i = 0; i != dist.size(); ++i) {
            if(mergeIndicator[i]) continue;
            for (auto j = i + 1; j != dist.size(); ++j) {
                if(mergeIndicator[j]) continue;
                if(KStest(dist[i],dist[j]))
                {
                    currentNumCluster--;
                    /*Merge Action, change mergeIndicator*/
                    for(int t=0;t!=dist[i].size();++t)
                    {
                        dist[i][t]+=dist[j][t];
                    }
                    mergeIndicator[j]=true;// j th cluster has been merged into i th cluster
                    for(int t=0;t!=clusterMemberShip[j].size();++t)
                    {
                        haplotypeCluster[site][t]=i;
                    }
                }
            }
            PrintDistributionAtSite(i,dist[i]);
        }
        std::cerr<<"finish of last round"<<std::endl;
    }
    return false;
}

bool PBWTWrapper::KStest(std::vector<int>& a, std::vector<int>& b) {

    int Dmax=0;
    int Dtmp=0;
    for(int i=0;i!=a.size();++i)
    {
        Dtmp=abs(a[i]-b[i]);
        if(Dtmp>Dmax) Dmax=Dtmp;
    }

    if(Dmax > 1.36*sqrt(double(a.back()+b.back())/(a.back()*b.back())))//1.36 is 0.05 significance parameter
        return false;//reject null hypo, they are different
    else
        return true;//accept null hypo, they are the same
}

int PBWTWrapper::setHaps(char **haps) {
    haplotype=haps;
    return 0;
}

int PBWTWrapper::PrintDistributionAtSite(int state,std::vector<int> &dist) {
    std::cerr<<"state:"<<state<<":\t";
    for(auto k:dist)
        std::cerr<<k<<"\t";
    std::cerr<<std::endl;
    return 0;
}
