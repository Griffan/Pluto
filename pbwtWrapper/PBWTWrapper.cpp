//
// Created by Fan Zhang on 7/20/15.
//

#include "PBWTWrapper.h"


PBWTWrapper::PBWTWrapper(const char **haplotype, int nhaps, int nsnps) {

    int version = 2;
    printf("Read %i SNPs %i haplotypes and %i individuals from PHASE format version %i\n", nsnps, nhaps, nhaps / 2.,
           version);
    PBWT *p = pbwtCreate(nhaps, nsnps);
    //p->chrom = strdup (fgetword(fp)) ; /* example 4th line is P followed by site positions */
    //p->sites = arrayCreate (4096, Site) ;
    //int i ; for (i = 0 ; i < p->N ; ++i) arrayp(p->sites,i,Site)->x = atoi (fgetword(fp)) ;

    int i,j;
    p->CompressedAllele = arrayCreate(4096 * 32, uchar);
    PbwtCursor *u = pbwtCursorCreate(p, TRUE, TRUE);
    for (i = 0; i < p->N; ++i) {
        for (j = 0; j < p->M; ++j)
            u->y[j] = haplotype[i][u->a[j]];
        pbwtCursorWriteForwards(u);
        if (nCheckPoint && !((i + 1) % nCheckPoint)) pbwtCheckPoint(u, p);
    }
    pbwtCursorToAFend(u, p);

    fprintf(stderr, "transform phase file");
    if (p->chrom) fprintf(stderr, " for chromosome %s", p->chrom);
    fprintf(stderr, ": M, N are\t%d\t%d; yz length is %ld\n", p->M, p->N, arrayMax(p->CompressedAllele));

    pbwtCursorDestroy(u);
    return p;

}

PBWTWrapper::PBWTWrapper(int nhaps, int nsnps) {
    pbwtCore = pbwtCreate(nhaps, nsnps);
    int i,j;
    pbwtCore->CompressedAllele = arrayCreate(4096 * 32, uchar);
    InitializeCursor(TRUE, TRUE);

}



int PBWTWrapper::InitializeCursor(BOOL isForwards, BOOL isStart) {

    forwardCursor = pbwtCursorCreate(pbwtCore, TRUE, TRUE);
    return 0;
}


int PBWTWrapper::InitializeReverseCursor(BOOL isForwards, BOOL isStart) {
    backwardCursor = pbwtCursorCreate(pbwtCoreReverse, FALSE, TRUE);
    return 0;
}


int PBWTWrapper::CursorForwards() {//so far only implemented for test purpose

    for(int k=0;k!=forwardCursor->N;++k)
        CursorForwardsTo(k,10);
    //copy end of a to PBWT
    pbwtCursorToAFend(forwardCursor,PBWT);
    //update crossover rate?
    return 0;
}

int PBWTWrapper::CursorBackwards() {
    return 0;
}

int PBWTWrapper::CursorForwardsTo(int k, int T) {/*T is the length that how far you look back*/
    int i, i0 = 0, ia, ib, na = 0, nb = 0, dmin, k;
    int group=0;

    /*coppy array d*/
    int * tmpD = new int [forwardCursor->M+1];
    memcpy (tmpD, forwardCursor->d, (forwardCursor->M+1) * sizeof(int));
    /*reprot haolotype cluster based on prefix, so current site not included*/

    int u = 0, v = 0;
    int p = k + 1;
    int q = k + 1;

    for (i = 0; i < forwardCursor->M; ++i) {

        if (forwardCursor->d[i] > p) p = forwardCursor->d[i];
        if (forwardCursor->d[i] > q) q = forwardCursor->d[i];

        if (tmpD[i] > k - T) {//if current sequence and last sequence have common sequence longer than T
            if (na && nb)        /* then there is something to report */
            {
                for (ia = i0; ia < i; ++ia)
//                    for (ib = ia + 1, dmin = 0; ib < i; ++ib) {//compare pairwsely from last >T sequence to current one
//                        if (tmpD[ib] > dmin) dmin = tmpD[ib];//update shortest common length
//                        if (forwardCursor->y[ib] != forwardCursor->y[ia]) {//sequence ib and ia end differently at k
//                            //(*report)(forwardCursor->a[ia], forwardCursor->a[ib], dmin, k);
//
//                        }
//                    }
                {
                    forwardCursor->haplotypeCluster[k][ia]=group;
                }
            }
            na = 0;
            nb = 0;
            i0 = i;
            group++;
        }
        if (forwardCursor->y[i] == 0) {
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
    forwardCursor->c=na;
    memcpy (forwardCursor->a + u, forwardCursor->b, v * sizeof(int));
    memcpy (forwardCursor->d + u, forwardCursor->e, v * sizeof(int));
    forwardCursor->d[0] = k + 2;
    forwardCursor->d[forwardCursor->M] = k + 2; /* sentinels */
    //pbwtCursorForwardsReadAD(forwardCursor, k);
   // updateCursorForwards();//

    return 0;
}

int PBWTWrapper::CursorBackwardsTo() {

    return 0;
}
