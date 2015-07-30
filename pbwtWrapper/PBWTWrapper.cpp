//
// Created by Fan Zhang on 7/20/15.
//

#include "PBWTWrapper.h"
#include "../SinglePhasing/libStatGen/general/Error.h"


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
            u->sortedY[j] = haplotype[i][u->a[j]];
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
    // TODO:initialize it with information already get from forward loop
    backwardCursor = pbwtCursorCreate(pbwtCore, FALSE, TRUE);
    return 0;
}


int PBWTWrapper::CursorForwards() {//so far only implemented for test purpose



    for(int k=0;k!=pbwtCore->N;++k) {
        //copy haplotypes into forwardCursor->y
        CopyHap(k,forwardCursor);
        CursorForwardsTo(k, 10);
    }
    //copy end of a to PBWT
    pbwtCursorToAFend(forwardCursor,pbwtCore);
    //update crossover rate?
    return 0;
}

int PBWTWrapper::CursorBackwards() {
	//in the following code, we assume forwardCursor is ready and just finished forward loop
    int i, j, M = pbwtCore->M;

  //  if (pbwtCore->aFend)
		//uF = pbwtCursorCreate(pbwtCore, TRUE, FALSE);
  //  else {
		//uF = pbwtCursorCreate(pbwtCore, TRUE, TRUE);
  //      for (i = 0; i < p->N; ++i)    /* first run forwards to the end */
  //          pbwtCursorForwardsRead(uF);
  //      pbwtCursorToAFend(uF, p);
  //      error("Please double check the completeness of PBWT structure, I can't find aFend!");
  //  }

    /* use p->aFend also to start the reverse cursor - this gives better performance */
    if (!pbwtCore->aRstart) pbwtCore->aRstart = new int [M];
    memcpy (pbwtCore->aRstart, forwardCursor->a, M * sizeof(int));// is Rstart the same as Fend and uF->a?
    pbwtCore->ReverseCompressedAllele = arrayReCreate (pbwtCore->ReverseCompressedAllele, arrayMax(pbwtCore->CompressedAllele), uchar);// I didn't actually use this array
    //PbwtCursor *uR = pbwtCursorCreate(pbwtCore, FALSE, TRUE); /* will pick up aRstart */
	reverseCursor = pbwtCursorCreate(pbwtCore, FALSE, TRUE); /* will pick up aRstart */

    //isolated from context
    for (i = pbwtCore->N; i--;) {
        //CopyHap(i,reverseCursor);
        CursorBackwardsTo(i,5);
    }
    //isolated from context


    /* save uR->a, which is the lexicographic order of the sequences */
    if (!pbwtCore->aRend) pbwtCore->aRend = myalloc (M, int);
    memcpy (pbwtCore->aRend, uR->a, M * sizeof(int));//the end when loop from back to the original first

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

int PBWTWrapper::CursorForwardsTo(int k, int T) {
/*T is the length that how far you look back
 *This function must be called along the sites, no skip permitted;
 *Mask the site you want to skip at the begining if you have to.
 */
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
                   haplotypeCluster[k][ia]=group;
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
    numCluster[k]=group;
    forwardCursor->c=na;
    memcpy (forwardCursor->a + u, forwardCursor->b, v * sizeof(int));
    memcpy (forwardCursor->d + u, forwardCursor->e, v * sizeof(int));
    forwardCursor->d[0] = k + 2;
    forwardCursor->d[forwardCursor->M] = k + 2; /* sentinels */
    //pbwtCursorForwardsReadAD(forwardCursor, k);
   // updateCursorForwards();//

    return 0;
}

int PBWTWrapper::ObtainHapFromSinglePhasing(char ** haps)
{
	haplotype = haps;
	pbwtCore->CompressedAllele = arrayCreate(4096 * 32, uchar);
	forwardCursor = pbwtCursorCreate(pbwtCore, TRUE, TRUE);
	for (int i = 0; i < pbwtCore->N; ++i)
	{
		for (int j = 0; j < pbwtCore->M; ++j) forwardCursor->sortedY[j] = haplotype[i][forwardCursor->a[j]];
		pbwtCursorWriteForwards(forwardCursor);
		if (nCheckPoint && !((i + 1) % nCheckPoint)) pbwtCheckPoint(forwardCursor, pbwtCore);
	}
	pbwtCursorToAFend(forwardCursor, pbwtCore);
	return 0;
}
int PBWTWrapper::CursorBackwardsTo(int k, int T=5) {
    int j;
	int M = pbwtCore->M;
    uchar *x = new uchar[M];

	//current status: forwardCursor's sortedY and a both stopped at the final site
	//

        //pbwtCursorReadBackwards(uF);
        for (j = 0; j < M; ++j) x[forwardCursor->a[j]] = forwardCursor->sortedY[j];//transform original order into forwardEnd order
		for (j = 0; j < M; ++j) reverseCursor->sortedY[j] = x[reverseCursor->a[j]];// I think uR->a is the same as uF->a
        //pbwtCursorWriteForwards(uR);
    delete [] x;

    return 0;
}

int PBWTWrapper::CopyHap(int k,PbwtCursor* Cursor) {//this function has the same effect as forward/backward read
	for (int i = 0; i != Cursor->M; ++i)
        Cursor->sortedY[i]= haplotype[k][i];
    return 0;
}
