//
// Created by Fan Zhang on 7/20/15.
//

#include <unordered_map>
#include <map>
#include <cmath>
#include "PBWTWrapper.h"
#include "iostream"
#include "pbwt/pbwt.h"

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

PBWTWrapper::PBWTWrapper(int nhaps, int nsnps):a(nsnps,std::vector<int>(nhaps,0)),alpha(a),d(a), sortedY(nsnps,std::vector<uchar>(nhaps,0)),haplotypeCluster(a),clusterAllele(nsnps,std::vector<uchar>()){
    N=nsnps;
    M=nhaps;//last two haps are slots for current individual need to be phased
    pbwtCore = pbwtCreate(nhaps, nsnps);
    //pbwtCore->CompressedAllele = arrayCreate(4096 * 32, uchar);

    forwardCursor = pbwtCursorCreate(pbwtCore, TRUE, TRUE);

    reverseCursor = pbwtCursorCreate(pbwtCore, FALSE, TRUE);

}


int PBWTWrapper::CursorForwards() {//so far only implemented for test purpose


    //PrintVector(forwardCursor->a,M,"end arrary aFend check 0");
    for (int k = 0; k != pbwtCore->N; ++k) {
        CursorForwardsTo(k, 200);
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
    int i, i0 = 0, ia, ib, na = 0, nb = 0, dmin;
    int group = 0;

    /*coppy array d*/
   // int *lastD = new int[forwardCursor->M + 1];
    uchar *lastY = new uchar[forwardCursor->M +1];
   // memcpy(lastD, forwardCursor->d, (forwardCursor->M + 1) * sizeof(int));
    memcpy(lastY, forwardCursor->sortedY, (forwardCursor->M + 1) * sizeof(uchar));
    //copy haplotypes into forwardCursor->y
    CopyHap(k, forwardCursor);

    if (k==pbwtCore->N-1)//deal with last columns
    {
        for (i = 0; i < forwardCursor->M; ++i) {
            haplotypeCluster[k][i] = 0;
        }

        clusterAllele[k].push_back(forwardCursor->sortedY[0]);
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
         if (k!=pbwtCore->N-1 and forwardCursor->d[i] > k - tmpT) {//if current sequence and last sequence have common sequence longer than T
            //if (na && nb)        /* then there is something to report */
            {
                for (ia = i0; ia < i; ++ia)
                {
                    haplotypeCluster[k][ia] = group;
                }

                clusterAllele[k].push_back(lastY[i0]);
            }
           // na = 0;
           // nb = 0;
            i0 = i;
            group++;
        }

        if (forwardCursor->sortedY[i] == 0) {
            forwardCursor->a[u] = forwardCursor->a[i];
            forwardCursor->d[u] = p;
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

    if( (k!=pbwtCore->N-1) and i0 < forwardCursor->M-1)
    {
        for (ia = i0; ia < forwardCursor->M; ++ia)
        {
            haplotypeCluster[k][ia] = group;
        }

        clusterAllele[k].push_back(lastY[i0]);
    }

    MergeCluster(k);//TODO:implement this function
    //numCluster[k] = group;
    //forwardCursor->c = na;
    //numZero[k]=na;
    memcpy(forwardCursor->a + u, forwardCursor->b, v * sizeof(int));
    memcpy(forwardCursor->d + u, forwardCursor->e, v * sizeof(int));
    //forwardCursor->d[0] = k + 2;
    forwardCursor->d[forwardCursor->M] = k + 2; /* sentinels */
    a[k].assign(forwardCursor->a,forwardCursor->a+forwardCursor->M);
    d[k].assign(forwardCursor->d,forwardCursor->d+forwardCursor->M);
    sortedY[k].assign(forwardCursor->sortedY,forwardCursor->sortedY+forwardCursor->M);
    //delete [] lastD;
    delete [] lastY;

    //pbwtCursorForwardsReadAD(forwardCursor, k);
    // updateCursorForwards();//
    //fprintf(stderr,"k:%d,T:%d\t",k,T);PrintVector(forwardCursor->d,forwardCursor->M,"after tmpD");

    return 0;
}
/*
int PBWTWrapper::CursorBackwards() {
    //in the following code, we assume forwardCursor is ready and just finished forward loop
    int i, j, M = pbwtCore->M;

    //PrintVector(pbwtCore->aFend,M,"end arrary aFend check 2");

    if (pbwtCore->aFend)//if aFend exists in pbwtCore, we directly use it
    {
        //fprintf(stderr,"create new cursor true false\n");
        forwardCursor = pbwtCursorCreate(pbwtCore, TRUE, FALSE);
    }
    else {//then let's calculate it from beginning
        forwardCursor = pbwtCursorCreate(pbwtCore, TRUE, TRUE);
        for (i = 0; i < pbwtCore->N; ++i)    //first run forwards to the end
            pbwtCursorForwardsRead(forwardCursor);
        pbwtCursorToAFend(forwardCursor, pbwtCore);
        //error("Please double check the completeness of PBWT structure, I can't find aFend!");
    }

    // use p->aFend also to start the reverse cursor - this gives better performance
    if (!pbwtCore->aRstart) pbwtCore->aRstart = new int[M];
    memcpy(pbwtCore->aRstart, forwardCursor->a, M * sizeof(int));// is Rstart the same as Fend and uF->a?

    //pbwtCore->ReverseCompressedAllele = arrayReCreate (pbwtCore->ReverseCompressedAllele, arrayMax(pbwtCore->CompressedAllele),uchar);// I didn't actually use this array
    reverseCursor = pbwtCursorCreate(pbwtCore, FALSE, TRUE); //will pick up aRstart

    //isolated from context
    for (i = pbwtCore->N; i--;) {

        CursorBackwardsTo(i, 5);
    }
    //isolated from context

    for (i = pbwtCore->N; i--;) {

        UpdateTransVector(i);
    }
    PrintSummary();
    //save uR->a, which is the lexicographic order of the sequences
    if (!pbwtCore->aRend) pbwtCore->aRend = myalloc (M, int);
    memcpy(pbwtCore->aRend, reverseCursor->a, M * sizeof(int));//the end when loop from back to the original first

    // fprintf(stderr, "built reverse PBWT - size %ld\n", arrayMax(p->ReverseCompressedAllele));

    if (isCheck)            // print out the reversed haplotypes
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
}*/

int PBWTWrapper::CursorBackwards() {
    for (int i = pbwtCore->N; i--;) {

        CursorBackwardsTo(i, 5);
    }
    return 0;
}
/*
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
}*/
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
        }
        else {
            reverseCursor->b[v] = reverseCursor->a[i];
            ++v;
        }
    }

    memcpy(reverseCursor->a + u, reverseCursor->b, v * sizeof(int));
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

    //fprintf(stderr,"site:%d\tprevSite:%d\tsite:%d\n",site,getNumStates(prevSite),getNumStates(site));

    transVector.insert(std::make_pair(prevSite,std::vector<std::vector<float> >(getNumStates(prevSite),std::vector<float>(getNumStates(site), 0.))));

	std::vector<float> marginal(getNumStates(prevSite), 0.000001);
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
	for (int i = 0; i != getNumStates(prevSite); ++i)
	{
		for (int j = 0; j != getNumStates(site); ++j) {
            transVector[prevSite][i][j] /= marginal[i];
            //fprintf(stderr,"i:%d to j:%d is %f\t",i,j,transVector[prevSite][i][j]);
        }
       //fprintf(stderr,"\n");
	}

    //fprintf(stderr,"finish %d and:prevStates:%d,States:%d\n",prevSite,getNumStates(prevSite),getNumStates(site));

	return 0;
}
/*
int PBWTWrapper::MergeCluster(int site) {
    int oldNumCluster = getNumStates(site);
    int numHaps = haplotypeCluster[site].size();
    std::vector<uchar> tmpAllele;
    std::vector<std::vector<int> > clusterMemberShip(oldNumCluster,std::vector<int>());
    std::vector<std::vector<int> > dist(oldNumCluster,std::vector<int>(numHaps,0));
    std::unordered_map<int, bool> mergeIndicator;
    std::vector<unsigned long> order(forwardCursor->M,0);
    for(int i=0;i!=reverseCursor->M;++i)
    {
        order[reverseCursor->a[i]]=i;//record where the ith sequence now is
    }
    for (int i = 0; i != numHaps; ++i) {
        mergeIndicator[haplotypeCluster[site][i]]=false;//initialize states' merge status

        clusterMemberShip[haplotypeCluster[site][i]].push_back(i);//put haps in the same state into same vector

        dist[haplotypeCluster[site][i]][order[i]]=1;//record rank distribution for each cluster
    }
    for(int i=0;i!= dist.size();++i)
    {
        for(int j(0),v(0);j!=dist[i].size();++j)
        {
            if(dist[i][j]==1) v++;
            dist[i][j]=v;
        }
    }
    if(0&&DEBUG)
    {
        std::cerr<<"\nenter debug section:"<<std::endl;
        for (auto i = 0; i != dist.size(); ++i) {
            //PrintDistributionAtSite(i,haplotypeCluster[i]);
            PrintDistributionAtSite(i,dist[i]);
        }
        std::cerr<<"exit debug section!"<<std::endl;
        //exit(0);
    }

    int currentNumCluster=oldNumCluster-1;
    std::unordered_map<int,int> stateOrder;//mapping oldState to newOrder
    int tmpOrder(0);

    while(currentNumCluster!=oldNumCluster) {
        oldNumCluster=currentNumCluster;
        tmpAllele.clear();//TODO:tmpAllele order
        stateOrder.clear();
        tmpOrder=0;
        for (auto i = 0; i != dist.size(); ++i) {
            if(mergeIndicator[i]) continue;
            else
            {
                tmpAllele.push_back(clusterAllele[site][i]);
                stateOrder[i]=tmpOrder;
                tmpOrder++;
            }

            for (auto j = i + 1; j != dist.size(); ++j) {
                if(mergeIndicator[j]) continue;
                if(KStest(dist[i],dist[j]))
                {
                    currentNumCluster--;
                    //Merge Action, change mergeIndicator
                    for(int t=0;t!=dist[i].size();++t)
                    {
                        dist[i][t]+=dist[j][t];
                    }
                    mergeIndicator[j]=true;// j th cluster has been merged into i th cluster

                    for(int t=0;t!=clusterMemberShip[j].size();++t)
                    {
                        haplotypeCluster[site][clusterMemberShip[j][t]]=i;
                        clusterMemberShip[i].push_back(clusterMemberShip[j][t]);
                    }
                    clusterMemberShip[j].clear();
                    if(DEBUG)std::cerr<<"Merge\t"<<j<<"\tinto\t"<<i<<std::endl;

                    //if(clusterAllele[site][i]!=clusterAllele[site][j]) die((char*)"alert: two states ready to be merged have different allele");

                    break;
                }
                else
                {
                    if(DEBUG) std::cerr<<"Cannot merge\t"<<j<<"\tinto\t"<<i<<std::endl;
                }
            }
            //PrintDistributionAtSite(i,dist[i]);
        }
        if(DEBUG)std::cerr<<"finish of last round"<<std::endl;
    }
    for (int k = 0; k < haplotypeCluster[site].size(); ++k) {
        haplotypeCluster[site][k]=stateOrder[haplotypeCluster[site][k]];
    }
    clusterAllele[site]=tmpAllele;//update merged cluster allele
    return false;
}
*/

int PBWTWrapper::MergeCluster(int site) {
    int oldNumCluster = getNumStates(site);
    int numHaps = haplotypeCluster[site].size();
    std::vector<uchar> tmpAllele;
    std::vector<std::vector<int> > clusterMemberShip(oldNumCluster,std::vector<int>());
    std::vector<std::vector<int> > dist(oldNumCluster,std::vector<int>(numHaps,0));
    std::unordered_map<int, bool> mergeIndicator;
    std::vector<unsigned long> order(forwardCursor->M,0);
    for(int i=0;i!=reverseCursor->M;++i)
    {
        order[reverseCursor->a[i]]=i;//record where the ith sequence now is
    }
    for (int i = 0; i != numHaps; ++i) {
        mergeIndicator[haplotypeCluster[site][i]]=false;//initialize states' merge status

        clusterMemberShip[haplotypeCluster[site][i]].push_back(i);//put haps in the same state into same vector

        dist[haplotypeCluster[site][i]][order[i]]=1;//record rank distribution for each cluster
    }
    for(int i=0;i!= dist.size();++i)
    {
        for(int j(0),v(0);j!=dist[i].size();++j)
        {
            if(dist[i][j]==1) v++;
            dist[i][j]=v;
        }
    }
    if(0&&DEBUG)
    {
        std::cerr<<"\nenter debug section:"<<std::endl;
        for (auto i = 0; i != dist.size(); ++i) {
            //PrintDistributionAtSite(i,haplotypeCluster[i]);
            PrintDistributionAtSite(i,dist[i]);
        }
        std::cerr<<"exit debug section!"<<std::endl;
        //exit(0);
    }

    int currentNumCluster=oldNumCluster-1;
    std::unordered_map<int,int> stateOrder;//mapping oldState to newOrder
    int tmpOrder(0);

    while(currentNumCluster!=oldNumCluster) {
        oldNumCluster=currentNumCluster;
        tmpAllele.clear();//TODO:tmpAllele order
        stateOrder.clear();
        tmpOrder=0;

        for (auto i = 0; i != dist.size(); ++i) {
            if(mergeIndicator[i])
            {
                continue;
            }
            else
            {
                tmpAllele.push_back(clusterAllele[site][i]);
                stateOrder[i]=tmpOrder;
                tmpOrder++;
            }

            for (auto j = i + 1; j != dist.size(); ++j) {
                if(mergeIndicator[j]) continue;
                if(KStest(dist[i],dist[j]))
                {
                    currentNumCluster--;
                    //Merge Action, change mergeIndicator
                    for(int t=0;t!=dist[i].size();++t)
                    {
                        dist[i][t]+=dist[j][t];
                    }
                    mergeIndicator[j]=true;// j th cluster has been merged into i th cluster

                    for(int t=0;t!=clusterMemberShip[j].size();++t)
                    {
                        haplotypeCluster[site][clusterMemberShip[j][t]]=i;
                        clusterMemberShip[i].push_back(clusterMemberShip[j][t]);
                    }
                    clusterMemberShip[j].clear();
                    if(DEBUG)std::cerr<<"Merge\t"<<j<<"\tinto\t"<<i<<std::endl;

                    //if(clusterAllele[site][i]!=clusterAllele[site][j]) die((char*)"alert: two states ready to be merged have different allele");

                    break;
                }
            }
            //PrintDistributionAtSite(i,dist[i]);
        }
        if(DEBUG)std::cerr<<"finish of last round"<<std::endl;
    }
    //adjust d array and a array
    MoveSegment(clusterMemberShip);

    for (int k = 0; k < haplotypeCluster[site].size(); ++k) {
        haplotypeCluster[site][k]=stateOrder[haplotypeCluster[site][k]];
    }
    clusterAllele[site]=tmpAllele;//update merged cluster allele
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

    if(Dmax > 1.36*std::sqrt(double(a.back()+b.back())/(a.back()*b.back())))//1.36 is 0.05 significance parameter
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
    printf("mean edges/node = %f\tmean count/mode = %f\n",(float)totalEdges/totalNodes,(float)(pbwtCore->M)*(pbwtCore->N)/totalNodes);

    return 0;
}

int PBWTWrapper::MoveSegment(std::vector<std::vector<int> >& MemberShip) {//fromEnd don't include
    std::vector<int> tmpD,tmpA;
    for (int i = 0; i <MemberShip.size() ; ++i) {
        if(MemberShip[i].size()>0) {
            tmpD.push_back(forwardCursor->d[MemberShip[i][0]]);
            tmpA.push_back(forwardCursor->a[MemberShip[i][0]]);
        }
        //int d=forwardCursor->d[MemberShip[i][0]];
        for (int j = 1; j <MemberShip[i].size() ; ++j) {
            tmpD.push_back(0);
            tmpA.push_back(forwardCursor->a[MemberShip[i][j]]);
        }
    }
    std::copy(tmpD.begin(),tmpD.end(),forwardCursor->d);
    std::copy(tmpA.begin(),tmpA.end(),forwardCursor->a);
    return 0;
}
