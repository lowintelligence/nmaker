#include "traversal.h"
#include <time.h>
#include <math.h>

__OffloadVar_Macro__
static Body *part;

__OffloadVar_Macro__
static Node *tree; // thread number

#if !BIGQUEUE
__OffloadFunc_Macro__
int getGridIndex(int* gridP,int* curIndex,Domain* dp){
	*curIndex=*gridP;
	*gridP=*gridP+1;
	while(*curIndex<=*(gridP+1)&&(dp->cuboid->tag[*curIndex]!=TAG_OUTERCORE&&dp->cuboid->tag[*curIndex]!=TAG_FRONTIERS)){
		*curIndex=*gridP;
		*gridP=*gridP+1;
	}
	if(*curIndex>*(gridP+1))
		return 1;
	else
		return 0;
		
}
#endif	

__OffloadFunc_Macro__
#if BIGQUEUE
int threaded(int nTeam, int nMaster, int nSlave, pthread_t* thread, pthread_mutex_t *mutex, Block* pth, int ts, int te, PPQ *pq, TWQ *tq, int* lower,int* upper, Domain*dp,GlobalParam*gp,void* process(void*))
#else
int threaded(pthread_t* thread, pthread_mutex_t *mutex, Block* pth, int ts, int te, PPQ *pq, TWQ *tq, int* gridP,int* curIndex,int*numpart, Domain*dp,GlobalParam*gp,void* process(void*))
#endif
{
	int i,teamid;
	printf("before threaded\n");
	pthread_attr_t attr;
	cpu_set_t cpuset;

	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	//	sleep(2);

	for (i=ts; i<te+1; i++)
//	int tend = ((ts+4)>=(te+1)) ? (te+1) : (ts+4); 
//	for (i=ts; i<tend; i++)
	{
#if BIGQUEUE
		pth[i].blockid = i/nTeam;
		pth[i].teamid=i%nTeam;
		pth[i].nTeam = nTeam;
		pth[i].nMaster = nMaster;
		pth[i].nSlave = nSlave;
		pth[i].bsize = nMaster+nSlave;
		pth[i].tall = nTeam*(nMaster+nSlave);
#else
		pth[i].blockid=0;
		pth[i].teamid=i;
		pth[i].nTeam = te+1;
		pth[i].nSlave = te+1;
#endif
		teamid=pth[i].teamid;
		pth[i].tid = i;
		pth[i].PQ_ppnode = &pq[teamid];
		pth[i].PQ_treewalk = &tq[teamid];
#if BIGQUEUE
		pth[i].first=lower[teamid];
		pth[i].last=upper[teamid];
		pth[i].lower=lower;
		pth[i].upper=upper;
#else
		pth[i].curIndex=&curIndex[teamid];
		pth[i].numpart=&numpart[teamid];
		pth[i].gridP=gridP;
#endif
		pth[i].dp=dp;
		pth[i].gp=gp;
		pth[i].thread=thread;
		pth[i].pthArr=pth;
		pth[i].P_PQ_ppnode=pq;
		pth[i].P_PQ_treewalk=tq;
		pth[i].mutex=mutex;
#ifdef __MIC__
//		CPU_ZERO(&cpuset);
//		int core_id=pth[i].blockid+teamid*(nMaster+nSlave)+dp->rank/2*pth[i].tall+1;
////		int core_id=((i+dp->rank/2*pth[i].tall)%60)*4+((i+dp->rank/2*pth[i].tall)/60)+1;
//		CPU_SET(core_id, &cpuset);
//		pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpuset);
#endif
#if BIGQUEUE
//		printf("[%d/%d] blockid=%d, teamid=%d, lower=%d, upper=%d, Ngrid=%d.\n", dp->rank, i, pth[i].blockid, pth[i].teamid, pth[i].first, pth[i].last, dp->cuboid->nSide[0]*dp->cuboid->nSide[1]*dp->cuboid->nSide[2]);
#else
//		printf("[%d/%d] blockid=%d, nTeam=%d, teamid=%d, Ngrid=%d.\n", dp->rank, i, pth[i].blockid, pth[i].nTeam, pth[i].teamid, dp->cuboid->nSide[0]*dp->cuboid->nSide[1]*dp->cuboid->nSide[2]);
#endif
		pthread_create(&thread[i], &attr, process, (void*)&pth[i]);
	}
	pthread_attr_destroy(&attr);
}

#if BIGQUEUE
__OffloadFunc_Macro__
void* compPP(void* param){
	Array3 pa, pb, pc;
	int i,j,k,n,m,l;
	Block* pth=(Block*)param;
	int blockid=pth->blockid;
	int teamid=pth->teamid;
	int bsize=pth->bsize;
	int tid=pth->tid;
	int tall=pth->tall;
	PPQ* PQ_ppnode=pth->PQ_ppnode;
	TWQ* PQ_treewalk=pth->PQ_treewalk;
	Domain* dp=pth->dp;
	GlobalParam* gp=pth->gp;
	pthread_t* thread=pth->thread;
	Block* pthArr=pth->pthArr;

	//	printf("\nQueue process finishing, total %d pp pairs.\n", PQ_ppnode->length);
	pa.x = (PRECTYPE*)memalign(ALIGNCNT, sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pa.y = (PRECTYPE*)memalign(ALIGNCNT, sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pa.z = (PRECTYPE*)memalign(ALIGNCNT, sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pb.x = (PRECTYPE*)memalign(ALIGNCNT, sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pb.y = (PRECTYPE*)memalign(ALIGNCNT, sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pb.z = (PRECTYPE*)memalign(ALIGNCNT, sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pc.x = (PRECTYPE*)memalign(ALIGNCNT, sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pc.y = (PRECTYPE*)memalign(ALIGNCNT, sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pc.z = (PRECTYPE*)memalign(ALIGNCNT, sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	
//	printf("before ppkernel, %ld\n", (1<<(MAX_MORTON_LEVEL))*sizeof(PRECTYPE)*MAX_PACKAGE_SIZE);
	double dpp = 0.0;
	usleep(2000000);
	while(PQ_ppnode->length>0||PQ_ppnode->tag==0)
	{
		if(PQ_ppnode->length>0)
		{
			dpp += ProcessQP_PPnode(PQ_ppnode, ppkernel, pa, pb, pc, &pth->mutex[teamid]);
		}
	}
//	PPelement *pe=PQ_ppnode->elements;
//	n=PQ_ppnode->length;
//	for(i=PQ_ppnode->length/2;i<n;i++)
//	{
//		int nA, nB, dq;
//
//		nA = tree[pe[i].TA].nPart;
//		packarray3(&part[tree[pe[i].TA].firstpart], nA, pa);
//		packarray3(NULL, nA, pc);
//
//		if (pe[i].mask == 0)
//		{
//			nB = tree[pe[i].TB].nPart;
//			packarray3(&part[tree[pe[i].TB].firstpart], nB, pb);
//		}
//		else
//		{
//			nB = 1;
//			pb.x[0] = tree[pe[i].TB].masscenter[0];
//			pb.y[0] = tree[pe[i].TB].masscenter[1];
//			pb.z[0] = tree[pe[i].TB].masscenter[2];
//	//		pB.x=pB.y=pB.z=NULL;	
//		}
//		ppkernel(pa, nA, pb, nB, EPS2, pc);
//
//		if (pe[i].mask == 1)
//		{	
//			pc.x[0] *= tree[pe[i].TB].mass;
//			pc.y[1] *= tree[pe[i].TB].mass;
//			pc.z[2] *= tree[pe[i].TB].mass;
//		}
//
//		pusharray3(&part[tree[pe[i].TA].firstpart], nA, pc);
//	}


//	printf("after ppkernel\n");	
//	printf("\nSlave PP processing finished, time = %f.\n", dpp);

	free(pa.x);
	free(pa.y);
	free(pa.z);
	free(pb.x);
	free(pb.y);
	free(pb.z);
	free(pc.x);
	free(pc.y);
	free(pc.z);

	//	destroy_queue_pp(&PQ_ppnode);
	//	destroy_queue_tw(&PQ_treewalk);

}
#endif

__OffloadFunc_Macro__
void* teamMaster(void* param)
{
	Array3 pa, pb, pc;
	int i,j,k,n,m,l,i1,i2,i3,cella;
	Block* pth=(Block*)param;
	int teamid=pth->teamid;
	int tid=pth->tid;
	PPQ* PQ_ppnode=pth->PQ_ppnode;
	TWQ* PQ_treewalk=pth->PQ_treewalk;
	Domain* dp=pth->dp;
	GlobalParam* gp=pth->gp;
	pthread_t* thread=pth->thread;
	Block* pthArr=pth->pthArr;
	PPQ* P_PQ_ppnode=pth->P_PQ_ppnode;
	TWQ* P_PQ_treewalk=pth->P_PQ_treewalk;
#if BIGQUEUE
	int blockid=pth->blockid;
	int bsize=pth->bsize;
	int tall=pth->tall;
	int first=pth->first;
	int last=pth->last;
	int* lower=pth->lower;
	int* upper=pth->upper;
#else
	int* curIndex=pth->curIndex;
	int* gridP=pth->gridP;
	int* numpart=pth->numpart;
	pthread_mutex_t* mutex=pth->mutex;
#endif

//	if(tid*4+4 < pth->nTeam)
//	{
//		threaded(thread, mutex, pthArr, tid*4+4, pth->nTeam-1, P_PQ_ppnode, P_PQ_treewalk, gridP, curIndex-teamid, numpart-teamid, dp, gp,teamMaster);
//	}
	int In = dp->cuboid->nSide[0];
	int Jn = dp->cuboid->nSide[1];
	int Kn = dp->cuboid->nSide[2];
	int cnt = 0;
	double open_angle = 0.3;
	int Ngrid=In*Jn*Kn;

	DomainTree *dtp = dp->domtree;
	int npart = dp->NumPart;
//	cpu_set_t cpuset;
//	pthread_getaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
	//    part = dp->Part;
	//    tree = dtp->tree;
//	printf("[%d/%d] first=%d, last=%d, Ngrid=%d, bind=%d, set=%d, part0key=%ld.\n",dp->rank, tid, first, last, In*Jn*Kn, blockid+teamid*pth->bsize+dp->rank*tall, cpuset, part[0].pos[0]);

	pa.x = (PRECTYPE*)memalign(ALIGNCNT, sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pa.y = (PRECTYPE*)memalign(ALIGNCNT, sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pa.z = (PRECTYPE*)memalign(ALIGNCNT, sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pb.x = (PRECTYPE*)memalign(ALIGNCNT, sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pb.y = (PRECTYPE*)memalign(ALIGNCNT, sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pb.z = (PRECTYPE*)memalign(ALIGNCNT, sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pc.x = (PRECTYPE*)memalign(ALIGNCNT, sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pc.y = (PRECTYPE*)memalign(ALIGNCNT, sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pc.z = (PRECTYPE*)memalign(ALIGNCNT, sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));

	double tmutex=0.0, tqueue=0.0, ds, dt, dsq, dtq;
#if BIGQUEUE
	for(i=first;i<last;i++)
	{
#else
	int finish=0;
//	int mc=0;
	while(1)
	{
		ds=dtime();
		pthread_mutex_lock(mutex);
		finish=getGridIndex(gridP,curIndex,dp);
		pthread_mutex_unlock(mutex);
		tmutex+=dtime()-ds;
		if(finish)
			break;
		cnt=0;
		i=*curIndex;
#endif
		dsq=dtime();
		i1=i/(Jn*Kn);
		if(i1==0||i1==In-1)
			continue;	
		i2=(i-i1*Jn*Kn)/Kn;
		if(i2==0||i2==Jn-1)
			continue;	
		i3=i-i1*Jn*Kn-i2*Kn;
		if(i3==0||i3==Kn-1)
			continue;
		cella=i;
		if(dtp->numparticles[cella]<=0 || part[tree[dtp->root_cell[cella]].firstpart].key<dp->thisdom->LowerHilbertKey || part[tree[dtp->root_cell[cella]].firstpart].key>=dp->thisdom->UpperHilbertKey)
		{
//			if(dtp->numparticles[cella]>0)
//				printf("[%d]: %d, %d, %d, grid %ld ejected.\n", dp->rank, i, j, k, part[tree[dtp->root_cell[cella]].firstpart].key); 
			continue;
		}
		for (n=i1-1; n<i1+2; n++)
		{
			for (m=i2-1; m<i2+2; m++)
			{
				for (l=i3-1; l<i3+2; l++)
				{
					int cellb = n*(Jn*Kn)+m*Kn+l;
					if (dtp->numparticles[cella]>0 && dtp->numparticles[cellb]>0)
					{
//						if(dtp->root_cell[cella]==0)
//							printf("enqueue: i:%d, j:%d, k:%d, n:%d, m:%d, l:%d, ta:%d, tb:%d, na=%d, nb=%d\n", i, j, k, n, m, l, dtp->root_cell[cella], dtp->root_cell[cellb], tree[dtp->root_cell[cella]].nPart, tree[dtp->root_cell[cellb]].nPart);
						EnqueueP_Cell(PQ_treewalk, dtp->root_cell[cella], dtp->root_cell[cellb], open_angle);
						//								printf("E");
						cnt++;
					}
					else
					{
#if BIGQUEUE
						if(dtp->numparticles[cella]>0)
							printf("[%d]: (%d, %d, %d), (%d, %d, %d) grid %ld ejected.\n", dp->rank, i, j, k, n, m, l, part[tree[dtp->root_cell[cella]].firstpart].key); 
#else
						int ii,jj,kk;
						ii=cellb/Jn/Kn;
						jj=(cellb-ii*Jn*Kn)/Kn;
						kk=cellb-ii*Jn*Kn-jj*Kn;
						printf("i,j,k=(%d,%d,%d)\n",ii,jj,kk);
#endif
					}
				}
			}
		}
//					if (cnt==1) break;
//				}
//				if (cnt==1) break;
//			}
//			if (cnt==1) break;
//		}
//		if (cnt==1) break;
#if BIGQUEUE
	}
#else
		*numpart+=cnt;
#endif
//	printf("\ncnt = %d\n", cnt);

#if BIGQUEUE
	PQ_ppnode->tag=0;
	threaded(pth->nTeam, pth->nMaster, pth->nSlave, thread, pth->mutex, pthArr, pth->nTeam+pth->nSlave*teamid, pth->nTeam+pth->nSlave*(teamid+1)-1, P_PQ_ppnode, P_PQ_treewalk, lower, upper, dp,gp,compPP);

	ds=dtime();
	pthread_mutex_lock(&pth->mutex[teamid]);
#endif
	while(PQ_treewalk->length>0)
	{
		ProcessQP_Cell(PQ_treewalk, PQ_ppnode, dtt_process_cell);
	}
#if BIGQUEUE
	pthread_mutex_unlock(&pth->mutex[teamid]);
	tmutex+=dtime()-ds;
//	threaded(pth->nTeam, pth->nMaster, pth->nSlave, thread, pth->mutex, pthArr, pth->nTeam+pth->nSlave*teamid, pth->nTeam+pth->nSlave*(teamid+1)-1, P_PQ_ppnode, P_PQ_treewalk, lower, upper, dp,gp,compPP);
	PQ_ppnode->tag=1;
  	printf("\nQueue process finishing, total %d pp pairs.\n", PQ_ppnode->length);
#endif
	dtq=dtime();
	tqueue+=dtq-dsq;
	
//	printf("before ppkernel, %ld\n", (1<<(MAX_MORTON_LEVEL))*sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*18);
	//	sleep(5);
	double dpp = 0.0;
	while(PQ_ppnode->length>0)
	{
#if BIGQUEUE
		dpp += ProcessQP_PPnode(PQ_ppnode, ppkernel, pa, pb, pc, &pth->mutex[teamid]);
#else
		dpp += ProcessQP_PPnode(PQ_ppnode, ppkernel, pa, pb, pc);
#endif
	}

//	PPelement *pe=PQ_ppnode->elements;
//	n=PQ_ppnode->length/2;
//	for(i=0;i<n;i++)
//	{
//		int nA, nB, dq;
//
//		nA = tree[pe[i].TA].nPart;
//		packarray3(&part[tree[pe[i].TA].firstpart], nA, pa);
//		packarray3(NULL, nA, pc);
//
//		if (pe[i].mask == 0)
//		{
//			nB = tree[pe[i].TB].nPart;
//			packarray3(&part[tree[pe[i].TB].firstpart], nB, pb);
//		}
//		else
//		{
//			nB = 1;
//			pb.x[0] = tree[pe[i].TB].masscenter[0];
//			pb.y[0] = tree[pe[i].TB].masscenter[1];
//			pb.z[0] = tree[pe[i].TB].masscenter[2];
//	//		pb.x=pb.y=pb.z=NULL;	
//		}
//		ppkernel(pa, nA, pb, nB, EPS2, pc);
//
//		if (pe[i].mask == 1)
//		{	
//			pc.x[0] *= tree[pe[i].TB].mass;
//			pc.y[1] *= tree[pe[i].TB].mass;
//			pc.z[2] *= tree[pe[i].TB].mass;
//		}
//
//		pusharray3(&part[tree[pe[i].TA].firstpart], nA, pc);
//	}
	
//	printf("after ppkernel\n");	
//	printf("\nPP processing finished, time = %f.\n", dpp);
#if !BIGQUEUE
	}

//	printf("rank=%d\n",dp->rank);
//	printf("\nQueue process finishing, total %d pp pairs.\n", *numpart);
//	printf("thmchen= %d\n",mc);
#endif
	printf("mutex lock time = %f, queue time = %f\n",tmutex, tqueue);

	free(pa.x);
	free(pa.y);
	free(pa.z);
	free(pb.x);
	free(pb.y);
	free(pb.z);
	free(pc.x);
	free(pc.y);
	free(pc.z);
}


// Cao! Use this critertian for cell/cell acceptance judging. 
__OffloadFunc_Macro__
int accepted_cell_to_cell(int TA, int TB, double theta/* Here theta is the open_angle  */)
{
	double delta, dr;
	int d;

	//# Strict substraction of R for correct acceptance:
	//	double Rmax = SQROOT3 * (pTA->bmax + pTB->bmax); 
	//#	Alternative subtraction of R:
	double Rmax = tree[TB].width;
	double Bmax = tree[TA].width;
	//	double Rmax = tree[TA].width + tree[TB].width;
	//	double Bmax = 0.0;

	// Tree node include test.
	if (tree[TA].level < tree[TB].level)
	{
		if (tree[TA].firstpart<=tree[TB].firstpart && tree[TA].firstpart+tree[TA].nPart>=tree[TB].firstpart)
		{
			return 0;
		}
	}
	else if (tree[TA].level > tree[TB].level)
	{
		if (tree[TB].firstpart<=tree[TA].firstpart && tree[TB].firstpart+tree[TB].nPart>=tree[TA].firstpart)
		{
			return 0;
		}
	}

	dr = 0.0;

	//	printf("accept.\n");
	//	sleep(10);
	for (d=0; d<DIM; d++)
	{
		delta = tree[TB].masscenter[d] - tree[TA].masscenter[d];
		dr += delta*delta;
	}

	dr = SQRT(dr) - Bmax;
	if ( Rmax < theta * dr )
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

// Cao! Walk the tree using Dual Tree Traversel algorithm.
__OffloadFunc_Macro__
int dtt_process_cell(int TA, int TB, double theta, PPQ *PQ_ppnode, TWQ *PQ_treewalk)
{
	//	assert ((tree[TA].level >=0) && (tree[TA].level <= MAX_MORTON_LEVEL));
	//	assert ((tree[TB].level >=0) && (tree[TB].level <= MAX_MORTON_LEVEL));

	//	printf("process. TA=%d, TB=%d, theta=%.2f.\n", TA, TB, theta);

	if(TA == TB)
	{
		if (tree[TA].childnum == 0)
		{
			/* Cao! Here we send the tree node, and the Queue      */
			/*      processing function should pack the particles  */
			/*      and do some further optimizatoins.             */
			/*      PQ_ppnode: node pp packs, !0 means accepted    */
			/*      PQ_particle: save the packed pps.              */
			EnqueueP_PPnode(PQ_ppnode, TA, TB, 0);
		}
		else
		{
			int i;
			for (i=0; i<tree[TA].childnum; i++)
			{
				EnqueueP_Cell(PQ_treewalk, tree[TA].firstchild+i, TB, theta);
			}
		}
	}
	else
	{
		if (accepted_cell_to_cell(TA, TB, theta))
		{
			/* Cao! See comments above. */
			EnqueueP_PPnode(PQ_ppnode, TA, TB, 1);
		}
		else
		{
			if ( tree[TA].childnum==0 && tree[TB].childnum==0 )
			{
				/* Cao! See comments above. */
				EnqueueP_PPnode(PQ_ppnode, TA, TB, 0);
			}
			else if ( tree[TB].childnum==0 || (tree[TA].childnum>0 && tree[TA].width>tree[TB].width) )
				//			else if ( tree[TB].childnum==0 || tree[TA].childnum>0 ) // Open A until reach the leaf.
			{
				int i;
				for (i=0; i<tree[TA].childnum; i++)
				{
					EnqueueP_Cell(PQ_treewalk, tree[TA].firstchild+i, TB, theta);
				}
			}
			else
			{
				int i;
				for (i=0; i<tree[TB].childnum; i++)
				{
					EnqueueP_Cell(PQ_treewalk, TA, tree[TB].firstchild+i, theta);
				}
			}
		}
	}
}

void dtt_traversal(Domain *dp, GlobalParam *gp)
{
	int i, j, k, m, n, l;
	int In, Jn, Kn;
	int npart;
#if BIGQUEUE
	int nTeam_MIC = 240/dp->NumDom;
	int nMaster_MIC = 1;
	int nSlave_MIC = 1;
	int nTeam_CPU = 15;
	int nMaster_CPU = 1;
	int nSlave_CPU = 1;
#else
	int nTeam_MIC = 120*2/dp->NumDom;//actually the num of slave threads,240 threads each mic(2)
	int nSlave_MIC = nTeam_MIC;
	int nTeam_CPU = 16*2/dp->NumDom;
//	int nTeam_CPU = 15;//pure CPU like ivx
	int nSlave_CPU = nTeam_CPU;
#endif

	DomainTree *dtp = dp->domtree;
	npart = dp->NumPart;
	int nnode = dtp->NumNode;
	int cnt = 0;
	In = (dp->cuboid)->nSide[0];
	Jn = (dp->cuboid)->nSide[1];
	Kn = (dp->cuboid)->nSide[2];
	int Ngrid=In*Jn*Kn;
	part = dp->Part;
	tree = dp->domtree->tree;

	int myid=dp->rank;
	int micid=myid%2;
#ifdef __INTEL_OFFLOAD
	// Preparing for array offload.
	int *count = (dp->cuboid)->count;
	printf("micid %d on CPU master.\n", micid);
	int *off_rc = dtp->root_cell;
	int *off_pi = dtp->partidxes;
	int *off_nump = dtp->numparticles;
	int off_numtree = dtp->NumTree;
	int np = dp->NumDom;
	DomainInfo *off_dl = dp->DomList;
	SubCuboid *off_sc = dp->cuboid;
	int *off_tag = dp->cuboid->tag;

#if !BIGQUEUE
	double split=4.2;//change to 4.5
	static int stepping_num=-1;

	int splitgrid,splitpart;
	int* accum;	
	int start=Jn*Kn+Kn+1;
	int end=(In-2)*Jn*Kn+(Jn-2)*Kn+Kn-2;
	accum=(int*)malloc(sizeof(int)*(In-2)*(Jn-2)*(Kn-2));

	//	int mc=0;
	for (n=start,m=0; n<=end; n++) {
		if(dp->cuboid->tag[n] == TAG_OUTERCORE || dp->cuboid->tag[n] == TAG_FRONTIERS){
			//	    mc++;
			if(m>0)
				accum[m] = accum[m-1] + count[n];
			else
				accum[m] = count[n];
		}
		else
			continue;
		m++;
	}
	//	printf("mc= %d\n",mc);
	for (l=start,n=0; l<=end; l++) {
		if(dp->cuboid->tag[l] == TAG_OUTERCORE || dp->cuboid->tag[l] == TAG_FRONTIERS){
			if(((n>0&&(double)accum[n-1]<(double)accum[m-1]/split)||n==0)&&(double)accum[n]>=(double)accum[m-1]/split){
				splitgrid=l+1;
				splitpart=tree[dp->domtree->root_cell[splitgrid]].firstpart;
				break;
			}

		}
		else
			continue;
		n++;
	}
	stepping_num++;//the first time is 0,others are !0
	//				splitgrid=start;
	splitpart=tree[dp->domtree->root_cell[splitgrid]].firstpart;
#endif
	double dsc, dtc;
	dsc = dtime();

#if BIGQUEUE
#pragma offload target(mic:micid) \
	in (gp[0:1]) in (off_dl[0:np]) in (off_sc[0:1]) \
	in (part[0:npart]) in (tree[0:nnode]) out(part[0:npart])\
	in (count[0:Ngrid]) in (off_tag[0:Ngrid]) \
	in (off_rc[0:Ngrid]) in (off_pi[0:Ngrid]) in (off_nump[0:Ngrid])
#else
#pragma offload target(mic:micid) \
	in (gp[0:1]) in (off_dl[0:np]) in (off_sc[0:1]) \
	in (part[0:npart]) in (tree[0:nnode]) \
	in (count[0:Ngrid]) in (off_tag[0:Ngrid]) \
	in (off_rc[0:Ngrid]) in (off_pi[0:Ngrid]) in (off_nump[0:Ngrid]) out(part[splitpart:npart-splitpart]) signal(part)
#endif
	{
		double dso, dto;
		dso=dtime();
#if BIGQUEUE
		int tnum=nTeam_MIC * (nMaster_MIC + nSlave_MIC);
#else
		int tnum=nSlave_MIC;
#endif
		// Processing offload array pointers.
		DomainTree dt;
		Domain d;
		dt.root_cell = off_rc;
		dt.partidxes = off_pi;
		dt.numparticles = off_nump;
		dt.NumTree = off_numtree;
		dt.NumNode = nnode;
		d.rank = myid;
		d.NumDom = np;
		d.NumPart = npart;
		d.DomList = off_dl;
		d.thisdom = off_dl + myid;
		d.domtree = &dt;
		d.cuboid = off_sc;
		d.cuboid->tag = off_tag;
		d.cuboid->count = count;
		d.cuboid->mesh = NULL;
		d.Part = part;
		d.progress = 0;
		d.buff = NULL;

//		cpu_set_t cpuset;
//		CPU_ZERO(&cpuset);
//		int core_id=d.rank*tnum;
//		CPU_SET(core_id, &cpuset);
//		sched_setaffinity(0, sizeof(cpu_set_t), &cpuset);

		printf("micid %d on MIC master.\n", micid);
		TWQ *P_PQ_treewalk;
		PPQ *P_PQ_ppnode;

		Block* pth;
		pthread_t* thread;
		pthread_mutex_t *mutex;

		pth = (Block*)malloc(sizeof(Block)*tnum);
		thread = (pthread_t*)malloc(sizeof(pthread_t)*tnum);
		P_PQ_ppnode=(PPQ*)malloc(sizeof(PPQ)*nTeam_MIC);
		P_PQ_treewalk=(TWQ*)malloc(sizeof(TWQ)*nTeam_MIC);

		initGlobal(part, tree);
		printf("Init queue ppnode.\n");
		init_queue_pp(P_PQ_ppnode, nTeam_MIC);
		printf("Init queue treewalk.\n");
		init_queue_tw(P_PQ_treewalk, nTeam_MIC);

		printf("real domain In-2:%d, Jn-2:%d, Kn-2:%d\n", In-2, Jn-2, Kn-2);
#if BIGQUEUE
		int *accum,*lower,*upper,*numpart;
		accum=(int*)malloc(sizeof(int)*Ngrid);
		accum[0]=0;
		for (n=1; n<Ngrid; n++) {
			if(d.cuboid->tag[n] != TAG_FRONTIERS && d.cuboid->tag[n] != TAG_OUTERCORE)
			{
				accum[n] = accum[n-1];
			}
			else
			{
				accum[n] = accum[n-1] + count[n];
			}
		}
		int accum_thread;
		accum_thread=accum[Ngrid-1]/nTeam_MIC;
		lower = (int*)malloc(sizeof(int)*nTeam_MIC);
		upper = (int*)malloc(sizeof(int)*nTeam_MIC);
		numpart = (int*)malloc(sizeof(int)*nTeam_MIC);
		mutex = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t)*nTeam_MIC);

		lower[0] = 0;
		for (m=1, n=0; n<Ngrid; n++) {
			if ( accum[n] > m*accum_thread ) {
				lower[m] = upper[m-1] = n;
				m++;
			}
		}
		upper[nTeam_MIC-1]=Ngrid;

		numpart[0] = accum[upper[0]-1];
		for (m=1; m<nTeam_MIC; m++) {
			numpart[m] = accum[upper[m]-1] - accum[lower[m]-1];
			printf("[%d]: Lower[%d]=%d, Upper[%d]=%d.\n", d.rank, m, lower[m], m, upper[m]);
		}
		for (m=0; m<nTeam_MIC; m++)
		{
			pthread_mutex_init(&mutex[m], NULL);
		}
#else
		int *curIndex,*numpart;
		int *gridRange;
		gridRange = (int*)malloc(sizeof(int)*2);
		curIndex = (int*)malloc(sizeof(int)*nTeam_MIC);
		numpart = (int*)malloc(sizeof(int)*nTeam_MIC);
		mutex = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t));
		gridRange[0]=splitgrid;
		gridRange[1]=(In-2)*Jn*Kn+(Jn-2)*Kn+Kn-2;
		for (n=0; n<nTeam_MIC; n++)
		{
			curIndex[n]=gridRange[0];
			numpart[n]=0;
		}
		pthread_mutex_init(mutex, NULL);
#endif
		printf("before teamMaster\n");	
#if BIGQUEUE
		threaded(nTeam_MIC, nMaster_MIC, nSlave_MIC, thread, mutex, pth, 0, nTeam_MIC-1, P_PQ_ppnode, P_PQ_treewalk, lower, upper, &d, gp,teamMaster);
#else
		threaded(thread, mutex, pth, 0, nTeam_MIC-1, P_PQ_ppnode, P_PQ_treewalk, gridRange, curIndex, numpart, &d, gp,teamMaster);
#endif
		printf("after teamMaster\n");	
		for (i=0; i<tnum; i++)
		{
			pthread_join(thread[i], NULL);
//			printf("[%d]: thread %d joined.\n", d.rank, i);
		}	
		dto = dtime();
		printf("[%d], offload threaded part time = %f.\n", d.rank, dto - dso);
	}
#if BIGQUEUE
#else
	int tnum=nSlave_CPU;
	TWQ *P_PQ_treewalk;
	PPQ *P_PQ_ppnode;

	Block* pth;
	pthread_t* thread;
	pthread_mutex_t *mutex;

	pth = (Block*)malloc(sizeof(Block)*tnum);
	thread = (pthread_t*)malloc(sizeof(pthread_t)*tnum);
	mutex = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t));
	P_PQ_ppnode=(PPQ*)malloc(sizeof(PPQ)*nTeam_CPU);
	P_PQ_treewalk=(TWQ*)malloc(sizeof(TWQ)*nTeam_CPU);

	initGlobal(part, tree);
	printf("Init queue ppnode.\n");
	init_queue_pp(P_PQ_ppnode, nTeam_CPU);
	printf("Init queue treewalk.\n");
	init_queue_tw(P_PQ_treewalk, nTeam_CPU);

	printf("real domain In-2:%d, Jn-2:%d, Kn-2:%d\n", In-2, Jn-2, Kn-2);
	//////////////////////////////////	
	int *curIndex,*numpart;
	int *gridRange;
	gridRange = (int*)malloc(sizeof(int)*2);
	curIndex = (int*)malloc(sizeof(int)*nTeam_CPU);
	numpart = (int*)malloc(sizeof(int)*nTeam_CPU);
	gridRange[0]=Jn*Kn+Kn+1;
	gridRange[1]=splitgrid-1;
	for (n=0; n<nTeam_CPU; n++)
	{
		curIndex[n]=gridRange[0];
		numpart[n]=0;
	}
	///////////////////////////////////
	pthread_mutex_init(mutex, NULL);
	printf("before teamMaster\n");	
	threaded(thread, mutex, pth, 0, nTeam_CPU-1, P_PQ_ppnode, P_PQ_treewalk, gridRange,curIndex,numpart, dp,gp,teamMaster);
	printf("after teamMaster\n");	
	for (i=0; i<tnum; i++)
	{
		pthread_join(thread[i], NULL);
	}
#pragma offload_wait target(mic:micid) wait(part)
#endif
	dtc = dtime();
	printf("[%d], offload and transfer time = %f.\n", dp->rank, dtc - dsc);
#else
#if BIGQUEUE
	int tnum=nTeam_CPU * (nMaster_CPU + nSlave_CPU);
#else
	int tnum=nSlave_CPU;
#endif
	TWQ *P_PQ_treewalk;
	PPQ *P_PQ_ppnode;

	Block* pth;
	pthread_t* thread;
	pthread_mutex_t *mutex;

	pth = (Block*)malloc(sizeof(Block)*tnum);
	thread = (pthread_t*)malloc(sizeof(pthread_t)*tnum);
	P_PQ_ppnode=(PPQ*)malloc(sizeof(PPQ)*nTeam_CPU);
	P_PQ_treewalk=(TWQ*)malloc(sizeof(TWQ)*nTeam_CPU);

	initGlobal(part, tree);
	printf("Init queue ppnode.\n");
	init_queue_pp(P_PQ_ppnode, nTeam_CPU);
	printf("Init queue treewalk.\n");
	init_queue_tw(P_PQ_treewalk, nTeam_CPU);

	printf("real domain In-2:%d, Jn-2:%d, Kn-2:%d\n", In-2, Jn-2, Kn-2);

#if BIGQUEUE
	int *accum,*lower,*upper,*numpart;
	int accum_thread;
	int* count=(dp->cuboid)->count;
	accum_thread=npart/nTeam_CPU;
	accum=(int*)malloc(sizeof(int)*Ngrid);
	accum[0]=count[0];
	for (n=1; n<Ngrid; n++) {
		accum[n] = accum[n-1] + count[n];
	}
	lower = (int*)malloc(sizeof(int)*nTeam_CPU);
	upper = (int*)malloc(sizeof(int)*nTeam_CPU);
	numpart = (int*)malloc(sizeof(int)*nTeam_CPU);
	mutex = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t)*nTeam_CPU);

	lower[0] = 0;
	for (m=1, n=0; n<Ngrid; n++) {
		if ( accum[n] > m*accum_thread ) {
			lower[m] = upper[m-1] = n;
			m++;
		}
	}
	upper[nTeam_CPU-1]=Ngrid;

	numpart[0] = accum[upper[0]-1];
	for (m=1; m<nTeam_CPU; m++) {
		numpart[m] = accum[upper[m]-1] - accum[lower[m]-1];
	}
	for (m=0; m<nTeam_CPU; m++)
	{
		pthread_mutex_init(&mutex[m], NULL);
	}
#else
	int *curIndex,*numpart;
	int *gridRange;
	gridRange = (int*)malloc(sizeof(int)*2);
	curIndex = (int*)malloc(sizeof(int)*nTeam_CPU);
	numpart = (int*)malloc(sizeof(int)*nTeam_CPU);
	mutex = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t));
	gridRange[0]=Jn*Kn+Kn+1;
	gridRange[1]=(In-2)*Jn*Kn+(Jn-2)*Kn+Kn-2;
	for (n=0; n<nTeam_CPU; n++)
	{
		curIndex[n]=gridRange[0];
		numpart[n]=0;
	}
	pthread_mutex_init(mutex, NULL);
#endif
	printf("before teamMaster\n");	
#if BIGQUEUE
	threaded(nTeam_CPU, nMaster_CPU, nSlave_CPU, thread, mutex, pth, 0, nTeam_CPU-1, P_PQ_ppnode, P_PQ_treewalk, lower, upper, dp,gp,teamMaster);
#else
	threaded(thread, mutex, pth, 0, nTeam_CPU-1, P_PQ_ppnode, P_PQ_treewalk, gridRange,curIndex,numpart, dp,gp,teamMaster);
#endif
	printf("after teamMaster\n");	
	for (i=0; i<tnum; i++)
	{
		pthread_join(thread[i], NULL);
	}
#endif
//	for (i=0; i<Ngrid; i++)
//	{
//		if (dp->cuboid->tag[i] == TAG_FRONTIERS)
//		{
//			int pidx = tree[dp->domtree->root_cell[i]].firstpart;
//			printf("[%d] idx=%d, acc.x=%f, acc.y=%f, acc.z=%f.\n", dp->rank, i, part[pidx].acc[0], part[pidx].acc[1], part[pidx].acc[2]); 
//			break;
//		}
//	}
//	for (i=0; i<Ngrid; i++)
//	{
//		if (dp->cuboid->tag[i] == TAG_OUTERCORE)
//		{
//			int pidx = tree[dp->domtree->root_cell[i]].firstpart;
//			printf("[%d] idx=%d, acc.x=%f, acc.y=%f, acc.z=%f.\n", dp->rank, i, part[pidx].acc[0], part[pidx].acc[1], part[pidx].acc[2]); 
//			break;
//		}
//	}
//	for (i=Ngrid-1; i>=0; i--)
//	{
//		if (dp->cuboid->tag[i] == TAG_FRONTIERS)
//		{
//			int pidx = tree[dp->domtree->root_cell[i]].firstpart;
//			printf("[%d] idx=%d, acc.x=%f, acc.y=%f, acc.z=%f.\n", dp->rank, i, part[pidx].acc[0], part[pidx].acc[1], part[pidx].acc[2]); 
//			break;
//		}
//	}
//	for (i=Ngrid-1; i>=0; i--)
//	{
//		if (dp->cuboid->tag[i] == TAG_OUTERCORE)
//		{
//			int pidx = tree[dp->domtree->root_cell[i]].firstpart;
//			printf("[%d] idx=%d, acc.x=%f, acc.y=%f, acc.z=%f.\n", dp->rank, i, part[pidx].acc[0], part[pidx].acc[1], part[pidx].acc[2]); 
//			break;
//		}
//	}
//	if (dp->rank==12)
//	{
//		FILE *fp;
//		fp = fopen("res12", "w+");
//
//		for (i=0; i<Ngrid; i++)
//		{
//			if (dp->cuboid->tag[i] == TAG_FRONTIERS || dp->cuboid->tag[i] == TAG_OUTERCORE)
//			{
//				int pidx = tree[dp->domtree->root_cell[i]].firstpart;
//				int pnum = dp->domtree->numparticles[i];
//				for (j=0; j<pnum; j++)
//				{
//					fprintf(fp, "[%d] idx=%d, px=%f, py=%f, pz=%f, ax=%f, ay=%f, az=%f.\n", dp->rank, i, part[pidx+j].pos[0], part[pidx+j].pos[1], part[pidx+j].pos[2], part[pidx+j].acc[0], part[pidx+j].acc[1], part[pidx+j].acc[2]); 
//				}
//			}
//		}
//		fclose(fp);
//	}
}
