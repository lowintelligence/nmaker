#include "traversal.h"
#include <time.h>
#include <math.h>

__OffloadVar_Macro__
static Body *part;

__OffloadVar_Macro__
static Node *tree; // thread number

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

__OffloadFunc_Macro__
int threaded(int nTeam, int nSlave, pthread_t* thread, pthread_mutex_t *mutex, pthread_barrier_t *bar, PPParameter *pppar, Block* pth, int ts, int te, PPQ *pq, TWQ *tq, int* gridP,int* curIndex,int*numpart, Domain*dp,GlobalParam*gp,void* process(void*))
{
	int i,teamid,core_id=0;
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
		pth[i].blockid=i/nTeam;
		pth[i].teamid=i%nTeam;
		pth[i].nTeam = nTeam;
		pth[i].nSlave = nSlave;

		teamid=pth[i].teamid;
		pth[i].tid = i;
		pth[i].PQ_ppnode = &pq[teamid];
		pth[i].PQ_treewalk = &tq[teamid];

		pth[i].curIndex=&curIndex[teamid];
		pth[i].numpart=&numpart[teamid];
		pth[i].gridP=gridP;
		pth[i].bar=&bar[teamid];
		pth[i].pppar=&pppar[teamid];

		pth[i].dp=dp;
		pth[i].gp=gp;
		pth[i].thread=thread;
		pth[i].pthArr=pth;
		pth[i].P_PQ_ppnode=pq;
		pth[i].P_PQ_treewalk=tq;
		pth[i].mutex=mutex;

#ifdef __MIC__
		CPU_ZERO(&cpuset);
		core_id=pth[i].blockid+teamid*nSlave+dp->rank/2*nTeam*nSlave+1;
//		int core_id=((i+dp->rank/2*pth[i].tall)%60)*4+((i+dp->rank/2*pth[i].tall)/60)+1;
		CPU_SET(core_id, &cpuset);
		pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpuset);
#endif
//		printf("[%d/%d] blockid=%d, nTeam=%d, teamid=%d, Ngrid=%d, core=%d.\n", dp->rank, i, pth[i].blockid, pth[i].nTeam, pth[i].teamid, dp->cuboid->nSide[0]*dp->cuboid->nSide[1]*dp->cuboid->nSide[2], core_id);
		pthread_create(&thread[i], &attr, process, (void*)&pth[i]);
	}
	pthread_attr_destroy(&attr);
}

__OffloadFunc_Macro__
void* compPP(void* param)
{
	int n;
	long flop=0;
	Block* pth=(Block*)param;
	int teamid=pth->teamid;

	pthread_barrier_wait(pth->bar);
	while(pth->pppar->finish==0)
	{
		if(pth->pppar->calcm)
		{
//			flop += pth->pppar->nA*pth->pppar->nB*23;
			ppmkernel(pth->pppar->pa, pth->pppar->nA, pth->pppar->pb, pth->pppar->mass, pth->pppar->nB, EPS2, pth->pppar->pc, pth->nSlave, pth->blockid);
		}
		else
		{
//			flop += pth->pppar->nA*pth->pppar->nB*22;
			ppkernel(pth->pppar->pa, pth->pppar->nA, pth->pppar->pb, pth->pppar->nB, EPS2, pth->pppar->pc, pth->nSlave, pth->blockid);
		}
		pthread_barrier_wait(pth->bar);
		pthread_barrier_wait(pth->bar);
	}
//	printf("[%d/%d] Flops = %ld.\n", pth->dp->rank, pth->teamid, flop);
}

__OffloadFunc_Macro__
void* teamMaster(void* param)
{
	Array3 pa, pb, pc;
	PRECTYPE *mass;
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

	int* curIndex=pth->curIndex;
	int* gridP=pth->gridP;
	int* numpart=pth->numpart;
	pthread_mutex_t* mutex=pth->mutex;

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
	mass = (PRECTYPE*)memalign(ALIGNCNT, sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));

	double tmutex=0.0, tqueue=0.0, ds, dt, dsq, dtq, dpp;

	pth->pppar->pa=pa;
	pth->pppar->pb=pb;
	pth->pppar->pc=pc;
	pth->pppar->mass=mass;
	pth->pppar->finish=0;
	int finish=0;
	dpp = 0.0;

//	int mc=0;
	while(1)
	{
		ds=dtime();
		pthread_mutex_lock(mutex);
		finish=getGridIndex(gridP,curIndex,dp);
		pthread_mutex_unlock(mutex);
		tmutex+=dtime()-ds;
		if(finish)
		{
			pth->pppar->finish=1;
			break;
		}
		cnt=0;
		i=*curIndex;

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
						int ii,jj,kk;
						ii=cellb/Jn/Kn;
						jj=(cellb-ii*Jn*Kn)/Kn;
						kk=cellb-ii*Jn*Kn-jj*Kn;
						printf("i,j,k=(%d,%d,%d)\n",ii,jj,kk);
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

		*numpart+=cnt;

//		printf("\ncnt = %d\n", cnt);

		while(PQ_treewalk->length>0)
		{
			ProcessQP_Cell(PQ_treewalk, pth, dtt_process_cell);
		}

		dtq=dtime();
		tqueue+=dtq-dsq;
		
//		printf("rank=%d\n",dp->rank);
//		printf("\nQueue process finishing, total %d pp pairs.\n", *numpart);
//		printf("thmchen= %d\n",mc);
	}
	printf("mutex lock time = %f, queue time = %f, pp memory time = %f \n",tmutex, tqueue, dpp);

	free(pa.x);
	free(pa.y);
	free(pa.z);
	free(pb.x);
	free(pb.y);
	free(pb.z);
	free(pc.x);
	free(pc.y);
	free(pc.z);
	free(mass);
	pthread_barrier_wait(pth->bar);
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
int dtt_process_cell(int TA, int TB, double theta, TWQ *PQ_treewalk, Block *pth)
{
	//	assert ((tree[TA].level >=0) && (tree[TA].level <= MAX_MORTON_LEVEL));
	//	assert ((tree[TB].level >=0) && (tree[TB].level <= MAX_MORTON_LEVEL));

	//	printf("process. TA=%d, TB=%d, theta=%.2f.\n", TA, TB, theta);
	Array3 pa = pth->pppar->pa;
	Array3 pb = pth->pppar->pb;
	Array3 pc = pth->pppar->pc;
	PRECTYPE *mass = pth->pppar->mass;

	int nA, nB, calcm;
	if(TA == TB)
	{
		if (tree[TA].childnum == 0)
		{
			/* Cao! Here we send the tree node, and the Queue      */
			/*      processing function should pack the particles  */
			/*      and do some further optimizatoins.             */
			/*      PQ_ppnode: node pp packs, !0 means accepted    */
			/*      PQ_particle: save the packed pps.              */
			nA = tree[TA].nPart;
			nB = nA;
			packarray3(&part[tree[TA].firstpart], nA, pa);
			packarray3(NULL, nA, pc);
			pth->pppar->nA = nA;
			pth->pppar->nB = nB;
			pth->pppar->calcm = 0;
			pthread_barrier_wait(pth->bar);
			ppkernel(pa, nA, pb, nB, EPS2, pc, pth->nSlave, pth->blockid);
			pthread_barrier_wait(pth->bar);
			pusharray3(&part[tree[TA].firstpart], nA, pc);
		}
		else
		{
			int i, j;
			for (i=0; i<tree[TA].childnum; i++)
			{
				int TA1 = tree[TA].firstchild+i;
				if (tree[TA1].childnum == 0)
				{
					calcm = 0;
					nA = tree[TA1].nPart;
					nB = 0;
					packarray3m(&part[tree[TA].firstpart], nA, pa, mass);
					packarray3(NULL, nA, pc);
					for (j=0; j<tree[TB].childnum; j++)
					{
						int TB1 = tree[TB].firstchild+j;
						if (accepted_cell_to_cell(TA1, TB1, theta))
						{
							calcm = 1;
							pb.x[nB] = tree[TB1].masscenter[0];
							pb.y[nB] = tree[TB1].masscenter[1];
							pb.z[nB] = tree[TB1].masscenter[2];
							mass[nB] = tree[TB1].mass;
							nB += 1;
						}
						else
						{
							if(tree[TB1].childnum == 0)
							{
								packarray3om(&part[tree[TB1].firstpart], nB, tree[TB1].nPart, pb, mass);
								nB += tree[TB1].nPart;
							}
							else
							{
								EnqueueP_Cell(PQ_treewalk, TA1, TB1, theta);
							}
						}
					}
					if (nB>0)
					{
						pth->pppar->nA = nA;
						pth->pppar->nB = nB;
						pth->pppar->calcm = calcm;
						pthread_barrier_wait(pth->bar);
						if (calcm == 0)
						{
							ppkernel(pa, nA, pb, nB, EPS2, pc, pth->nSlave, pth->blockid);
						}
						else
						{	
							ppmkernel(pa, nA, pb, mass, nB, EPS2, pc, pth->nSlave, pth->blockid);
						}
						pthread_barrier_wait(pth->bar);
						pusharray3(&part[tree[TA1].firstpart], nA, pc);
					}
				}
				else
				{
					EnqueueP_Cell(PQ_treewalk, TA1, TB, theta);
				}
			}
		}
	}
	else
	{
		calcm = 0;
		nA = tree[TA].nPart;
		nB = 0;
		packarray3m(&part[tree[TA].firstpart], nA, pa, mass);
		packarray3(NULL, nA, pc);
		if (accepted_cell_to_cell(TA, TB, theta))
		{
			calcm = 1;
			pb.x[nB] = tree[TB].masscenter[0];
			pb.y[nB] = tree[TB].masscenter[1];
			pb.z[nB] = tree[TB].masscenter[2];
			mass[nB] = tree[TB].mass;
			nB = 1;
		}
		else
		{
			if ( tree[TA].childnum==0 && tree[TB].childnum==0 )
			{
				packarray3(&part[tree[TB].firstpart], nB, pb);
				nB = tree[TB].nPart;
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
					int TB1 = tree[TB].firstchild+i;
					if (accepted_cell_to_cell(TA, TB1, theta))
					{
						calcm = 1;
						pb.x[nB] = tree[TB1].masscenter[0];
						pb.y[nB] = tree[TB1].masscenter[1];
						pb.z[nB] = tree[TB1].masscenter[2];
						mass[nB] = tree[TB1].mass;
						nB += 1;
					}
					else
					{
						if(tree[TB1].childnum == 0)
						{
							packarray3om(&part[tree[TB1].firstpart], nB, tree[TB1].nPart, pb, mass);
							nB += tree[TB1].nPart;
						}
						else
						{
							EnqueueP_Cell(PQ_treewalk, TA, TB1, theta);
						}
					}
				}
			}
		}
		if (nB>0)
		{
			pth->pppar->nA = nA;
			pth->pppar->nB = nB;
			pth->pppar->calcm = calcm;
			pthread_barrier_wait(pth->bar);
			if (calcm == 0)
			{
				ppkernel(pa, nA, pb, nB, EPS2, pc, pth->nSlave, pth->blockid);
			}
			else
			{	
				ppmkernel(pa, nA, pb, mass, nB, EPS2, pc, pth->nSlave, pth->blockid);
			}
			pthread_barrier_wait(pth->bar);
			pusharray3(&part[tree[TA].firstpart], nA, pc);
		}
	}
}

void dtt_traversal(Domain *dp, GlobalParam *gp)
{
	int i, j, k, m, n, l;
	int In, Jn, Kn;
	int npart;

	int nTeam_MIC = 120/dp->NumDom;//actually the num of slave threads,240 threads each mic(2)
	int nSlave_MIC = 4;
//	int nTeam_CPU = 16/dp->NumDom;
	int nTeam_CPU = 15;//pure CPU like ivx
	int nSlave_CPU = 2;

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

	double split=7.1;//change to 4.5
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

	double dsc, dtc;
	dsc = dtime();

#pragma offload target(mic:micid) \
	in (gp[0:1]) in (off_dl[0:np]) in (off_sc[0:1]) \
	in (part[0:npart]) in (tree[0:nnode]) \
	in (count[0:Ngrid]) in (off_tag[0:Ngrid]) \
	in (off_rc[0:Ngrid]) in (off_pi[0:Ngrid]) in (off_nump[0:Ngrid]) out(part[splitpart:npart-splitpart]) signal(part)
	{
		double dso, dto;
		dso=dtime();

		int tnum=nTeam_MIC * nSlave_MIC;

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
		pthread_barrier_t *bar;
		PPParameter *pppar;

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

		int *curIndex,*numpart;
		int *gridRange;
		gridRange = (int*)malloc(sizeof(int)*2);
		curIndex = (int*)malloc(sizeof(int)*nTeam_MIC);
		numpart = (int*)malloc(sizeof(int)*nTeam_MIC);
		mutex = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t));
		bar = (pthread_barrier_t*)malloc(sizeof(pthread_barrier_t)*nTeam_MIC);
		pppar = (PPParameter*)malloc(sizeof(PPParameter)*nTeam_MIC);

		gridRange[0]=splitgrid;
		gridRange[1]=(In-2)*Jn*Kn+(Jn-2)*Kn+Kn-2;
		for (n=0; n<nTeam_MIC; n++)
		{
			curIndex[n]=gridRange[0];
			numpart[n]=0;
			pthread_barrier_init(&bar[n], NULL, nSlave_MIC);
			pppar[n].finish=0;
		}
		pthread_mutex_init(mutex, NULL);

		printf("before teamMaster\n");	

		threaded(nTeam_MIC, nSlave_MIC, thread, mutex, bar, pppar, pth, 0, nTeam_MIC-1, P_PQ_ppnode, P_PQ_treewalk, gridRange, curIndex, numpart, &d, gp,teamMaster);
		if (nSlave_MIC > 1)
			threaded(nTeam_MIC, nSlave_MIC, thread, mutex, bar, pppar, pth, nTeam_MIC, tnum-1, P_PQ_ppnode, P_PQ_treewalk, gridRange, curIndex, numpart, &d, gp,compPP);

		printf("after teamMaster\n");	
		for (i=0; i<tnum; i++)
		{
			pthread_join(thread[i], NULL);
//			printf("[%d]: thread %d joined.\n", d.rank, i);
		}	
		dto = dtime();
		printf("[%d], offload threaded part time = %f.\n", d.rank, dto - dso);
	}

	int tnum=nTeam_CPU * nSlave_CPU;

	TWQ *P_PQ_treewalk;
	PPQ *P_PQ_ppnode;

	Block* pth;
	pthread_t* thread;
	pthread_mutex_t *mutex;
	pthread_barrier_t *bar;
	PPParameter *pppar;

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
	bar = (pthread_barrier_t*)malloc(sizeof(pthread_barrier_t)*nTeam_CPU);
	pppar = (PPParameter*)malloc(sizeof(PPParameter)*nTeam_CPU);
	gridRange[0]=Jn*Kn+Kn+1;
	gridRange[1]=splitgrid-1;
	for (n=0; n<nTeam_CPU; n++)
	{
		curIndex[n]=gridRange[0];
		numpart[n]=0;
		pthread_barrier_init(&bar[n], NULL, nSlave_CPU);
		pppar[n].finish=0;
	}
	///////////////////////////////////
	pthread_mutex_init(mutex, NULL);
	printf("before teamMaster\n");	
	threaded(nTeam_CPU, nSlave_CPU, thread, mutex, bar, pppar, pth, 0, nTeam_CPU-1, P_PQ_ppnode, P_PQ_treewalk, gridRange, curIndex, numpart, dp, gp,teamMaster);
	if (nSlave_CPU > 1)
		threaded(nTeam_CPU, nSlave_CPU, thread, mutex, bar, pppar, pth, nTeam_CPU, tnum-1, P_PQ_ppnode, P_PQ_treewalk, gridRange, curIndex, numpart, dp, gp,compPP);
	printf("after teamMaster\n");	
	for (i=0; i<tnum; i++)
	{
		pthread_join(thread[i], NULL);
	}
#pragma offload_wait target(mic:micid) wait(part)
	dtc = dtime();
	printf("[%d], offload and transfer time = %f.\n", dp->rank, dtc - dsc);
#else // No offload
	int tnum=nTeam_CPU * nSlave_CPU;

	TWQ *P_PQ_treewalk;
	PPQ *P_PQ_ppnode;

	Block* pth;
	pthread_t* thread;
	pthread_mutex_t *mutex;
	pthread_barrier_t *bar;
	PPParameter *pppar;

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

	int *curIndex,*numpart;
	int *gridRange;
	gridRange = (int*)malloc(sizeof(int)*2);
	curIndex = (int*)malloc(sizeof(int)*nTeam_CPU);
	numpart = (int*)malloc(sizeof(int)*nTeam_CPU);
	mutex = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t));
	bar = (pthread_barrier_t*)malloc(sizeof(pthread_barrier_t)*nTeam_CPU);
	pppar = (PPParameter*)malloc(sizeof(PPParameter)*nTeam_CPU);
	gridRange[0]=Jn*Kn+Kn+1;
	gridRange[1]=(In-2)*Jn*Kn+(Jn-2)*Kn+Kn-2;
	for (n=0; n<nTeam_CPU; n++)
	{
		curIndex[n]=gridRange[0];
		numpart[n]=0;
		pthread_barrier_init(&bar[n], NULL, nSlave_CPU);
		pppar[n].finish=0;
	}
	pthread_mutex_init(mutex, NULL);

	printf("before teamMaster\n");	

	threaded(nTeam_CPU, nSlave_CPU, thread, mutex, bar, pppar, pth, 0, nTeam_CPU-1, P_PQ_ppnode, P_PQ_treewalk, gridRange, curIndex, numpart, dp, gp,teamMaster);
	if (nSlave_CPU > 1)
		threaded(nTeam_CPU, nSlave_CPU, thread, mutex, bar, pppar, pth, nTeam_CPU, tnum-1, P_PQ_ppnode, P_PQ_treewalk, gridRange, curIndex, numpart, dp, gp,compPP);

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
