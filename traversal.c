#include "traversal.h"
#include <time.h>
#include <math.h>
#include <pthread.h>

__OffloadVar_Macro__
static Body *part;

__OffloadVar_Macro__
static Node *tree; // thread number

	__OffloadFunc_Macro__
int threaded(int nTeam, int nMaster, int nSlave, pthread_t* thread, pthread_mutex_t *mutex, Block* pth, int ts, int te, PPQ *pq, TWQ *tq, int* lower,int* upper, Domain*dp,GlobalParam*gp,void* process(void*))
{
	int i,teamid;
	printf("before threaded\n");
	pthread_attr_t attr;

	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	//	sleep(2);	
	for (i=ts; i<te+1; i++)
	{
		pth[i].blockid = i/nTeam;
		pth[i].teamid=i%nTeam;
		teamid=pth[i].teamid;
		pth[i].nTeam = nTeam;
		pth[i].nMaster = nMaster;
		pth[i].nSlave = nSlave;
		pth[i].bsize = nMaster+nSlave;
		pth[i].tid = i;
		pth[i].tall = nTeam*(nMaster+nSlave);
		pth[i].PQ_ppnode = &pq[teamid];
		pth[i].PQ_treewalk = &tq[teamid];
		pth[i].first=lower[teamid];
		pth[i].last=upper[teamid];

		pth[i].dp=dp;
		pth[i].gp=gp;
		pth[i].thread=thread;
		pth[i].pthArr=pth;
		pth[i].lower=lower;
		pth[i].upper=upper;
		pth[i].P_PQ_ppnode=pq;
		pth[i].P_PQ_treewalk=tq;
		pth[i].mutex=mutex;	
		printf("[%d/%d] blockid=%d, teamid=%d, lower=%d, upper=%d, Ngrid=%d.\n", dp->rank, i, pth[i].blockid, pth[i].teamid, pth[i].first, pth[i].last, dp->cuboid->nSide[0]*dp->cuboid->nSide[1]*dp->cuboid->nSide[2]);
		pthread_create(&thread[i], &attr, process, (void*)&pth[i]);
	}
	pthread_attr_destroy(&attr);
}

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
	pa.x = (PRECTYPE*)malloc(sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pa.y = (PRECTYPE*)malloc(sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pa.z = (PRECTYPE*)malloc(sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pb.x = (PRECTYPE*)malloc(sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pb.y = (PRECTYPE*)malloc(sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pb.z = (PRECTYPE*)malloc(sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pc.x = (PRECTYPE*)malloc(sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pc.y = (PRECTYPE*)malloc(sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pc.z = (PRECTYPE*)malloc(sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	
	printf("before ppkernel, %ld\n", (1<<(MAX_MORTON_LEVEL))*sizeof(PRECTYPE)*MAX_PACKAGE_SIZE);
	sleep(2);
	int end=0;	
	while(PQ_ppnode->length>0||PQ_ppnode->tag==0)
	{
		if(PQ_ppnode->length>0)
		{
			ProcessQP_PPnode(PQ_ppnode, ppkernel, pa, pb, pc, &pth->mutex[teamid]);
		}
	}
	printf("after ppkernel\n");	
	printf("\nPP processing finished.\n");

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

__OffloadFunc_Macro__
void* teamMaster(void* param)
{
	Array3 pa, pb, pc;
	int i,j,k,n,m,l,i1,i2,i3,cella;
	Block* pth=(Block*)param;
	int blockid=pth->blockid;
	int teamid=pth->teamid;
	int bsize=pth->bsize;
	int tid=pth->tid;
	int tall=pth->tall;
	PPQ* PQ_ppnode=pth->PQ_ppnode;
	TWQ* PQ_treewalk=pth->PQ_treewalk;
	int first=pth->first;
	int last=pth->last;
	Domain* dp=pth->dp;
	GlobalParam* gp=pth->gp;
	pthread_t* thread=pth->thread;
	Block* pthArr=pth->pthArr;
	int* lower=pth->lower;
	int* upper=pth->upper;
	PPQ* P_PQ_ppnode=pth->P_PQ_ppnode;
	TWQ* P_PQ_treewalk=pth->P_PQ_treewalk;

	int In = dp->cuboid->nSide[0];
	int Jn = dp->cuboid->nSide[1];
	int Kn = dp->cuboid->nSide[2];
	int cnt = 0;
	double open_angle = 0.3;

	DomainTree *dtp = dp->domtree;
	int npart = dp->NumPart;
	//    part = dp->Part;
	//    tree = dtp->tree;
	printf("[%d/%d] first=%d, last=%d, Ngrid=%d, part0key=%ld.\n",dp->rank, tid, first, last, In*Jn*Kn, part[0].pos[0]);
	for(i=first;i<last;i++)
	{
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
		for (n=i1-1; n<i1+2; n++)
		{
			for (m=i2-1; m<i2+2; m++)
			{
				for (l=i3-1; l<i3+2; l++)
				{
					int cellb = n*(Jn*Kn)+m*Kn+l;
					if (dtp->numparticles[cella]>0 && dtp->numparticles[cellb]>0)
					{
						if(dtp->root_cell[cella]==0)
							printf("enqueue: i:%d, j:%d, k:%d, n:%d, m:%d, l:%d, ta:%d, tb:%d, na=%d, nb=%d\n", i, j, k, n, m, l, dtp->root_cell[cella], dtp->root_cell[cellb], tree[dtp->root_cell[cella]].nPart, tree[dtp->root_cell[cellb]].nPart);
						EnqueueP_Cell(PQ_treewalk, dtp->root_cell[cella], dtp->root_cell[cellb], open_angle);
						//								printf("E");
						cnt++;
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
	}
	printf("\ncnt = %d\n", cnt);

	PQ_ppnode->tag=0;
	threaded(pth->nTeam, pth->nMaster, pth->nSlave, thread, pth->mutex, pthArr, pth->nTeam+pth->nSlave*teamid, pth->nTeam+pth->nSlave*(teamid+1)-1, P_PQ_ppnode, P_PQ_treewalk, lower, upper, dp,gp,compPP);

	pthread_mutex_lock(&pth->mutex[teamid]);
	while(PQ_treewalk->length>0)
	{
		ProcessQP_Cell(PQ_treewalk, PQ_ppnode, dtt_process_cell);
	}
	pthread_mutex_unlock(&pth->mutex[teamid]);
	PQ_ppnode->tag=1;
	printf("\nQueue process finishing, total %d pp pairs.\n", PQ_ppnode->length);

	pa.x = (PRECTYPE*)malloc(sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pa.y = (PRECTYPE*)malloc(sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pa.z = (PRECTYPE*)malloc(sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pb.x = (PRECTYPE*)malloc(sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pb.y = (PRECTYPE*)malloc(sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pb.z = (PRECTYPE*)malloc(sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pc.x = (PRECTYPE*)malloc(sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pc.y = (PRECTYPE*)malloc(sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	pc.z = (PRECTYPE*)malloc(sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*(1<<(MAX_MORTON_LEVEL*2)));
	
	printf("before ppkernel, %ld\n", (1<<(MAX_MORTON_LEVEL))*sizeof(PRECTYPE)*MAX_PACKAGE_SIZE*18);
	//	sleep(5);
	while(PQ_ppnode->length>0)
	{
		ProcessQP_PPnode(PQ_ppnode, ppkernel, pa, pb, pc, &pth->mutex[teamid]);
	}
	printf("after ppkernel\n");	
	printf("\nPP processing finished.\n");

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

	int nTeam_MIC = 60;
	int nMaster_MIC = 1;
	int nSlave_MIC = 1;
	int nTeam_CPU = 16/dp->NumDom;
	int nMaster_CPU = 1;
	int nSlave_CPU = 1;

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

#pragma offload target(mic:micid) \
	in (gp[0:1]) in (off_dl[0:np]) in (off_sc[0:1]) \
	in (part[0:npart]) in (tree[0:nnode]) \
	in (count[0:Ngrid]) in (off_tag[0:Ngrid]) \
	in (off_rc[0:Ngrid]) in (off_pi[0:Ngrid]) in (off_nump[0:Ngrid])
	{
		int tnum=nTeam_MIC * (nMaster_MIC + nSlave_MIC);
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

		printf("micid %d on MIC master.\n", micid);
		TWQ *P_PQ_treewalk;
		PPQ *P_PQ_ppnode;

		Block* pth;
		pthread_t* thread;
		pthread_mutex_t *mutex;

		pth = (Block*)malloc(sizeof(Block)*tnum);
		thread = (pthread_t*)malloc(sizeof(pthread_t)*tnum);
		mutex = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t)*nTeam_MIC);
		P_PQ_ppnode=(PPQ*)malloc(sizeof(PPQ)*nTeam_MIC);
		P_PQ_treewalk=(TWQ*)malloc(sizeof(TWQ)*nTeam_MIC);

		initGlobal(part, tree);
		printf("Init queue ppnode.\n");
		init_queue_pp(P_PQ_ppnode, nTeam_MIC);
		printf("Init queue treewalk.\n");
		init_queue_tw(P_PQ_treewalk, nTeam_MIC);

		printf("real domain In-2:%d, Jn-2:%d, Kn-2:%d\n", In-2, Jn-2, Kn-2);
		int *accum,*lower,*upper,*numpart;
		int accum_thread;
		accum_thread=npart/nTeam_MIC;
		accum=(int*)malloc(sizeof(int)*Ngrid);
		accum[0]=count[0];
		for (n=1; n<Ngrid; n++) {
			accum[n] = accum[n-1] + count[n];
		}
		lower = (int*)malloc(sizeof(int)*nTeam_MIC);
		upper = (int*)malloc(sizeof(int)*nTeam_MIC);
		numpart = (int*)malloc(sizeof(int)*nTeam_MIC);

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
		printf("before teamMaster\n");	
		threaded(nTeam_MIC, nMaster_MIC, nSlave_MIC, thread, mutex, pth, 0, nTeam_MIC-1, P_PQ_ppnode, P_PQ_treewalk, lower, upper, &d, gp,teamMaster);
		printf("after teamMaster\n");	
		for (i=0; i<tnum; i++)
		{
			pthread_join(thread[i], NULL);
//			printf("[%d]: thread %d joined.\n", d.rank, i);
		}	
	}
#else
	int tnum=nTeam_CPU * (nMaster_CPU + nSlave_CPU);
	TWQ *P_PQ_treewalk;
	PPQ *P_PQ_ppnode;

	Block* pth;
	pthread_t* thread;
	pthread_mutex_t *mutex;

	pth = (Block*)malloc(sizeof(Block)*tnum);
	thread = (pthread_t*)malloc(sizeof(pthread_t)*tnum);
	mutex = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t)*nTeam_CPU);
	P_PQ_ppnode=(PPQ*)malloc(sizeof(PPQ)*nTeam_CPU);
	P_PQ_treewalk=(TWQ*)malloc(sizeof(TWQ)*nTeam_CPU);

	initGlobal(part, tree);
	printf("Init queue ppnode.\n");
	init_queue_pp(P_PQ_ppnode, nTeam_CPU);
	printf("Init queue treewalk.\n");
	init_queue_tw(P_PQ_treewalk, nTeam_CPU);

	printf("real domain In-2:%d, Jn-2:%d, Kn-2:%d\n", In-2, Jn-2, Kn-2);
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
	printf("before teamMaster\n");	
	threaded(nTeam_CPU, nMaster_CPU, nSlave_CPU, thread, mutex, pth, 0, nTeam_CPU-1, P_PQ_ppnode, P_PQ_treewalk, lower, upper, dp,gp,teamMaster);
	printf("after teamMaster\n");	
	for (i=0; i<tnum; i++)
	{
		pthread_join(thread[i], NULL);
	}
#endif
}
