#define _GNU_SOURCE
#include <sched.h>
#include <pthread.h>
#include <math.h>

#include "global.h"
#include "domain.h"
#include "queue.h"
#include "ppkernel.h"
#include "traversal.h"

__OffloadVar_Macro__
static CalcBody *part;

__OffloadVar_Macro__
static Body *fpart;

__OffloadVar_Macro__
static Node *tree; // thread number

__OffloadFunc_Macro__
int getGridIndex(int *gridP, int *curIndex, Domain *dp)
{
	*curIndex = gridP[0];
	gridP[0]++;
	while(*curIndex<gridP[1] && (dp->tag[*curIndex]!=TAG_OUTERCORE && dp->tag[*curIndex]!=TAG_FRONTIERS))
	{
		*curIndex = gridP[0];
		gridP[0]++;
	}
	if(*curIndex >= gridP[1])
		return 1;
	else
		return 0;
}

__OffloadFunc_Macro__
void* compPP(void* param)
{
	DttBlock* pth=(DttBlock*)param;

	pthread_barrier_wait(pth->bar);
	while(pth->pppar->finish==0)
	{
		if(pth->pppar->calcm)
		{
			ppmkernel(pth->pppar->pa, pth->pppar->nA, pth->pppar->pb, pth->pppar->mass, pth->pppar->nB, pth->constants, pth->pppar->pc, pth->nSlave, pth->blockid);
		}
		else
		{
			ppkernel(pth->pppar->pa, pth->pppar->nA, pth->pppar->pb, pth->pppar->nB, pth->constants, pth->pppar->pc, pth->nSlave, pth->blockid);
		}
		pthread_barrier_wait(pth->bar);
		pthread_barrier_wait(pth->bar);
	}
}

__OffloadFunc_Macro__
void* teamMaster(void* param)
{
	Array3 pa, pb, pc;
	Real *mass;

	int n, m, l, i1, i2, i3, cella;
	DttBlock *pth=(DttBlock*)param;

	Domain *dp = pth->dp;
	Constants *constants = pth->constants;
	TWQ* pq_tw = pth->pq_tw;

	int *curIndex = pth->curIndex;
	int *gridP = pth->gridP;
	pthread_mutex_t *mutex = pth->mutex;

	int In = dp->nSide[0];
	int Jn = dp->nSide[1];
	int Kn = dp->nSide[2];

	int cnt = 0;

	pa.x = (Real*)xmemalign(sizeof(Real)*constants->MAX_PACKAGE_SIZE*(1<<(constants->MAX_TREE_LEVEL*2)), 8301);
	pa.y = (Real*)xmemalign(sizeof(Real)*constants->MAX_PACKAGE_SIZE*(1<<(constants->MAX_TREE_LEVEL*2)), 8302);
	pa.z = (Real*)xmemalign(sizeof(Real)*constants->MAX_PACKAGE_SIZE*(1<<(constants->MAX_TREE_LEVEL*2)), 8303);
	pb.x = (Real*)xmemalign(sizeof(Real)*constants->MAX_PACKAGE_SIZE*(1<<(constants->MAX_TREE_LEVEL*2)), 8304);
	pb.y = (Real*)xmemalign(sizeof(Real)*constants->MAX_PACKAGE_SIZE*(1<<(constants->MAX_TREE_LEVEL*2)), 8305);
	pb.z = (Real*)xmemalign(sizeof(Real)*constants->MAX_PACKAGE_SIZE*(1<<(constants->MAX_TREE_LEVEL*2)), 8306);
	pc.x = (Real*)xmemalign(sizeof(Real)*constants->MAX_PACKAGE_SIZE*(1<<(constants->MAX_TREE_LEVEL*2)), 8307);
	pc.y = (Real*)xmemalign(sizeof(Real)*constants->MAX_PACKAGE_SIZE*(1<<(constants->MAX_TREE_LEVEL*2)), 8308);
	pc.z = (Real*)xmemalign(sizeof(Real)*constants->MAX_PACKAGE_SIZE*(1<<(constants->MAX_TREE_LEVEL*2)), 8309);
	mass = (Real*)xmemalign(sizeof(Real)*constants->MAX_PACKAGE_SIZE*(1<<(constants->MAX_TREE_LEVEL*2)), 8310);

	double tmutex=0.0, tqueue=0.0, tstamp;

	pth->pppar->pa=pa;
	pth->pppar->pb=pb;
	pth->pppar->pc=pc;
	pth->pppar->mass=mass;
	pth->pppar->finish=0;

	int finish=0;
	double dpp = 0.0;

	while(1)
	{
		tstamp = dtime();

		pthread_mutex_lock(mutex);
		finish = getGridIndex(gridP, curIndex, dp);
		pthread_mutex_unlock(mutex);

		tmutex += dtime() - tstamp;

		if(finish)
		{
			pth->pppar->finish=1;
			break;
		}

		cella = *curIndex;
		DBG_INFO(DBG_CALC, "[%d-%d] Calculate the %d/%d mesh.\n", dp->rank, pth->teamid, cella, pth->gridP[1]);

		tstamp = dtime();

		i1 = cella/(Jn*Kn);
		i2 = (cella-i1*Jn*Kn)/Kn;
		i3 = cella%Kn;
		if(i1==0 || i1==In-1 || i2==0 || i2==Jn-1 || i3==0 || i3==Kn-1)
		{
			continue;
		}

		if(dp->count[cella]<=0 || dp->pidx[cella]>=dp->npart)
		{
			continue;
		}
		for (n=i1-1; n<i1+2; n++)
		{
			for (m=i2-1; m<i2+2; m++)
			{
				for (l=i3-1; l<i3+2; l++)
				{
					int cellb = n*(Jn*Kn)+m*Kn+l;
					if (dp->count[cella]>0 && dp->count[cellb]>0)
					{
						enqueue_tw(pq_tw, dp->mroot[cella], dp->mroot[cellb]);
						cnt++;
					}
					else
					{
						DBG_INFO(DBG_LOGIC, "[%d] i, j, k = (%d, %d, %d)\n", dp->rank, n, m, l);
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

		while(pq_tw->length>0)
		{
			process_queue_tw(pq_tw, (void*)pth, dtt_process_cell);
		}

		tqueue += dtime() - tstamp;
	}

	DBG_INFO(DBG_LOGIC, "[%d-%d] Local tree processing counted %d pairs.\n", dp->rank, pth->teamid, cnt);
	DBG_INFO(DBG_TIME, "[%d-%d] mutex lock time = %f, queue & pp time = %f\n", dp->rank, pth->teamid, tmutex, tqueue);

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

__OffloadFunc_Macro__
int dtt_threaded(Domain *dp, Constants *constants, int gs, int ge, int nTeam, int nSlave, long *counter)
{
		int i;
		TWQ *twqueues;

		DttBlock *pth;
		pthread_t *thread;

		int tnum = nTeam*nSlave;

		pth = (DttBlock*)xmalloc(sizeof(DttBlock)*tnum, 8012);
		thread = (pthread_t*)xmalloc(sizeof(pthread_t)*tnum, 8013);
		twqueues=(TWQ*)xmalloc(sizeof(TWQ)*nTeam, 8014);

		initGlobal(part, tree);
		init_queues_tw(twqueues, nTeam);

		int gridRange[2];
		gridRange[0]=gs;
		gridRange[1]=ge;

		int *curIndex;
		pthread_mutex_t mutex;
		pthread_barrier_t *bar;
		PPParameter *pppar;

		curIndex = (int*)xmalloc(sizeof(int)*nTeam, 8015);
		bar = (pthread_barrier_t*)xmalloc(sizeof(pthread_barrier_t)*nTeam, 8016);
		pppar = (PPParameter*)xmalloc(sizeof(PPParameter)*nTeam, 8017);

		for (i=0; i<nTeam; i++)
		{
			curIndex[i]=gridRange[0];
			pthread_barrier_init(&bar[i], NULL, nSlave);
			pppar[i].finish=0;
		}
		pthread_mutex_init(&mutex, NULL);

		int teamid, core_id;
		pthread_attr_t attr;
		cpu_set_t cpuset;

		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

		for (i=0; i<tnum; i++)
		{
			teamid = i%nTeam;

			pth[i].blockid = i/nTeam;
			pth[i].teamid = teamid;
			pth[i].nSlave = nSlave;

			pth[i].counter = 0;

			pth[i].pq_tw = &twqueues[teamid];
			pth[i].dp = dp;
			pth[i].constants = constants;
			pth[i].pppar = &pppar[teamid];

			pth[i].mutex = &mutex;
			pth[i].bar = &bar[teamid];

			pth[i].gridP = gridRange;
			pth[i].curIndex = &curIndex[teamid];

#ifdef __MIC__
			CPU_ZERO(&cpuset);
			core_id = pth[i].blockid+teamid*nSlave+(dp->rank%constants->NP_PER_NODE)/constants->NMIC_PER_NODE*tnum+1;
			CPU_SET(core_id, &cpuset);
			pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpuset);
#endif
			DBG_INFO(DBG_LOGIC, "[%d-%d] blockid=%d, teamid=%d, Ngrid=%d, core=%d.\n", dp->rank, i, pth[i].blockid, teamid, dp->nSide[0]*dp->nSide[1]*dp->nSide[2], core_id);
			pthread_create(&thread[i], &attr, ((i<nTeam) ? teamMaster : compPP), (void*)&pth[i]);
		}
		pthread_attr_destroy(&attr);

		for (i=0; i<tnum; i++)
		{
			pthread_join(thread[i], NULL);
		}

		for (i=0; i<nTeam; i++)
		{
			*counter += pth[i].counter;
		}	

		free(pth);
		free(thread);
		destroy_queues_tw(twqueues, nTeam);
		free(twqueues);
		free(curIndex);
		free(bar);
		free(pppar);
}


// Use this critertian for cell/cell acceptance judging. 
__OffloadFunc_Macro__
int accepted_cell_to_cell(Int TA, Int TB, double theta/* Here theta is the open_angle  */)
{
	double delta, dr;
	int d;

	//# Strict substraction of R for correct acceptance:
	//	double Rmax = SQROOT3 * (pTA->bmax + pTB->bmax); 
	//#	Alternative subtraction of R:
	//	double Rmax = tree[TB].width;
	//	double Bmax = tree[TA].width;
	//	double Rmax = tree[TA].width + tree[TB].width;
	//	double Bmax = 0.0;
	double Rmax = tree[TA].rmax + tree[TB].rmax;
	double Bmax = 0.0;

	// Tree node include test.
	if (TA == TB)
	{
//		DBG_INFO(DBG_TMP, "[Decline] %f/%f, %f/%f, %f/%f, Rmax=%f, dr=%f, nA=%d, nB=%d.\n", tree[TB].massCenter[0], tree[TA].massCenter[0], tree[TB].massCenter[1], tree[TA].massCenter[1], tree[TB].massCenter[2], tree[TA].massCenter[2], Rmax, dr, tree[TA].nPart, tree[TB].nPart);
		return 0;
	}
	if (tree[TA].level < tree[TB].level)
	{
		if (tree[TA].firstPart<=tree[TB].firstPart && tree[TA].firstPart+tree[TA].nPart>tree[TB].firstPart)
		{
//			DBG_INFO(DBG_TMP, "[Decline] %f/%f, %f/%f, %f/%f, Rmax=%f, dr=%f, nA=%d, nB=%d.\n", tree[TB].massCenter[0], tree[TA].massCenter[0], tree[TB].massCenter[1], tree[TA].massCenter[1], tree[TB].massCenter[2], tree[TA].massCenter[2], Rmax, dr, tree[TA].nPart, tree[TB].nPart);
			return 0;
		}
	}
	else if (tree[TA].level > tree[TB].level)
	{
		if (tree[TB].firstPart<=tree[TA].firstPart && tree[TB].firstPart+tree[TB].nPart>tree[TA].firstPart)
		{
//			DBG_INFO(DBG_TMP, "[Decline] %f/%f, %f/%f, %f/%f, Rmax=%f, dr=%f, nA=%d, nB=%d.\n", tree[TB].massCenter[0], tree[TA].massCenter[0], tree[TB].massCenter[1], tree[TA].massCenter[1], tree[TB].massCenter[2], tree[TA].massCenter[2], Rmax, dr, tree[TA].nPart, tree[TB].nPart);
			return 0;
		}
	}

	dr = 0.0;

	for (d=0; d<DIM; d++)
	{
		delta = tree[TB].massCenter[d] - tree[TA].massCenter[d];
		dr += delta*delta;
	}

	dr = sqrt(dr) - Bmax;
	if ( Rmax < theta * dr )
	{
//		DBG_INFO(DBG_TMP, "[Accept] %f/%f, %f/%f, %f/%f, Rmax=%f, dr=%f, nA=%d, nB=%d.\n", tree[TB].massCenter[0], tree[TA].massCenter[0], tree[TB].massCenter[1], tree[TA].massCenter[1], tree[TB].massCenter[2], tree[TA].massCenter[2], Rmax, dr, tree[TA].nPart, tree[TB].nPart);
		return 1;
	}
	else
	{
//		DBG_INFO(DBG_TMP, "[Decline] %f/%f, %f/%f, %f/%f, Rmax=%f, dr=%f, nA=%d, nB=%d.\n", tree[TB].massCenter[0], tree[TA].massCenter[0], tree[TB].massCenter[1], tree[TA].massCenter[1], tree[TB].massCenter[2], tree[TA].massCenter[2], Rmax, dr, tree[TA].nPart, tree[TB].nPart);
		return 0;
	}
}

// Walk the tree using Dual Tree Traversel algorithm.
__OffloadFunc_Macro__
int dtt_process_cell(Int TA, Int TB, TWQ *pq_tw, void *param)
{
	DttBlock *pth = (DttBlock*) param;
	//	assert ((tree[TA].level >=0) && (tree[TA].level <= pth->constants->MAX_TREE_LEVEL));
	//	assert ((tree[TB].level >=0) && (tree[TB].level <= pth->constants->MAX_TREE_LEVEL));

	//	printf("process. TA=%d, TB=%d, theta=%.2f.\n", TA, TB, theta);
	//
	Array3 pa = pth->pppar->pa;
	Array3 pb = pth->pppar->pb;
	Array3 pc = pth->pppar->pc;
	Real *mass = pth->pppar->mass;

	double theta = pth->constants->OPEN_ANGLE;
	Real partmass = pth->constants->PART_MASS;

	int minlevel = pth->constants->MIN_TREE_LEVEL;
	int maxlevel = pth->constants->MAX_TREE_LEVEL;

	Int *pT = (Int*)xmalloc(maxlevel*sizeof(Int), 8801);

#ifdef __MIC__
	Int ppsize = pth->constants->MIC_PP_THRESHOLD;
#endif

	int nA, nB, calcm;

	if(TA == TB)
	{
#ifdef __MIC__
		if (tree[TA].nChildren == 0 || tree[TA].nPart <= ppsize)
#else
		if (tree[TA].nChildren == 0)
#endif
		{
			/*      Here we send the tree node, and the Queue      */
			/*      processing function should pack the particles  */
			/*      and do some further optimizatoins.             */
			/*      PQ_particle: save the packed pps.              */
			nA = tree[TA].nPart;
			nB = nA;
			packarray3(&part[tree[TA].firstPart], nA, pa);
			packarray3(&part[tree[TA].firstPart], nB, pb);
			packarray3(NULL, nA, pc);
			pth->pppar->nA = nA;
			pth->pppar->nB = nB;
			pth->pppar->calcm = 0;

//			DBG_INFO(DBG_TMP, "[Accum] TA==TB/Dec  TA=%ld, nA=%d, nB=%d, calcm=%d\n", TA, nA, nB, calcm); 

			pth->counter += nA*(nB*22+2);
			pthread_barrier_wait(pth->bar);
			ppkernel(pa, nA, pb, nB, pth->constants, pc, pth->nSlave, pth->blockid);
			pthread_barrier_wait(pth->bar);
			pusharray3m(&part[tree[TA].firstPart], nA, pc, partmass);
		}
		else
		{
			int i, j;
			for (i=0; i<tree[TA].nChildren; i++)
			{
				int TA1 = tree[TA].firstChild+i;
				if (tree[TA1].nChildren == 0)
				{
					calcm = 0;
					nA = tree[TA1].nPart;
					nB = 0;
					for (j=0; j<tree[TB].nChildren; j++)
					{
						int TB1 = tree[TB].firstChild+j;
						if (accepted_cell_to_cell(TA1, TB1, theta))
						{
							calcm = 1;
							pb.x[nB] = tree[TB1].massCenter[0];
							pb.y[nB] = tree[TB1].massCenter[1];
							pb.z[nB] = tree[TB1].massCenter[2];
							mass[nB] = tree[TB1].mass;
//							DBG_INFO(DBG_TMP, "[AcceptInfo]  %f  %f  %f  %f  %f  %f\n", pb.x[nB], pb.y[nB], pb.z[nB], mass[nB], tree[TB1].mass,(Real)tree[TB1].nPart*partmass);
//							int xx;
//							for (xx=0; xx<tree[TB1].nPart; xx++)
//							{
//								pb.x[nB+xx] = tree[TB1].massCenter[0];
//								pb.y[nB+xx] = tree[TB1].massCenter[1];
//								pb.z[nB+xx] = tree[TB1].massCenter[2];
//								mass[nB+xx] = partmass;
//							}
//							nB += tree[TB1].nPart;
							nB += 1;
//							DBG_INFO(DBG_TMP, "[Accum] TA==TB/TB1_Acc  TA1=%ld, nA=%d, nB=%d, calcm=%d\n", TA1, nA, nB, calcm); 
						}
						else
						{
							if(tree[TB1].nChildren == 0)
							{
								packarray3om(&part[tree[TB1].firstPart], nB, tree[TB1].nPart, pb, mass);
								nB += tree[TB1].nPart;
//								DBG_INFO(DBG_TMP, "[Accum] TA==TB/TB1_Dec  TA1=%ld, nA=%d, nB=%d, calcm=%d\n", TA1, nA, nB, calcm); 
							}
							else
							{
								enqueue_tw(pq_tw, TA1, TB1);
							}
						}
					}
					if (nB>0)
					{
						packarray3(&part[tree[TA1].firstPart], nA, pa);
						packarray3(NULL, nA, pc);
						pth->pppar->nA = nA;
						pth->pppar->nB = nB;
						pth->pppar->calcm = calcm;
						if (calcm == 0)
						{
							pthread_barrier_wait(pth->bar);
							pth->counter += nA*(nB*22+2);
							ppkernel(pa, nA, pb, nB, pth->constants, pc, pth->nSlave, pth->blockid);
							pthread_barrier_wait(pth->bar);
							pusharray3m(&part[tree[TA1].firstPart], nA, pc, partmass);
						}
						else
						{	
							pthread_barrier_wait(pth->bar);
							pth->counter += nA*(nB*23+1);
							ppmkernel(pa, nA, pb, mass, nB, pth->constants, pc, pth->nSlave, pth->blockid);
							pthread_barrier_wait(pth->bar);
							pusharray3(&part[tree[TA1].firstPart], nA, pc);
						}
					}
				}
				else
				{
					enqueue_tw(pq_tw, TA1, TB);
				}
			}
		}
	}
	else
	{
		calcm = 0;
		nA = tree[TA].nPart;
		nB = 0;
		if (accepted_cell_to_cell(TA, TB, theta))
		{
			calcm = 1;
			pb.x[nB] = tree[TB].massCenter[0];
			pb.y[nB] = tree[TB].massCenter[1];
			pb.z[nB] = tree[TB].massCenter[2];
			mass[nB] = tree[TB].mass;
//			DBG_INFO(DBG_TMP, "[AcceptInfo]  %f  %f  %f  %f  %f  %f\n", pb.x[nB], pb.y[nB], pb.z[nB], mass[nB], tree[TB].mass, (Real)tree[TB].nPart*partmass);
//			int xx;
//			for (xx=0; xx<tree[TB].nPart; xx++)
//			{
//				pb.x[nB+xx] = tree[TB].massCenter[0];
//				pb.y[nB+xx] = tree[TB].massCenter[1];
//				pb.z[nB+xx] = tree[TB].massCenter[2];
//				mass[nB+xx] = partmass;
//			}
//			nB += tree[TB].nPart;
			nB = 1;
//			DBG_INFO(DBG_TMP, "[Accum] TB_Acc  TA=%ld, nA=%d, nB=%d, calcm=%d\n", TA, nA, nB, calcm); 
		}
		else
		{
#ifdef __MIC__
			if ( (tree[TA].nChildren==0 || tree[TA].nPart<=ppsize) && 
				(tree[TB].nChildren==0 || tree[TB].nPart<=ppsize) )
#else
			if ( tree[TA].nChildren==0 && tree[TB].nChildren==0 )
#endif
			{
				packarray3(&part[tree[TB].firstPart], nB, pb);
				nB = tree[TB].nPart;
//				DBG_INFO(DBG_TMP, "[Accum] TB_Dec  TA=%ld, nA=%d, nB=%d, calcm=%d\n", TA, nA, nB, calcm); 
			}
#ifdef __MIC__
			else if ((tree[TB].nChildren==0 || tree[TB].nPart<=ppsize) || 
				((tree[TA].nChildren>0 && tree[TA].nPart>ppsize) && tree[TA].width>tree[TB].width))
//			else if (((tree[TA].nChildren>0 && tree[TA].nPart>ppsize) && tree[TA].width>tree[TB].width))
#else
			else if (tree[TB].nChildren==0 || (tree[TA].nChildren>0 && tree[TA].width>=tree[TB].width))
//			else if ( tree[TA].nChildren>0 && tree[TA].width>tree[TB].width )
#endif
//			else if ( tree[TB].nChildren==0 || tree[TA].nChildren>0 ) // Open A until reach the leaf.
			{
				int i;
				for (i=0; i<tree[TA].nChildren; i++)
				{
					enqueue_tw(pq_tw, tree[TA].firstChild+i, TB);
				}
			}
#ifdef __MIC__
			else if ( tree[TA].nChildren==0 || tree[TA].nPart<=ppsize )
#else
			else if ( tree[TA].nChildren==0 )
#endif
			{
				Int TB1, TBp;
				int level = tree[TB].level;
				pT[level] = TB;

				TBp = TB;
				level++;
				pT[level] = tree[TB].firstChild;
				TB1 = pT[level];
				while (level>tree[TB].level && TB1<tree[TBp].firstChild+tree[TBp].nChildren)
				{
//					DBG_INFO(DBG_TMP, "[Search] level=%d, TB=%d,%d, TB1=%d,%d,%d, TBp=%d,%d,%d,%d.\n", level, TB, tree[TB].level, TB1, tree[TB1].level, tree[TB1].nPart, TBp, tree[TBp].level, tree[TBp].firstChild, tree[TBp].nChildren);
					if (accepted_cell_to_cell(TA, TB1, theta))
					{
						calcm = 1;
						pb.x[nB] = tree[TB1].massCenter[0];
						pb.y[nB] = tree[TB1].massCenter[1];
						pb.z[nB] = tree[TB1].massCenter[2];
						mass[nB] = tree[TB1].mass;
//						DBG_INFO(DBG_TMP, "[AcceptInfo]  %f  %f  %f  %f  %f  %f\n", pb.x[nB], pb.y[nB], pb.z[nB], mass[nB], tree[TB1].mass, (Real)tree[TB1].nPart*partmass);
//						int xx;
//						for (xx=0; xx<tree[TB1].nPart; xx++)
//						{
//							pb.x[nB+xx] = tree[TB1].massCenter[0];
//							pb.y[nB+xx] = tree[TB1].massCenter[1];
//							pb.z[nB+xx] = tree[TB1].massCenter[2];
//							mass[nB+xx] = partmass;
//						}
//						nB += tree[TB1].nPart;
						nB += 1;
//						DBG_INFO(DBG_TMP, "[Accum] TA/TB1_Acc  TA=%ld, nA=%d, nB=%d, calcm=%d\n", TA, nA, nB, calcm); 

						TB1++;
						while(level>tree[TB].level && TB1>=tree[TBp].firstChild+tree[TBp].nChildren)
						{
							level--;
							TB1 = pT[level]+1;
							if(level>tree[TB].level)
							{
								TBp = pT[level-1];
							}
//							DBG_INFO(DBG_TMP, "[Search up1] level=%d, TB=%ld,%d, TB1=%ld, %d, TBp=%ld,%d,%ld,%d.\n", level, TB, tree[TB].level, TB1, tree[TB1].level, TBp, tree[TBp].level, tree[TBp].firstChild, tree[TBp].nChildren);
						}
						if (level>tree[TB].level && TB1<tree[TBp].firstChild+tree[TBp].nChildren)
						{
							pT[level] = TB1;
						}
//						DBG_INFO(DBG_TMP, "[Search ver1] level=%d, TB=%ld,%d, TB1=%ld, %d, TBp=%ld,%d,%ld,%d.\n", level, TB, tree[TB].level, TB1, tree[TB1].level, TBp, tree[TBp].level, tree[TBp].firstChild, tree[TBp].nChildren);
					}
#ifdef __MIC__
					else if (tree[TB1].nChildren == 0 || tree[TB1].nChildren <= ppsize)
#else
					else if (tree[TB1].nChildren == 0)
#endif
					{
//						DBG_INFO(DBG_TMP, "[Search t21] level=%d, nB=%d, TB1=%d,%d,%d, TBp=%d,%d,%d,%d.\n", level, nB, TB1, tree[TB1].level, tree[TB1].nPart, TBp, tree[TBp].level, tree[TBp].firstChild, tree[TBp].nChildren);
						packarray3om(&part[tree[TB1].firstPart], nB, tree[TB1].nPart, pb, mass);
						nB += tree[TB1].nPart;
//						DBG_INFO(DBG_TMP, "[Search t22] level=%d, nB=%d, TB1=%d,%d,%d, TBp=%d,%d,%d,%d.\n", level, nB, TB1, tree[TB1].level, tree[TB1].nPart, TBp, tree[TBp].level, tree[TBp].firstChild, tree[TBp].nChildren);
//						DBG_INFO(DBG_TMP, "[Accum] TA/TB1_Dec  TA=%ld, nA=%d, nB=%d, calcm=%d\n", TA, nA, nB, calcm); 

						TB1++;
						while(level>tree[TB].level && TB1>=tree[TBp].firstChild+tree[TBp].nChildren)
						{
							level--;
							TB1 = pT[level]+1;
							if(level>tree[TB].level)
							{
								TBp = pT[level-1];
							}
//							DBG_INFO(DBG_TMP, "[Search up2] level=%d, TB=%ld,%d, TB1=%ld, %d, TBp=%ld,%d,%ld,%d.\n", level, TB, tree[TB].level, TB1, tree[TB1].level, TBp, tree[TBp].level, tree[TBp].firstChild, tree[TBp].nChildren);
						}
						if (level>tree[TB].level && TB1<tree[TBp].firstChild+tree[TBp].nChildren)
						{
							pT[level] = TB1;
						}
//						DBG_INFO(DBG_TMP, "[Search ver2] level=%d, TB=%ld,%d, TB1=%ld, %d, TBp=%ld,%d,%ld,%d.\n", level, TB, tree[TB].level, TB1, tree[TB1].level, TBp, tree[TBp].level, tree[TBp].firstChild, tree[TBp].nChildren);
					}
					else
					{
						TBp = TB1;
						level++;
						pT[level] = tree[TB1].firstChild;
						TB1 = pT[level];
//						DBG_INFO(DBG_TMP, "[Search down] level=%d, TB=%ld,%d, TB1=%ld, %d, TBp=%ld,%d,%ld,%d.\n", level, TB, tree[TB].level, TB1, tree[TB1].level, TBp, tree[TBp].level, tree[TBp].firstChild, tree[TBp].nChildren);
					}
				}
			}
			else
			{
				int i;
				for (i=0; i<tree[TB].nChildren; i++)
				{
					int TB1 = tree[TB].firstChild+i;
					if (accepted_cell_to_cell(TA, TB1, theta))
					{
						calcm = 1;
						pb.x[nB] = tree[TB1].massCenter[0];
						pb.y[nB] = tree[TB1].massCenter[1];
						pb.z[nB] = tree[TB1].massCenter[2];
						mass[nB] = tree[TB1].mass;
//						DBG_INFO(DBG_TMP, "[AcceptInfo]  %f  %f  %f  %f  %f  %f\n", pb.x[nB], pb.y[nB], pb.z[nB], mass[nB], tree[TB1].mass, (Real)tree[TB1].nPart*partmass);
//						int xx;
//						for (xx=0; xx<tree[TB1].nPart; xx++)
//						{
//							pb.x[nB+xx] = tree[TB1].massCenter[0];
//							pb.y[nB+xx] = tree[TB1].massCenter[1];
//							pb.z[nB+xx] = tree[TB1].massCenter[2];
//							mass[nB+xx] = partmass;
//						}
//						nB += tree[TB1].nPart;
						nB += 1;
//						DBG_INFO(DBG_TMP, "[Accum] TB1/TB1_Acc  TA=%ld, nA=%d, nB=%d, calcm=%d\n", TA, nA, nB, calcm); 
					}
					else
					{
#ifdef __MIC__
//						if((tree[TB1].nChildren == 0 || tree[TB1].nPart <= ppsize))
						if(tree[TB1].nChildren == 0 || tree[TB1].nPart <= ppsize)
#else
//						if(tree[TB1].nChildren == 0)
						if(tree[TB1].nChildren == 0 && tree[TA].level > minlevel)
#endif
						{
							packarray3om(&part[tree[TB1].firstPart], nB, tree[TB1].nPart, pb, mass);
							nB += tree[TB1].nPart;
//							DBG_INFO(DBG_TMP, "[Accum] TB1/TB1_Dec  TA=%ld, nA=%d, nB=%d, calcm=%d\n", TA, nA, nB, calcm); 
						}
						else
						{
							int j;
							for (j=0; j<tree[TA].nChildren; j++)
							{
								int TA1 = tree[TA].firstChild+j;
								enqueue_tw(pq_tw, TA1, TB1);
							}
						}
					}
				}
			}
		}
		if (nB>0)
		{
			packarray3(&part[tree[TA].firstPart], nA, pa);
			packarray3(NULL, nA, pc);
			pth->pppar->nA = nA;
			pth->pppar->nB = nB;
			pth->pppar->calcm = calcm;
			if (calcm == 0)
			{
				pthread_barrier_wait(pth->bar);
				pth->counter += nA*(nB*22+2);
				ppkernel(pa, nA, pb, nB, pth->constants, pc, pth->nSlave, pth->blockid);
				pthread_barrier_wait(pth->bar);
				pusharray3m(&part[tree[TA].firstPart], nA, pc, partmass);
			}
			else
			{	
				pthread_barrier_wait(pth->bar);
				pth->counter += nA*(nB*23+1);
				ppmkernel(pa, nA, pb, mass, nB, pth->constants, pc, pth->nSlave, pth->blockid);
				pthread_barrier_wait(pth->bar);
				pusharray3(&part[tree[TA].firstPart], nA, pc);
			}
		}
	}
	free(pT);
}

void tree_traversal(Domain *dp, Constants *constants)
{
	Int i, j, k, m, n, l;

	int npnode = constants->NP_PER_NODE;
	int nmicpn = constants->NMIC_PER_NODE;
	int nTeam_MIC = constants->NTEAM_MIC;
	int nSlave_MIC = constants->NTHREAD_MIC*nmicpn/npnode/nTeam_MIC;
	int nTeam_CPU = constants->NTEAM_CPU;
	int nSlave_CPU = constants->NTHREAD_CPU/nTeam_CPU;

	Int npart = dp->npart;
	Int allpart = npart + dp->npart_adj;
	Int nnode = dp->nnode;

	long counter = 0;

	int In = dp->nSide[0];
	int Jn = dp->nSide[1];
	int Kn = dp->nSide[2];
	int Ngrid = In*Jn*Kn;
	int gstart = (Jn+1)*Kn+1; // (1, 1, 1)
	int gend = (In-2)*Jn*Kn+(Jn-2)*Kn+Kn-1; // (In-2, Jn-2, Kn-2) + 1

	fpart = dp->part;
	tree = dp->tree;
	part = (CalcBody*)xmemalign(sizeof(CalcBody)*allpart, 8001);

	for (i=0; i<allpart; i++)
	{
		part[i].pos[0] = fpart[i].pos[0];
		part[i].pos[1] = fpart[i].pos[1];
		part[i].pos[2] = fpart[i].pos[2];
		part[i].acc[0] = (Real) 0.0;
		part[i].acc[1] = (Real) 0.0;
		part[i].acc[2] = (Real) 0.0;
		part[i].mass = fpart[i].mass;
	}

	int myid=dp->rank;
	int nproc = dp->nproc;
	int micid=myid%npnode%nmicpn;

#ifdef __INTEL_OFFLOAD
	// Preparing for array offload.
	int *count = dp->count;
	int *mroot = dp->mroot;
	int *pidx = dp->pidx;
	int *tag = dp->tag;

	Int splitnum = (Int)((double)npart*constants->CPU_RATIO);

	Int splitpart;
	int splitgrid;

	long counter_mic = 0;

	for (n=gstart; n<gend; n++)
	{
		if((dp->tag[n] == TAG_OUTERCORE || dp->tag[n] == TAG_FRONTIERS) && (dp->pidx[n] >= splitnum))
		{
			break;
		}
	}
	splitgrid = n;
	if (n==gend)
	{
		splitpart = npart;
	}
	else
	{
		splitpart = dp->pidx[n];
	}

	double tstamp, tstampm, tstampc;
	tstamp = dtime();

#pragma offload target(mic:micid) \
	in (constants[0:1]) in (dp[0:1]) in (part[0:allpart]) in (tree[0:nnode]) \
	in (count[0:Ngrid]) in (tag[0:Ngrid]) in (mroot[0:Ngrid]) in (pidx[0:Ngrid]) \
	out(part[splitpart:npart-splitpart]) signal(part)
	{
		tstampm = dtime();

		// Processing offload array pointers.
		dp->count = count;
		dp->mroot = mroot;
		dp->pidx = pidx;
		dp->tag = tag;

		DBG_INFO(DBG_TMP, "[%d] Treewakling: start walking %d-%d with %d/%d Teams on MIC%d.\n", dp->rank, splitgrid, gend, nTeam_MIC, nSlave_MIC, micid);
		dtt_threaded(dp, constants, splitgrid, gend, nTeam_MIC, nSlave_MIC, &counter_mic);

		tstampm = dtime() - tstampm;
		DBG_INFO(DBG_TIME, "[%d] offload threaded part time = %f.\n", dp->rank, tstampm);
	}

	DBG_INFO(DBG_TMP, "[%d] Treewakling: start walking %d-%d with %d/%d Teams.\n", dp->rank, gstart, splitgrid, nTeam_CPU, nSlave_CPU);
	dtt_threaded(dp, constants, gstart, splitgrid, nTeam_CPU, nSlave_CPU, &counter);
	tstampc = dtime() - tstamp;
#pragma offload_wait target(mic:micid) wait(part)
	counter += counter_mic;

	DBG_INFO(DBG_TIME, "[%d] Treewalking: offload and transfer time = %f, ratio = %.5f, cnt = %ld.\n", dp->rank, dtime()-tstamp, constants->CPU_RATIO, counter);
	if(constants->DYNAMIC_LOAD && (1.15*tstampc < tstampm || 0.9*tstampc > tstampm))
	{
		double vratio = ((tstampm / tstampc) - 1) * 0.7;
		int vgrid = (int)((double) (splitgrid - gstart) * vratio);
		if (vgrid != 0)
		{
			constants->CPU_RATIO += constants->CPU_RATIO*vratio;
		}
	}
#else // __INTEL_OFFLOAD No offload
	double tstamp;
	tstamp = dtime();

	DBG_INFO(DBG_TMP, "[%d] Treewakling: start walking %d-%d with %d/%d Teams.\n", dp->rank, gstart, gend, nTeam_CPU, nSlave_CPU);
	dtt_threaded(dp, constants, gstart, gend, nTeam_CPU, nSlave_CPU, &counter);

	DBG_INFO(DBG_TIME, "[%d] Treewalking time = %f, cnt = %ld.\n", dp->rank, dtime()-tstamp, counter);
#endif

	for (i=0; i<npart; i++)
	{
#ifdef NMK_NAIVE_GRAVITY
		fpart[i].acc[0] += part[i].acc[0];
		fpart[i].acc[1] += part[i].acc[1];
		fpart[i].acc[2] += part[i].acc[2];
#else
		fpart[i].acc[0] += constants->GRAV_CONST*part[i].acc[0];
		fpart[i].acc[1] += constants->GRAV_CONST*part[i].acc[1];
		fpart[i].acc[2] += constants->GRAV_CONST*part[i].acc[2];
#endif
	}

	free(part);

#ifdef NMK_VERIFY
	char file[256];
	sprintf(file, "pt.acc.%d", dp->rank);
	FILE *fp = fopen(file, "w+");
	for (i=0; i<npart; i++)
	{
		fprintf(fp, "%d  %f  %f  %f\n", fpart[i].id, fpart[i].acc[0], fpart[i].acc[1], fpart[i].acc[2]);
//			fprintf(fp, "%d  %f  %f  %f\n", fpart[i].id, part[i].acc[0]*constants->GRAV_CONST, part[i].acc[1]*constants->GRAV_CONST, part[i].acc[2]*constants->GRAV_CONST);
	}
	fclose(fp);
	MPI_Barrier(MPI_COMM_WORLD);
	exit(0);
#endif
}
