#ifndef _TRAVERSAL_H_
#define _TRAVERSAL_H_

#include "global.h"
#include "domain.h"
#include "queue.h"
#include "ppkernel.h"

#define SQROOT3 1.73205080705688773

typedef struct _taskinfo
{
	int gridid;
	int npart;
} GridTask;

typedef struct _dttblock
{
	int blockid;//id in one block
	int teamid;// id of team
	int nSlave;

	long counter;

	TWQ* pq_tw;
	Domain* dp;
	Constants* constants;
	PPParameter *pppar;

	pthread_mutex_t *mutex;
	pthread_barrier_t *bar;

    int* gridP;
	int* curIndex;
	GridTask* gtask;
	int maxpart; 
} DttBlock;

__OffloadFunc_Macro__
int accepted_cell_to_cell(Int TA, Int TB, double theta);

__OffloadFunc_Macro__
int dtt_process_cell(Int TA, Int TB, TWQ* pq_tw, void *param);

void tree_traversal(Domain *dp, Constants *constants);
int checkaccpm(int rank, Body *part, Int n);
int checkpart(int rank, Body *part, Int n, Real err);

#endif /* _TRAVERSAL_H_ */
