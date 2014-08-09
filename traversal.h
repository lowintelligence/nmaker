#ifndef _TRAVERSAL_H_
#define _TRAVERSAL_H_

#include "data.h"
#include "offload.h"
#include "domain.h"
#include "subtree.h"
#include "parameter.h"
#include "subcuboid.h"
#include "ppkernel.h"
#include "queue.h"

#define SQROOT3 1.73205080705688773

#define NTEAM   16
#define MASTER  1
#define NSLAVE  2

__OffloadVar_Macro__
typedef struct _block
{
	int blockid;//id in one block
	int bsize;//size of block
	int tid;// thread id in all
	int tall;//total num of threads

	int teamid; //id of teams
	PPQ* PQ_ppnode;
	TWQ* PQ_treewalk;
	PPQ* PQ_ptnode; //to be implemented
	int first;
	int last;

	Domain* dp;
	GlobalParam* gp;

	pthread_t* thread;
	struct _block* pth;
	int* lower;
	int* upper;
	PPQ* P_PQ_ppnode;
	TWQ* P_PQ_treewalk;
	
} Block;

// Qiao's functions.
//int accepted_cell_branes(int one, int ns, Node *tree, Body *part);
//void naive_walk_cell(int add, int level, int body, Node *tree, Body *part);
//void inner_traversal(Domain *dp, GlobalParam *gp, int nThread);

// Cao!'s functions.
__DefOffloadFunc__
int accepted_cell_to_cell(int TA, int TB, double theta);

__DefOffloadFunc__
int dtt_process_cell(int TA, int TB, double theta, PPQ* PQ_ppnode, TWQ* PQ_treewalk);

__DefOffloadFunc__
void dtt_traversal(Domain *dp, GlobalParam *gp);

#endif /* TRAVERSAL_H */
