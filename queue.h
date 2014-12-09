/*
 * =====================================================================================
 *
 *       Filename:  queue.h
 *
 *    Description:  Queue interface and scheduler header file
 *
 *        Version:  1.0
 *        Created:  07/24/2014 02:26:25 PM
 *       Revision:  none
 *       Compiler:  icc
 *
 *         Author:  Cao Zongyan (), zycao@sccas.cn
 *        Company:  Supercomputing Center, CNIC, CAS, China
 *
 * =====================================================================================
 */


#ifndef _QUEUE_H_
#define _QUEUE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "data.h"
#include "dtime.h"
#include "domain.h"
#include "subtree.h"
#include "parameter.h"
#include "subcuboid.h"
#include "ppkernel.h"
#include "offload.h"
#include <pthread.h>

#define QUEUE_BLOCK_SIZE 2097152

typedef struct
{
	int TA;
	int TB;
	int mask;
} PPelement;

typedef struct
{
	int size;
	int length;
	int head;
	int tail;
	PPelement* elements;
	int tag;//the queue is not to be added
} PPQ;

__OffloadFunc_Macro__
int packarray3(SBody* pp, int n, Array3 pa);

__OffloadFunc_Macro__
int packarray3m(SBody* pp, int n, Array3 pa, PRECTYPE *mass);

__OffloadFunc_Macro__
int packarray3o(SBody* pp, int offset, int n, Array3 pa);

__OffloadFunc_Macro__
int packarray3om(SBody* pp, int offset, int n, Array3 pa, PRECTYPE *mass);

__OffloadFunc_Macro__
int pusharray3(SBody* pp, int n, Array3 pa);

__OffloadFunc_Macro__
int init_queue_pp(PPQ *pq, int n);

__OffloadFunc_Macro__
int enqueue_pp(PPQ *pq, int ta, int tb, int m);

__OffloadFunc_Macro__
int dequeue_pp(PPQ *pq, PPelement *peout);

__OffloadFunc_Macro__
int destroy_queue_pp(PPQ *pq);

__OffloadFunc_Macro__
int EnqueueP_PPnode(PPQ *PQ_ppnode, int TA, int TB, int mask);

__OffloadFunc_Macro__
double ProcessQP_PPnode(PPQ *PQ_ppnode, int process(Array3, int, Array3, int, PRECTYPE, Array3), Array3 pA, Array3 pB, Array3 pC);

typedef struct
{
	int TA;
	int TB;
	double theta;
} TWelement;

typedef struct
{
	int size;
	int length;
	int head;
	int tail;
	TWelement* elements;
} TWQ;

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
	int nTeam;
	int nMaster;
	int nSlave;

	Domain* dp;
	GlobalParam* gp;

	pthread_t* thread;
	pthread_mutex_t *mutex;
	pthread_barrier_t *bar;
	PPParameter *pppar;
	struct _block* pthArr;
	int* lower;
	int* upper;
	PPQ* P_PQ_ppnode;
	TWQ* P_PQ_treewalk;

    int* gridP;
	int* numpart;
	int* curIndex;	
} Block;

__OffloadFunc_Macro__
int init_queue_tw(TWQ *pq, int n);

__OffloadFunc_Macro__
int enqueue_tw(TWQ *pq, int ta, int tb, double theta);

__OffloadFunc_Macro__
int dequeue_tw(TWQ *pq, TWelement *peout);

__OffloadFunc_Macro__
int destroy_queue_tw(TWQ *pq, int n);

__OffloadFunc_Macro__
int EnqueueP_Cell(TWQ *PQ_treewalk, int TA, int TB, double theta);

__OffloadFunc_Macro__
int ProcessQP_Cell(TWQ *PQ_treewalk, Block* pth, int process(int, int, double, TWQ*, Block*));

__OffloadFunc_Macro__
void initGlobal(SBody* p, Node* t);
#endif
