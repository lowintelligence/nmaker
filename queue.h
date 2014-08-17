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
#include "domain.h"
#include "subtree.h"
#include "parameter.h"
#include "subcuboid.h"
#include "ppkernel.h"
#include "offload.h"
#include <pthread.h>

#define NTEAM   60
#define MASTER  1
#define NSLAVE  1

#define QUEUE_BLOCK_SIZE 16777216

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
int init_queue_pp(PPQ *pq);

__OffloadFunc_Macro__
int enqueue_pp(PPQ *pq, int ta, int tb, int m);

__OffloadFunc_Macro__
int dequeue_pp(PPQ *pq, PPelement *peout);

__OffloadFunc_Macro__
int destroy_queue_pp(PPQ *pq);

__OffloadFunc_Macro__
int EnqueueP_PPnode(PPQ *PQ_ppnode, int TA, int TB, int mask);

__OffloadFunc_Macro__
int ProcessQP_PPnode(PPQ *PQ_ppnode, int process(Array3, int, Array3, int, PRECTYPE, Array3), Array3 pA, Array3 pB, Array3 pC, pthread_mutex_t *mutex);

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

	Domain* dp;
	GlobalParam* gp;

	pthread_t* thread;
	pthread_mutex_t *mutex;
	struct _block* pthArr;
	int* lower;
	int* upper;
	PPQ* P_PQ_ppnode;
	TWQ* P_PQ_treewalk;
	
} Block;

__OffloadFunc_Macro__
int init_queue_tw(TWQ *pq);

__OffloadFunc_Macro__
int enqueue_tw(TWQ *pq, int ta, int tb, double theta);

__OffloadFunc_Macro__
int dequeue_tw(TWQ *pq, TWelement *peout);

__OffloadFunc_Macro__
int destroy_queue_tw(TWQ *pq);

__OffloadFunc_Macro__
int EnqueueP_Cell(TWQ *PQ_treewalk, int TA, int TB, double theta);

__OffloadFunc_Macro__
int ProcessQP_Cell(TWQ *PQ_treewalk, PPQ *PQ_ppnode, int process(int, int, double, PPQ*, TWQ*));

__OffloadFunc_Macro__
void initGlobal(Body* p, Node* t);
#endif
