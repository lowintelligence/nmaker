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
#include <pthread.h>

#define NTEAM   12 
#define MASTER  1
#define NSLAVE  1

#define QUEUE_BLOCK_SIZE 33554432

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

int init_queue_pp(PPQ *pq);

int enqueue_pp(PPQ *pq, int ta, int tb, int m);

int dequeue_pp(PPQ *pq, PPelement *peout);

int destroy_queue_pp(PPQ *pq);

int EnqueueP_PPnode(PPQ *PQ_ppnode, int TA, int TB, int mask);
int ProcessQP_PPnode(PPQ *PQ_ppnode, int process(Array3, int, Array3, int, PRECTYPE, Array3), Array3 pA, Array3 pB, Array3 pC);

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
	struct _block* pthArr;
	int* lower;
	int* upper;
	PPQ* P_PQ_ppnode;
	TWQ* P_PQ_treewalk;
	
} Block;

int init_queue_tw(TWQ *pq);

int enqueue_tw(TWQ *pq, int ta, int tb, double theta);

int dequeue_tw(TWQ *pq, TWelement *peout);

int destroy_queue_tw(TWQ *pq);

int EnqueueP_Cell(TWQ *PQ_treewalk, int TA, int TB, double theta);
int ProcessQP_Cell(TWQ *PQ_treewalk, PPQ *PQ_ppnode, int process(int, int, double, PPQ*, TWQ*));

void initGlobal(Body* p, Node* t);
#endif


