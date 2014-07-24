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

#define QUEUE_BLOCK_SIZE 1048576

typedef struct
{
	Array3 A;
	Array3 B;
	Array3 C;
	int nA;
	int nB;	
} PPelement;

typedef struct
{
	int size;
	int length;
	int head;
	int tail;
	PPelement* elements;
} PPQ;

int init_queue_pp(PPQ *pq);

int enqueue_pp(PPQ *pq, Array3 pa, Array3 pb, Array3 pc, int na, int nb);

int dequeue_pp(PPQ *pq, PPelement *peout);

int destroy_queue_pp(PPQ *pq);

int EnqueueP_PPnode(PPQ *PQ_ppnode, int TA, int TB, int mask);
int ProcessQP_PPnode(PPQ *PQ_ppnode, int process(Array3, int, Array3, int, PRECTYPE, Array3));

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

int init_queue_tw(TWQ *pq);

int enqueue_tw(TWQ *pq, int ta, int tb, double theta);

int dequeue_tw(TWQ *pq, TWelement *peout);

int destroy_queue_tw(TWQ *pq);

int EnqueueP_Cell(TWQ *PQ_treewalk, int TA, int TB, double theta);
int ProcessQP_Cell(TWQ *PQ_treewalk, PPQ *PQ_ppnode, int process(int, int, double, PPQ*, TWQ*));

void initGlobal(Body* p, Node* t);
#endif
