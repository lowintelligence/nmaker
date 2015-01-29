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

#include "global.h"

#define QUEUE_BLOCK_SIZE 1048576

typedef struct
{
	Int TA;
	Int TB;
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

__OffloadFunc_Macro__
int init_queues_tw(TWQ *pq, int n);

__OffloadFunc_Macro__
int enqueue_tw(TWQ *pq, Int ta, Int tb);

__OffloadFunc_Macro__
int dequeue_tw(TWQ *pq, TWelement *peout);

__OffloadFunc_Macro__
int destroy_queues_tw(TWQ *pq, int n);

__OffloadFunc_Macro__
int process_queue_tw(TWQ *PQ_treewalk, void* param, int process(Int, Int, TWQ*, void*));

__OffloadFunc_Macro__
void initGlobal(CalcBody* p, Node* t);
#endif /* _QUEUE_H_ */
