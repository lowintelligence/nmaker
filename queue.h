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

int init_queue_pp(PPQ *pq)
{
	pq->size = QUEUE_BLOCK_SIZE;
	pq->length = 0;
	pq->head = 0;
	pq->tail = 1;
	pq->elements = (PPelement*)malloc(sizeof(PPelement)*pq->size);
	return 0;
}

int enqueue_pp(PPQ *pq, Array3 pa, Array3 pb, Array3 pc, int na, int nb)
{
	PPelement *pe;

	if (pq->length == pq->size-256)
	{
		pq->size += QUEUE_BLOCK_SIZE;
		pq->elements = (PPelement*)realloc(pq->elements, pq->size*sizeof(PPelement));
		if (pq->tail < pq->head)
		{
			pe = pq->elements;
			memcpy((void*)pe+pq->size-QUEUE_BLOCK_SIZE, (void*)pe, sizeof(PPelement)*pq->tail);
		}
	}
	pe[pq->tail].A = pa;
	pe[pq->tail].B = pb;
	pe[pq->tail].C = pc;
	pe[pq->tail].nA = na;
	pe[pq->tail].nB = nb;
	pq->tail = (pq->tail+1)%pq->size;
	pq->length++;
	return 0;
}

int dequeue_pp(PPQ *pq, PPelement *peout)
{
	PPelement *pe = pq->elements+pq->head;

	if (pq->length>0)
	{
		peout->A = pe->A;
		peout->B = pe->B;
		peout->C = pe->C;
		peout->nA = pe->nA;
		peout->nA = pe->nB;
		pq->head = (pq->head+1)%pq->head;
		pq->length--;
		return 0;
	}
	else
	{
		return -1;
	}
}

int destroy_queue_pp(PPQ *pq)
{
	if (pq->size > 0)
	{
		free(pq->elements);
		pq->elements = NULL;
		pq->head = pq->tail = -1;
		pq->size = 0;
		pq->length = 0;
	}
}

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

int init_queue_tw(TWQ *pq)
{
	pq->size = QUEUE_BLOCK_SIZE;
	pq->length = 0;
	pq->head = 0;
	pq->tail = 1;
	pq->elements = (TWelement*)malloc(sizeof(TWelement)*pq->size);
	return 0;
}

int enqueue_tw(TWQ *pq, int ta, int tb, double theta)
{
	TWelement *pe;

	if (pq->length == pq->size-256)
	{
		pq->size += QUEUE_BLOCK_SIZE;
		pq->elements = (TWelement*)realloc(pq->elements, pq->size*sizeof(TWelement));
		if (pq->tail < pq->head)
		{
			pe = pq->elements;
			memcpy((void*)pe+pq->size-QUEUE_BLOCK_SIZE, (void*)pe, sizeof(TWelement)*pq->tail);
		}
	}
	pe[pq->tail].TA = ta;
	pe[pq->tail].TB = tb;
	pe[pq->tail].theta = theta;
	pq->tail = (pq->tail+1)%pq->size;
	pq->length++;
	return 0;
}

int dequeue_tw(TWQ *pq, TWelement *peout)
{
	TWelement *pe = pq->elements+pq->head;

	if (pq->length>0)
	{
		peout->TA = pe->TA;
		peout->TA = pe->TB;
		peout->theta = pe->theta;
		pq->head = (pq->head+1)%pq->head;
		pq->length--;
		return 0;
	}
	else
	{
		return -1;
	}
}

int destroy_queue_tw(TWQ *pq)
{
	if (pq->size > 0)
	{
		free(pq->elements);
		pq->elements = NULL;
		pq->head = pq->tail = -1;
		pq->size = 0;
		pq->length = 0;
	}
}

int EnqueueP_Cell(TWQ *PQ_treewalk, int TA, int TB, double theta);
int ProcessQP_Cell(TWQ *PQ_treewalk, PPQ *PQ_ppnode, int process(int, int, double, PPQ*, TWQ*));

#endif
