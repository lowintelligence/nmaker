/*
 * =====================================================================================
 *
 *       Filename:  queue.c
 *
 *    Description:  Queue interface and scheduler
 *
 *        Version:  1.0
 *        Created:  07/24/2014 02:24:46 PM
 *       Revision:  none
 *       Compiler:  icc
 *
 *         Author:  Cao Zongyan (), zycao@sccas.cn
 *        Company:  Supercomputing Center, CNIC, CAS, China
 *
 * =====================================================================================
 */
#include "global.h"
#include "queue.h"

__OffloadVar_Macro__
static CalcBody* part;

__OffloadVar_Macro__
static Node* tree;

void initGlobal(CalcBody* p, Node* t)
{
	part = p;
	tree = t;
}

int init_queues_tw(TWQ *p_pq, int n)
{
	int i;
	for(i=0;i<n;i++)
	{
		p_pq[i].size = QUEUE_BLOCK_SIZE;
		p_pq[i].length = 0;
		p_pq[i].head = 0;
		p_pq[i].tail = 0;
		p_pq[i].elements = (TWelement*)xmalloc(sizeof(TWelement)*p_pq[i].size, 9101);
	}
	return 0;
}

int enqueue_tw(TWQ *pq, Int ta, Int tb)
{
	TWelement *pe=pq->elements;

	if (pq->length == pq->size-256)
	{
		DBG_INFO(DBG_MEMORY, "Enlarge Q0 s=%d, l=%d, TA=%d, nPa=%d, TB=%d, nPb=%d, hd=%d, tl=%d.\n", pq->size, pq->length, ta, tree[ta].nPart, tb, tree[tb].nPart, pq->head, pq->tail);
		pq->size += QUEUE_BLOCK_SIZE;
		pq->elements = (TWelement*)xrealloc(pq->elements, pq->size*sizeof(TWelement), 9102);
		pe = pq->elements;
		if (pq->tail < pq->head)
		{
			if (pq->tail>QUEUE_BLOCK_SIZE)
			{
				memcpy((void*)(pe+pq->size-QUEUE_BLOCK_SIZE), (void*)pe, sizeof(TWelement)*QUEUE_BLOCK_SIZE);
				memcpy((void*)pe, (void*)(pe+QUEUE_BLOCK_SIZE), sizeof(TWelement)*(pq->tail-QUEUE_BLOCK_SIZE));
			}
			else
			{
				memcpy((void*)(pe+pq->size-QUEUE_BLOCK_SIZE), (void*)pe, sizeof(TWelement)*pq->tail);
			}
			pq->tail=(pq->head+pq->length)%pq->size;
		}
		DBG_INFO(DBG_MEMORY, "Enlarge Q1 s=%d, l=%d, TA=%d, nPa=%d, TB=%d, nPb=%d, hd=%d, tl=%d.\n", pq->size, pq->length, ta, tree[ta].nPart, tb, tree[tb].nPart, pq->head, pq->tail);
	}
	pe[pq->tail].TA = ta;
	pe[pq->tail].TB = tb;
	pq->tail = (pq->tail+1)%pq->size;
//	printf("Enqueue TA=%d, nPa=%d, TB=%d, nPb=%d.\n", ta, tree[ta].nPart, tb, tree[tb].nPart);
	pq->length++;
	return 0;
}

int dequeue_tw(TWQ *pq, TWelement *peout)
{
	TWelement *pe = pq->elements+pq->head;

//	printf("w head=%d, tail=%d \n", pq->head, pq->tail);
	if (pq->length>0)
	{
		peout->TA = pe->TA;
		peout->TB = pe->TB;
		pq->head = (pq->head+1)%pq->size;
		pq->length--;
//		printf("w head=%d, tail=%d \n", pq->head, pq->tail);
//		printf("Dequeue TA=%d, nPa=%d, TB=%d, nPb=%d.\n", peout->TA, tree[peout->TA].nPart, peout->TB, tree[peout->TB].nPart);
		return 0;
	}
	else
	{
		return -1;
	}
}

int destroy_queues_tw(TWQ *pq, int n)
{
	int i;
	for (i=0; i<n; i++)
	{
	if (pq[i].size > 0)
	{
		free(pq[i].elements);
		pq[i].elements = NULL;
		pq[i].head = pq[i].tail = -1;
		pq[i].size = 0;
		pq[i].length = 0;
	}
	}
}

int process_queue_tw(TWQ *PQ_treewalk, void *param, int process(Int, Int, TWQ*, void*))
{
	TWelement el;
	dequeue_tw(PQ_treewalk, &el);
//	printf("w");
	process(el.TA, el.TB, PQ_treewalk, param);
}
