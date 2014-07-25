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
#include "queue.h"

static Body* part;
static Node* tree;

void initGlobal(Body* p, Node* t)
{
	part = p;
	tree = t;
}

int init_queue_pp(PPQ *pq)
{
	pq->size = QUEUE_BLOCK_SIZE;
	pq->length = 0;
	pq->head = 0;
	pq->tail = 0;
	pq->elements = (PPelement*)malloc(sizeof(PPelement)*pq->size);
	return 0;
}

int enqueue_pp(PPQ *pq, int ta, int tb, int m)
{
	PPelement *pe=pq->elements;

	if (pq->length == pq->size-256)
	{
		printf("Enlarge size Q1=%d, TA=%d, TB=%d, head=%d, tail=%d.\n", pq->size, ta, tb, pq->head, pq->tail);
		pq->size += QUEUE_BLOCK_SIZE;
		pq->elements = (PPelement*)realloc(pq->elements, pq->size*sizeof(PPelement));
		pe = pq->elements;
		if (pq->tail < pq->head)
		{
			memcpy((void*)pe+pq->size-QUEUE_BLOCK_SIZE, (void*)pe, sizeof(PPelement)*pq->tail);
			pq->tail=pq->head+pq->length;
		}
		printf("Enlarge size Q2=%d, TA=%d, TB=%d, head=%d, tail=%d.\n", pq->size, ta, tb, pq->head, pq->tail);
	}
	pe[pq->tail].TA = ta;
	pe[pq->tail].TB = tb;
	pe[pq->tail].mask = m;
	pq->tail = (pq->tail+1)%pq->size;
	pq->length++;
	return 0;
}

int dequeue_pp(PPQ *pq, PPelement *peout)
{
	PPelement *pe = pq->elements+pq->head;

	if (pq->length>0)
	{
		peout->TA = pe->TA;
		peout->TB = pe->TB;
		peout->mask = pe->mask;
		pq->head = (pq->head+1)%pq->size;
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


int init_queue_tw(TWQ *pq)
{
	pq->size = QUEUE_BLOCK_SIZE;
	pq->length = 0;
	pq->head = 0;
	pq->tail = 0;
	pq->elements = (TWelement*)malloc(sizeof(TWelement)*pq->size);
	return 0;
}

int enqueue_tw(TWQ *pq, int ta, int tb, double theta)
{
	TWelement *pe=pq->elements;

	if (pq->length == pq->size-256)
	{
		printf("Enlarge size0=%d, TA=%d, nPa=%d, TB=%d, nPb=%d, head=%d, tail=%d.\n", pq->size, ta, tree[ta].nPart, tb, tree[tb].nPart, pq->head, pq->tail);
		pq->size += QUEUE_BLOCK_SIZE;
		pq->elements = (TWelement*)realloc(pq->elements, pq->size*sizeof(TWelement));
		pe = pq->elements;
		if (pq->tail < pq->head)
		{
			memcpy((void*)pe+pq->size-QUEUE_BLOCK_SIZE, (void*)pe, sizeof(TWelement)*pq->tail);
			pq->tail = pq->head+pq->length;
		}
		printf("Enlarge size1=%d, TA=%d, nPa=%d, TB=%d, nPb=%d, head=%d, tail=%d.\n", pq->size, ta, tree[ta].nPart, tb, tree[tb].nPart, pq->head, pq->tail);
	}
	pe[pq->tail].TA = ta;
	pe[pq->tail].TB = tb;
	pe[pq->tail].theta = theta;
	pq->tail = (pq->tail+1)%pq->size;
//	if(pq->length==0) printf("Enqueue TA=%d, nPa=%d, TB=%d, nPb=%d.\n", ta, tree[ta].nPart, tb, tree[tb].nPart);
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
		peout->theta = pe->theta;
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


packarray3(Body* pp, int n, Array3 pa)
{
	int i;
//	pa.x = pa.y = pa.z = NULL;
	if(pp)
	{
		for (i=0;i<n;i++)
		{
			pa.x[i] = pp[i].pos[0];
			pa.y[i] = pp[i].pos[1];
			pa.z[i] = pp[i].pos[2];
		}
	}
	else
	{
		for (i=0;i<n;i++)
		{
			pa.x[i] = 0;
			pa.y[i] = 0;
			pa.z[i] = 0;
		}
	}
}

int EnqueueP_PPnode(PPQ *PQ_ppnode, int TA, int TB, int mask)
{
	enqueue_pp(PQ_ppnode, TA, TB, mask);
}

int ProcessQP_PPnode(PPQ *PQ_ppnode, int process(Array3, int, Array3, int, PRECTYPE, Array3), Array3 pA, Array3 pB, Array3 pC)
{
	PPelement el;
	int nA, nB;

	dequeue_pp(PQ_ppnode, &el);
//	printf("length=%d, TA=%d, nA=%d, TB=%d, nB=%d, %.2f\n", PQ_ppnode->length, el.TA, tree[el.TA].nPart, el.TB, tree[el.TB].nPart, pC.x[0]);
//	sleep(3);

	nA = tree[el.TA].nPart;
	packarray3(&part[tree[el.TA].firstpart], nA, pA);
	packarray3(NULL, nA, pC);

	if (el.mask == 0)
	{
		nB = tree[el.TB].nPart;
		packarray3(&part[tree[el.TB].firstpart], nB, pB);
	}
	else
	{
		nB = 1;
		pB.x[0] = tree[el.TB].masscenter[0];
		pB.y[0] = tree[el.TB].masscenter[1];
		pB.z[0] = tree[el.TB].masscenter[2];
//		pB.x=pB.y=pB.z=NULL;	
	}
	process(pA, nA, pB, nB, EPS2, pC);
}

int EnqueueP_Cell(TWQ *PQ_treewalk, int TA, int TB, double theta)
{
	enqueue_tw(PQ_treewalk, TA, TB, theta);
//	printf("e");
}

int ProcessQP_Cell(TWQ *PQ_treewalk, PPQ *PQ_ppnode, int process(int, int, double, PPQ*, TWQ*))
{
	TWelement el;
	dequeue_tw(PQ_treewalk, &el);
//	printf("w");
	process(el.TA, el.TB, el.theta, PQ_ppnode, PQ_treewalk);
}
