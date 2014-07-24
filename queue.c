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

Array3 packarray3(Body* pp, int n)
{
	int i;
	Array3 pa;
	pa.x = (PRECTYPE*)malloc(sizeof(PRECTYPE)*n);
	pa.y = (PRECTYPE*)malloc(sizeof(PRECTYPE)*n);
	pa.z = (PRECTYPE*)malloc(sizeof(PRECTYPE)*n);
	for (i=1; i<n; i*=4)
	{
		printf("z");
	}
	return pa;
}

int EnqueueP_PPnode(PPQ *PQ_ppnode, int TA, int TB, int mask)
{
	Array3 pA, pB, pC;
	int nA, nB;

	nA = tree[TA].nPart;
	pA = packarray3(&part[tree[TA].firstpart], nA);
	pC = packarray3(NULL, nA);

	if (mask == 0)
	{
		nB = tree[TB].nPart;
		pB = packarray3(&part[tree[TB].firstpart], nB);
	}
	else
	{
		nB = 1;
		pB.x = (PRECTYPE*)malloc(sizeof(PRECTYPE));
		pB.x[0] = tree[TB].masscenter[0];
		pB.y[0] = tree[TB].masscenter[1];
		pB.z[0] = tree[TB].masscenter[2];	
	}
	enqueue_pp(PQ_ppnode, pA, pB, pC, nA, nB);
}

int ProcessQP_PPnode(PPQ *PQ_ppnode, int process(Array3, int, Array3, int, PRECTYPE, Array3))
{
	PPelement el;
	dequeue_pp(PQ_ppnode, &el);
	printf("p");
	process(el.A, el.nA, el.B, el.nB, EPS2, el.C);
	free(el.A.x);
	free(el.A.y);
	free(el.A.z);
	free(el.B.x);
	free(el.B.y);
	free(el.B.z);
	free(el.C.x);
	free(el.C.y);
	free(el.C.z);
}

int EnqueueP_Cell(TWQ *PQ_treewalk, int TA, int TB, double theta)
{
	enqueue_tw(PQ_treewalk, TA, TB, theta);
	printf("e");
}

int ProcessQP_Cell(TWQ *PQ_treewalk, PPQ *PQ_ppnode, int process(int, int, double, PPQ*, TWQ*))
{
	TWelement el;
	dequeue_tw(PQ_treewalk, &el);
	process(el.TA, el.TB, el.theta, PQ_ppnode, PQ_treewalk);
	printf("w");
}
