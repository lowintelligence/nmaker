/*
 * =====================================================================================
 *
 *       Filename:  offload.h
 *
 *    Description:  Offload controller head file
 *
 *        Version:  1.0
 *        Created:  07/24/2014 02:23:34 PM
 *       Revision:  none
 *       Compiler:  icc
 *
 *         Author:  Cao Zongyan (), zycao@sccas.cn
 *        Company:  Supercomputing Center, CNIC, CAS, China
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>

#define NTEAM   16
#define MASTER  1
#define NSLAVE  2
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


