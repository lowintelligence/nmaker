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

#ifdef __INTEL_OFFLOAD
__declspec(target(mic))
#endif
typedef struct _block
{
	int blockid;
	int bsize;
	int tid;
	int tall;
	double *A;
	double *B;
	pthread_barrier_t *bbar;
	pthread_barrier_t *bglobal;
} Block;


