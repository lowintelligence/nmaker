/*
 * =====================================================================================
 *
 *       Filename:  ppkernel.h
 *
 *    Description:  Header files of ppkernel.c
 *
 *        Version:  1.0
 *        Created:  07/19/2014 04:00:10 PM
 *       Revision:  none
 *       Compiler:  icc
 *
 *         Author:  Cao Zongyan (), zycao@sccas.cn
 *        Company:  Supercomputing Center, CNIC, CAS, China
 *
 * =====================================================================================
 */

#define ALIGNCNT 64
#define __single_prec

#ifdef __single_prec
    #define PRECTYPE float
    #define SQRT sqrtf
    #define INVSQRT invsqrtf
    #define N_CACHE 1024
    #define CLCNT 64
    #define UNROLL 4
#else
    #define PRECTYPE double
    #define SQRT sqrt
    #define INVSQRT invsqrt
    #define N_CACHE 512
    #define CLCNT 16
    #define UNROLL 2
#endif

#ifdef _OPENMP
	#define __MULTI_THREAD_
#endif

typedef struct {
	PRECTYPE *x;
	PRECTYPE *y;
	PRECTYPE *z;
} array3;

typedef struct {
	array3 pos;
	array3 acc;
} PPPack;

// To be finished later if needed.
//int get_global_tnum();
//int get_global_tid();
//int get_local_tnum();
//int get_local_tid();

int get_block_tnum(int bid);
int get_block_tid(int bid);

int ppkernel(PPPack A, int la, PPPack B, int lb, PRECTYPE eps2);
