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

#ifndef _PPKERNEL_H_
#define _PPKERNEL_H_

#include "offload.h"

#define ALIGNCNT 64
#define __single_prec
#define EPS2 0.00026

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
} Array3;

typedef struct {
    Array3 pos;
    Array3 acc;
	PRECTYPE *mass;
} PPPack;

typedef struct {
    Array3 pa;
    Array3 pb;
    Array3 pc;
	PRECTYPE *mass;
	int nA;
	int nB;
	int finish;
	int calcm;
} PPParameter;
// To be finished later if needed.
//int get_global_tnum();
//int get_global_tid();
//int get_local_tnum();
//int get_local_tid();

__OffloadFunc_Macro__
int get_block_tnum(int bid);

__OffloadFunc_Macro__
int get_block_tid(int bid);

#ifdef __MULTI_THREAD_
__OffloadFunc_Macro__
int ppkernel(Array3 A, int la, Array3 B, int lb, PRECTYPE eps2, Array3 C);

__OffloadFunc_Macro__
int ppmkernel(Array3 A, int la, Array3 B, PRECTYPE *Bm, int lb, PRECTYPE eps2, Array3 C);
#else
__OffloadFunc_Macro__
int ppkernel(Array3 A, int la, Array3 B, int lb, PRECTYPE eps2, Array3 C, int tnum, int tid);

__OffloadFunc_Macro__
int ppmkernel(Array3 A, int la, Array3 B, PRECTYPE *Bm, int lb, PRECTYPE eps2, Array3 C, int tnum, int tid);
#endif

#endif
