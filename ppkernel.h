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

#ifdef _OPENMP
	#define __MULTI_THREAD_
#endif

#ifndef __MULTI_THREAD_
#include "global.h"
#else
#include "proto.h"
#include "nmk_utils.h"
#endif

#define ALIGNCNT 64

#ifdef NMK_SINGLE_PREC
    #define N_CACHE 1024
    #define CLCNT 16
    #define UNROLL 4
#else
    #define N_CACHE 512
    #define CLCNT 16
    #define UNROLL 2
#endif

typedef struct {
	Real *x;
	Real *y;
	Real *z;
} Array3;

typedef struct {
    Array3 pos;
    Array3 acc;
	Real *mass;
} PPPack;

typedef struct {
    Array3 pa;
    Array3 pb;
    Array3 pc;
	Real *mass;
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

__OffloadFunc_Macro__
int packarray3(CalcBody* pp, int n, Array3 pa);

__OffloadFunc_Macro__
int packarray3m(CalcBody* pp, int n, Array3 pa, Real *mass);

__OffloadFunc_Macro__
int packarray3o(CalcBody* pp, int offset, int n, Array3 pa);

__OffloadFunc_Macro__
int packarray3om(CalcBody* pp, int offset, int n, Array3 pa, Real *mass);

__OffloadFunc_Macro__
int pusharray3(CalcBody* pp, int n, Array3 pa);

__OffloadFunc_Macro__
int pusharray3m(CalcBody* pp, int n, Array3 pa, Real m);

#ifdef __MULTI_THREAD_
__OffloadFunc_Macro__
int ppkernel(Array3 A, int la, Array3 B, int lb, Real eps2, Array3 C);

__OffloadFunc_Macro__
int ppmkernel(Array3 A, int la, Array3 B, Real *Bm, int lb, Real eps2, Array3 C);
#else // __MULTI_THREAD_
__OffloadFunc_Macro__
int ppkernel(Array3 A, int la, Array3 B, int lb, Constants *constants, Array3 C, int tnum, int tid);

__OffloadFunc_Macro__
int ppmkernel(Array3 A, int la, Array3 B, Real *Bm, int lb, Constants *constants, Array3 C, int tnum, int tid);
#endif // __MULTI_THREAD_

#endif
