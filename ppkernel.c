#include "ppkernel.h"

#ifdef __MULTI_THREAD_
#include "proto.h"
#include "nmk_utils.h"
#include <omp.h>
#define NMK_NAIVE_GRAVITY
#else
#include "global.h"
#endif

#include <math.h>
#include <assert.h>

#define RSQRTPI 0.564189583547756286948

int packarray3(CalcBody* pp, int n, Array3 pa)
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

int packarray3m(CalcBody* pp, int n, Array3 pa, Real *mass)
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
			mass[i] = pp[i].mass;
		}
	}
	else
	{
		for (i=0;i<n;i++)
		{
			pa.x[i] = 0;
			pa.y[i] = 0;
			pa.z[i] = 0;
			mass[i] = 0;
		}
	}
}

int packarray3o(CalcBody* pp, int offset, int n, Array3 pa)
{
	int i;
//	pa.x = pa.y = pa.z = NULL;
	if(pp)
	{
		for (i=0;i<n;i++)
		{
			pa.x[i+offset] = pp[i].pos[0];
			pa.y[i+offset] = pp[i].pos[1];
			pa.z[i+offset] = pp[i].pos[2];
		}
	}
	else
	{
		for (i=0;i<n;i++)
		{
			pa.x[i+offset] = 0;
			pa.y[i+offset] = 0;
			pa.z[i+offset] = 0;
		}
	}
}

int packarray3om(CalcBody* pp, int offset, int n, Array3 pa, Real *mass)
{
	int i;
//	pa.x = pa.y = pa.z = NULL;
	if(pp)
	{
		for (i=0;i<n;i++)
		{
			pa.x[i+offset] = pp[i].pos[0];
			pa.y[i+offset] = pp[i].pos[1];
			pa.z[i+offset] = pp[i].pos[2];
			mass[i+offset] = pp[i].mass;
		}
	}
	else
	{
		for (i=0;i<n;i++)
		{
			pa.x[i+offset] = (Real) 0;
			pa.y[i+offset] = (Real) 0;
			pa.z[i+offset] = (Real) 0;
			mass[i+offset] = (Real) 0;
		}
	}
}

int pusharray3(CalcBody* pp, int n, Array3 pa)
{
	int i;
//	pa.x = pa.y = pa.z = NULL;
	assert(pp);
		
	for (i=0;i<n;i++)
	{
		pp[i].acc[0] += pa.x[i];
		pp[i].acc[1] += pa.y[i];
		pp[i].acc[2] += pa.z[i];
	}
}

int pusharray3m(CalcBody* pp, int n, Array3 pa, Real m)
{
	int i;
//	pa.x = pa.y = pa.z = NULL;
	assert(pp);
		
	for (i=0;i<n;i++)
	{
		pp[i].acc[0] += pa.x[i]*m;
		pp[i].acc[1] += pa.y[i]*m;
		pp[i].acc[2] += pa.z[i]*m;
	}
}

int get_block_tnum(int bid)
{
#ifdef __MULTI_THREAD_
	return omp_get_num_threads();
#else
	return 1;
#endif
}

int get_block_tid(int bid)
{
#ifdef __MULTI_THREAD_
	return omp_get_thread_num();
#else
	return 0;
#endif
}

#define __FUNC_CONCAT(x, y) x##y
#define __PP_SUB_SQ(a, b, c) (a*a + b*b + c*c)
#define __PP_CUBIC(a) (a * a * a)

#ifdef NMK_SINGLE_PREC
#define	__PP_RSQRT(a) ((Real) 1.0) / SQRT(__PP_CUBIC(a))
#else
#define __PP_RSQRT(a)	INVSQRT(__PP_CUBIC(a))
#endif // NMK_SINGLE_PREC

#define _PP_SUB(_d, _dnum, _didx) \
	dx##_dnum = _d.x[_didx] - x2; \
	dy##_dnum = _d.y[_didx] - y2; \
	dz##_dnum = _d.z[_didx] - z2;

#define _PP_SUB_REV(_d, _dnum, _didx) \
	dx##_dnum = x2 - _d.x[_didx]; \
	dy##_dnum = y2 - _d.y[_didx]; \
	dz##_dnum = z2 - _d.z[_didx];

#define PP_SUB_U1(_d, _didx) \
	_PP_SUB(_d, 1, _didx)

#define PP_SUB_U2(_d, _didx) \
	_PP_SUB(_d, 1, _didx) \
	_PP_SUB(_d, 2, _didx+N_CACHE/UNROLL)

#define PP_SUB_U4(_d, _didx) \
	_PP_SUB(_d, 1, _didx) \
	_PP_SUB(_d, 2, _didx+N_CACHE/UNROLL) \
	_PP_SUB(_d, 3, _didx+2*N_CACHE/UNROLL) \
	_PP_SUB(_d, 4, _didx+3*N_CACHE/UNROLL)

#define PP_SUB_REV(_d, _didx) \
	_PP_SUB_REV(_d, 1, _didx)

#define PP_SUB(_u) __FUNC_CONCAT(PP_SUB_U, _u)
#define PP_SUB_NOUNROLL PP_SUB(1)
#define PP_SUB_UNROLL PP_SUB(UNROLL)

#define _PP_DISTANCE_N(_dnum) \
	dr2##_dnum = eps2 + __PP_SUB_SQ(dx##_dnum, dy##_dnum, dz##_dnum);

#define _PP_DISTANCE_F(_dnum) \
	dd##_dnum = __PP_SUB_SQ(dx##_dnum, dy##_dnum, dz##_dnum); \
	dr##_dnum = SQRT(dd##_dnum); \
	dg##_dnum = ERFC(dr##_dnum*inv2rs) + dr##_dnum*invpirs*EXP(dd##_dnum*invrs2); \
	dr2##_dnum = eps2+dd##_dnum;
//	dd##_dnum = __PP_SUB_SQ(dx##_dnum, dy##_dnum, dz##_dnum); \
//	dr2##_dnum = (Real) 1.0/SQRT(eps2+dd##_dnum); \
//	dr##_dnum = (Real) 1.0/dr2##_dnum; \
//	dg##_dnum = ERFC(dr##_dnum*inv2rs) + dr##_dnum*invpirs*EXP(dd##_dnum*invrs2); 

#define _PP_DISTANCE_T(_dnum) \
	dd##_dnum = __PP_SUB_SQ(dx##_dnum, dy##_dnum, dz##_dnum); \
	dr##_dnum = SQRT(dd##_dnum); 

#define _PP_BM_0(_dm, _midx, _exp) _exp
#define _PP_BM_1(_dm, _midx, _exp) _dm[_midx] * (_exp)

#define _PP_GRAV_N(_dnum, _m, _dm, _midx) \
	dr3##_dnum = _PP_BM_##_m(_dm, _midx, __PP_RSQRT(dr2##_dnum));

#define _PP_GRAV_F(_dnum, _m, _dm, _midx) \
	dr3##_dnum = _PP_BM_##_m(_dm, _midx, dg##_dnum * __PP_RSQRT(dr2##_dnum));

#define _PP_GRAV_T(_dnum, _m, _dm, _midx) \
	idx##_dnum = (int)(rad2idx * dr##_dnum); \
	dr3##_dnum = (idx##_dnum<TABLEN) ? _PP_BM_##_m(_dm, _midx, val[idx##_dnum] + ( dr##_dnum - dlt * idx##_dnum ) * slp[idx##_dnum]) : 0 ;

#ifdef NMK_NAIVE_GRAVITY
	#define _PP_DISTANCE _PP_DISTANCE_N
	#define _PP_GRAV _PP_GRAV_N
#else // NMK_NAIVE_GRAVITY
#ifdef NMK_PP_TAB
	#define _PP_DISTANCE _PP_DISTANCE_T
	#define _PP_GRAV _PP_GRAV_T
#else // NMK_PP_TAB
	#define _PP_DISTANCE _PP_DISTANCE_F
	#define _PP_GRAV _PP_GRAV_F
#endif // NMK_PP_TAB
#endif // NMK_NAIVE_GRAVITY

#define PP_COE_U1(_m, _dm, _midx) \
	_PP_DISTANCE(1) \
	_PP_GRAV(1, _m, _dm, _midx)

#define PP_COE_U2(_m, _dm, _midx) \
	_PP_DISTANCE(1) \
	_PP_DISTANCE(2) \
	_PP_GRAV(1, _m, _dm, _midx) \
	_PP_GRAV(2, _m, _dm, _midx+N_CACHE/UNROLL)

#define PP_COE_U4(_m, _dm, _midx) \
	_PP_DISTANCE(1) \
	_PP_DISTANCE(2) \
	_PP_DISTANCE(3) \
	_PP_DISTANCE(4) \
	_PP_GRAV(1, _m, _dm, _midx) \
	_PP_GRAV(2, _m, _dm, _midx+N_CACHE/UNROLL) \
	_PP_GRAV(3, _m, _dm, _midx+2*N_CACHE/UNROLL) \
	_PP_GRAV(4, _m, _dm, _midx+3*N_CACHE/UNROLL)
	
#define PP_COE(_u) __FUNC_CONCAT(PP_COE_U, _u)
#define PP_COE_NOUNROLL PP_COE(1)
#define PP_COE_UNROLL PP_COE(UNROLL)

#define _PP_ACC_U1(_v) _v##1*dr31
#define _PP_ACC_U2(_v) _v##1*dr31 + _v##2*dr32
#define _PP_ACC_U4(_v) _v##1*dr31 + _v##2*dr32 + _v##3*dr33 + _v##4*dr34

#define _PP_ACC_REV(_d, _didx) \
	_d.x[_didx] += _PP_ACC_U1(dx); \
	_d.y[_didx] += _PP_ACC_U1(dy); \
	_d.z[_didx] += _PP_ACC_U1(dz); 

#define PP_ACC_REV _PP_ACC_REV

#define _PP_ACC(_u) __FUNC_CONCAT(_PP_ACC_U, _u)

#define PP_ACC_NOUNROLL \
	ax2 += _PP_ACC(1)(dx); \
	ay2 += _PP_ACC(1)(dy); \
	az2 += _PP_ACC(1)(dz); 

#define PP_ACC_UNROLL \
	ax2 += _PP_ACC(UNROLL)(dx); \
	ay2 += _PP_ACC(UNROLL)(dy); \
	az2 += _PP_ACC(UNROLL)(dz); 




#ifdef __MULTI_THREAD_
int ppkernel(Array3 A, int la, Array3 B, int lb, Real eps2, Array3 C)
#else
int ppkernel(Array3 A, int la, Array3 B, int lb, Constants *constants, Array3 C, int tnum, int tid)
#endif
{
    int n;

#ifndef __MULTI_THREAD_
#ifdef NMK_PP_TAB
	int idx1, idx2, idx3, idx4;
	Real dlt  = constants->delta;
	Real *val = constants->value;
	Real *slp = constants->slope;
	Real rad2idx = TABLEN/(constants->CUTOFF_SCALE);
#else
	Real eps2 = constants->EPS2;
	Real invrs = (Real) 1.0 / constants->SPLIT_SCALE;
	Real invrs2 = (Real) -0.25 * invrs * invrs;
	Real invpirs = (Real) RSQRTPI * invrs;
	Real inv2rs = (Real) 0.5 * invrs;
	Real rc = constants->CUTOFF_SCALE;
#endif // NMK_PP_TAB
#endif

    /* Origin version of ppkernel 
    printf("Starting Compute\n");
    
    tstart = dtime();
    
    for ( n=0; n<la; n++ ) {
        for ( m=0; m<lb; m++) {
            dx1 = B.x[m] -  A.x[n];
            dy1 = B.y[m] -  A.y[n];
            dz1 = B.z[m] -  A.z[n];
            
            dr21 = eps2 + dx*dx + dy*dy + dz*dz;
            dr31 = SQRT(dr2)*dr2;
            
            C.x[n] += dx1/dr31;
            C.y[n] += dy1/dr31;
            C.z[n] += dz1/dr31;
        }
    }

    tstop = dtime();
    
    ttime = tstop - tstart;
    printf("dtime = %lf \n", ttime);
    */


#ifdef __MULTI_THREAD_
    double tstart, tstop, ttime;
    tstart = dtime();

#pragma omp parallel
#endif
	{
		int j, k, m, nb, mb, mt;

		Real x2, y2, z2, ax2, ay2, az2;
		Real dx1, dy1, dz1, dd1, dg1, dr1, dr21, dr31;
		Real dx2, dy2, dz2, dd2, dg2, dr2, dr22, dr32;
		Real dx3, dy3, dz3, dd3, dg3, dr3, dr23, dr33;
		Real dx4, dy4, dz4, dd4, dg4, dr4, dr24, dr34;

#ifdef __MULTI_THREAD_
		int tid, tnum;
		tid = get_block_tid(0);
		tnum = get_block_tnum(0);
#endif
//		printf ("%ld, %d, %ld, %d, %ld, %d, %d.\n", A.x, la, B.x, lb, C.x, tnum, tid);
//		return (0);

#if 1 // Set 0 to debug	
		nb = (la+CLCNT-1)/CLCNT;
		mb = (lb+N_CACHE-1)/N_CACHE;

		if( (nb >= (tnum<<1)) && (lb > CLCNT) ) // The num of A and B are all large enough
		{
//			int nbs = tid*(nb/tnum)+((nb%tnum)+tid)/tnum*((nb+tid)%tnum);
//			int nbt = nbs + (nb+tid)/tnum;
//			int nbs = tid*(nb/tnum)+(nb%tnum-(tnum+(nb%tnum)-tid)/tnum*(nb%tnum-tid));
//			int nbt = nbs + (nb+tnum-1-tid)/tnum;
			int pernb = nb/tnum;
			int modnb = nb%tnum;
			int tpass = modnb-tid;
			int tex = (tnum+tpass-1)/tnum;
			int nbs = tid*pernb+(modnb-tex*tpass);
			int nbt = nbs + pernb + tex;
//			printf ("Thread[%d]: nbs = %d, nbt = %d.\n", tid, nbs, nbt);

			for ( n=nbs; n<nbt; n++ ) {

				int nt = (n==nb-1) ? la : (n+1)*CLCNT;

				for (m=0; m<mb; m++) {
#ifdef __MULTI_THREAD_
					mt = ((m+tid)%mb)*N_CACHE + (((m+tid)%mb==mb-1 && lb%N_CACHE) ? lb-(mb-1)*N_CACHE : N_CACHE/UNROLL);
#else
					mt = m*N_CACHE + ((m==mb-1) ? lb-(mb-1)*N_CACHE : N_CACHE/UNROLL);
#endif

					if ( mt == lb ) // No unrolling at the last B block
					{
//						printf ("Thread[%d]: ms = %d, mt = %d.\n", tid, ((m+tid)%mb)*N_CACHE, mt);
//#pragma prefetch A.x:1:CLCNT
//#pragma prefetch A.y:1:CLCNT
//#pragma prefetch A.z:1:CLCNT
						for (j=n*CLCNT; j<nt; j++)
						{
							x2=A.x[j];
							y2=A.y[j];
							z2=A.z[j];

							ax2=0;
							ay2=0;
							az2=0;

#pragma ivdep
#pragma vector aligned
#ifdef __MULTI_THREAD_
							for (k=((m+tid)%mb)*N_CACHE; k<mt; k++)
#else
							for (k=m*N_CACHE; k<mt; k++)
#endif
							{
								PP_SUB_NOUNROLL(B, k) ;

								PP_COE_NOUNROLL(0, Bm, k) ;

								PP_ACC_NOUNROLL ;
							}

							C.x[j] += ax2;
							C.y[j] += ay2;
							C.z[j] += az2;
						}
					}
					else // Unrolling the loop
					{
//#pragma prefetch A.x:1:CLCNT
//#pragma prefetch A.y:1:CLCNT
//#pragma prefetch A.z:1:CLCNT
						for (j=n*CLCNT; j<nt; j++)
						{
							x2=A.x[j];
							y2=A.y[j];
							z2=A.z[j];

							ax2=0;
							ay2=0;
							az2=0;

#pragma ivdep
#pragma vector aligned
#ifdef __MULTI_THREAD_
							for (k=((m+tid)%mb)*N_CACHE; k<mt; k++)
#else
							for (k=m*N_CACHE; k<mt; k++)
#endif
							{
								PP_SUB_UNROLL(B, k) ;

								PP_COE_UNROLL(0, Bm, k) ;

								PP_ACC_UNROLL ;
							}

							C.x[j] += ax2;
							C.y[j] += ay2;
							C.z[j] += az2;
						}
					}
				}
			}
		}
		else if ( lb > CLCNT )
		{
//			int ns = tid*(la/tnum)+((la%tnum)+tid)/tnum*((la+tid)%tnum);
//			int nt = ns + (la+tid)/tnum;
//			int ns = tid*(la/tnum)+(la%tnum-(tnum+(la%tnum)-tid)/tnum*(la%tnum-tid));
//			int nt = ns + (la+tnum-1-tid)/tnum;
			int pern = la/tnum;
			int modn = la%tnum;
			int tpass = modn-tid;
			int tex = (tnum+tpass-1)/tnum;
			int ns = tid*pern+(modn-tex*tpass);
			int nt = ns + pern + tex;
//			printf ("Thread[%d]: ns = %d, nt = %d.\n", tid, ns, nt);

			for (m=0; m<mb; m++)
			{
#ifdef __MULTI_THREAD_
				mt = ((m+tid)%mb)*N_CACHE + (((m+tid)%mb==mb-1 && lb%N_CACHE) ? lb-(mb-1)*N_CACHE : N_CACHE/UNROLL);
#else
				mt = m*N_CACHE + ((m==mb-1) ? lb-(mb-1)*N_CACHE : N_CACHE/UNROLL);
#endif

				for (j=ns; j<nt; j++)
				{
					x2=A.x[j];
					y2=A.y[j];
					z2=A.z[j];

					ax2=0;
					ay2=0;
					az2=0;

					if ( mt == lb ) // No unrolling at the last B block
					{
#pragma vector aligned
#pragma ivdep
#ifdef __MULTI_THREAD_
						for (k=((m+tid)%mb)*N_CACHE; k<mt; k++)
#else
						for (k=m*N_CACHE; k<mt; k++)
#endif
						{
							PP_SUB_NOUNROLL(B, k) ;

							PP_COE_NOUNROLL(0, Bm, k) ;

							PP_ACC_NOUNROLL ;
						}
					}
					else // Unrolling
					{
#pragma vector aligned
#pragma ivdep
#ifdef __MULTI_THREAD_
						for (k=((m+tid)%mb)*N_CACHE; k<mt; k++)
#else
						for (k=m*N_CACHE; k<mt; k++)
#endif
						{
							PP_SUB_UNROLL(B, k) ;

							PP_COE_UNROLL(0, Bm, k) ;

							PP_ACC_UNROLL ;
						}
					}

					C.x[j] += ax2;
					C.y[j] += ay2;
					C.z[j] += az2;
				}
			}
		}
		else
		{
//			int ns = tid*(la/tnum)+((la%tnum)+tid)/tnum*((la+tid)%tnum);
//			int nt = ns + (la+tid)/tnum;
//			int ns = tid*(la/tnum)+(la%tnum-(tnum+(la%tnum)-tid)/tnum*(la%tnum-tid));
//			int nt = ns + (la+tnum-1-tid)/tnum;
			int pern = la/tnum;
			int modn = la%tnum;
			int tpass = modn-tid;
			int tex = (tnum+tpass-1)/tnum;
			int ns = tid*pern+(modn-tex*tpass);
			int nt = ns + pern + tex;
//			printf ("Thread[%d]: ns = %d, nt = %d.\n", tid, ns, nt);

			for (m=0; m<lb; m++)
			{
				x2 = B.x[m];
				y2 = B.y[m];
				z2 = B.z[m];
#pragma ivdep
				for (k=ns; k<nt; k++)
				{
					PP_SUB_REV(A, k) ;

					PP_COE_UNROLL(0, Bm, m) ;

					PP_ACC_REV(C, k) ;
				}
			}
		}

#else // #if 1 Debug
		// Original ppkernel without mass.
		int pern = la/tnum;
		int modn = la%tnum;
		int tpass = modn-tid;
		int tex = (tnum+tpass-1)/tnum;
		int ns = tid*pern+(modn-tex*tpass);
		int nt = ns + pern + tex;
		for (j=ns; j<nt; j++)
		{
			x2=A.x[j];
			y2=A.y[j];
			z2=A.z[j];

			ax2=0;
			ay2=0;
			az2=0;

			for (k=0; k<lb; k++)
			{
				PP_SUB_NOUNROLL(B, k) ;
				
				PP_COE_NOUNROLL(0, Bm, k) ;

				PP_ACC_NOUNROLL ;
			}

			C.x[j] += ax2;
			C.y[j] += ay2;
			C.z[j] += az2;
		}
#endif // #if 1
	}
#ifdef __MULTI_THREAD_
    tstop = dtime();
    ttime = tstop - tstart;
    printf("dtime = %lf \n", ttime);
#endif
    return 0;
}

#ifdef __MULTI_THREAD_
int ppmkernel(Array3 A, int la, Array3 B, Real *Bm, int lb, Real eps2, Array3 C)
#else
int ppmkernel(Array3 A, int la, Array3 B, Real *Bm, int lb, Constants *constants, Array3 C, int tnum, int tid)
#endif
{
    int n;

#ifndef __MULTI_THREAD_
#ifdef NMK_PP_TAB
	int idx1, idx2, idx3, idx4;
	Real dlt  = constants->delta;
	Real *val = constants->value;
	Real *slp = constants->slope;
	Real rad2idx = TABLEN/(constants->CUTOFF_SCALE);
#else
	Real eps2 = constants->EPS2;
	Real invrs = (Real) 1.0 / constants->SPLIT_SCALE;
	Real invrs2 = (Real) -0.25 * invrs * invrs;
	Real invpirs = (Real) RSQRTPI * invrs;
	Real inv2rs = (Real) 0.5 * invrs;
	Real rc = constants->CUTOFF_SCALE;
#endif // NMK_PP_TAB
#endif

    /* Origin version of ppkernel with mass
    printf("Starting Compute\n");
    
    tstart = dtime();
    
    for ( n=0; n<la; n++ ) {
        for ( m=0; m<lb; m++) {
            dx1 = B.x[m] -  A.x[n];
            dy1 = B.y[m] -  A.y[n];
            dz1 = B.z[m] -  A.z[n];
            
            dr21 = eps2 + dx*dx + dy*dy + dz*dz;
            dr31 = Bm[m] * SQRT(dr2)*dr2;
            
            C.x[n] += dx1/dr31;
            C.y[n] += dy1/dr31;
            C.z[n] += dz1/dr31;
        }
    }

    tstop = dtime();
    
    ttime = tstop - tstart;
    printf("dtime = %lf \n", ttime);
    */


#ifdef __MULTI_THREAD_
    double tstart, tstop, ttime;
    tstart = dtime();

#pragma omp parallel
#endif
	{
		int j, k, m, nb, mb, mt;

		Real x2, y2, z2, ax2, ay2, az2;
		Real dx1, dy1, dz1, dd1, dg1, dr1, dr21, dr31;
		Real dx2, dy2, dz2, dd2, dg2, dr2, dr22, dr32;
		Real dx3, dy3, dz3, dd3, dg3, dr3, dr23, dr33;
		Real dx4, dy4, dz4, dd4, dg4, dr4, dr24, dr34;

#ifdef __MULTI_THREAD_
		int tid, tnum;
		tid = get_block_tid(0);
		tnum = get_block_tnum(0);
#endif

#if 1 // Set 0 to debug	
		nb = (la+CLCNT-1)/CLCNT;
		mb = (lb+N_CACHE-1)/N_CACHE;

		if( (nb >= (tnum<<1)) && (lb > CLCNT) ) // The num of A and B are all large enough
		{
//			int nbs = tid*(nb/tnum)+((nb%tnum)+tid)/tnum*((nb+tid)%tnum);
//			int nbt = nbs + (nb+tid)/tnum;
//			int nbs = tid*(nb/tnum)+(nb%tnum-(tnum+(nb%tnum)-tid)/tnum*(nb%tnum-tid));
//			int nbt = nbs + (nb+tnum-1-tid)/tnum;
			int pernb = nb/tnum;
			int modnb = nb%tnum;
			int tpass = modnb-tid;
			int tex = (tnum+tpass-1)/tnum;
			int nbs = tid*pernb+(modnb-tex*tpass);
			int nbt = nbs + pernb + tex;
//			printf ("Thread[%d]: nbs = %d, nbt = %d.\n", tid, nbs, nbt);

			for ( n=nbs; n<nbt; n++ ) {

				int nt = (n==nb-1) ? la : (n+1)*CLCNT;

				for (m=0; m<mb; m++) {
#ifdef __MULTI_THREAD_
					mt = ((m+tid)%mb)*N_CACHE + (((m+tid)%mb==mb-1 && lb%N_CACHE) ? lb-(mb-1)*N_CACHE : N_CACHE/UNROLL);
#else
					mt = m*N_CACHE + ((m==mb-1) ? lb-(mb-1)*N_CACHE : N_CACHE/UNROLL);
#endif

					if ( mt == lb ) // No unrolling at the last B block
					{
//						printf ("Thread[%d]: ms = %d, mt = %d.\n", tid, ((m+tid)%mb)*N_CACHE, mt);
//#pragma prefetch A.x:1:CLCNT
//#pragma prefetch A.y:1:CLCNT
//#pragma prefetch A.z:1:CLCNT
						for (j=n*CLCNT; j<nt; j++)
						{
							x2=A.x[j];
							y2=A.y[j];
							z2=A.z[j];

							ax2=0;
							ay2=0;
							az2=0;

#pragma ivdep
#pragma vector aligned
#ifdef __MULTI_THREAD_
							for (k=((m+tid)%mb)*N_CACHE; k<mt; k++)
#else
							for (k=m*N_CACHE; k<mt; k++)
#endif
							{
								PP_SUB_NOUNROLL(B, k) ;

								PP_COE_NOUNROLL(1, Bm, k) ;

								PP_ACC_NOUNROLL ;
							}

							C.x[j] += ax2;
							C.y[j] += ay2;
							C.z[j] += az2;
						}
					}
					else // Unrolling the loop
					{
//#pragma prefetch A.x:1:CLCNT
//#pragma prefetch A.y:1:CLCNT
//#pragma prefetch A.z:1:CLCNT
						for (j=n*CLCNT; j<nt; j++)
						{
							x2=A.x[j];
							y2=A.y[j];
							z2=A.z[j];

							ax2=0;
							ay2=0;
							az2=0;

#pragma ivdep
#pragma vector aligned
#ifdef __MULTI_THREAD_
							for (k=((m+tid)%mb)*N_CACHE; k<mt; k++)
#else
							for (k=m*N_CACHE; k<mt; k++)
#endif
							{
								PP_SUB_UNROLL(B, k) ;

								PP_COE_UNROLL(1, Bm, k) ;

								PP_ACC_UNROLL ;
							}

							C.x[j] += ax2;
							C.y[j] += ay2;
							C.z[j] += az2;
						}
					}
				}
			}
		}
		else if ( lb > CLCNT )
		{
//			int ns = tid*(la/tnum)+((la%tnum)+tid)/tnum*((la+tid)%tnum);
//			int nt = ns + (la+tid)/tnum;
//			int ns = tid*(la/tnum)+(la%tnum-(tnum+(la%tnum)-tid)/tnum*(la%tnum-tid));
//			int nt = ns + (la+tnum-1-tid)/tnum;
			int pern = la/tnum;
			int modn = la%tnum;
			int tpass = modn-tid;
			int tex = (tnum+tpass-1)/tnum;
			int ns = tid*pern+(modn-tex*tpass);
			int nt = ns + pern + tex;
//			printf ("Thread[%d]: ns = %d, nt = %d.\n", tid, ns, nt);

			for (m=0; m<mb; m++)
			{
#ifdef __MULTI_THREAD_
				mt = ((m+tid)%mb)*N_CACHE + (((m+tid)%mb==mb-1 && lb%N_CACHE) ? lb-(mb-1)*N_CACHE : N_CACHE/UNROLL);
#else
				mt = m*N_CACHE + ((m==mb-1) ? lb-(mb-1)*N_CACHE : N_CACHE/UNROLL);
#endif

				for (j=ns; j<nt; j++)
				{
					x2=A.x[j];
					y2=A.y[j];
					z2=A.z[j];

					ax2=0;
					ay2=0;
					az2=0;

					if ( mt == lb ) // No unrolling at the last B block
					{
#pragma vector aligned
#pragma ivdep
#ifdef __MULTI_THREAD_
						for (k=((m+tid)%mb)*N_CACHE; k<mt; k++)
#else
						for (k=m*N_CACHE; k<mt; k++)
#endif
						{
							PP_SUB_NOUNROLL(B, k) ;

							PP_COE_NOUNROLL(1, Bm, k) ;

							PP_ACC_NOUNROLL ;
						}
					}
					else // Unrolling
					{
#pragma vector aligned
#pragma ivdep
#ifdef __MULTI_THREAD_
						for (k=((m+tid)%mb)*N_CACHE; k<mt; k++)
#else
						for (k=m*N_CACHE; k<mt; k++)
#endif
						{
							PP_SUB_UNROLL(B, k) ;

							PP_COE_UNROLL(1, Bm, k) ;

							PP_ACC_UNROLL ;
						}
					}

					C.x[j] += ax2;
					C.y[j] += ay2;
					C.z[j] += az2;
				}
			}
		}
		else
		{
//			int ns = tid*(la/tnum)+((la%tnum)+tid)/tnum*((la+tid)%tnum);
//			int nt = ns + (la+tid)/tnum;
//			int ns = tid*(la/tnum)+(la%tnum-(tnum+(la%tnum)-tid)/tnum*(la%tnum-tid));
//			int nt = ns + (la+tnum-1-tid)/tnum;
			int pern = la/tnum;
			int modn = la%tnum;
			int tpass = modn-tid;
			int tex = (tnum+tpass-1)/tnum;
			int ns = tid*pern+(modn-tex*tpass);
			int nt = ns + pern + tex;
//			printf ("Thread[%d]: ns = %d, nt = %d.\n", tid, ns, nt);

			for (m=0; m<lb; m++)
			{
				x2 = B.x[m];
				y2 = B.y[m];
				z2 = B.z[m];
#pragma ivdep
				for (k=ns; k<nt; k++)
				{
					PP_SUB_REV(A, k) ;

					PP_COE_UNROLL(1, Bm, m) ;

					PP_ACC_REV(C, k) ;
				}
			}
		}

#else // #if 1 Debug
		// Original ppkernel without mass.
		int pern = la/tnum;
		int modn = la%tnum;
		int tpass = modn-tid;
		int tex = (tnum+tpass-1)/tnum;
		int ns = tid*pern+(modn-tex*tpass);
		int nt = ns + pern + tex;
		for (j=ns; j<nt; j++)
		{
			x2=A.x[j];
			y2=A.y[j];
			z2=A.z[j];

			ax2=0;
			ay2=0;
			az2=0;

			for (k=0; k<lb; k++)
			{
				PP_SUB_NOUNROLL(B, k) ;
				
				PP_COE_NOUNROLL(1, Bm, k) ;

				PP_ACC_NOUNROLL ;
			}

			C.x[j] += ax2;
			C.y[j] += ay2;
			C.z[j] += az2;
		}
#endif // #if 1
	}

#ifdef __MULTI_THREAD_
    tstop = dtime();
    ttime = tstop - tstart;
    printf("dtime = %lf \n", ttime);
#endif

    return 0;
}


