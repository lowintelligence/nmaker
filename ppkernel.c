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

#define RSQRTPI 1.772453850905516027298

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

#ifdef __MULTI_THREAD_
int ppkernel(Array3 A, int la, Array3 B, int lb, Real eps2, Array3 C)
#else
int ppkernel(Array3 A, int la, Array3 B, int lb, Constants *constants, Array3 C, int tnum, int tid)
#endif
{
    int n;

#ifndef __MULTI_THREAD_
	Real eps2 = constants->EPS2;
	Real invrs = (Real) 1.0 / constants->SPLIT_SCALE;
	Real invrs2 = (Real) -0.5 * invrs * invrs;
	Real invpirs = (Real) RSQRTPI * invrs;
	Real inv2rs = (Real) 0.5 * invrs;
	Real rc = constants->CUTOFF_SCALE;
	Real G = constants->GRAV_CONST;
#endif

#ifdef TABLEN
    int m;
    int idx;
    Real dlt  = constants->delta;
    Real *val = constants->value;
    Real *slp = constants->slope;
    Real dx1, dy1, dz1, dr21;
    Real rad2idx = TABLEN/(constants->CUTOFF_SCALE);


    for ( n=0; n<la; n++ ) {
        for ( m=0; m<lb; m++) {
            dx1 = B.x[m] -  A.x[n];
            dy1 = B.y[m] -  A.y[n];
            dz1 = B.z[m] -  A.z[n];
            
            dr21 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
            idx = (int)(rad2idx*dr21);
            if (idx<TABLEN) {
                dr21 = val[idx] + ( dr21 - dlt * idx ) * slp[idx];
            }
            else {
                dr21 = 0.0;                
            }                     
            C.x[n] += dx1*dr21;
            C.y[n] += dy1*dr21;
            C.z[n] += dz1*dr21;

        }
    }
    return 0;
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
		int j, k, m, nb, mb, nt, mt;

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

#if 1 // Set 0 to Debug	
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

				nt = (n==nb-1) ? la : (n+1)*CLCNT;

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
								dx1 = B.x[k] - x2;
								dy1 = B.y[k] - y2;
								dz1 = B.z[k] - z2;

#ifdef NMK_NAIVE_GRAVITY
								dr21 = eps2 + dx1*dx1 + dy1*dy1 + dz1*dz1;
#ifdef NMK_SINGLE_PREC
								dr31 = ((Real) 1.0) / SQRT(dr21 * dr21 * dr21);
#else
								dr31 = INVSQRT(dr21 * dr21 * dr21);
#endif

#else // NMK_NAIVE_GRAVITY
								dd1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
								dr1 = SQRT(dd1);
								dg1 = ERFC(dr1*inv2rs) + dr1*invpirs*EXP(dd1*invrs2);
								dr21 = eps2+dd1;
#ifdef NMK_SINGLE_PREC
								dr31 = dg1 * ((Real) 1.0) / SQRT(dr21 * dr21 * dr21);
#else
								dr31 = dg1 * INVSQRT(dr21 * dr21 * dr21);
#endif

#endif // NMK_NAIVE_GRAVITY

								ax2 += dx1*dr31 ;
								ay2 += dy1*dr31 ;
								az2 += dz1*dr31 ;
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
								dx1 = B.x[k] - x2;
								dy1 = B.y[k] - y2;
								dz1 = B.z[k] - z2;
#if (UNROLL > 1)
								dx2 = B.x[k+N_CACHE/UNROLL] - x2;
								dy2 = B.y[k+N_CACHE/UNROLL] - y2;
								dz2 = B.z[k+N_CACHE/UNROLL] - z2;
#if (UNROLL >= 4)
								dx3 = B.x[k+2*N_CACHE/UNROLL] - x2;
								dy3 = B.y[k+2*N_CACHE/UNROLL] - y2;
								dz3 = B.z[k+2*N_CACHE/UNROLL] - z2;

								dx4 = B.x[k+3*N_CACHE/UNROLL] - x2;
								dy4 = B.y[k+3*N_CACHE/UNROLL] - y2;
								dz4 = B.z[k+3*N_CACHE/UNROLL] - z2;
#endif				
#endif

#ifdef NMK_NAIVE_GRAVITY
								dr21 = eps2 + dx1*dx1 + dy1*dy1 + dz1*dz1;
#if (UNROLL > 1)                                     
								dr22 = eps2 + dx2*dx2 + dy2*dy2 + dz2*dz2;
			//					dr21 = invsqrtf(eps2 + dx1*dx1 + dy1*dy1 + dz1*dz1);
			//					dr22 = invsqrtf(eps2 + dx2*dx2 + dy2*dy2 + dz2*dz2);
#if (UNROLL >= 4)
								dr23 = eps2 + dx3*dx3 + dy3*dy3 + dz3*dz3;
								dr24 = eps2 + dx4*dx4 + dy4*dy4 + dz4*dz4;
			//					dr23 = invsqrtf(eps2 + dx3*dx3 + dy3*dy3 + dz3*dz3);
			//					dr24 = invsqrtf(eps2 + dx4*dx4 + dy4*dy4 + dz4*dz4);
#endif
#endif
#ifdef NMK_SINGLE_PREC
								dr31 = ((Real) 1.0) / SQRT(dr21 * dr21 * dr21);
#if (UNROLL > 1)                                                  
								dr32 = ((Real) 1.0) / SQRT(dr22 * dr22 * dr22);
#if (UNROLL >= 4)                                                 
								dr33 = ((Real) 1.0) / SQRT(dr23 * dr23 * dr23);
																			  
								dr34 = ((Real) 1.0) / SQRT(dr24 * dr24 * dr24);
#endif
#endif
#else
								dr31 = INVSQRT(dr21 * dr21 * dr21);
#if (UNROLL > 1)
								dr32 = INVSQRT(dr22 * dr22 * dr22);
#if (UNROLL >= 4)
								dr33 = INVSQRT(dr23 * dr23 * dr23);

								dr34 = INVSQRT(dr24 * dr24 * dr24);
#endif
#endif
#endif

#else // NMK_NAIVE_GRAVITY
								dd1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
								dr1 = SQRT(dd1);
								dg1 = ERFC(dr1*inv2rs) + dr1*invpirs*EXP(dd1*invrs2);
								dr21 = eps2+dd1;
#if (UNROLL > 1)
								dd2 = dx2*dx2 + dy2*dy2 + dz2*dz2;
								dr2 = SQRT(dd2);
								dg2 = ERFC(dr2*inv2rs) + dr2*invpirs*EXP(dd2*invrs2);
								dr22 = eps2+dd2;
#if (UNROLL >=4)
								dd3 = dx3*dx3 + dy3*dy3 + dz3*dz3;
								dr3 = SQRT(dd3);
								dg3 = ERFC(dr3*inv2rs) + dr3*invpirs*EXP(dd3*invrs2);
								dr23 = eps2+dd3;

								dd4 = dx4*dx4 + dy4*dy4 + dz4*dz4;
								dr4 = SQRT(dd4);
								dg4 = ERFC(dr4*inv2rs) + dr4*invpirs*EXP(dd4*invrs2);
								dr24 = eps2+dd4;
#endif
#endif

#ifdef NMK_SINGLE_PREC
								dr31 = dg1 * ((Real) 1.0) / SQRT(dr21 * dr21 * dr21);
#if (UNROLL > 1)                                                        
								dr32 = dg2 * ((Real) 1.0) / SQRT(dr22 * dr22 * dr22);
#if (UNROLL >= 4)                                                       
								dr33 = dg3 * ((Real) 1.0) / SQRT(dr23 * dr23 * dr23);
									          									  
								dr34 = dg4 * ((Real) 1.0) / SQRT(dr24 * dr24 * dr24);
#endif
#endif
#else // NMK_SINGLE_PREC
								dr31 = dg1 * INVSQRT(dr21 * dr21 * dr21);
#if (UNROLL > 1)                             
								dr32 = dg2 * INVSQRT(dr22 * dr22 * dr22);
#if (UNROLL >= 4)                            
								dr33 = dg3 * INVSQRT(dr23 * dr23 * dr23);
                                             
								dr34 = dg4 * INVSQRT(dr24 * dr24 * dr24);
#endif
#endif
#endif // NMK_SINGLE_PREC

#endif // NMK_NAIVE_GRAVITY

#if (UNROLL == 4)
								ax2 += dx1*dr31 + dx2*dr32 + dx3*dr33 + dx4*dr34;
								ay2 += dy1*dr31 + dy2*dr32 + dy3*dr33 + dy4*dr34;
								az2 += dz1*dr31 + dz2*dr32 + dz3*dr33 + dz4*dr34;
#endif
#if (UNROLL == 2) 
								ax2 += dx1*dr31 + dx2*dr32;
								ay2 += dy1*dr31 + dy2*dr32;
								az2 += dz1*dr31 + dz2*dr32;
#endif
#if (UNROLL == 1)
								ax2 += dx1*dr31 ;
								ay2 += dy1*dr31 ;
								az2 += dz1*dr31 ;
#endif
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
							dx1 = B.x[k] - x2;
							dy1 = B.y[k] - y2;
							dz1 = B.z[k] - z2;

#ifdef NMK_NAIVE_GRAVITY
							dr21 = eps2 + dx1*dx1 + dy1*dy1 + dz1*dz1;
#ifdef NMK_SINGLE_PREC
							dr31 = ((Real) 1.0) / SQRT(dr21 * dr21 * dr21);
#else
							dr31 = INVSQRT(dr21 * dr21 * dr21);
#endif

#else // NMK_NAIVE_GRAVITY
							dd1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
							dr1 = SQRT(dd1);
							dg1 = ERFC(dr1*inv2rs) + dr1*invpirs*EXP(dd1*invrs2);
							dr21 = eps2+dd1;
#ifdef NMK_SINGLE_PREC
							dr31 = dg1 * ((Real) 1.0) / SQRT(dr21 * dr21 * dr21);
#else
							dr31 = dg1 * INVSQRT(dr21 * dr21 * dr21);
#endif

#endif // NMK_NAIVE_GRAVITY

							ax2 += dx1*dr31 ;
							ay2 += dy1*dr31 ;
							az2 += dz1*dr31 ;
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
							dx1 = B.x[k] - x2;
							dy1 = B.y[k] - y2;
							dz1 = B.z[k] - z2;
#if (UNROLL > 1)
							dx2 = B.x[k+N_CACHE/UNROLL] - x2;
							dy2 = B.y[k+N_CACHE/UNROLL] - y2;
							dz2 = B.z[k+N_CACHE/UNROLL] - z2;
#if (UNROLL >= 4)
							dx3 = B.x[k+2*N_CACHE/UNROLL] - x2;
							dy3 = B.y[k+2*N_CACHE/UNROLL] - y2;
							dz3 = B.z[k+2*N_CACHE/UNROLL] - z2;

							dx4 = B.x[k+3*N_CACHE/UNROLL] - x2;
							dy4 = B.y[k+3*N_CACHE/UNROLL] - y2;
							dz4 = B.z[k+3*N_CACHE/UNROLL] - z2;
#endif				
#endif

#ifdef NMK_NAIVE_GRAVITY
							dr21 = eps2 + dx1*dx1 + dy1*dy1 + dz1*dz1;
#if (UNROLL > 1)                                     
							dr22 = eps2 + dx2*dx2 + dy2*dy2 + dz2*dz2;
		//					dr21 = invsqrtf(eps2 + dx1*dx1 + dy1*dy1 + dz1*dz1);
		//					dr22 = invsqrtf(eps2 + dx2*dx2 + dy2*dy2 + dz2*dz2);
#if (UNROLL >= 4)
							dr23 = eps2 + dx3*dx3 + dy3*dy3 + dz3*dz3;
							dr24 = eps2 + dx4*dx4 + dy4*dy4 + dz4*dz4;
		//					dr23 = invsqrtf(eps2 + dx3*dx3 + dy3*dy3 + dz3*dz3);
		//					dr24 = invsqrtf(eps2 + dx4*dx4 + dy4*dy4 + dz4*dz4);
#endif
#endif
#ifdef NMK_SINGLE_PREC
							dr31 = ((Real) 1.0) / SQRT(dr21 * dr21 * dr21);
#if (UNROLL > 1)                                                  
							dr32 = ((Real) 1.0) / SQRT(dr22 * dr22 * dr22);
#if (UNROLL >= 4)                                                 
							dr33 = ((Real) 1.0) / SQRT(dr23 * dr23 * dr23);
																		  
							dr34 = ((Real) 1.0) / SQRT(dr24 * dr24 * dr24);
#endif
#endif
#else
							dr31 = INVSQRT(dr21 * dr21 * dr21);
#if (UNROLL > 1)
							dr32 = INVSQRT(dr22 * dr22 * dr22);
#if (UNROLL >= 4)
							dr33 = INVSQRT(dr23 * dr23 * dr23);

							dr34 = INVSQRT(dr24 * dr24 * dr24);
#endif
#endif
#endif

#else // NMK_NAIVE_GRAVITY
							dd1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
							dr1 = SQRT(dd1);
							dg1 = ERFC(dr1*inv2rs) + dr1*invpirs*EXP(dd1*invrs2);
							dr21 = eps2+dd1;
#if (UNROLL > 1)
							dd2 = dx2*dx2 + dy2*dy2 + dz2*dz2;
							dr2 = SQRT(dd2);
							dg2 = ERFC(dr2*inv2rs) + dr2*invpirs*EXP(dd2*invrs2);
							dr22 = eps2+dd2;
#if (UNROLL >=4)
							dd3 = dx3*dx3 + dy3*dy3 + dz3*dz3;
							dr3 = SQRT(dd3);
							dg3 = ERFC(dr3*inv2rs) + dr3*invpirs*EXP(dd3*invrs2);
							dr23 = eps2+dd3;

							dd4 = dx4*dx4 + dy4*dy4 + dz4*dz4;
							dr4 = SQRT(dd4);
							dg4 = ERFC(dr4*inv2rs) + dr4*invpirs*EXP(dd4*invrs2);
							dr24 = eps2+dd4;
#endif
#endif

#ifdef NMK_SINGLE_PREC
							dr31 = dg1 * ((Real) 1.0) / SQRT(dr21 * dr21 * dr21);
#if (UNROLL > 1)                                                        
							dr32 = dg2 * ((Real) 1.0) / SQRT(dr22 * dr22 * dr22);
#if (UNROLL >= 4)                                                       
							dr33 = dg3 * ((Real) 1.0) / SQRT(dr23 * dr23 * dr23);

							dr34 = dg4 * ((Real) 1.0) / SQRT(dr24 * dr24 * dr24);
#endif
#endif
#else // NMK_SINGLE_PREC
							dr31 = dg1 * INVSQRT(dr21 * dr21 * dr21);
#if (UNROLL > 1)                             
							dr32 = dg2 * INVSQRT(dr22 * dr22 * dr22);
#if (UNROLL >= 4)                            
							dr33 = dg3 * INVSQRT(dr23 * dr23 * dr23);

							dr34 = dg4 * INVSQRT(dr24 * dr24 * dr24);
#endif
#endif
#endif // NMK_SINGLE_PREC

#endif // NMK_NAIVE_GRAVITY

#if (UNROLL == 4)
							ax2 += dx1*dr31 + dx2*dr32 + dx3*dr33 + dx4*dr34;
							ay2 += dy1*dr31 + dy2*dr32 + dy3*dr33 + dy4*dr34;
							az2 += dz1*dr31 + dz2*dr32 + dz3*dr33 + dz4*dr34;
#endif
#if (UNROLL == 2) 
							ax2 += dx1*dr31 + dx2*dr32;
							ay2 += dy1*dr31 + dy2*dr32;
							az2 += dz1*dr31 + dz2*dr32;
#endif
#if (UNROLL == 1)
							ax2 += dx1*dr31 ;
							ay2 += dy1*dr31 ;
							az2 += dz1*dr31 ;
#endif
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
					dx1 = x2 - A.x[k];
					dy1 = y2 - A.y[k];
					dz1 = z2 - A.z[k];

#ifdef NMK_NAIVE_GRAVITY
					dr21 = eps2 + dx1*dx1 + dy1*dy1 + dz1*dz1;
#ifdef NMK_SINGLE_PREC
					dr31 = ((Real) 1.0) / SQRT(dr21 * dr21 * dr21);
#else
					dr31 = INVSQRT(dr21 * dr21 * dr21);
#endif

#else // NMK_NAIVE_GRAVITY
					dd1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
					dr1 = SQRT(dd1);
					dg1 = ERFC(dr1*inv2rs) + dr1*invpirs*EXP(dd1*invrs2);
					dr21 = eps2+dd1;
#ifdef NMK_SINGLE_PREC
					dr31 = dg1 * ((Real) 1.0) / SQRT(dr21 * dr21 * dr21);
#else
					dr31 = dg1 * INVSQRT(dr21 * dr21 * dr21);
#endif

#endif // NMK_NAIVE_GRAVITY

					C.x[k] += dx1*dr31 ;
					C.y[k] += dy1*dr31 ;
					C.z[k] += dz1*dr31 ;
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
		nt = ns + pern + tex;
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
				dx1 = B.x[k] - x2;
				dy1 = B.y[k] - y2;
				dz1 = B.z[k] - z2;

#ifdef NMK_NAIVE_GRAVITY
				dr21 = eps2 + dx1*dx1 + dy1*dy1 + dz1*dz1;
#ifdef NMK_SINGLE_PREC
				dr31 = ((Real) 1.0) / SQRT(dr21 * dr21 * dr21);
#else
				dr31 = INVSQRT(dr21 * dr21 * dr21);
#endif

#else // NMK_NAIVE_GRAVITY
					dd1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
					dr1 = SQRT(dd1);
					dg1 = ERFC(dr1*inv2rs) + dr1*invpirs*EXP(dd1*invrs2);
					dr21 = eps2+dd1;
#ifdef NMK_SINGLE_PREC
					dr31 = dg1 * ((Real) 1.0) / SQRT(dr21 * dr21 * dr21);
#else
					dr31 = dg1 * INVSQRT(dr21 * dr21 * dr21);
#endif
#endif // NMK_NAIVE_GRAVITY

				ax2 += dx1*dr31 ;
				ay2 += dy1*dr31 ;
				az2 += dz1*dr31 ;
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
	Real eps2 = constants->EPS2;
	Real invrs = (Real) 1.0 / constants->SPLIT_SCALE;
	Real invrs2 = (Real) -0.5 * invrs * invrs;
	Real invpirs = (Real) RSQRTPI * invrs;
	Real inv2rs = (Real) 0.5 * invrs;
	Real rc = constants->CUTOFF_SCALE;
	Real G = constants->GRAV_CONST;
#endif

#ifdef TABLEN
    int m;
    int idx;
    Real dlt  = constants->delta;
    Real *val = constants->value;
    Real *slp = constants->slope;
    Real dx1, dy1, dz1, dr21;
    Real rad2idx = TABLEN/(constants->CUTOFF_SCALE);

    for ( n=0; n<la; n++ ) {
        for ( m=0; m<lb; m++) {
            dx1 = B.x[m] -  A.x[n];
            dy1 = B.y[m] -  A.y[n];
            dz1 = B.z[m] -  A.z[n];
            
            dr21 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
            idx = (int)(rad2idx*dr21);
            if (idx<TABLEN) {
                dr21 = val[idx] + ( dr21 - dlt * idx ) * slp[idx];
            }
            else {
                dr21 = 0.0;                
            }                     

            C.x[n] += dx1*dr21*Bm[m];
            C.y[n] += dy1*dr21*Bm[m];
            C.z[n] += dz1*dr21*Bm[m];

        }
    }
    return 0;
#endif


#ifdef __MULTI_THREAD_
    double tstart, tstop, ttime;
    tstart = dtime();

#pragma omp parallel
#endif
	{
		int j, k, m, nb, mb, nt, mt;

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

#if 1
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

			for ( n=nbs; n<nbt; n++ )
			{
				nt = (n==nb-1) ? la : (n+1)*CLCNT;

				for (m=0; m<mb; m++)
				{
#ifdef __MULTI_THREAD_
					mt = ((m+tid)%mb)*N_CACHE + (((m+tid)%mb==mb-1 && lb%N_CACHE) ? lb-(mb-1)*N_CACHE : N_CACHE/UNROLL);
#else
					mt = m*N_CACHE + ((m==mb-1) ? lb-(mb-1)*N_CACHE : N_CACHE/UNROLL);
#endif

					if ( mt == lb ) // No unrolling at the last B block
					{
//						printf ("Thread[%d]: ms = %d, mt = %d.\n", tid, (mb-1)*N_CACHE, mt);
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

#pragma vector aligned
#pragma ivdep
#ifdef __MULTI_THREAD_
							for (k=((m+tid)%mb)*N_CACHE; k<mt; k++)
#else
							for (k=m*N_CACHE; k<mt; k++)
#endif
							{
								dx1 = B.x[k] - x2;
								dy1 = B.y[k] - y2;
								dz1 = B.z[k] - z2;

#ifdef NMK_NAIVE_GRAVITY
								dr21 = eps2 + dx1*dx1 + dy1*dy1 + dz1*dz1;
#ifdef NMK_SINGLE_PREC
								dr31 = Bm[k] * (((Real) 1.0) / SQRT(dr21 * dr21 * dr21));
#else
								dr31 = Bm[k] * INVSQRT(dr21 * dr21 * dr21);
#endif

#else // NMK_NAIVE_GRAVITY
								dd1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
								dr1 = SQRT(dd1);
								dg1 = ERFC(dr1*inv2rs) + dr1*invpirs*EXP(dd1*invrs2);
								dr21 = eps2+dd1;
#ifdef NMK_SINGLE_PREC
								dr31 = Bm[k] * dg1 * ((Real) 1.0) / SQRT(dr21 * dr21 * dr21);
#else
								dr31 = Bm[k] * dg1 * INVSQRT(dr21 * dr21 * dr21);
#endif

#endif // NMK_NAIVE_GRAVITY

								ax2 += dx1*dr31 ;
								ay2 += dy1*dr31 ;
								az2 += dz1*dr31 ;
							}

							C.x[j] += ax2;
							C.y[j] += ay2;
							C.z[j] += az2;
						}
					}
					else // Unrolling the loop
					{
//						printf ("Thread[%d]: ms = %d, mt = %d.\n", tid, m*N_CACHE, mt);
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
								dx1 = B.x[k] - x2;
								dy1 = B.y[k] - y2;
								dz1 = B.z[k] - z2;
#if (UNROLL > 1)
								dx2 = B.x[k+N_CACHE/UNROLL] - x2;
								dy2 = B.y[k+N_CACHE/UNROLL] - y2;
								dz2 = B.z[k+N_CACHE/UNROLL] - z2;
#if (UNROLL >= 4)
								dx3 = B.x[k+2*N_CACHE/UNROLL] - x2;
								dy3 = B.y[k+2*N_CACHE/UNROLL] - y2;
								dz3 = B.z[k+2*N_CACHE/UNROLL] - z2;

								dx4 = B.x[k+3*N_CACHE/UNROLL] - x2;
								dy4 = B.y[k+3*N_CACHE/UNROLL] - y2;
								dz4 = B.z[k+3*N_CACHE/UNROLL] - z2;
#endif				
#endif

#ifdef NMK_NAIVE_GRAVITY
								dr21 = eps2 + dx1*dx1 + dy1*dy1 + dz1*dz1;
#if (UNROLL > 1)                                     
								dr22 = eps2 + dx2*dx2 + dy2*dy2 + dz2*dz2;
			//					dr21 = invsqrtf(eps2 + dx1*dx1 + dy1*dy1 + dz1*dz1);
			//					dr22 = invsqrtf(eps2 + dx2*dx2 + dy2*dy2 + dz2*dz2);
#if (UNROLL >= 4)
								dr23 = eps2 + dx3*dx3 + dy3*dy3 + dz3*dz3;
								dr24 = eps2 + dx4*dx4 + dy4*dy4 + dz4*dz4;
			//					dr23 = invsqrtf(eps2 + dx3*dx3 + dy3*dy3 + dz3*dz3);
			//					dr24 = invsqrtf(eps2 + dx4*dx4 + dy4*dy4 + dz4*dz4);
#endif
#endif
#ifdef NMK_SINGLE_PREC
								dr31 = Bm[k] * (((Real) 1.0) / SQRT(dr21 * dr21 * dr21));
#if (UNROLL > 1)                                                  
								dr32 = Bm[k+N_CACHE/UNROLL] * (((Real) 1.0) / SQRT(dr22 * dr22 * dr22));
#if (UNROLL >= 4)                                                 
								dr33 = Bm[k+2*N_CACHE/UNROLL] * (((Real) 1.0) / SQRT(dr23 * dr23 * dr23));
																			  
								dr34 = Bm[k+3*N_CACHE/UNROLL] * (((Real) 1.0) / SQRT(dr24 * dr24 * dr24));
#endif
#endif
#else
								dr31 = Bm[k] * INVSQRT(dr21 * dr21 * dr21);
#if (UNROLL > 1)
								dr32 = Bm[k+N_CACHE/UNROLL] * INVSQRT(dr22 * dr22 * dr22);
#if (UNROLL >= 4)
								dr33 = Bm[k+2*N_CACHE/UNROLL] * INVSQRT(dr23 * dr23 * dr23);

								dr34 = Bm[k+3*N_CACHE/UNROLL] * INVSQRT(dr24 * dr24 * dr24);
#endif
#endif
#endif

#else // NMK_NAIVE_GRAVITY
								dd1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
								dr1 = SQRT(dd1);
								dg1 = ERFC(dr1*inv2rs) + dr1*invpirs*EXP(dd1*invrs2);
								dr21 = eps2+dd1;
#if (UNROLL > 1)
								dd2 = dx2*dx2 + dy2*dy2 + dz2*dz2;
								dr2 = SQRT(dd2);
								dg2 = ERFC(dr2*inv2rs) + dr2*invpirs*EXP(dd2*invrs2);
								dr22 = eps2+dd2;
#if (UNROLL >=4)
								dd3 = dx3*dx3 + dy3*dy3 + dz3*dz3;
								dr3 = SQRT(dd3);
								dg3 = ERFC(dr3*inv2rs) + dr3*invpirs*EXP(dd3*invrs2);
								dr23 = eps2+dd3;

								dd4 = dx4*dx4 + dy4*dy4 + dz4*dz4;
								dr4 = SQRT(dd4);
								dg4 = ERFC(dr4*inv2rs) + dr4*invpirs*EXP(dd4*invrs2);
								dr24 = eps2+dd4;
#endif
#endif

#ifdef NMK_SINGLE_PREC
								dr31 = Bm[k] * dg1 * ((Real) 1.0) / SQRT(dr21 * dr21 * dr21);
#if (UNROLL > 1)                                                        
								dr32 = Bm[k+N_CACHE/UNROLL] * dg2 * ((Real) 1.0) / SQRT(dr22 * dr22 * dr22);
#if (UNROLL >= 4)                                                                             
								dr33 = Bm[k+2*N_CACHE/UNROLL] * dg3 * ((Real) 1.0) / SQRT(dr23 * dr23 * dr23);
									                                									  
								dr34 = Bm[k+3*N_CACHE/UNROLL] * dg4 * ((Real) 1.0) / SQRT(dr24 * dr24 * dr24);
#endif
#endif
#else // NMK_SINGLE_PREC
								dr31 = Bm[k] * dg1 * INVSQRT(dr21 * dr21 * dr21);
#if (UNROLL > 1)                             
								dr32 = Bm[k+N_CACHE/UNROLL] * dg2 * INVSQRT(dr22 * dr22 * dr22);
#if (UNROLL >= 4)                                                  
								dr33 = Bm[k+2*N_CACHE/UNROLL] * dg3 * INVSQRT(dr23 * dr23 * dr23);
                                                                   
								dr34 = Bm[k+3*N_CACHE/UNROLL] * dg4 * INVSQRT(dr24 * dr24 * dr24);
#endif
#endif
#endif // NMK_SINGLE_PREC

#endif // NMK_NAIVE_GRAVITY

#if (UNROLL == 4)
								ax2 += dx1*dr31 + dx2*dr32 + dx3*dr33 + dx4*dr34;
								ay2 += dy1*dr31 + dy2*dr32 + dy3*dr33 + dy4*dr34;
								az2 += dz1*dr31 + dz2*dr32 + dz3*dr33 + dz4*dr34;
#endif
#if (UNROLL == 2) 
								ax2 += dx1*dr31 + dx2*dr32;
								ay2 += dy1*dr31 + dy2*dr32;
								az2 += dz1*dr31 + dz2*dr32;
#endif
#if (UNROLL == 1)
								ax2 += dx1*dr31 ;
								ay2 += dy1*dr31 ;
								az2 += dz1*dr31 ;
#endif
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
							dx1 = B.x[k] - x2;
							dy1 = B.y[k] - y2;
							dz1 = B.z[k] - z2;

#ifdef NMK_NAIVE_GRAVITY
							dr21 = eps2 + dx1*dx1 + dy1*dy1 + dz1*dz1;
#ifdef NMK_SINGLE_PREC
							dr31 = Bm[k] * (((Real) 1.0) / SQRT(dr21 * dr21 * dr21));
#else
							dr31 = Bm[k] * (INVSQRT(dr21 * dr21 * dr21));
#endif

#else // NMK_NAIVE_GRAVITY
							dd1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
							dr1 = SQRT(dd1);
							dg1 = ERFC(dr1*inv2rs) + dr1*invpirs*EXP(dd1*invrs2);
							dr21 = eps2+dd1;
#ifdef NMK_SINGLE_PREC
							dr31 = Bm[k] * dg1 * ((Real) 1.0) / SQRT(dr21 * dr21 * dr21);
#else
							dr31 = Bm[k] * dg1 * INVSQRT(dr21 * dr21 * dr21);
#endif

#endif // NMK_NAIVE_GRAVITY

							ax2 += dx1*dr31 ;
							ay2 += dy1*dr31 ;
							az2 += dz1*dr31 ;
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
							dx1 = B.x[k] - x2;
							dy1 = B.y[k] - y2;
							dz1 = B.z[k] - z2;
#if (UNROLL > 1)
							dx2 = B.x[k+N_CACHE/UNROLL] - x2;
							dy2 = B.y[k+N_CACHE/UNROLL] - y2;
							dz2 = B.z[k+N_CACHE/UNROLL] - z2;
#if (UNROLL >= 4)
							dx3 = B.x[k+2*N_CACHE/UNROLL] - x2;
							dy3 = B.y[k+2*N_CACHE/UNROLL] - y2;
							dz3 = B.z[k+2*N_CACHE/UNROLL] - z2;

							dx4 = B.x[k+3*N_CACHE/UNROLL] - x2;
							dy4 = B.y[k+3*N_CACHE/UNROLL] - y2;
							dz4 = B.z[k+3*N_CACHE/UNROLL] - z2;
#endif				
#endif

#ifdef NMK_NAIVE_GRAVITY
							dr21 = eps2 + dx1*dx1 + dy1*dy1 + dz1*dz1;
#if (UNROLL > 1)                                     
							dr22 = eps2 + dx2*dx2 + dy2*dy2 + dz2*dz2;
		//					dr21 = invsqrtf(eps2 + dx1*dx1 + dy1*dy1 + dz1*dz1);
		//					dr22 = invsqrtf(eps2 + dx2*dx2 + dy2*dy2 + dz2*dz2);
#if (UNROLL >= 4)
							dr23 = eps2 + dx3*dx3 + dy3*dy3 + dz3*dz3;
							dr24 = eps2 + dx4*dx4 + dy4*dy4 + dz4*dz4;
		//					dr23 = invsqrtf(eps2 + dx3*dx3 + dy3*dy3 + dz3*dz3);
		//					dr24 = invsqrtf(eps2 + dx4*dx4 + dy4*dy4 + dz4*dz4);
#endif
#endif
#ifdef NMK_SINGLE_PREC
							dr31 = Bm[k] * (((Real) 1.0) / SQRT(dr21 * dr21 * dr21));
#if (UNROLL > 1)                                                  
							dr32 = Bm[k] * (((Real) 1.0) / SQRT(dr22 * dr22 * dr22));
#if (UNROLL >= 4)                                                 
							dr33 = Bm[k] * (((Real) 1.0) / SQRT(dr23 * dr23 * dr23));
																		  
							dr34 = Bm[k] * (((Real) 1.0) / SQRT(dr24 * dr24 * dr24));
#endif
#endif
#else
							dr31 = Bm[k] * INVSQRT(dr21 * dr21 * dr21);
#if (UNROLL > 1)                           
							dr32 = Bm[k] * INVSQRT(dr22 * dr22 * dr22);
#if (UNROLL >= 4)                          
							dr33 = Bm[k] * INVSQRT(dr23 * dr23 * dr23);
                                    	   
							dr34 = Bm[k] * INVSQRT(dr24 * dr24 * dr24);
#endif
#endif
#endif

#else // NMK_NAIVE_GRAVITY
							dd1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
							dr1 = SQRT(dd1);
							dg1 = ERFC(dr1*inv2rs) + dr1*invpirs*EXP(dd1*invrs2);
							dr21 = eps2+dd1;
#if (UNROLL > 1)
							dd2 = dx2*dx2 + dy2*dy2 + dz2*dz2;
							dr2 = SQRT(dd2);
							dg2 = ERFC(dr2*inv2rs) + dr2*invpirs*EXP(dd2*invrs2);
							dr22 = eps2+dd2;
#if (UNROLL >=4)
							dd3 = dx3*dx3 + dy3*dy3 + dz3*dz3;
							dr3 = SQRT(dd3);
							dg3 = ERFC(dr3*inv2rs) + dr3*invpirs*EXP(dd3*invrs2);
							dr23 = eps2+dd3;

							dd4 = dx4*dx4 + dy4*dy4 + dz4*dz4;
							dr4 = SQRT(dd4);
							dg4 = ERFC(dr4*inv2rs) + dr4*invpirs*EXP(dd4*invrs2);
							dr24 = eps2+dd4;
#endif
#endif

#ifdef NMK_SINGLE_PREC
							dr31 = Bm[k] * dg1 * ((Real) 1.0) / SQRT(dr21 * dr21 * dr21);
#if (UNROLL > 1)                                                        
							dr32 = Bm[k+N_CACHE/UNROLL] * dg2 * ((Real) 1.0) / SQRT(dr22 * dr22 * dr22);
#if (UNROLL >= 4)                                                                             
							dr33 = Bm[k+2*N_CACHE/UNROLL] * dg3 * ((Real) 1.0) / SQRT(dr23 * dr23 * dr23);

							dr34 = Bm[k+3*N_CACHE/UNROLL] * dg4 * ((Real) 1.0) / SQRT(dr24 * dr24 * dr24);
#endif
#endif
#else // NMK_SINGLE_PREC
							dr31 = Bm[k] * dg1 * INVSQRT(dr21 * dr21 * dr21);
#if (UNROLL > 1)                             
							dr32 = Bm[k+N_CACHE/UNROLL] * dg2 * INVSQRT(dr22 * dr22 * dr22);
#if (UNROLL >= 4)                                                  
							dr33 = Bm[k+2*N_CACHE/UNROLL] * dg3 * INVSQRT(dr23 * dr23 * dr23);

							dr34 = Bm[k+3*N_CACHE/UNROLL] * dg4 * INVSQRT(dr24 * dr24 * dr24);
#endif
#endif
#endif // NMK_SINGLE_PREC

#endif // NMK_NAIVE_GRAVITY

#if (UNROLL == 4)
							ax2 += dx1*dr31 + dx2*dr32 + dx3*dr33 + dx4*dr34;
							ay2 += dy1*dr31 + dy2*dr32 + dy3*dr33 + dy4*dr34;
							az2 += dz1*dr31 + dz2*dr32 + dz3*dr33 + dz4*dr34;
#endif
#if (UNROLL == 2) 
							ax2 += dx1*dr31 + dx2*dr32;
							ay2 += dy1*dr31 + dy2*dr32;
							az2 += dz1*dr31 + dz2*dr32;
#endif
#if (UNROLL == 1)
							ax2 += dx1*dr31 ;
							ay2 += dy1*dr31 ;
							az2 += dz1*dr31 ;
#endif
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
					dx1 = x2 - A.x[k];
					dy1 = y2 - A.y[k];
					dz1 = z2 - A.z[k];

#ifdef NMK_NAIVE_GRAVITY
					dr21 = eps2 + dx1*dx1 + dy1*dy1 + dz1*dz1;
#ifdef NMK_SINGLE_PREC
					dr31 = Bm[m] * (((Real) 1.0) / SQRT(dr21 * dr21 * dr21));
#else
					dr31 = Bm[m] * INVSQRT(dr21 * dr21 * dr21);
#endif

#else // NMK_NAIVE_GRAVITY
					dd1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
					dr1 = SQRT(dd1);
					dg1 = ERFC(dr1*inv2rs) + dr1*invpirs*EXP(dd1*invrs2);
					dr21 = eps2+dd1;
#ifdef NMK_SINGLE_PREC
					dr31 = Bm[m] * dg1 * ((Real) 1.0) / SQRT(dr21 * dr21 * dr21);
#else
					dr31 = Bm[m] * dg1 * INVSQRT(dr21 * dr21 * dr21);
#endif

#endif // NMK_NAIVE_GRAVITY

					C.x[k] += dx1*dr31 ;
					C.y[k] += dy1*dr31 ;
					C.z[k] += dz1*dr31 ;
				}
			}
		}
#else
		// Original ppkernel with mass.
		int pern = la/tnum;
		int modn = la%tnum;
		int tpass = modn-tid;
		int tex = (tnum+tpass-1)/tnum;
		int ns = tid*pern+(modn-tex*tpass);
		nt = ns + pern + tex;
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
				dx1 = B.x[k] - x2;
				dy1 = B.y[k] - y2;
				dz1 = B.z[k] - z2;

#ifdef NMK_NAIVE_GRAVITY
				dr21 = eps2 + dx1*dx1 + dy1*dy1 + dz1*dz1;
#ifdef NMK_SINGLE_PREC
				dr31 = Bm[k] * (((Real) 1.0) / SQRT(dr21 * dr21 * dr21));
#else
				dr31 = Bm[k] * (INVSQRT(dr21 * dr21 * dr21));
#endif

#else // NMK_NAIVE_GRAVITY
					dd1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
					dr1 = SQRT(dd1);
					dg1 = ERFC(dr1*inv2rs) + dr1*invpirs*EXP(dd1*invrs2);
					dr21 = eps2+dd1;
#ifdef NMK_SINGLE_PREC
					dr31 = Bm[m] * dg1 * ((Real) 1.0) / SQRT(dr21 * dr21 * dr21);
#else
					dr31 = Bm[m] * dg1 * INVSQRT(dr21 * dr21 * dr21);
#endif
#endif // NMK_NAIVE_GRAVITY

				ax2 += dx1*dr31 ;
				ay2 += dy1*dr31 ;
				az2 += dz1*dr31 ;
			}
			C.x[j] += ax2;
			C.y[j] += ay2;
			C.z[j] += az2;
		}
#endif
	}

#ifdef __MULTI_THREAD_
    tstop = dtime();
    ttime = tstop - tstart;
    printf("dtime = %lf \n", ttime);
#endif

    return 0;
}


