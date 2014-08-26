#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <malloc.h>
#ifdef __MULTI_THREAD_
#include <omp.h>
#endif

#include "ppkernel.h"
//#include "dtime.h"

__OffloadFunc_Macro__
double dtime()
{
	double tseconds = 0.0;
	struct timeval mytime;	
	gettimeofday(&mytime, (struct timezone*)0);
	tseconds = (double)(mytime.tv_sec + mytime.tv_usec*1.0e-6);
	return (tseconds);
}

int get_block_tnum(int bid)
{
	return omp_get_num_threads();
}

int get_block_tid(int bid)
{
	return omp_get_thread_num();
}

#ifdef __MULTI_THREAD_
int ppkernel(Array3 A, int la, Array3 B, int lb, PRECTYPE eps2, Array3 C)
#else
int ppkernel(Array3 A, int la, Array3 B, int lb, PRECTYPE eps2, Array3 C, int tnum, int tid)
#endif
{
    int n;
    double tstart, tstop, ttime;

    /*
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

//    tstart = dtime();

//#pragma omp parallel
	{
		int j, k, m, nb, mb, nt, mt;

		PRECTYPE x2, y2, z2, ax2, ay2, az2;
		PRECTYPE dx1, dy1, dz1, dr21, dr31;
		PRECTYPE dx2, dy2, dz2, dr22, dr32;
		PRECTYPE dx3, dy3, dz3, dr23, dr33;
		PRECTYPE dx4, dy4, dz4, dr24, dr34;

#ifdef __MULTI_THREAD_
		int tid, tnum;
		tid = get_block_tid(0);
		tnum = get_block_tnum(0);
#endif
//		printf ("%ld, %d, %ld, %d, %ld, %d, %d.\n", A.x, la, B.x, lb, C.x, tnum, tid);
//		return (0);
	
		nb = (la+CLCNT-1)/CLCNT;
		mb = (lb+N_CACHE-1)/N_CACHE;

		if( (nb >= (tnum<<1)) && (lb > CLCNT) ) // The num of A and B are all large enough
		{
//			int nbs = tid*(nb/tnum)+((nb%tnum)+tid)/tnum*((nb+tid)%tnum);
//			int nbt = nbs + (nb+tid)/tnum;
			int nbs = tid*(nb/tnum)+(tnum+(nb%tnum)-tid)/tnum*(tid%tnum);
			int nbt = nbs + (nb+tnum-1-tid)/tnum;
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

								dr21 = eps2 + dx1*dx1 + dy1*dy1 + dz1*dz1;
#ifdef __single_prec
								dr31 = ((PRECTYPE) 1.0) / SQRT(dr21 * dr21 * dr21);
#else
								dr31 = INVSQRT(dr21 * dr21 * dr21);
#endif

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
#ifdef __single_prec
								dr31 = ((PRECTYPE) 1.0) / SQRT(dr21 * dr21 * dr21);
#if (UNROLL > 1)                                                  
								dr32 = ((PRECTYPE) 1.0) / SQRT(dr22 * dr22 * dr22);
#if (UNROLL >= 4)                                                 
								dr33 = ((PRECTYPE) 1.0) / SQRT(dr23 * dr23 * dr23);
																			  
								dr34 = ((PRECTYPE) 1.0) / SQRT(dr24 * dr24 * dr24);
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
			int ns = tid*(la/tnum)+(tnum+(la%tnum)-tid)/tnum*(tid%tnum);
			int nt = ns + (la+tnum-1-tid)/tnum;

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

							dr21 = eps2 + dx1*dx1 + dy1*dy1 + dz1*dz1;
#ifdef __single_prec
							dr31 = ((PRECTYPE) 1.0) / SQRT(dr21 * dr21 * dr21);
#else
							dr31 = INVSQRT(dr21 * dr21 * dr21);
#endif

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
#ifdef __single_prec
							dr31 = ((PRECTYPE) 1.0) / SQRT(dr21 * dr21 * dr21);
#if (UNROLL > 1)                                                  
							dr32 = ((PRECTYPE) 1.0) / SQRT(dr22 * dr22 * dr22);
#if (UNROLL >= 4)                                                 
							dr33 = ((PRECTYPE) 1.0) / SQRT(dr23 * dr23 * dr23);
																		  
							dr34 = ((PRECTYPE) 1.0) / SQRT(dr24 * dr24 * dr24);
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
			int ns = tid*(la/tnum)+(tnum+(la%tnum)-tid)/tnum*(tid%tnum);
			int nt = ns + (la+tnum-1-tid)/tnum;
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

					dr21 = eps2 + dx1*dx1 + dy1*dy1 + dz1*dz1;
#ifdef __single_prec
					dr31 = ((PRECTYPE) 1.0) / SQRT(dr21 * dr21 * dr21);
#else
					dr31 = INVSQRT(dr21 * dr21 * dr21);
#endif

					C.x[k] += dx1*dr31 ;
					C.y[k] += dy1*dr31 ;
					C.z[k] += dz1*dr31 ;
				}
			}
		}
	}

//    tstop = dtime();
//    ttime = tstop - tstart;
//    printf("dtime = %lf \n", ttime);

    return 0;
}

#ifdef __MULTI_THREAD_
int ppmkernel(Array3 A, int la, Array3 B, PRECTYPE *Bm, int lb, PRECTYPE eps2, Array3 C)
#else
int ppmkernel(Array3 A, int la, Array3 B, PRECTYPE *Bm, int lb, PRECTYPE eps2, Array3 C, int tnum, int tid)
#endif
{
    int n;
    double tstart, tstop, ttime;

    /*
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

//    tstart = dtime();

//#pragma omp parallel
	{
		int j, k, m, nb, mb, nt, mt;

		PRECTYPE x2, y2, z2, ax2, ay2, az2;
		PRECTYPE dx1, dy1, dz1, dr21, dr31;
		PRECTYPE dx2, dy2, dz2, dr22, dr32;
		PRECTYPE dx3, dy3, dz3, dr23, dr33;
		PRECTYPE dx4, dy4, dz4, dr24, dr34;

#ifdef __MULTI_THREAD_
		int tid, tnum;
		tid = get_block_tid(0);
		tnum = get_block_tnum(0);
#endif
	
		nb = (la+CLCNT-1)/CLCNT;
		mb = (lb+N_CACHE-1)/N_CACHE;

		if( (nb >= (tnum<<1)) && (lb > CLCNT) ) // The num of A and B are all large enough
		{
//			int nbs = tid*(nb/tnum)+((nb%tnum)+tid)/tnum*((nb+tid)%tnum);
//			int nbt = nbs + (nb+tid)/tnum;
			int nbs = tid*(nb/tnum)+(tnum+(nb%tnum)-tid)/tnum*(tid%tnum);
			int nbt = nbs + (nb+tnum-1-tid)/tnum;
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

								dr21 = eps2 + dx1*dx1 + dy1*dy1 + dz1*dz1;
#ifdef __single_prec
								dr31 = Bm[k] * (((PRECTYPE) 1.0) / SQRT(dr21 * dr21 * dr21));
#else
								dr31 = Bm[k] * INVSQRT(dr21 * dr21 * dr21);
#endif

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
#ifdef __single_prec
								dr31 = Bm[k] * (((PRECTYPE) 1.0) / SQRT(dr21 * dr21 * dr21));
#if (UNROLL > 1)                                                  
								dr32 = Bm[k+N_CACHE/UNROLL] * (((PRECTYPE) 1.0) / SQRT(dr22 * dr22 * dr22));
#if (UNROLL >= 4)                                                 
								dr33 = Bm[k+2*N_CACHE/UNROLL] * (((PRECTYPE) 1.0) / SQRT(dr23 * dr23 * dr23));
																			  
								dr34 = Bm[k+3*N_CACHE/UNROLL] * (((PRECTYPE) 1.0) / SQRT(dr24 * dr24 * dr24));
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
			int ns = tid*(la/tnum)+(tnum+(la%tnum)-tid)/tnum*(tid%tnum);
			int nt = ns + (la+tnum-1-tid)/tnum;

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

							dr21 = eps2 + dx1*dx1 + dy1*dy1 + dz1*dz1;
#ifdef __single_prec
							dr31 = Bm[k] * (((PRECTYPE) 1.0) / SQRT(dr21 * dr21 * dr21));
#else
							dr31 = Bm[k] * (INVSQRT(dr21 * dr21 * dr21));
#endif

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
#ifdef __single_prec
							dr31 = Bm[k] * (((PRECTYPE) 1.0) / SQRT(dr21 * dr21 * dr21));
#if (UNROLL > 1)                                                  
							dr32 = Bm[k] * (((PRECTYPE) 1.0) / SQRT(dr22 * dr22 * dr22));
#if (UNROLL >= 4)                                                 
							dr33 = Bm[k] * (((PRECTYPE) 1.0) / SQRT(dr23 * dr23 * dr23));
																		  
							dr34 = Bm[k] * (((PRECTYPE) 1.0) / SQRT(dr24 * dr24 * dr24));
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
			int ns = tid*(la/tnum)+(tnum+(la%tnum)-tid)/tnum*(tid%tnum);
			int nt = ns + (la+tnum-1-tid)/tnum;
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

					dr21 = eps2 + dx1*dx1 + dy1*dy1 + dz1*dz1;
#ifdef __single_prec
					dr31 = Bm[k] * (((PRECTYPE) 1.0) / SQRT(dr21 * dr21 * dr21));
#else
					dr31 = Bm[k] * INVSQRT(dr21 * dr21 * dr21);
#endif

					C.x[k] += dx1*dr31 ;
					C.y[k] += dy1*dr31 ;
					C.z[k] += dz1*dr31 ;
				}
			}
		}
	}

//    tstop = dtime();
//    ttime = tstop - tstart;
//    printf("dtime = %lf \n", ttime);

    return 0;
}


