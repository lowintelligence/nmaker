#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <malloc.h>
#include <omp.h>

#include "ppkernel.h"

double dtime()
{
	double tseconds = 0.0;
	struct timeval mytime;	
	gettimeofday(&mytime, (struct timezone*)0);
	tseconds = (double)(mytime.tv_sec + mytime.tv_usec*1.0e-6);
	return (tseconds);
}

int ppkernel(PPPack *A, int la, PPPack *B, int lb, PRECTYPE eps2)
{
    int n;
    double tstart, tstop, ttime;

    /*
    printf("Starting Compute\n");
    
    tstart = dtime();
    
    for ( n=0; n<la; n++ ) {
        for ( m=0; m<lb; m++) {
            dx1 = B->pos.x[m] -  A->pos.x[n];
            dy1 = B->pos.y[m] -  A->pos.y[n];
            dz1 = B->pos.z[m] -  A->pos.z[n];
            
            dr21 = eps2 + dx*dx + dy*dy + dz*dz;
            dr31 = SQRT(dr2)*dr2;
            
            A->acc.x[n] += dx1/dr31;
            A->acc.y[n] += dy1/dr31;
            A->acc.z[n] += dz1/dr31;
        }
    }

    tstop = dtime();
    
    ttime = tstop - tstart;
    printf("dtime = %lf \n", ttime);
    */
    
    tstart = dtime();

#pragma omp parallel
	{
		int j, k, m, nb, mb, nt, mt, tid, tnum;

		PRECTYPE px2, py2, pz2, ax2, ay2, az2;
		PRECTYPE dx1, dy1, dz1, dr21, dr31;
		PRECTYPE dx2, dy2, dz2, dr22, dr32;
		PRECTYPE dx3, dy3, dz3, dr23, dr33;
		PRECTYPE dx4, dy4, dz4, dr24, dr34;

		tid = omp_get_thread_num();
		tnum = omp_get_num_threads();
	
		nb = (la+L2CNT-1)/L2CNT;
		mb = (lb+L1CNT-1)/L1CNT;

		if( nb >= (tnum<<1) )
		{
#pragma omp for schedule(static)
			for ( n=0; n<nb; n++ ) {

				nt = (n==nb-1) ? la : (n+1)*L2CNT;

				for (m=0; m<mb; m++) {
#ifdef _OPENMP
					mt = ((m+tid)%mb)*L1CNT + (((m+tid)%mb==mb-1 && lb%L1CNT) ? lb-(mb-1)*L1CNT : L1CNT/UNROLL);
#else
					mt = m*L1CNT + ((m==mb-1) ? (lb-m*L1CNT)/UNROLL : L1CNT/UNROLL);
#endif

					if ( mt == lb ) // No unrolling at the last B block
					{
//#pragma prefetch A->pos.x:1:L2CNT
//#pragma prefetch A->pos.y:1:L2CNT
//#pragma prefetch A->pos.z:1:L2CNT
						for (j=n*L2CNT; j<nt; j++)
						{
							px2=A->pos.x[j];
							py2=A->pos.y[j];
							pz2=A->pos.z[j];

							ax2=0;
							ay2=0;
							az2=0;

#pragma ivdep
#pragma vector aligned
#ifdef _OPENMP
							for (k=((m+tid)%mb)*L1CNT; k<mt; k++)
#else
							for (k=m*L1CNT; k<mt; k++)
#endif
							{
								dx1 = B->pos.x[k] - px2;
								dy1 = B->pos.y[k] - py2;
								dz1 = B->pos.z[k] - pz2;

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

							A->acc.x[j] += ax2;
							A->acc.y[j] += ay2;
							A->acc.z[j] += az2;
						}
					}
					else // Unrolling the loop
					{
//#pragma prefetch A->pos.x:1:L2CNT
//#pragma prefetch A->pos.y:1:L2CNT
//#pragma prefetch A->pos.z:1:L2CNT
						for (j=n*L2CNT; j<nt; j++)
						{
							px2=A->pos.x[j];
							py2=A->pos.y[j];
							pz2=A->pos.z[j];

							ax2=0;
							ay2=0;
							az2=0;

#pragma ivdep
#pragma vector aligned
#ifdef _OPENMP
							for (k=((m+tid)%mb)*L1CNT; k<mt; k++)
#else
							for (k=m*L1CNT; k<mt; k++)
#endif
							{
								dx1 = B->pos.x[k] - px2;
								dy1 = B->pos.y[k] - py2;
								dz1 = B->pos.z[k] - pz2;
#if (UNROLL > 1)
								dx2 = B->pos.x[k+L1CNT/UNROLL] - px2;
								dy2 = B->pos.y[k+L1CNT/UNROLL] - py2;
								dz2 = B->pos.z[k+L1CNT/UNROLL] - pz2;
#if (UNROLL >= 4)
								dx3 = B->pos.x[k+2*L1CNT/UNROLL] - px2;
								dy3 = B->pos.y[k+2*L1CNT/UNROLL] - py2;
								dz3 = B->pos.z[k+2*L1CNT/UNROLL] - pz2;

								dx4 = B->pos.x[k+3*L1CNT/UNROLL] - px2;
								dy4 = B->pos.y[k+3*L1CNT/UNROLL] - py2;
								dz4 = B->pos.z[k+3*L1CNT/UNROLL] - pz2;
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

							A->acc.x[j] += ax2;
							A->acc.y[j] += ay2;
							A->acc.z[j] += az2;
						}
					}
				}
			}
		}
		else
		{
#pragma omp for schedule(static)
			for ( n=0; n<la; n++ ) {

				px2=A->pos.x[n];
				py2=A->pos.y[n];
				pz2=A->pos.z[n];

				ax2=0;
				ay2=0;
				az2=0;

				for (m=0; m<mb; m++)
				{
#ifdef _OPENMP
					mt = ((m+tid)%mb)*L1CNT + (((m+tid)%mb==mb-1 && lb%L1CNT) ? lb-(mb-1)*L1CNT : L1CNT/UNROLL);
#else
					mt = m*L1CNT + ((m==mb-1) ? (lb-m*L1CNT)/UNROLL : L1CNT/UNROLL);
#endif

					if ( mt == lb ) // No unrolling at the last B block
					{
#pragma ivdep
#pragma vector aligned
#ifdef _OPENMP
						for (k=((m+tid)%mb)*L1CNT; k<mt; k++)
#else
						for (k=m*L1CNT; k<mt; k++)
#endif
						{
							dx1 = B->pos.x[k] - px2;
							dy1 = B->pos.y[k] - py2;
							dz1 = B->pos.z[k] - pz2;

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
#pragma ivdep
#pragma vector aligned
#ifdef _OPENMP
						for (k=((m+tid)%mb)*L1CNT; k<mt; k++)
#else
						for (k=m*L1CNT; k<mt; k++)
#endif
						{
							dx1 = B->pos.x[k] - px2;
							dy1 = B->pos.y[k] - py2;
							dz1 = B->pos.z[k] - pz2;
#if (UNROLL > 1)
							dx2 = B->pos.x[k+L1CNT/UNROLL] - px2;
							dy2 = B->pos.y[k+L1CNT/UNROLL] - py2;
							dz2 = B->pos.z[k+L1CNT/UNROLL] - pz2;
#if (UNROLL >= 4)
							dx3 = B->pos.x[k+2*L1CNT/UNROLL] - px2;
							dy3 = B->pos.y[k+2*L1CNT/UNROLL] - py2;
							dz3 = B->pos.z[k+2*L1CNT/UNROLL] - pz2;

							dx4 = B->pos.x[k+3*L1CNT/UNROLL] - px2;
							dy4 = B->pos.y[k+3*L1CNT/UNROLL] - py2;
							dz4 = B->pos.z[k+3*L1CNT/UNROLL] - pz2;
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
				}

				A->acc.x[n] += ax2;
				A->acc.y[n] += ay2;
				A->acc.z[n] += az2;
			}
		}
	}

    tstop = dtime();
    ttime = tstop - tstart;
    //printf("dtime = %lf \n", ttime);

    return 0;
}


