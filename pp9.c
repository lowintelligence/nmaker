#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <malloc.h>
#include <omp.h>

#include "offload.h"
#include "ppkernel.h"

#ifdef __PAPI
#include <papiStdEventDefs.h>
#include <papi.h>
#endif

__OffloadVar_Macro__
PPPack part;

__OffloadFunc_Macro__
double dtime();

float ran3(long *idum);

void initialize_particle(int nPart, double box)
{
    int n, m, d;
    long seed = 16435234;
    part.pos.x = (PRECTYPE*)memalign(ALIGNCNT, sizeof(PRECTYPE)*nPart);
    part.pos.y = (PRECTYPE*)memalign(ALIGNCNT, sizeof(PRECTYPE)*nPart);
    part.pos.z = (PRECTYPE*)memalign(ALIGNCNT, sizeof(PRECTYPE)*nPart);
    part.acc.x = (PRECTYPE*)memalign(ALIGNCNT, sizeof(PRECTYPE)*nPart);
    part.acc.y = (PRECTYPE*)memalign(ALIGNCNT, sizeof(PRECTYPE)*nPart);
    part.acc.z = (PRECTYPE*)memalign(ALIGNCNT, sizeof(PRECTYPE)*nPart);

    for (n=0; n<nPart; n++) {
         part.pos.x[n] = box*ran3(&seed);
         part.pos.y[n] = box*ran3(&seed);
         part.pos.z[n] = box*ran3(&seed);
//         part.pos.x[n] = (PRECTYPE) 0.0;
//         part.pos.y[n] = (PRECTYPE) 0.0;
//         part.pos.z[n] = (PRECTYPE) n;
         part.acc.x[n] = (PRECTYPE) 0.0;
         part.acc.y[n] = (PRECTYPE) 0.0;
         part.acc.z[n] = (PRECTYPE) 0.0;
    }
}


int main(int argc, char *argv[])
{
    int n, nPart, lb, la;
    double tstart, tstop, ttime;
    double gflops = 0.0;
    float a = 1.1;
//	PRECTYPE eps = 1.0/100.0/40.0;
	PRECTYPE eps = 0.00026;
#ifdef __PAPI
	int EventSet = PAPI_NULL;
	long long value;
#endif
	    
    nPart = 1048576;
	if (argc>1)
	{
		nPart=atoi(argv[1]);
	}
	lb = la = nPart;
	if (argc>2)
	{
		lb=atoi(argv[2]);
	}
	if (argc>3)
	{
		la=atoi(argv[3]);
	}
    /*
    printf("Initializing\r\n");
    initialize_particle(nPart, 1.0);

    printf("Starting Compute\n");
    
    tstart = dtime();
    
    for ( n=0; n<nPart-1; n++ ) {
        for ( m=n+1; m<nPart; m++) {
            dx = part[m].pos[0] -  part[n].pos[0];
            dy = part[m].pos[1] -  part[n].pos[1];
            dz = part[m].pos[2] -  part[n].pos[2];
            
            dr2 = dx*dx + dy*dy + dz*dz;
            dr3 = sqrt(dr2)*dr2;
            
            part[n].acc[0] += dx/dr3;
            part[n].acc[1] += dy/dr3;
            part[n].acc[2] += dz/dr3;
            
            part[m].acc[0] -= dx/dr3;
            part[m].acc[1] -= dy/dr3;
            part[m].acc[2] -= dz/dr3;
        }
    }

    tstop = dtime();
    
    free(part);
    
    ttime = tstop - tstart;
    printf("dtime = %lf \n", ttime);
    */
    initialize_particle(nPart, 1.0);
 //   omp_set_num_threads(32);
  //  kmp_set_defaults("KMP_AFFINITY=compact");
//    kmp_set_defaults("KMP_AFFINITY=scatter");
#ifdef __INTEL_OFFLOAD
	int nCPU;
	int nMIC0;
	int nMIC1;

//	nCPU = 0;
//	nMIC0 = nPart / 2;
	nCPU = (int) ((double)nPart / (double)CLCNT / 6.72) * CLCNT; 
	nMIC0 = (int) ((double)nCPU * 2.86) / CLCNT * CLCNT; 
	nMIC1 = nPart - nCPU - nMIC0;
	printf ("nCPU=%d, nMIC0=%d, nMIC1=%d.\n", nCPU, nMIC0, nMIC1);

	PRECTYPE *px = part.pos.x;
	PRECTYPE *py = part.pos.y;
	PRECTYPE *pz = part.pos.z;
	PRECTYPE *ax = part.acc.x;
	PRECTYPE *ay = part.acc.y;
	PRECTYPE *az = part.acc.z;
	PRECTYPE *ax0 = part.acc.x + nCPU;
	PRECTYPE *ay0 = part.acc.y + nCPU;
	PRECTYPE *az0 = part.acc.z + nCPU;
	PRECTYPE *ax1 = part.acc.x + nCPU + nMIC0;
	PRECTYPE *ay1 = part.acc.y + nCPU + nMIC0;
	PRECTYPE *az1 = part.acc.z + nCPU + nMIC0;
#pragma offload target(mic:0) in (px[0:nPart]:alloc_if(1) free_if(0)) in (py[0:nPart]:alloc_if(1) free_if(0)) in (pz[0:nPart]:alloc_if(1) free_if(0)) out (ax0[0:nMIC0]:alloc_if(1) free_if(0)) out (ay0[0:nMIC0]:alloc_if(1) free_if(0)) out (az0[0:nMIC0]:alloc_if(1) free_if(0)) signal(nMIC0)
{
	int n;
	Array3 PA, PB, AC;
	PB.x = px;
	PB.y = py;
	PB.z = pz;
	PA.x = px + nCPU;
	PA.y = py + nCPU;
	PA.z = pz + nCPU;
	AC.x = ax0;
	AC.y = ay0;
	AC.z = az0;
    kmp_set_defaults("KMP_AFFINITY=compact");
#pragma omp parallel for schedule(static)
    for (n=0; n<nMIC0; n++) {
         AC.x[n] = 0;
         AC.y[n] = 0;
         AC.z[n] = 0;
    }
}
#pragma offload target(mic:1) in (px[0:nPart]:alloc_if(1) free_if(0)) in (py[0:nPart]:alloc_if(1) free_if(0)) in (pz[0:nPart]:alloc_if(1) free_if(0)) out (ax1[0:nMIC1]:alloc_if(1) free_if(0)) out (ay1[0:nMIC1]:alloc_if(1) free_if(0)) out (az1[0:nMIC1]:alloc_if(1) free_if(0)) signal(nMIC1)
{
	int n;
	Array3 PA, PB, AC;
	PB.x = px;
	PB.y = py;
	PB.z = pz;
	PA.x = px + nCPU + nMIC0;
	PA.y = py + nCPU + nMIC0;
	PA.z = pz + nCPU + nMIC0;
	AC.x = ax1;
	AC.y = ay1;
	AC.z = az1;
    kmp_set_defaults("KMP_AFFINITY=compact");
#pragma omp parallel for schedule(static)
    for (n=0; n<nMIC1; n++) {
         AC.x[n] = 0;
         AC.y[n] = 0;
         AC.z[n] = 0;
    }
}

    kmp_set_defaults("KMP_AFFINITY=compact");
#pragma omp parallel for schedule(static)
    for (n=0; n<nCPU; n++) {
         part.acc.x[n] = 0;
         part.acc.y[n] = 0;
         part.acc.z[n] = 0;
    }
#pragma offload_wait target(mic:0) wait(nMIC0)
#pragma offload_wait target(mic:1) wait(nMIC1)
//	printf ("nCPU=%d, nMIC0=%d, nMIC1=%d.\n", nCPU, nMIC0, nMIC1);
#else
#pragma omp parallel for schedule(static)
    for (n=0; n<la; n++) {
         part.acc.x[n] = 0;
         part.acc.y[n] = 0;
         part.acc.z[n] = 0;
    }
#endif   
    tstart = dtime();

#ifdef __PAPI
	PAPI_library_init(PAPI_VER_CURRENT);
	PAPI_create_eventset(&EventSet);
#ifdef __single_prec
	PAPI_add_event(EventSet, PAPI_SP_OPS);
#else
	PAPI_add_event(EventSet, PAPI_DP_OPS);
#endif
	PAPI_start(EventSet);
#endif

#ifdef __INTEL_OFFLOAD
#pragma offload target(mic:0) in (px[0:nPart]:alloc_if(0)) in (py[0:nPart]:alloc_if(0)) in (pz[0:nPart]:alloc_if(0)) out (ax0[0:nMIC0]:alloc_if(0)) out (ay0[0:nMIC0]:alloc_if(0)) out (az0[0:nMIC0]:alloc_if(0)) signal(nMIC0)
{
//	printf ("MIC0: la=%d, nCPU=%d, nMIC0=%d, nMIC1=%d.\n", la, nCPU, nMIC0, nMIC1);
	Array3 PA, PB, AC;
	PB.x = px;
	PB.y = py;
	PB.z = pz;
	PA.x = px + nCPU;
	PA.y = py + nCPU;
	PA.z = pz + nCPU;
	AC.x = ax0;
	AC.y = ay0;
	AC.z = az0;
	ppkernel(PA, nMIC0, PB, lb, eps, AC);
}
#pragma offload target(mic:1) in (px[0:nPart]:alloc_if(0)) in (py[0:nPart]:alloc_if(0)) in (pz[0:nPart]:alloc_if(0)) out (ax1[0:nMIC1]:alloc_if(0)) out (ay1[0:nMIC1]:alloc_if(0)) out (az1[0:nMIC1]:alloc_if(0)) signal(nMIC1)
{
//	printf ("MIC1: la=%d, nCPU=%d, nMIC0=%d, nMIC1=%d.\n", la, nCPU, nMIC0, nMIC1);
	Array3 PA, PB, AC;
	PB.x = px;
	PB.y = py;
	PB.z = pz;
	PA.x = px + nCPU + nMIC0;
	PA.y = py + nCPU + nMIC0;
	PA.z = pz + nCPU + nMIC0;
	AC.x = ax1;
	AC.y = ay1;
	AC.z = az1;
	ppkernel(PA, nMIC1, PB, lb, eps, AC);
}
//	printf ("la=%d, nCPU=%d, nMIC0=%d, nMIC1=%d.\n", la, nCPU, nMIC0, nMIC1);
	ppkernel(part.pos, nCPU, part.pos, lb, eps, part.acc);
#pragma offload_wait target(mic:0) wait(nMIC0)
#pragma offload_wait target(mic:1) wait(nMIC1)
#else
	ppkernel(part.pos, la, part.pos, lb, eps, part.acc);
#endif

#ifdef __PAPI
	PAPI_stop(EventSet, &value);
#endif

    tstop = dtime();
    ttime = tstop - tstart;
    printf("dtime = %lf \n", ttime);
#ifdef __PAPI
	printf("PAPI Counter: total %lld FLOPS, speed %f Mflops.\n", value, (double)(value/1048576) / ttime);
	PAPI_shutdown();
#endif
	printf("acc.x.y.z = %lf,%lf,%lf \n", part.acc.x[la/2-1],part.acc.y[la/2-1],part.acc.z[la/2-1]);
	printf("acc.x.y.z = %lf,%lf,%lf \n", part.acc.x[la/2],part.acc.y[la/2],part.acc.z[la/2]);
 /*
  
    gflops = (double)(1.0e-9 * LOOP_COUNT * MAXFLOPS_ITERS * FLOPSPERCALC);
    if ((ttime)>0.0)
    {
        printf("GFlops = %10.3lf, Secs = %10.3lf, GFlops per sec = %10.3lf\r\n", 
            gflops, ttime, gflops/ttime);
        
    }
   */
    return 0;
}



#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran3(long *idum)
{
    static int inext,inextp;
    static long ma[56];
    static int iff=0;
    long mj,mk;
    int i,ii,k;
    
    if (*idum < 0 || iff == 0) {
        iff=1;
        mj=MSEED-(*idum < 0 ? -*idum : *idum);
        mj %= MBIG;
        ma[55]=mj;
        mk=1;
        for (i=1; i<=54; i++) {
            ii=(21*i) % 55;
            ma[ii]=mk;
            mk=mj-mk;
            if (mk < MZ) mk += MBIG;
            mj=ma[ii];
        }
        for (k=1; k<=4; k++)
            for (i=1; i<=55; i++) {
                ma[i] -= ma[1+(i+30) % 55];
                if (ma[i] < MZ) ma[i] += MBIG;
            }
        inext=0;
        inextp=31;
        *idum=1;
    }
    if (++inext == 56) inext=1;
    if (++inextp == 56) inextp=1;
    mj=ma[inext]-ma[inextp];
    if (mj < MZ) mj += MBIG;
    ma[inext]=mj;
    return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC


