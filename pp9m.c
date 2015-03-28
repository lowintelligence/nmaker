#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <malloc.h>
#include <omp.h>

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
    part.pos.x = (Real*)memalign(ALIGNCNT, sizeof(Real)*nPart);
    part.pos.y = (Real*)memalign(ALIGNCNT, sizeof(Real)*nPart);
    part.pos.z = (Real*)memalign(ALIGNCNT, sizeof(Real)*nPart);
    part.acc.x = (Real*)memalign(ALIGNCNT, sizeof(Real)*nPart);
    part.acc.y = (Real*)memalign(ALIGNCNT, sizeof(Real)*nPart);
    part.acc.z = (Real*)memalign(ALIGNCNT, sizeof(Real)*nPart);
    part.mass  = (Real*)memalign(ALIGNCNT, sizeof(Real)*nPart);
	
#ifdef NMK_VERIFY
	srand48(0);
#endif
    for (n=0; n<nPart; n++) {
#ifdef NMK_VERIFY
         part.pos.x[n] = drand48()*24998.0+87501.0;
         part.pos.y[n] = drand48()*24998.0+87501.0;
         part.pos.z[n] = drand48()*24998.0+87501.0;
#else
	     part.pos.x[n] = box*ran3(&seed);
         part.pos.y[n] = box*ran3(&seed);
         part.pos.z[n] = box*ran3(&seed);
#endif
         part.acc.x[n] = (Real) 0.0;
         part.acc.y[n] = (Real) 0.0;
         part.acc.z[n] = (Real) 0.0;
    }
    for (n=0; n<nPart; n++) {
#ifdef NMK_VERIFY
		 part.mass[n] = (3.0*(0.25)*0.01)/1080886.296;
	     part.mass[n] *= pow(200000.0,3.0)/(nPart);
#else
         part.mass[n] = box * 0.25;
#endif
	}
#if (defined NMK_VERIFY) && (defined NMK_GEN_VER_FILE)
	FILE *fp = fopen("pp.ori", "w+");
	for (n=0; n<nPart; n++)
	{
		fprintf(fp, "%d  %f  %f  %f\n", n+1, part.pos.x[n], part.pos.y[n], part.pos.z[n]);
	}
	fclose(fp);
#endif
}


int main(int argc, char *argv[])
{
    int n, m, nPart, lb, la;
    double tstart, tstop, ttime;
    double gflops = 0.0;
    float a = 1.1;
#ifdef NMK_VERIFY	
	FILE *fp;
#endif
//	Real eps = 1.0/100.0/40.0;
	Real eps = 0.00026;
	Real dx, dy, dz, dr2, dr3;
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
    
    printf("Initializing\r\n");
#ifdef NMK_VERIFY
    initialize_particle(nPart, 200000.0);
#else
    initialize_particle(nPart, 1.0);
#endif

#if (defined NMK_VERIFY) && (defined NMK_GEN_VER_FILE)
    printf("Starting Compute\n");
    
    tstart = dtime();
    
#pragma omp parallel for schedule(static) private(dx, dy, dz, dr2, dr3)
    for ( n=0; n<la; n++ ) {
		Real accx = 0.0;
		Real accy = 0.0;
		Real accz = 0.0;
        for ( m=0; m<lb; m++) {
            dx = part.pos.x[m] -  part.pos.x[n];
            dy = part.pos.y[m] -  part.pos.y[n];
            dz = part.pos.z[m] -  part.pos.z[n];
            
            dr2 = eps + dx*dx + dy*dy + dz*dz;
            dr3 = (Real) 1.0 / sqrt(dr2*dr2*dr2);
            
            accx += part.mass[m]*dx*dr3;
            accy += part.mass[m]*dy*dr3;
            accz += part.mass[m]*dz*dr3;
        }
		part.acc.x[n] = accx;
		part.acc.y[n] = accy;
		part.acc.z[n] = accz;
    }

    tstop = dtime();
    
    ttime = tstop - tstart;
    printf("dtime = %lf \n", ttime);
    
	fp = fopen("ppn.acc", "w+");
	double G1=43007.1;
	for (n=0; n<nPart; n++)
	{
		fprintf(fp, "%d  %f  %f  %f\n", n+1, part.acc.x[n]*G1, part.acc.y[n]*G1, part.acc.z[n]*G1);
		part.acc.x[n] = 0;
		part.acc.y[n] = 0;
		part.acc.z[n] = 0;
	}
	fclose(fp);
#endif

 //   omp_set_num_threads(32);
  //  kmp_set_defaults("KMP_AFFINITY=compact");
//    kmp_set_defaults("KMP_AFFINITY=scatter");

#ifdef __INTEL_OFFLOAD
	int nCPU;
	int nMIC0;
	int nMIC1;

//	nCPU = 0;
//	nMIC0 = nPart / 2;
	nCPU = (int) ((double)(nPart + CLCNT - 1) / (double)CLCNT / 6.76) * CLCNT; 
	nMIC0 = ((int) ((double)nCPU * 2.88) + CLCNT - 1) / CLCNT  * CLCNT; 
	nMIC1 = nPart - nCPU - nMIC0;
	printf ("nCPU=%d, nMIC0=%d, nMIC1=%d.\n", nCPU, nMIC0, nMIC1);

	Real *px = part.pos.x;
	Real *py = part.pos.y;
	Real *pz = part.pos.z;
	Real *mass = part.mass;
	Real *ax = part.acc.x;
	Real *ay = part.acc.y;
	Real *az = part.acc.z;
	Real *ax0 = part.acc.x + nCPU;
	Real *ay0 = part.acc.y + nCPU;
	Real *az0 = part.acc.z + nCPU;
	Real *ax1 = part.acc.x + nCPU + nMIC0;
	Real *ay1 = part.acc.y + nCPU + nMIC0;
	Real *az1 = part.acc.z + nCPU + nMIC0;
#pragma offload target(mic:0) in (px[0:nPart]:alloc_if(1) free_if(0)) in (py[0:nPart]:alloc_if(1) free_if(0)) in (pz[0:nPart]:alloc_if(1) free_if(0)) in (mass[0:nPart]:alloc_if(1) free_if(0)) out (ax0[0:nMIC0]:alloc_if(1) free_if(0)) out (ay0[0:nMIC0]:alloc_if(1) free_if(0)) out (az0[0:nMIC0]:alloc_if(1) free_if(0)) signal(nMIC0)
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
#pragma offload target(mic:1) in (px[0:nPart]:alloc_if(1) free_if(0)) in (py[0:nPart]:alloc_if(1) free_if(0)) in (pz[0:nPart]:alloc_if(1) free_if(0)) in (mass[0:nPart]:alloc_if(1) free_if(0)) out (ax1[0:nMIC1]:alloc_if(1) free_if(0)) out (ay1[0:nMIC1]:alloc_if(1) free_if(0)) out (az1[0:nMIC1]:alloc_if(1) free_if(0)) signal(nMIC1)
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
#pragma offload target(mic:0) in (px[0:nPart]:alloc_if(0)) in (py[0:nPart]:alloc_if(0)) in (pz[0:nPart]:alloc_if(0)) in (mass[0:nPart]:alloc_if(0)) out (ax0[0:nMIC0]:alloc_if(0)) out (ay0[0:nMIC0]:alloc_if(0)) out (az0[0:nMIC0]:alloc_if(0)) signal(nMIC0)
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
	ppmkernel(PA, nMIC0, PB, mass, lb, eps, AC);
}
#pragma offload target(mic:1) in (px[0:nPart]:alloc_if(0)) in (py[0:nPart]:alloc_if(0)) in (pz[0:nPart]:alloc_if(0)) in (mass[0:nPart]:alloc_if(0)) out (ax1[0:nMIC1]:alloc_if(0)) out (ay1[0:nMIC1]:alloc_if(0)) out (az1[0:nMIC1]:alloc_if(0)) signal(nMIC1)
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
	ppmkernel(PA, nMIC1, PB, mass, lb, eps, AC);
}
//	printf ("la=%d, nCPU=%d, nMIC0=%d, nMIC1=%d.\n", la, nCPU, nMIC0, nMIC1);
	ppmkernel(part.pos, nCPU, part.pos, mass, lb, eps, part.acc);
#pragma offload_wait target(mic:0) wait(nMIC0)
#pragma offload_wait target(mic:1) wait(nMIC1)
#else
	ppmkernel(part.pos, la, part.pos, part.mass, lb, eps, part.acc);
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
#ifdef NMK_VERIFY
	fp = fopen("pp.acc", "w+");
	double G=43007.1;
	for (n=0; n<nPart; n++)
	{
		fprintf(fp, "%d  %f  %f  %f\n", n+1, part.acc.x[n]*G, part.acc.y[n]*G, part.acc.z[n]*G);
	}
	fclose(fp);
#endif
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


