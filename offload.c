/*
 * =====================================================================================
 *
 *       Filename:  off.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07/22/2014 07:16:20 PM
 *       Revision:  none
 *       Compiler:  gcc
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

#define N 512

double dtime()
{
	double tseconds = 0.0;
	struct timeval mytime;	
	gettimeofday(&mytime, (struct timezone*)0);
	tseconds = (double)(mytime.tv_sec + mytime.tv_usec*1.0e-6);
	return (tseconds );
}

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

#ifdef __INTEL_OFFLOAD
__attribute__ ((target(mic)))
#endif
void* sayhello(void *pB)
{
	Block B = *((Block*)pB);
	B.B[N+B.tid] = B.A[N+B.tid] + 1.0;
#ifdef __MIC__
	printf("++Hello in thread %d of block %d bid %d on MIC.\n", B.tid, B.tid/B.bsize, B.blockid);
#else
	printf("==Hello in thread %d of block %d bid %d on CPU.\n", B.tid, B.tid/B.bsize, B.blockid);
#endif

//	pthread_barrier_wait(B.bglobal);
}

#ifdef __INTEL_OFFLOAD
__attribute__ ((target(mic)))
#endif
int threaded(int tnum, int ts, int te, double *A, double *B)
{
	Block *pth;
	pthread_t *thread;

	pth = (Block*)malloc(sizeof(Block)*tnum);
	thread = (pthread_t*)malloc(sizeof(pthread_t)*tnum);
	pthread_attr_t attr;
	pthread_barrier_t bbar, bglobal;

	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	pthread_barrier_init(&bglobal, NULL, te-ts+2);

	int i;
	for (i=ts; i<te+1; i++)
	{
		pth[i].blockid = i%4;
		pth[i].bsize = 4;
		pth[i].tid = i;
		pth[i].tall = tnum;
		A[N+i] = (double)i*2.0;
		pth[i].A = A;
		pth[i].B = B;
		pth[i].bglobal = &bglobal;
		pthread_create(&thread[i], &attr, sayhello, (void*)&pth[i]);
	}
	pthread_attr_destroy(&attr);

#ifdef __MIC__
	printf("**Hello in main thread on MIC.\n");
#else
	printf("--Hello in main thread on CPU.\n");
#endif

//	pthread_barrier_wait(&bglobal);

	for (i=ts; i<te+1; i++)
	{
		pthread_join(thread[i], NULL);
	}

}

int main(int argc, char **argv)
{
	int tnum, haha1, haha2;
	double ts, te;
	double A[N*N*N], B[N*N*N];
	double *pa=A, *pb=B;

	tnum = 16;

	if(argc>1)
	{
		tnum=atoi(argv[1]);
		tnum=tnum-tnum%4;
	}
	int micno;
	micno = 0;
	if(argc>2)
	{
		micno=atoi(argv[2]);
	}

#ifdef __INTEL_OFFLOAD
#pragma offload target (mic:micno) in (pa:length(N*N*N) free_if(0)) in (pb:length(N*N*N) free_if(0)) signal(haha1)
{	
	threaded(tnum, 2, 15, pa, pb);
}
//#pragma offload target (mic:1) in (pa:length(N*N*N) free_if(0)) in (pb:length(N*N*N) free_if(0)) signal(haha2)
//{	
//	threaded(tnum, 9, 15, pa, pb);
//}
#endif
	threaded(tnum, 0, 1, pa, pb);
#ifdef __INTEL_OFFLOAD
#pragma offload_wait target(mic:micno) wait(haha1)
//#pragma offload_wait target(mic:1) wait(haha2)
#endif
	ts = dtime();
	int i;
	for (i=0; i<10; i++)
	{
#ifdef __INTEL_OFFLOAD
#pragma offload target (mic:micno) in (pa:length(N*N*N) free_if(0)) in (pb:length(N*N*N) free_if(0)) signal(haha1)
{	
	threaded(tnum, 2, 15, pa, pb);
}
//#pragma offload target (mic:1) in (pa:length(N*N*N) free_if(0)) in (pb:length(N*N*N) free_if(0)) signal(haha2)
//{	
//	threaded(tnum, 9, 15, pa, pb);
//}
#endif
	threaded(tnum, 0, 1, pa, pb);
#ifdef __INTEL_OFFLOAD
#pragma offload_wait target(mic:micno) wait(haha1)
//#pragma offload_wait target(mic:1) wait(haha2)
#endif
	}
	te = dtime();
	printf("Total time: %.2f, speed %.3f MB/s.\n", (te-ts)/10, ((double)40*N*N*N*sizeof(double)/1048576)/(te-ts));

	return 0;
}
