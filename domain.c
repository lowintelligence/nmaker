#include "domain.h"
#include "parameter.h"
#include "stdlib.h"

/* for certain simultion contained 256^3 particles */

void load_particle_into_domain(Domain* dp, int myid, int numprocs) {
    int n;
    double box;
    int  n_start = (int)(myid*256*256*256)/numprocs;
    int  n_count = (int)(256*256*256)/numprocs;

    if (myid == numprocs-1)
        n_count = 256*256*256 - n_start;
    char fname[80]="/home/qwang/data/z0_256_100mpch/snapshot256_100mpch_z0.dat";

    dp->NumPart = n_count;
    dp->Part = (Body*) malloc(sizeof(Body)*(n_count) );

    read_Particle(dp->Part , fname, n_start, n_count);
//   printf(" [%d/%d] %15ld %15ld \n", myid, numprocs, n_start, n_count);

    dp->NumDom = numprocs;
    dp->DomList = (DomainInfo*)malloc(sizeof(DomainInfo)*(dp->NumDom) );

//    printf("%lf %lf %lf\n", dp->Part[0].pos[0], dp->Part[0].pos[1], dp->Part[0].pos[2]);
//    printf("%lf %lf %lf\n", dp->Part[n_count-1].pos[0], dp->Part[n_count-1].pos[1], dp->Part[n_count-1].pos[2]);

    for (n=0; n<numprocs; n++) {
        dp->DomList[n].NumPart = (int)(256*256*256)/numprocs;
    }
    dp->DomList[numprocs-1].NumPart = 256*256*256-(int)((numprocs-1)*256*256*256)/numprocs;
    dp->rank = myid;
    dp->thisdom = &(dp->DomList[myid]);
    /*
       printf("[%d] ", myid);
       for (n=0; n<numprocs; n++) {
           printf("%d ",dp->DomList[n].NumPart );
       }
       printf("[%d] \n", myid);
       */
}

void init_particle_into_domain(Domain* dp, GlobalParam* gp,int myid, int numprocs) {
    int n;
    double box;
	int nsize=1<<gp->NumBits;
	unsigned long int nparl=1<<gp->PartBits;
    unsigned long int n_start = myid*nparl*nparl*nparl/numprocs;
    int  n_count = (int)(nparl/numprocs*nparl*nparl);

    if (myid == numprocs-1)
        n_count = nparl*nparl*nparl - n_start;
//    char fname[80]="/home/qwang/data/z0_256_100mpch/snapshot256_100mpch_z0.dat";

    dp->NumPart = n_count;
    dp->Part = (Body*) malloc(sizeof(Body)*(n_count) );
    srand48(8888*myid);
	double velocity=5000.0;
    for (n=0; n<n_count; n++)
    {
		dp->Part[n].pos[0]=drand48()*gp->BoxSize;
		dp->Part[n].pos[1]=drand48()*gp->BoxSize;
		dp->Part[n].pos[2]=drand48()*gp->BoxSize;
		if(dp->Part[n].pos[0]>99999.97f)
			dp->Part[n].pos[0]=99999.97f;
		if(dp->Part[n].pos[1]>99999.97f)
			dp->Part[n].pos[1]=99999.97f;
		if(dp->Part[n].pos[2]>99999.97f)
			dp->Part[n].pos[2]=99999.97f;
//		dp->Part[n].vel[0]=drand48()*velocity;
//		dp->Part[n].vel[1]=drand48()*velocity;
//		dp->Part[n].vel[2]=drand48()*velocity;
		dp->Part[n].ID=n_start+n+1;
		dp->Part[n].mass=1.0;
    }
//    read_Particle(dp->Part , fname, n_start, n_count);
	
   printf(" [%d/%d] %15ld %15ld \n", myid, numprocs, n_start, n_count);

    dp->NumDom = numprocs;
    dp->DomList = (DomainInfo*)malloc(sizeof(DomainInfo)*(dp->NumDom) );

//    printf("%lf %lf %lf\n", dp->Part[0].pos[0], dp->Part[0].pos[1], dp->Part[0].pos[2]);
//    printf("%lf %lf %lf\n", dp->Part[n_count-1].pos[0], dp->Part[n_count-1].pos[1], dp->Part[n_count-1].pos[2]);

    for (n=0; n<numprocs; n++) {
        dp->DomList[n].NumPart = n_count;
    }
    dp->rank = myid;
    dp->thisdom = &(dp->DomList[myid]);
    /*
       printf("[%d] ", myid);
       for (n=0; n<numprocs; n++) {
           printf("%d ",dp->DomList[n].NumPart );
       }
       printf("[%d] \n", myid);
       */
}

void free_domain(Domain* dp) {
    free(dp->Part);
    free(dp->DomList);
//    free_subcuboid(dp->cuboid);
}

