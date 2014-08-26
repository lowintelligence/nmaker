#include "parameter.h"
#include <assert.h>
#include <stdio.h>
#include <mpi.h>

void load_parameters(GlobalParam* gp, const char ParamFileName[]) {
    // to do
}

void setup_parameters_numproc(GlobalParam* gp);



void setup_parameters_numproc(GlobalParam* gp) {
    int numprocs;
 	unsigned long int parl;
    assert(gp != NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    gp->BoxSize = 100000.0;
    gp->OmegaM0 = 0.3;
    gp->OmegaX0 = 0.7;
    gp->Hubble0 = 0.7;

    /* system parameters*/
    gp->NumSocket = numprocs;
    gp->NumCpuPerSocket = 1;
    gp->NumMicPerSocket = 0;
    gp->NumThreadPerSocket = (int)(32/numprocs);
    if (gp->NumThreadPerSocket <1 )
        gp->NumThreadPerSocket = 1;


    gp->MemoBytesPerSocket = (long)(33664344064/numprocs);

    /* runtime parameters */
//    gp->NumBits = 13;
//    gp->NumBits = 9;
//      gp->NumBits = 8;
    gp->NumBits = 4;



    gp->PartBits = gp->NumBits + 4 ;
    gp->NumGridPerSide = 1<<gp->NumBits;//128;
    parl =(1<<gp->PartBits);//2048*2048*2048;
    gp->TotNumPart =parl*parl*parl;//2048*2048*2048;
	printf("NumPart=%ld,%d \n",gp->TotNumPart,gp->PartBits);
    gp->OpenAngle = 0.3;
    gp->CutRadius = 5.6;
    gp->SoftenLength = 0.001;
    sprintf(gp->WorkPath,"/home/qwang/") ;
    sprintf(gp->SnapName,"snapshot.dat");
    gp->NumFile = 1;

    /* status parameters */
    gp->redshift = 0.0;
}

void setup_parameters_two_sockets(GlobalParam* gp) {
    assert(gp != NULL);

    gp->BoxSize = 100000.0;
    gp->OmegaM0 = 0.3;
    gp->OmegaX0 = 0.7;
    gp->Hubble0 = 0.7;

    /* system parameters*/
    gp->NumSocket = 2;
    gp->NumCpuPerSocket = 1;
    gp->NumMicPerSocket = 0;
    gp->NumThreadPerSocket = 4;
    gp->MemoBytesPerSocket = 16832172032L;

    /* runtime parameters */
    gp->NumBits = 13;
    gp->NumGridPerSide = 256;
    gp->TotNumPart =256*256*256;
    gp->OpenAngle = 0.3;
    gp->CutRadius = 5.6;
    gp->SoftenLength = 0.001;
    sprintf(gp->WorkPath,"/home/qwang/") ;
    sprintf(gp->SnapName,"snapshot.dat");
    gp->NumFile = 1;

    /* status parameters */
    gp->redshift = 0.0;
}


void setup_parameters_one_socket(GlobalParam* gp) {
    assert(gp != NULL);

    gp->BoxSize = 100000.0;
    gp->OmegaM0 = 0.3;
    gp->OmegaX0 = 0.7;
    gp->Hubble0 = 0.7;

    /* system parameters*/
    gp->NumSocket = 1;
    gp->NumCpuPerSocket = 1;
    gp->NumMicPerSocket = 0;
    gp->NumThreadPerSocket = 8;
    gp->MemoBytesPerSocket = 33664344064L;

    /* runtime parameters */
    gp->NumBits = 13;
    gp->NumGridPerSide = 256;
    gp->TotNumPart =256*256*256;
    gp->OpenAngle = 0.3;
    gp->CutRadius = 5.6;
    gp->SoftenLength = 0.001;
    sprintf(gp->WorkPath,"/home/qwang/") ;
    sprintf(gp->SnapName,"snapshot.dat");
    gp->NumFile = 1;

    /* status parameters */
    gp->redshift = 0.0;
}

void setup_parameters(GlobalParam* gp) {
    setup_parameters_numproc(gp) ;
//    setup_parameters_two_sockets(gp);
//    setup_parameters_one_socket(gp);
}

void print_parameters(GlobalParam* gp) {
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (myid == 0) {
        printf("\n * * * * *  global parameters * * * * * \n");
        printf("boxsize = %.2f OmegaM = %.2f OmegaX = %.2f hubble = %.2f\n",gp->BoxSize, gp->OmegaM0 ,gp->OmegaX0,gp->Hubble0 );
        printf("NumSocket = %d NumThreadPerSock = %d\n ", gp->NumSocket, gp->NumThreadPerSocket);
        printf(" NumParticles %d\n", gp->TotNumPart);
        printf(" * * * * *  global parameters * * * * * \n\n");
        fflush(stdout);
    }
}

