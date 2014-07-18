#include "partition.h"


typedef struct {
    float pos[3];
    float w;
} Sample;


Sample* gather_total_sample(Domain* dp, int *total_sample)
{
    int n, total, start, end, myid, nParticle ;
    int nDomain = dp->NumDom;
    int sample_gap = 16;//1.0;
    Sample *sample;
    Sample *sendBuffer;
    int *nSample;
    int *dispSample;
    code *codeSplit;

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    nSample    = (int*)malloc( (int) sizeof(int)*nDomain );
    dispSample = (int*)malloc( (int) sizeof(int)*nDomain );

    total = 0;
    for (n=0; n<nDomain; n++) {
        nSample[n] = dp->DomList[n].NumPart/sample_gap;
        total  += nSample[n];
    }
    *total_sample = total;

    dispSample[0] = 0;
    for (n=1; n<nDomain; n++)
        dispSample[n] = nSample[n-1] + dispSample[n-1];

    sendBuffer = (Sample*)malloc( (int)sizeof(Sample)*nSample[myid] );
    sample = (Sample*)malloc( (int)sizeof(Sample)*total );

    n=0;
    nParticle = nSample[myid];

    while( n < nParticle ) {
        sendBuffer[n].pos[0] = dp->Part[n].pos[0];
        sendBuffer[n].pos[1] = dp->Part[n].pos[1];
        sendBuffer[n].pos[2] = dp->Part[n].pos[2];
        sendBuffer[n].w = 1.0f;
        n+=sample_gap;
    }//Cao! n+=nblock and sample_ratio=0.1
    MPI_Datatype mpi_sample_type;
    MPI_Type_contiguous(sizeof(Sample), MPI_CHAR, &mpi_sample_type );
    MPI_Type_commit(&mpi_sample_type);

    MPI_Barrier(MPI_COMM_WORLD );
    MPI_Allgatherv( sendBuffer, nSample[myid], mpi_sample_type, sample,
                    nSample, dispSample, mpi_sample_type, MPI_COMM_WORLD );

    MPI_Type_free(&mpi_sample_type);

    /*
        n=0;
        printf("[%d] n=%10d, %f %f %f\n", myid, n, sendBuffer[n].pos[0], sendBuffer[n].pos[1], sendBuffer[n].pos[2]);
        n=nParticle-1;
        printf("[%d] n=%10d, %f %f %f\n", myid, n, sendBuffer[n].pos[0], sendBuffer[n].pos[1], sendBuffer[n].pos[2]);
        if ( 0==myid) {
            printf(" after all gatherc \n");
            n=0;
            printf("n=%10d, %f %f %f\n",  n, sample[n].pos[0], sample[n].pos[1], sample[n].pos[2]);
            n += nSample[0]-1;
            printf("n=%10d, %f %f %f\n",  n, sample[n].pos[0], sample[n].pos[1], sample[n].pos[2]);
            n++;
            printf("n=%10d, %f %f %f\n", n, sample[n].pos[0], sample[n].pos[1], sample[n].pos[2]);
            n += nSample[1]-1;
            printf("n=%10d, %f %f %f\n",  n, sample[n].pos[0], sample[n].pos[1], sample[n].pos[2]);
            n++;
            printf("n=%10d, %f %f %f\n", n, sample[n].pos[0], sample[n].pos[1], sample[n].pos[2]);
            n += nSample[2]-1;
            printf("n=%10d, %f %f %f\n",  n, sample[n].pos[0], sample[n].pos[1], sample[n].pos[2]);
        }
    */

    free(sendBuffer);
    free(nSample);
    free(dispSample);

    return sample;
}



typedef struct {
    code key;
    float weight;
} CodeSort;

int ascend(const void *a ,const void * b) {
    return ( (*(CodeSort*)a).key > (*(CodeSort*)b).key ) ? 1 : 0;
}

typedef struct {
    int rank;
    int nPart;
    int nBits;
    Sample *part;
    double duration;
    double boxsize;
    CodeSort *buff;
} TaskSortData;


void* partial_sort_key(void* data) {
    int n, x, y, z;
    double real2int;
    TaskSortData *p = (TaskSortData*)data;
    double BOXSIZE = p->boxsize;
    int NBITS = p->nBits;

    clock_t time_start, time_end;
    time_start = clock();
#ifndef SILENCE
    printf("task[%d]:%s\n", p->rank, __func__);
#endif
    real2int = (double)(1<<NBITS)/BOXSIZE;
    long int seed = 1243981;
    for (n=0; n< p->nPart; n++) {
        x = (int)(p->part[n].pos[0]*real2int);
        y = (int)(p->part[n].pos[1]*real2int);
        z = (int)(p->part[n].pos[2]*real2int);
        p->buff[n].key = coding_3d_filling_curve(x, y, z, NBITS);
        p->buff[n].weight = p->part[n].w;

    }
    qsort(p->buff, p->nPart, sizeof(CodeSort), ascend);

    time_end = clock();
    p->duration = (double)(time_end-time_start)/CLOCKS_PER_SEC;
    return NULL;
}

void determine_partition_code(Sample* sample, int NumPart, int NumThread, double box, int nbits, Domain *dp)
{
    int n, i, id;

    int *curr, *lend;
    pthread_t *task;
    TaskSortData *sortdata;
    CodeSort *codebuff, *codelist;
    code min, CODE_MAXIMUM, CODE_MINIMUM;

    CODE_MINIMUM = 0L;
    CODE_MAXIMUM = (1L<<(nbits*3));

    assert(NumThread >0 || NumPart>NumThread) ;

    task = (pthread_t*) malloc((int)sizeof(pthread_t)*NumThread);
    sortdata =(TaskSortData*) malloc((int)sizeof(TaskSortData)*NumThread);
    codelist = (CodeSort*)malloc((int)sizeof(CodeSort)*NumPart);

    for (i=0; i<NumThread; i++) {
        sortdata[i].rank = i;
        sortdata[i].nPart = (i+NumPart)/NumThread;
        if (0==i) {
            sortdata[i].part = sample;
            sortdata[i].buff = codelist;
        }
        else {
            sortdata[i].part = sortdata[i-1].part + sortdata[i-1].nPart ;
            sortdata[i].buff = sortdata[i-1].buff + sortdata[i-1].nPart ;
        }
        sortdata[i].duration = 0.0;
        sortdata[i].nBits = nbits;
        sortdata[i].boxsize = box;
    }
    int myid, numprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
//   MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

#ifndef SILENCE
    if ( 0 == myid)
        printf("npart_sample= %d  in %d Domain(s) \n",  sortdata[i-1].nPart, dp->NumDom);
#endif

    for (i=0; i<NumThread; i++) {

        pthread_create(&task[i], NULL, partial_sort_key, &sortdata[i]);
    }

    for (i=0; i<NumThread; i++) {
        pthread_join(task[i], NULL);
        sortdata[i].duration /= NumThread;

#ifndef SILENCE
        printf(" ->time[%d] = %.3f \n", i, sortdata[i].duration);
#endif
    }

    /*
    for (n=0; n<NumPart; n++)
        printf("%d %ld\n", n, codelist[n].key);
    */

    /* merge sort */
    if ( 1 < NumThread) {
        curr = (int*)malloc((int)sizeof(int)*NumThread);
        lend = (int*)malloc((int)sizeof(int)*NumThread);
        codebuff = (CodeSort*)malloc((int)sizeof(CodeSort)*NumPart);
        curr[0] = 0;
        lend[NumThread-1] = NumPart;
        for (i=1; i<NumThread; i++) {
            curr[i] = curr[i-1] + sortdata[i-1].nPart ;
            lend[i-1] = curr[i];
        }

        for (i=0; i<NumThread; i++) {
            printf("curr[%d] = %d , lend[%d] = %d\n", i, curr[i], i ,lend[i]);
        }
        n = 0;
        while (n<NumPart) {
            min = 0xFFFFFFFFFFFFFFFL;
            for (i=0; i<NumThread; i++)
                if (curr[i]<lend[i] && codelist[curr[i]].key<min ) {
                    min=codelist[curr[i]].key;
                    id = i;
                }
            codebuff[n].key = min;
            codebuff[n].weight = codelist[curr[id]].weight;

            curr[id]++;
            n++;
        }
        for (n=0; n<NumPart; n++) {
            codelist[n].key = codebuff[n].key;
            codelist[n].weight = codebuff[n].weight;
        }

        for (n=1; n<NumPart; n++) {
            if (codelist[n].key < codelist[n-1].key) {
                printf("error %s\n", __func__);
                exit(0);
            }
        }

        free(curr);
        free(lend);
        free(codebuff);
    }
    /* merge sort */

    for (n=1; n<NumPart; n++) {
        codelist[n].weight += codelist[n-1].weight;
    }
    int iDom, nGroup = dp->NumDom; //= nDomain;
    float next_weight, total_weight = codelist[NumPart-1].weight;
    next_weight = total_weight/nGroup;
    dp->DomList[0].LowerHilbertKey = CODE_MINIMUM;
    dp->DomList[(dp->NumDom)-1].UpperHilbertKey = CODE_MAXIMUM;
    iDom = 0;
    for (n=0; n<NumPart; n++) {
        if ( codelist[n].weight > next_weight ) {
            next_weight += (total_weight+iDom)/nGroup;
            dp->DomList[iDom].UpperHilbertKey = codelist[n].key;
            dp->DomList[iDom+1].LowerHilbertKey = codelist[n].key;
            iDom++;
        }
    }

#ifndef SILENCE
    if ( 0 == myid) {
        for (iDom = 0; iDom < nGroup; iDom ++) {
            printf(" split - Domain[%d] %lld %lld\n", iDom, dp->DomList[iDom].LowerHilbertKey ,dp->DomList[iDom].UpperHilbertKey);
        }
    }
#endif
    /*
       FILE *fd = fopen("ratio.dat","w");
       for (iDom = 0; iDom < nGroup; iDom ++) {
           double radio = (double)(dp->DomList[iDom].UpperHilbertKey - dp->DomList[iDom].LowerHilbertKey)/CODE_MAXIMUM;
           fprintf(fd,"%d %e\n", iDom, radio);
       }
       fclose(fd);
       */
    free(codelist);
    free(task);
    free(sortdata);

    return;
}

/********************** coding this domain *********************/
typedef struct {
    int rank;
    int nPart; // total particle in this domain
    int nBits;
    double boxsize;
    Body* part; //
} TaskCoding;

void* partial_coding_domain(void* param)
{
    TaskCoding *p = (TaskCoding*)param;

    int n, x, y, z;
    int nbits = p->nBits;
    double real2int;

    real2int = (double)(1<<nbits)/(p->boxsize);
    for (n=0; n< (p->nPart); n++) {
        x = (int)(p->part[n].pos[0]*real2int);
        y = (int)(p->part[n].pos[1]*real2int);
        z = (int)(p->part[n].pos[2]*real2int);
        p->part[n].key = coding_3d_filling_curve(x, y, z, nbits);
    }
    return NULL;
}

void coding_domain_particle(Body* part, int NumPart, int NumThread, double box, int nbits)
{
    int myid, i, cnt;
    pthread_t *task;
    TaskCoding *taskcoding;
    task = (pthread_t*) malloc((int)sizeof(pthread_t)*NumThread);
    taskcoding = (TaskCoding*)malloc((int)sizeof(TaskCoding)*NumThread);

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);


    for (cnt=0, i=0; i<NumThread; i++) {
        taskcoding[i].rank    = i;
        taskcoding[i].nPart   = (i+NumPart)/NumThread;
        taskcoding[i].part    = &part[cnt];
        taskcoding[i].boxsize = box;
        taskcoding[i].nBits   = nbits;
        cnt += taskcoding[i].nPart;
    }

    for (i=0; i<NumThread; i++) {

        pthread_create(&task[i], NULL, partial_coding_domain, &taskcoding[i]);
    }

    for (i=0; i<NumThread; i++) {
        pthread_join(task[i], NULL);
    }
    free(task);
    free(taskcoding);
}
/**********************  coding this domain *********************/

/********************** coding & tagging this domain *********************/
typedef struct {
    int rank;
    int nPart; // total particle in this domain
    int nBits; // for coding
    int nLevel;  // for tagging
    double boxsize;
    Body *part;
    code *split; // must be a power of 2
} TaskCodeTag;


/* tagging the destination of particle */
void* partial_coding_tagging_domain(void* param)
{
    TaskCodeTag *p = (TaskCodeTag*)param;
    code *split = (p->split);
    int n, x, y, z;
    int tag, l;
    int level = 1<< (p->nLevel);
    int nbits = p->nBits;
    double real2int;
    code key;

    level >>= 1; // 4 ... 2 ... 1 in case of nLevel = 3
    printf("func:%s: level = %d\n" , __func__, level);
    real2int = (double)(1<<nbits)/(p->boxsize);
    for (n=0; n< (p->nPart); n++) {
        x = (int)(p->part[n].pos[0]*real2int);
        y = (int)(p->part[n].pos[1]*real2int);
        z = (int)(p->part[n].pos[2]*real2int);
        key = coding_3d_filling_curve(x, y, z, nbits);

        l = tag = level ; // 4
        while ( l >>= 1 ) {
            if (key < split[tag]) {
                tag -= l;
            }
            else {
                tag += l;
            }
        }
        if (key < split[tag]) {
            tag -= 1;
        }
        p->part[n].key = key;
        p->part[n].tag = tag;
    }

    return NULL;
}

void coding_tagging_domain_particle(Body* part, int NumPart, int NumThread, double box, int nbits, int nSplit, Domain *dp)
{
    int myid, numprocs, i, cnt;
    int nLevel, splitLength;
    code *split;

    pthread_t *task;
    TaskCodeTag *taskcodetag;
    task = (pthread_t*) malloc((int)sizeof(pthread_t)*NumThread);
    taskcodetag = (TaskCodeTag*)malloc((int)sizeof(TaskCodeTag)*NumThread);

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);



    for (nLevel = 0; nSplit > (1<<nLevel); nLevel++) ;

    splitLength = 1<<nLevel;

    printf(" myid = %d/%d, nLevel = %d, splitLength =%d, nPart = %ld\n",
           myid, numprocs, nLevel, splitLength, dp->NumPart);
    split = (code*) malloc( (int) sizeof(code) * splitLength );
    for (i=0; i<splitLength; i++) {
        if (i<nSplit)
            split[i] = dp->DomList[i].LowerHilbertKey;
        else
            split[i] =  dp->DomList[nSplit-1].UpperHilbertKey;

        //      printf("%d %lld\n", i, split[i]);
    }

    for (cnt=0, i=0; i<NumThread; i++) {
        taskcodetag[i].rank    = i;
        taskcodetag[i].nPart   = (i+NumPart)/NumThread;
        taskcodetag[i].part    = &part[cnt];
        taskcodetag[i].boxsize = box;
        taskcodetag[i].nBits   = nbits;
        taskcodetag[i].nLevel  = nLevel;
        taskcodetag[i].split   = split;
        cnt += taskcodetag[i].nPart;
    }

    for (i=0; i<NumThread; i++) {
        pthread_create(&task[i], NULL, partial_coding_tagging_domain, &taskcodetag[i]);
    }

    for (i=0; i<NumThread; i++) {
        pthread_join(task[i], NULL);
    }
    free(split);
    free(task);
    free(taskcodetag);
}
/********************** coding & tagging this domain *********************/



/********************** partition *********************/
typedef struct {
    int rank;
    int *CountPart; // particle counting in the buffer to domain[i]
    int DomNumPart; // total particle in this domain
    int ExpRank;
    Body* ExpPart; // export those particles to ExpRank
    Body* DomPart; //
    code LowerBound;
    code UpperBound;
} TaskPartition;

void partition_sorted_particles(Domain *dp, int NumPart, int NumThread, double box, int nbits)
{
    int myid, i;
    pthread_t *task;
    TaskPartition *taskpartition;
    taskpartition=(TaskPartition*)malloc((int)sizeof(TaskPartition)*NumThread);

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

}
/********************** partition *********************/

/********************** partition *********************/
void partition_domains(Domain *dp, int NumThread)
{
    int i, n, idx, disp;
    int myid, numprocs;
    int nDomain = dp->NumDom;
    int nPart = dp->NumPart;
    Body *part = dp->Part;
    Body *sendbuff;
    int *sendcount;
    int *recvcount;
    int *senddispl;
    int *recvdispl;
    int totrecv, totsend;

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    /*** compute displacement ***/

    sendcount = (int*)malloc(sizeof(int)*nDomain);
    recvcount = (int*)malloc(sizeof(int)*nDomain);
    senddispl = (int*)malloc(sizeof(int)*nDomain);
    recvdispl = (int*)malloc(sizeof(int)*nDomain);

    for (i=0; i<nDomain; i++ ) {
        sendcount[i] = 0;
        recvcount[i] = 0;
        senddispl[i] = 0;
        recvdispl[i] = 0;
    }

    for (n=0; n<nPart; n++) {
        idx = part[n].tag;
        sendcount[idx] ++;
    }
    for (i=1, senddispl[0] = 0; i<nDomain; i++ ) {
        senddispl[i] = senddispl[i-1] + sendcount[i-1];
    }

    //  for (i=0; i<nDomain; i++ ){
    //      printf("[%d]sendcount.%d = %d\n", myid, i, sendcount[i]);
    //      printf("[%d]senddispl.%d = %d\n", myid, i, senddispl[i]);
    //  }
    /*** build buffer by tag ***/
    sendbuff = (Body*) malloc(sizeof(Body)*nPart);

    for (n=0; n<nPart; n++) {
        idx = part[n].tag;
        disp = senddispl[idx];
        sendbuff[disp] = part[n];
        senddispl[idx]++;
    }
    MPI_Barrier( MPI_COMM_WORLD );

    /*    printf("[%d] nDom = %d, NumPart = %ld\n", myid, dp->NumDom, dp->NumPart); */

    /* alltoall for recv_cont */

    MPI_Alltoall(sendcount, 1, MPI_INT, recvcount, 1, MPI_INT, MPI_COMM_WORLD);


    for (i=1, senddispl[0] = 0; i<nDomain; i++ ) {
        senddispl[i] = senddispl[i-1] + sendcount[i-1];
    }
    for (i=1, recvdispl[0] = 0; i<nDomain; i++ ) {
        recvdispl[i] = recvdispl[i-1] + recvcount[i-1];
    }
    totrecv = recvdispl[nDomain-1] + recvcount[nDomain-1] ;
    totsend = senddispl[nDomain-1] + sendcount[nDomain-1] ;

    for (n=0; n<nPart; n++) {
        if (sendbuff[n].ID ==0 || part[n].ID == 0)
            // if (sendbuff[n].ID ==0 )
        {
            printf("error - %d %f\n", n, sendbuff[n].pos[0]);
            exit(0);
        }
    }

//    for (i=0; i<nDomain; i++ ){
//        printf("[%d]recvcount.%d = %d\n", myid, i, recvcount[i]);
//        printf("[%d]recvdispl.%d = %d\n", myid, i, recvdispl[i]);
//    }


    printf(" [%d] #part = %d --> %d\n", myid, totsend, totrecv);
    /* memory fragment */
    free(dp->Part);
    dp->Part = (Body*)malloc(sizeof(Body)*totrecv);
    dp->NumPart = totrecv;
    part = dp->Part;

    MPI_Barrier( MPI_COMM_WORLD );
    /* alltoallv */

    MPI_Datatype mpi_body;
    MPI_Type_contiguous(sizeof(Body), MPI_CHAR, &mpi_body );
    MPI_Type_commit(&mpi_body);

    //  MPI_Barrier( MPI_COMM_WORLD );

    int mpi_err =
        MPI_Alltoallv( sendbuff, sendcount, senddispl, mpi_body,
                       part, recvcount, recvdispl, mpi_body, MPI_COMM_WORLD );

    if (mpi_err != 0 ) {
        printf("error = %d \n", mpi_err);
        exit(0);
    }

    MPI_Type_free(&mpi_body);


    /* change the dp->nPart */

    code maxkey, minkey;
    maxkey = dp->DomList[myid].UpperHilbertKey;
    minkey = dp->DomList[myid].LowerHilbertKey;

    for (n=0; n<totrecv; n++) {

        if (part[n].key < minkey ||  part[n].key > maxkey )
        {
            printf("error %d ,key = %ld\n", n, part[n].key);
            exit(0);
        }
    }

    free(sendbuff);
    free(sendcount);
    free(recvcount);
    free(senddispl);
    free(recvdispl);
}
/********************** partition *********************/

void decompose_domain(Domain* dp, GlobalParam *gp)
{
    int myrank, numprocs, NumPart;
    int totSample;
    int n, nSplit;
    Sample *sample;
    int nDomain = dp->NumDom;
    NumPart = dp->NumPart;
    double boxSize = gp->BoxSize;
    int nBits = gp->NumBits;
    int nThread = gp->NumThreadPerSocket;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    sample = gather_total_sample(dp, &totSample);

#ifndef SILENCE
    if ( 0 == myrank )
        printf(" total_sample = %d nBits = %d\n", totSample, nBits);
#endif

    determine_partition_code(sample, totSample, nThread, boxSize, nBits, dp);
    free(sample);

//   coding_domain_particle(dp->Part, NumPart, nThread, boxSize, nBits);

    coding_tagging_domain_particle(dp->Part, NumPart, nThread, boxSize, nBits, dp->NumDom,dp);

//   partition_sorted_particles(dp, NumPart, nThread, boxSize, nBits);
    partition_domains(dp, nThread);
}

/**************** adjust partition ****************/

/* somtimes you can adjust domain boundary slightly, so we have another performance  */

/**************** adjust partition ****************/

void deconstruct_domain() {

}




