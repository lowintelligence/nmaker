#include "global.h"
#include "partition.h"
#include <assert.h>

int increase(const void* p1,const void* p2) {
    if ( *((Peano*)p1) > *((Peano*)p2) )
        return 1;
    else
        return 0;
}

static int snap_index = 0;
static int snap_count = 0;


void partition_system(Constants *constants, Status *status, System *sys)
{
    int r, n, s;
    int nsample, npart, *nsample_rank, nsample_tot;
    int rank, Nproc;
    const int SAMPLE_STAMP = 100000;
    const int RANK_ROOT = 0;
    double accept;
    Int *npart_rank;
    Peano *sample, *sample_total;
    Peano *lower, *upper, PEANO_MAX;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &Nproc);

	/* I don't know what usage will be with PEANO_BITS, 
	 * just disable it, forct it equal to MESH_BITS. */
	assert(constants->PEANO_BITS == constants->MESH_BITS);

    npart_rank = (Int*)xmalloc(Nproc*sizeof(Int), 3001);
    for (r=0; r<Nproc; r++)
        npart_rank[r] = status->num_part_rank[r];

    npart = sys->num_part;
    Body *part = sys->part;

    srand48(rank);
//    accept = (double)(status->frac_sample)*(status->weight_rank[rank])/(double)npart;
    accept = 0.1;
    Real boxsize = constants->BOX_SIZE;
    int d, x, y, z;
    const int NBITS = constants->PEANO_BITS;
    double real2int = (double)(1<<NBITS)/boxsize;
    for (n=0, s=0; n<npart; n++) {
        for (d=0; d<3; d++ ) {
            if ( part[n].pos[d] >= boxsize)
                part[n].pos[d] -= boxsize;
            if ( part[n].pos[d] <0.0 )
                part[n].pos[d] += boxsize;
        }

        x = (int)(part[n].pos[0]*real2int);
        y = (int)(part[n].pos[1]*real2int);
        z = (int)(part[n].pos[2]*real2int);
        part[n].key = encoding_peano(x, y, z, NBITS);

        assert(part[n].pos[0]>=0.0);
        assert(part[n].pos[0]<boxsize);
        assert(part[n].pos[1]>=0.0);
        assert(part[n].pos[1]<boxsize);
        assert(part[n].pos[2]>=0.0);
        assert(part[n].pos[2]<boxsize);

        if ( drand48() < accept*part[n].nstep ) {
            part[n].tag = SAMPLE_STAMP;
            s++;
        }
    }
//    for (n=0; n<10; n++)
    //      printf(" [%d]  %d  %f  %d %lu\n", rank, n, part[n].pos[0], NBITS, part[n].key);


    nsample = s;
    sample = (Peano*)xmalloc(nsample*sizeof(Peano), 3002);
    for (n=0, s=0; n<npart; n++)
        if ( SAMPLE_STAMP == part[n].tag ) {
            sample[s++] = part[n].key;
            part[n].tag = 0;
        }


//    printf("[%d] nsample = %d\n", rank, nsample);

//    qsort(sample, nsample, sizeof(Peano), increase);
//if (0==rank)


    nsample_rank = (int*)xmalloc(Nproc*sizeof(int), 3003);
    MPI_Allgather(&nsample, 1, MPI_INT, nsample_rank, 1, MPI_INT, MPI_COMM_WORLD);
    int *disp_rank = (int*)xmalloc(Nproc*sizeof(int), 3004);
    int *recvpeano = (int*)xmalloc(Nproc*sizeof(int), 3005);
    for (r=0, nsample_tot = 0; r<Nproc; r++) {
        nsample_tot += nsample_rank[r];
        disp_rank[r] = (nsample_tot - nsample_rank[r])*sizeof(Peano);
        recvpeano[r] = nsample_rank[r]*sizeof(Peano);
    }
    sample_total = (Peano*)xmalloc(nsample_tot*sizeof(Peano), 3006);
//printf( "[%d] sample_total = %d\n", rank, nsample_tot);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgatherv(sample, nsample*sizeof(Peano), MPI_CHAR, sample_total,
                   recvpeano, disp_rank, MPI_CHAR, MPI_COMM_WORLD);

    free(recvpeano);
    free(disp_rank);

    qsort(sample_total, nsample_tot, sizeof(Peano), increase);
// merger sort //


    lower = (Peano*)xmalloc(Nproc*sizeof(Peano), 3007);
    upper = (Peano*)xmalloc(Nproc*sizeof(Peano), 3008);
    PEANO_MAX = ((Peano)1<<(3*NBITS));

    int cut;
    double d_cut = 1.0/Nproc;
    lower[0] = (Peano)0;
    for (r=1; r<Nproc; r++) {
        cut = (int)( r*d_cut*nsample_tot );
//    printf("[%d], cut = %d\n", rank, cut);
        upper[r-1] = lower[r] = sample_total[cut];
//       if (RANK_ROOT == rank)
        //         printf(" %d %lu\n", cut, sample_total[cut] );
    }
    upper[Nproc-1] = PEANO_MAX;
//printf("[%d] lower = %lu, upper = %lu\n", rank, lower[rank], upper[rank]);
    free(sample_total);
    free(nsample_rank);

    int *sendcount, *recvcount;
    Peano p_lower, p_upper;
    p_lower = lower[rank];
    p_upper = upper[rank];

    sendcount = (int*)xmalloc(Nproc*sizeof(int),3009);
    recvcount = (int*)xmalloc(Nproc*sizeof(int),3010);
    for (r=0; r<Nproc; r++)
        sendcount[r] = recvcount[r] = 0;

    for (n=0; n<npart; n++)
	{
        if ( p_lower <= part[n].key && part[n].key < p_upper )
		{
            sendcount[rank]++;
            part[n].tag = (Int)rank;
        }
        else
		{
            for (r=0; r<Nproc; r++)
			{
                if (lower[r] > part[n].key) break;
			}
            r--;
            sendcount[r]++;
            part[n].tag = (Int)(r);
        }
    }

//   for (r=0; r<Nproc; r++)
    //      printf(" sendcount[%d->%d] = %d \n", rank, r, sendcount[r]);

    MPI_Alltoall(sendcount, 1, MPI_INT, recvcount, 1, MPI_INT, MPI_COMM_WORLD);

//    for (r=0; r<Nproc; r++)
    //      printf(" recvcount[%d->%d] = %d \n", r, rank, recvcount[r]);

    int  *senddisp, *recvdisp;
    Body *sendbuff, *recvbuff;
    int  dest, source, flag;
    MPI_Status mpi_status;

    senddisp = (int*)xmalloc(Nproc*sizeof(int), 3011);
    recvdisp = (int*)xmalloc(Nproc*sizeof(int), 3012);
    senddisp[0] = recvdisp[0] = 0;
    for (r=1; r<Nproc; r++) {
        senddisp[r] = sendcount[r-1] + senddisp[r-1];
        recvdisp[r] = recvcount[r-1] + recvdisp[r-1];
    }
//    for (r=0; r<Nproc; r++)
    //      printf("[%d] before {%d} senddisp[r] %d\n", rank, r, senddisp[r]);

    sendbuff = (Body*)xmalloc(npart*sizeof(Body), 3013);
    for (n=0; n<npart; n++) {
        r = (int)(part[n].tag);
        if(r>=Nproc) {
            printf(" -   r = %d [%d] \n", r, rank);
            MPI_Abort(MPI_COMM_WORLD, 0);
        }
        sendbuff[senddisp[r]] = part[n];
        senddisp[r]++;
    }
//    for (r=0; r<Nproc; r++)
//      printf("[%d] after  {%d} senddisp[r] %d\n",rank,  r, senddisp[r]);

    senddisp[0] = recvdisp[0] =  0;
    int total_recv = recvcount[0];
    for (r=1; r<Nproc; r++) {
        senddisp[r] = sendcount[r-1] + senddisp[r-1];
        recvdisp[r] = recvcount[r-1] + recvdisp[r-1];
        total_recv += recvcount[r];
    }
    free(part);
    sys->part = (Body*)xmalloc(total_recv*sizeof(Body),3014);
    recvbuff = sys->part;

    flag = 10;
    /*
        for (r=0; r<Nproc; r++) {
            dest  = (rank+r)%Nproc;
            source = (rank-r+Nproc)%Nproc;
            MPI_Sendrecv(sendbuff+senddisp[dest], sendcount[dest]*sizeof(Body),
                MPI_CHAR, dest, flag, recvbuff+recvdisp[source],
                recvcount[source]*sizeof(Body), MPI_CHAR, source, flag,
                MPI_COMM_WORLD, &mpi_status);
        }
    */
//////////////////////////
    for (r=0; r<Nproc; r++) {
        sendcount[r] *= sizeof(Body);
        recvcount[r] *= sizeof(Body);
        senddisp[r]  *= sizeof(Body);
        recvdisp[r]  *= sizeof(Body);
    }
    MPI_Alltoallv(sendbuff, sendcount, senddisp,
                  MPI_CHAR,  recvbuff,
                  recvcount, recvdisp, MPI_CHAR,
                  MPI_COMM_WORLD);
/////////////////////////////////
    sys->num_part = (Int)total_recv;

    MPI_Allgather(&(sys->num_part), sizeof(Int), MPI_CHAR, status->num_part_rank, sizeof(Int), MPI_CHAR, MPI_COMM_WORLD);

    double mean_part = (double)(constants->TOTAL_NUM_PART)/Nproc;
    for (r=0; r<Nproc; r++) {
        status->lower_rank[r]  = lower[r];
        status->upper_rank[r]  = upper[r];
        status->weight_rank[r] = 1.0;
//        status->weight_rank[r] = 1.0+0.1*rank;
    }
    free(sendbuff);
    MPI_Barrier(MPI_COMM_WORLD);
//   sleep(1.0*rank);
    LOG_INFO(LOG_VB1, "[%d] npart = %lu  lower = %lu, upper = %lu\n", rank, sys->num_part, lower[rank], upper[rank]);

    /*
        FILE *fd;
        char temp[50];
        sprintf(temp,"snap_%d.txt", rank);
        fd = fopen(temp,"w");

        for (n=0; n<(sys->num_part); n++) {
            fprintf(fd, "%f %f %f\n", sys->part[n].pos[0], sys->part[n].pos[1], sys->part[n].pos[2] );
        }
        fclose(fd);
    */
    for (n=0; n<sys->num_part; n++) {
        if (sys->part[n].tag != rank) {
            LOG_INFO(LOG_FATAL, " partition error !\n");
            system_exit(0);
        }
    }

    /*
    if ( 0== rank && 0 == snap_count % 10 ) {

        FILE *fd;
        char fname[30];
        sprintf(fname, "sample%d.txt", snap_index);
        fd = fopen(fname, "w");
        for (n=0, s=0; n<npart; n++)
            if ( 0 == ( (int)(part[n].id)% 30000 ) )
                fprintf(fd, "%f %f %f\n", part[n].pos[0], part[n].pos[1], part[n].pos[2]);


        fclose(fd);

        sprintf(fname, "slice%d.txt", snap_index);
        fd = fopen(fname, "w");
        for (n=0, s=0; n<npart; n++)
            if (  part[n].pos[0] < 5000.0 && part[n].pos[1] < 100000.0 && part[n].pos[2] < 100000.0)
                fprintf(fd, "%f %f\n", part[n].pos[1], part[n].pos[2]);


        fclose(fd);


        snap_index ++;
    }
    snap_count ++;
    */


    free(recvdisp);
    free(senddisp);

    free(recvcount);
    free(sendcount);

    free(lower);
    free(upper);

    free(sample);
    free(npart_rank);
}



