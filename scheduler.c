#include "scheduler.h"
#include <assert.h>

void broadcast_frontiers(Domain *dp, GlobalParam *gp)
{
    int i, j, k, n, nb, cntpart, start, idx, dest, In, Jn, Kn, index_length;
    int x, y, z;
    int s, nSock = gp->NumSocket;
    int *tag  = dp->cuboid->tag;
    int *count = dp->cuboid->count;
    int npart;
    int *numlock, *numpart;
    int *numsend, *numrecv;
    BuffData **buffsend, **buffrecv;
    Body *part = dp->Part;
    Body *adj_part;

    numsend = (int*) malloc(sizeof(int)*nSock);
    numrecv = (int*) malloc(sizeof(int)*nSock);
    numlock = (int*) malloc(sizeof(int)*nSock);
    numpart = (int*) malloc(sizeof(int)*nSock);

    buffsend = (BuffData**) malloc(sizeof(BuffData*)*nSock);
    buffrecv = (BuffData**) malloc(sizeof(BuffData*)*nSock);

    for (s=0; s<nSock; s++) {
        numpart[s] = 0;
        numsend[s] = 0;
        numrecv[s] = 0;
        numlock[s] = 1;
        buffsend[s] = NULL;
        buffrecv[s] = NULL;
    }

    In = dp->cuboid->nSide[0];
    Jn = dp->cuboid->nSide[1];
    Kn = dp->cuboid->nSide[2];

    int count_frontiers = 0;
    for (i=0; i<In; i++) {
        for (j=0; j<Jn; j++) {
            for (k=0; k<Kn; k++) {
                idx = (i*Jn+j)*Kn+k;
                if ( tag[idx] == TAG_FRONTIERS )
                    count_frontiers += count[idx];
            }
        }
    }
    printf(" [%d] frontiers = %d\n", dp->rank, count_frontiers);
    /*
        for (i=0; i<In; i++)
            for (j=0; j<Jn; j++)
                for (k=0; k<Kn; k++) {
                    idx = (i*Jn + j )*Kn + k;
                    if (TAG_FRONTIERS == tag[idx]) {
                        for (x=i-1; x<=i+1; x++)
                            for (y=j-1; y<=j+1; y++)
                                for (z=k-1; z<=k+1; z++)  {
                                    dest = tag[(x*Jn+y)*Kn+z];
                                    if ( dest  >= 0 && numlock[dest] ) {
                                        numsend[dest] += count[idx];
                                        numlock[dest] = 0;
                                    }
                                }
                        for (x=i-1; x<=i+1; x++)
                            for (y=j-1; y<=j+1; y++)
                                for (z=k-1; z<=k+1; z++)
                                    if ( (dest=tag[(x*Jn+y)*Kn+z])  >= 0 )
                                        numlock[dest] = 1;

                    }
                }

        */


    int cnt, active[27];
    for (i=0; i<In; i++)
        for (j=0; j<Jn; j++)
            for (k=0; k<Kn; k++) {
                idx = (i*Jn + j )*Kn + k;
                if (TAG_FRONTIERS == tag[idx]) {
                    cnt = 0;
                    for (x=i-1; x<=i+1; x++)
                        for (y=j-1; y<=j+1; y++)
                            for (z=k-1; z<=k+1; z++)  {
                                dest = tag[(x*Jn+y)*Kn+z];
                                if ( dest  >= 0 && numlock[dest] ) {
                                    numsend[dest] += count[idx];
                                    numlock[dest] = 0;

                                    active[cnt++] = dest;

                                }
                            }

                    while ( 0< cnt-- )
                        numlock[ active[cnt] ] = 1;

                    assert(cnt==-1);


                }
            }

    for (s=0; s<nSock; s++) {
        printf(" send %10d part(s) from[%d]to[%d]\n", numsend[s], dp->rank, s);
    }

    MPI_Alltoall(numsend, 1, MPI_INT, numrecv, 1, MPI_INT, TREE_COMM_WORLD);

    for (s=0; s<nSock; s++) {
//        printf(" recv %10d part(s) from[%d]to[%d]\n", numrecv[s], s, dp->rank);
    }

    int recv_tot_part, send_tot_part;

    send_tot_part = recv_tot_part = 0;
    for (s=0; s<nSock; s++) {
        recv_tot_part += numrecv[s];
        send_tot_part += numsend[s];
    }
    printf(" [%2d] tot_send = %d  tot_recv = %d domain_part = %d\n", dp->rank, send_tot_part, recv_tot_part, dp->NumPart);

    for (s=0; s<nSock; s++) {
        if ( 0<numsend[s] )
            buffsend[s] = (BuffData*)malloc(sizeof(BuffData) * numsend[s] );
        if ( 0<numrecv[s] )
            buffrecv[s] = (BuffData*)malloc(sizeof(BuffData) * numrecv[s] );
    }

    int nbits = gp->NumBits;
    double boxsize = gp->BoxSize;
    double nbox_boxsize = (1<<nbits)/boxsize;
    double minidx[3], npad;

    npad = 1;
    minidx[0] = dp->cuboid->minimum[0];
    minidx[1] = dp->cuboid->minimum[1];
    minidx[2] = dp->cuboid->minimum[2];
    npart = dp->NumPart;
//printf("%f %f %f\n", part[npart-1].pos[0], part[npart-1].pos[1], part[npart-1].pos[2]);
    MPI_Barrier(TREE_COMM_WORLD);

    for (n=0; n<npart; n++) {
        if (TAG_FRONTIERS == part[n].tag) {
            i = (int)(part[n].pos[0]*nbox_boxsize) - minidx[0] + npad;
            j = (int)(part[n].pos[1]*nbox_boxsize) - minidx[1] + npad;
            k = (int)(part[n].pos[2]*nbox_boxsize) - minidx[2] + npad;
            //printf("%d %d %d\n", x, y, z);
            cnt =0 ;
            for (x=i-1; x<=i+1; x++)
                for (y=j-1; y<=j+1; y++)
                    for (z=k-1; z<=k+1; z++) {
                        dest = tag[(x*Jn+y)*Kn+z];
                        if ( dest >= 0 && numlock[dest] ) {
                            //numsend[dest]+= count[idx];
                            nb = numpart[dest];

                            buffsend[dest][nb].pos[0] = part[n].pos[0];
                            buffsend[dest][nb].pos[1] = part[n].pos[1];
                            buffsend[dest][nb].pos[2] = part[n].pos[2];

                            numpart[dest]++;
                            numlock[dest] = 0;
                            active[cnt++] = dest;

                        }
                    }

            while ( 0< cnt-- )
                numlock[ active[cnt] ] = 1;
        }
    }

    for (s=0; s<nSock; s++) {
        assert(numlock[s]);
        assert(numpart[s] == numsend[s]);
    }

    int data_btye = sizeof(BuffData);
    int rank, recv, send;
    MPI_Status *status;
    MPI_Request *req_send, *req_recv;

    req_send = (MPI_Request*)malloc(sizeof(MPI_Request)*nSock);
    req_recv = (MPI_Request*)malloc(sizeof(MPI_Request)*nSock);
    status = (MPI_Status*)malloc(sizeof(MPI_Status)*nSock);

    MPI_Datatype mpi_buffdata_type;
    MPI_Type_contiguous(sizeof(BuffData), MPI_CHAR, &mpi_buffdata_type );
    MPI_Type_commit(&mpi_buffdata_type);

    rank = dp->rank;

    MPI_Barrier(TREE_COMM_WORLD);

    for (s=0; s<nSock; s++) {
        send = (rank+s)%nSock;
        recv = (rank-s+nSock)%nSock;
        printf("%d:  %d ->[%d] :%d\n",s,recv, rank, numrecv[recv]);

        MPI_Sendrecv(buffsend[send], numsend[send], mpi_buffdata_type, send, 0, buffrecv[recv], numrecv[recv], mpi_buffdata_type, recv, 0, TREE_COMM_WORLD, &status[s]);
    }

    printf(" tot_adjoining_part = %d \n",recv_tot_part );
    cnt = 0;
    adj_part = (Body*)malloc(sizeof(Body)*recv_tot_part);
    for (s=0; s<nSock; s++) {
        for (n=0; n<numrecv[s]; n++) {
            adj_part[cnt].pos[0] = buffrecv[s][n].pos[0];
            adj_part[cnt].pos[1] = buffrecv[s][n].pos[1];
            adj_part[cnt].pos[2] = buffrecv[s][n].pos[2];
            cnt++;
        }
    }

    /*
      char fname[200] ;
      sprintf(fname, "test_adj_%d.dat", rank);
           FILE *fd = fopen(fname,"w");
           for (n=0; n<recv_tot_part; n++){
               fprintf(fd, "%f %f %f\n", adj_part[n].pos[0], adj_part[n].pos[1], adj_part[n].pos[2]);
           }
           fclose(fd);
      */



    MPI_Barrier(TREE_COMM_WORLD);

    free(req_send);
    free(req_recv);
    free(status);
    MPI_Type_free(&mpi_buffdata_type);

    free(numpart);
    free(numlock);
    free(numsend);
    free(numrecv);
    for (s=0; s<nSock; s++) {
        if ( NULL != buffsend[s] )
            free(buffsend[s]);
        if ( NULL != buffrecv[s] )
            free(buffrecv[s]);
    }
    free(buffsend);
    free(buffrecv);
}

