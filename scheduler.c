#include "scheduler.h"
#include <assert.h>

#define TMPNUM 2097152

int inbox(int x, int y, int z, int *minidx, int *maxidx)
{
	int in=0;
	if (x>=minidx[0]-1 && x<=maxidx[0]+1)
	{
		if (y>=minidx[1]-1 && y<=maxidx[1]+1)
		{
			if (z>=minidx[2]-1 && z<=maxidx[2]+1)
			{
				in = 1;
			}
		}
	}
//	if (in==1) printf("1");
	return in;
}

int getidx(int x, int y, int z, int Jn, int Kn, int *minidx, int npad)
{
	int idx;
	idx = ((x-minidx[0]+npad)*Jn + (y-minidx[1]+npad))*Kn + (z-minidx[2]+npad);
	return idx;
}

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
	Body **buffsend, **buffrecv;
    Body *part = dp->Part;
    Body *adj_part;

    numsend = (int*) malloc(sizeof(int)*nSock);
    numrecv = (int*) malloc(sizeof(int)*nSock);
    numlock = (int*) malloc(sizeof(int)*nSock);
    numpart = (int*) malloc(sizeof(int)*nSock);

    buffsend = (Body**) malloc(sizeof(Body*)*nSock);
    buffrecv = (Body**) malloc(sizeof(Body*)*nSock);

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
            buffsend[s] = (Body*)malloc(sizeof(Body) * numsend[s] );
        if ( 0<numrecv[s] )
            buffrecv[s] = (Body*)malloc(sizeof(Body) * numrecv[s] );
    }

    int nbits = gp->NumBits;
	int Gside = 1<<nbits;
	int Ngrid = In*Jn*Kn;
    double boxsize = gp->BoxSize;
    double nbox_boxsize = (1<<nbits)/boxsize;
    int minidx[3], maxidx[3], npad;

    npart = dp->NumPart;
	npad = 1;
    minidx[0] = dp->cuboid->minimum[0];
    minidx[1] = dp->cuboid->minimum[1];
    minidx[2] = dp->cuboid->minimum[2];
    maxidx[0] = dp->cuboid->maximum[0];
    maxidx[1] = dp->cuboid->maximum[1];
    maxidx[2] = dp->cuboid->maximum[2];
	
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

                            buffsend[dest][nb] = part[n];

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

    int data_btye = sizeof(Body);
    int rank, recv, send;
    MPI_Status *status;
    MPI_Request *req_send, *req_recv;

    req_send = (MPI_Request*)malloc(sizeof(MPI_Request)*nSock);
    req_recv = (MPI_Request*)malloc(sizeof(MPI_Request)*nSock);
    status = (MPI_Status*)malloc(sizeof(MPI_Status)*nSock);

    MPI_Datatype mpi_buffdata_type;
    MPI_Type_contiguous(sizeof(Body), MPI_CHAR, &mpi_buffdata_type );
    MPI_Type_commit(&mpi_buffdata_type);

    rank = dp->rank;

    MPI_Barrier(TREE_COMM_WORLD);

    for (s=0; s<nSock; s++) {
        send = (rank+s)%nSock;
        recv = (rank-s+nSock)%nSock;
        printf("%d:  %d ->[%d] :%d\n",s,recv, rank, numrecv[recv]);

        MPI_Sendrecv(buffsend[send], numsend[send], mpi_buffdata_type, send, 0, buffrecv[recv], numrecv[recv], mpi_buffdata_type, recv, 0, TREE_COMM_WORLD, &status[s]);
    }

    printf("[%d] tot_adjoining_part = %d \n",dp->rank, recv_tot_part );
    cnt = 0;
	part = (Body*)realloc((void*)dp->Part,sizeof(Body)*(npart+recv_tot_part+TMPNUM));
	dp->Part = part;
//	dp->NumPart += recv_tot_part;
    adj_part = part + npart;

    for (s=0; s<nSock; s++) {
        for (n=0; n<numrecv[s]; n++) {
            adj_part[cnt] = buffrecv[s][n];

			x = (int)(adj_part[cnt].pos[0]*nbox_boxsize);
			y = (int)(adj_part[cnt].pos[1]*nbox_boxsize);
			z = (int)(adj_part[cnt].pos[2]*nbox_boxsize);
			int idx=getidx(x, y, z, Jn, Kn, minidx, npad);
//			printf("[%d]  box=%f, idx=%d, x=%d, y=%d, z=%d, px=%f, py=%f, pz=%f.\n", dp->rank, 1/nbox_boxsize, idx, x, y, z, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
			int tmpidx;
			if (x==Gside-1)
			{
				if (minidx[0]==0)
				{
					tmpidx = getidx(x-Gside, y, z, Jn, Kn, minidx, npad);
					if (inbox(x-Gside,y,z,minidx,maxidx) && dp->cuboid->tag[tmpidx]>=0)
					{
						adj_part[cnt].pos[0] -= boxsize;
						if(tmpidx>=0 && tmpidx<Ngrid) 
						{
							adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
							adj_part[cnt].group = tmpidx;
							dp->cuboid->count[tmpidx]++;
							if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
								printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
						}
						else
						{
							printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
						}
						cnt++;
						adj_part[cnt] = buffrecv[s][n];
					}
					if (y==Gside-1)
					{
						if (minidx[1]==0)
						{
							tmpidx = getidx(x-Gside, y-Gside, z, Jn, Kn, minidx, npad);
							if ((inbox(x-Gside,y-Gside,z,minidx,maxidx)) && dp->cuboid->tag[tmpidx]>=0)
							{
								adj_part[cnt].pos[0] -= boxsize;
								adj_part[cnt].pos[1] -= boxsize;
								if(tmpidx>=0 && tmpidx<Ngrid) 
								{
									adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
									adj_part[cnt].group = tmpidx;
									dp->cuboid->count[tmpidx]++;
									if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
										printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
								}
								else
								{
									printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
								}
								cnt++;
								adj_part[cnt] = buffrecv[s][n];
							}
							if (z==Gside-1)
							{
								tmpidx = getidx(x-Gside, y-Gside, z-Gside, Jn, Kn, minidx, npad);
								if (minidx[2]==0 && (inbox(x-Gside,y-Gside,z-Gside,minidx,maxidx)) && dp->cuboid->tag[tmpidx]>=0)
								{
									adj_part[cnt].pos[0] -= boxsize;
									adj_part[cnt].pos[1] -= boxsize;
									adj_part[cnt].pos[2] -= boxsize;
									tmpidx = idx-Gside*Jn*Kn-Gside*Kn-Gside;
									if(tmpidx>=0 && tmpidx<Ngrid) 
									{
										adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
										adj_part[cnt].group = tmpidx;
										dp->cuboid->count[tmpidx]++;
										if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
											printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
									}
									else
									{
										printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
									}
									cnt++;
									adj_part[cnt] = buffrecv[s][n];
								}
							}
							if(z==0)
							{
								tmpidx = getidx(x-Gside, y-Gside, z+Gside, Jn, Kn, minidx, npad);
								if (maxidx[2]==Gside-1 && (inbox(x-Gside,y-Gside,z+Gside,minidx,maxidx)) && dp->cuboid->tag[tmpidx]>=0)
								{
									adj_part[cnt].pos[0] -= boxsize;
									adj_part[cnt].pos[1] -= boxsize;
									adj_part[cnt].pos[2] += boxsize;
									if(tmpidx>=0 && tmpidx<Ngrid) 
									{
										adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
										adj_part[cnt].group = tmpidx;
										dp->cuboid->count[tmpidx]++;
										if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
											printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
									}
									else
									{
										printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
									}
									cnt++;
									adj_part[cnt] = buffrecv[s][n];
								}
							}
						}
					}
					if (y==0)
					{
						if (maxidx[1]==Gside-1)
						{
							tmpidx = getidx(x-Gside, y+Gside, z, Jn, Kn, minidx, npad);
							if ((inbox(x-Gside,y+Gside,z,minidx,maxidx)) && dp->cuboid->tag[tmpidx]>=0)
							{
								adj_part[cnt].pos[0] -= boxsize;
								adj_part[cnt].pos[1] += boxsize;
								if(tmpidx>=0 && tmpidx<Ngrid) 
								{
									adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
									adj_part[cnt].group = tmpidx;
									dp->cuboid->count[tmpidx]++;
									if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
										printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
								}
								else
								{
									printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
								}
								cnt++;
								adj_part[cnt] = buffrecv[s][n];
							}
							if (z==Gside-1)
							{
								tmpidx = getidx(x-Gside, y+Gside, z-Gside, Jn, Kn, minidx, npad);
								if (minidx[2]==0 && (inbox(x-Gside,y+Gside,z-Gside,minidx,maxidx)) && dp->cuboid->tag[tmpidx]>=0)
								{
									adj_part[cnt].pos[0] -= boxsize;
									adj_part[cnt].pos[1] += boxsize;
									adj_part[cnt].pos[2] -= boxsize;
									if(tmpidx>=0 && tmpidx<Ngrid) 
									{
										adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
										adj_part[cnt].group = tmpidx;
										dp->cuboid->count[tmpidx]++;
										if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
											printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
									}
									else
									{
										printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
									}
									cnt++;
									adj_part[cnt] = buffrecv[s][n];
								}
							}
							if(z==0)
							{
								tmpidx = getidx(x-Gside, y+Gside, z+Gside, Jn, Kn, minidx, npad);
								if (maxidx[2]==Gside-1 && (inbox(x-Gside,y+Gside,z+Gside,minidx,maxidx)) && dp->cuboid->tag[tmpidx]>=0)
								{
									adj_part[cnt].pos[0] -= boxsize;
									adj_part[cnt].pos[1] += boxsize;
									adj_part[cnt].pos[2] += boxsize;
									if(tmpidx>=0 && tmpidx<Ngrid) 
									{
										adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
										adj_part[cnt].group = tmpidx;
										dp->cuboid->count[tmpidx]++;
										if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
											printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
									}
									else
									{
										printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
									}
									cnt++;
									adj_part[cnt] = buffrecv[s][n];
								}
							}
						}
					}
					if (z==Gside-1)
					{
						tmpidx = getidx(x-Gside, y, z-Gside, Jn, Kn, minidx, npad);
						if (minidx[2]==0 && (inbox(x-Gside,y,z-Gside,minidx,maxidx)) && dp->cuboid->tag[tmpidx]>=0)
						{
							adj_part[cnt].pos[0] -= boxsize;
							adj_part[cnt].pos[2] -= boxsize;
							if(tmpidx>=0 && tmpidx<Ngrid) 
							{
								adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
								adj_part[cnt].group = tmpidx;
								dp->cuboid->count[tmpidx]++;
								if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
									printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
							}
							else
							{
								printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
							}
							cnt++;
							adj_part[cnt] = buffrecv[s][n];
						}
					}
					if(z==0)
					{
						tmpidx = getidx(x-Gside, y, z+Gside, Jn, Kn, minidx, npad);
						if (maxidx[2]==Gside-1 && (inbox(x-Gside,y,z+Gside,minidx,maxidx)) && dp->cuboid->tag[tmpidx]>=0)
						{
							adj_part[cnt].pos[0] -= boxsize;
							adj_part[cnt].pos[2] += boxsize;
							tmpidx = idx-Gside*Jn*Kn+Gside;
							if(tmpidx>=0 && tmpidx<Ngrid) 
							{
								adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
								adj_part[cnt].group = tmpidx;
								dp->cuboid->count[tmpidx]++;
								if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
									printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
							}
							else
							{
								printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
							}
							cnt++;
							adj_part[cnt] = buffrecv[s][n];
						}
					}
				}
			}
			if (x==0)
			{
				if (maxidx[0]==Gside-1)
				{
					tmpidx = getidx(x+Gside, y, z, Jn, Kn, minidx, npad);
					if ((inbox(x+Gside,y,z,minidx,maxidx)) && dp->cuboid->tag[tmpidx]>=0)
					{
						adj_part[cnt].pos[0] += boxsize;
						if(tmpidx>=0 && tmpidx<Ngrid) 
						{
							adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
							adj_part[cnt].group = tmpidx;
							dp->cuboid->count[tmpidx]++;
							if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
								printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
						}
						else
						{
							printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
						}
						cnt++;
						adj_part[cnt] = buffrecv[s][n];
					}
					if (y==Gside-1)
					{
						if (minidx[1]==0)
						{
							tmpidx = getidx(x+Gside, y-Gside, z, Jn, Kn, minidx, npad);
							if ((inbox(x+Gside,y-Gside,z,minidx,maxidx)) && dp->cuboid->tag[tmpidx]>=0)
							{
								adj_part[cnt].pos[0] += boxsize;
								adj_part[cnt].pos[1] -= boxsize;
								if(tmpidx>=0 && tmpidx<Ngrid) 
								{
									adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
									adj_part[cnt].group = tmpidx;
									dp->cuboid->count[tmpidx]++;
									if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
										printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
								}
								else
								{
									printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
								}
								cnt++;
								adj_part[cnt] = buffrecv[s][n];
							}
							if (z==Gside-1)
							{
								tmpidx = getidx(x+Gside, y-Gside, z-Gside, Jn, Kn, minidx, npad);
								if (minidx[2]==0 && (inbox(x+Gside,y-Gside,z-Gside,minidx,maxidx)) && dp->cuboid->tag[tmpidx]>=0)
								{
									adj_part[cnt].pos[0] += boxsize;
									adj_part[cnt].pos[1] -= boxsize;
									adj_part[cnt].pos[2] -= boxsize;
									if(tmpidx>=0 && tmpidx<Ngrid) 
									{
										adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
										adj_part[cnt].group = tmpidx;
										dp->cuboid->count[tmpidx]++;
										if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
											printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
									}
									else
									{
										printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
									}
									cnt++;
									adj_part[cnt] = buffrecv[s][n];
								}
							}
							if(z==0)
							{
								tmpidx = getidx(x+Gside, y-Gside, z+Gside, Jn, Kn, minidx, npad);
								if (maxidx[2]==Gside-1 && (inbox(x+Gside,y-Gside,z+Gside,minidx,maxidx)) && dp->cuboid->tag[tmpidx]>=0)
								{
									adj_part[cnt].pos[0] += boxsize;
									adj_part[cnt].pos[1] -= boxsize;
									adj_part[cnt].pos[2] += boxsize;
									if(tmpidx>=0 && tmpidx<Ngrid) 
									{
										adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
										adj_part[cnt].group = tmpidx;
										dp->cuboid->count[tmpidx]++;
										if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
											printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
									}
									else
									{
										printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
									}
									cnt++;
									adj_part[cnt] = buffrecv[s][n];
								}
							}
						}
					}
					if (y==0)
					{
						if (maxidx[1]==Gside-1)
						{
							tmpidx = getidx(x+Gside, y+Gside, z, Jn, Kn, minidx, npad);
							if ((inbox(x+Gside,y+Gside,z,minidx,maxidx)) && dp->cuboid->tag[tmpidx]>=0)
							{
								adj_part[cnt].pos[0] += boxsize;
								adj_part[cnt].pos[1] += boxsize;
								if(tmpidx>=0 && tmpidx<Ngrid) 
								{
									adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
									adj_part[cnt].group = tmpidx;
									dp->cuboid->count[tmpidx]++;
									if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
										printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
								}
								else
								{
									printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
								}
								cnt++;
								adj_part[cnt] = buffrecv[s][n];
							}
							if (z==Gside-1)
							{
								tmpidx = getidx(x+Gside, y+Gside, z-Gside, Jn, Kn, minidx, npad);
								if (minidx[2]==0 && (inbox(x+Gside,y+Gside,z-Gside,minidx,maxidx)) && dp->cuboid->tag[tmpidx]>=0)
								{
									adj_part[cnt].pos[0] += boxsize;
									adj_part[cnt].pos[1] += boxsize;
									adj_part[cnt].pos[2] -= boxsize;
									if(tmpidx>=0 && tmpidx<Ngrid) 
									{
										adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
										adj_part[cnt].group = tmpidx;
										dp->cuboid->count[tmpidx]++;
										if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
											printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
									}
									else
									{
										printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
									}
									cnt++;
									adj_part[cnt] = buffrecv[s][n];
								}
							}
							if(z==0)
							{
								tmpidx = getidx(x+Gside, y+Gside, z+Gside, Jn, Kn, minidx, npad);
								if (maxidx[2]==Gside-1 && (inbox(x+Gside,y+Gside,z+Gside,minidx,maxidx)) && dp->cuboid->tag[tmpidx]>=0)
								{
									adj_part[cnt].pos[0] += boxsize;
									adj_part[cnt].pos[1] += boxsize;
									adj_part[cnt].pos[2] += boxsize;
									tmpidx = idx+Gside*Jn*Kn+Gside*Kn+Gside;
									if(tmpidx>=0 && tmpidx<Ngrid) 
									{
										adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
										adj_part[cnt].group = tmpidx;
										dp->cuboid->count[tmpidx]++;
										if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
											printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
									}
									else
									{
										printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
									}
									cnt++;
									adj_part[cnt] = buffrecv[s][n];
								}
							}
						}
					}
					if (z==Gside-1)
					{
						tmpidx = getidx(x+Gside, y, z-Gside, Jn, Kn, minidx, npad);
						if (minidx[2]==0 && (inbox(x+Gside,y,z-Gside,minidx,maxidx)) && dp->cuboid->tag[tmpidx]>=0)
						{
							adj_part[cnt].pos[0] += boxsize;
							adj_part[cnt].pos[2] -= boxsize;
							if(tmpidx>=0 && tmpidx<Ngrid) 
							{
								adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
								adj_part[cnt].group = tmpidx;
								dp->cuboid->count[tmpidx]++;
								if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
									printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
							}
							else
							{
								printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
							}
							cnt++;
							adj_part[cnt] = buffrecv[s][n];
						}
					}
					if(z==0)
					{
						tmpidx = getidx(x+Gside, y, z+Gside, Jn, Kn, minidx, npad);
						if (maxidx[2]==Gside-1 && (inbox(x+Gside,y,z+Gside,minidx,maxidx)) && dp->cuboid->tag[tmpidx]>=0)
						{
							adj_part[cnt].pos[0] += boxsize;
							adj_part[cnt].pos[2] += boxsize;
							if(tmpidx>=0 && tmpidx<Ngrid) 
							{
								adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
								adj_part[cnt].group = tmpidx;
								dp->cuboid->count[tmpidx]++;
								if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
									printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
							}
							else
							{
								printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
							}
							cnt++;
							adj_part[cnt] = buffrecv[s][n];
						}
					}
				}
			}
			if (y==Gside-1)
			{
				if (minidx[1]==0)
				{
					tmpidx = getidx(x, y-Gside, z, Jn, Kn, minidx, npad);
					if ((inbox(x,y-Gside,z,minidx,maxidx)) && dp->cuboid->tag[tmpidx]>=0)
					{
						adj_part[cnt].pos[1] -= boxsize;
						if(tmpidx>=0 && tmpidx<Ngrid) 
						{
							adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
							adj_part[cnt].group = tmpidx;
							dp->cuboid->count[tmpidx]++;
							if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
								printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
						}
						else
						{
							printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
						}
						cnt++;
						adj_part[cnt] = buffrecv[s][n];
					}
					if (z==Gside-1)
					{
						tmpidx = getidx(x, y-Gside, z-Gside, Jn, Kn, minidx, npad);
						if (minidx[2]==0 && (inbox(x,y-Gside,z-Gside,minidx,maxidx)) && dp->cuboid->tag[tmpidx]>=0)
						{
							adj_part[cnt].pos[1] -= boxsize;
							adj_part[cnt].pos[2] -= boxsize;
							if(tmpidx>=0 && tmpidx<Ngrid) 
							{
								adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
								adj_part[cnt].group = tmpidx;
								dp->cuboid->count[tmpidx]++;
								if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
									printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
							}
							else
							{
								printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
							}
							cnt++;
							adj_part[cnt] = buffrecv[s][n];
						}
					}
					if(z==0)
					{
						tmpidx = getidx(x, y-Gside, z+Gside, Jn, Kn, minidx, npad);
						if (maxidx[2]==Gside-1 && (inbox(x,y-Gside,z+Gside,minidx,maxidx)) && dp->cuboid->tag[tmpidx]>=0)
						{
							adj_part[cnt].pos[1] -= boxsize;
							adj_part[cnt].pos[2] += boxsize;
							if(tmpidx>=0 && tmpidx<Ngrid) 
							{
								adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
								adj_part[cnt].group = tmpidx;
								dp->cuboid->count[tmpidx]++;
								if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
									printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
							}
							else
							{
								printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
							}
							cnt++;
							adj_part[cnt] = buffrecv[s][n];
						}
					}
				}
			}
			if (y==0)
			{
				if (maxidx[1]==Gside-1)
				{
					tmpidx = getidx(x, y+Gside, z, Jn, Kn, minidx, npad);
					if ((inbox(x,y+Gside,z,minidx,maxidx)) && dp->cuboid->tag[tmpidx]>=0)
					{
						adj_part[cnt].pos[1] += boxsize;
						if(tmpidx>=0 && tmpidx<Ngrid) 
						{
							adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
							adj_part[cnt].group = tmpidx;
							dp->cuboid->count[tmpidx]++;
							if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
								printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
						}
						else
						{
							printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
						}
						cnt++;
						adj_part[cnt] = buffrecv[s][n];
					}
					if (z==Gside-1)
					{
						tmpidx = getidx(x, y+Gside, z-Gside, Jn, Kn, minidx, npad);
						if (minidx[2]==0 && (inbox(x,y+Gside,z-Gside,minidx,maxidx)) && dp->cuboid->tag[tmpidx]>=0)
						{
							adj_part[cnt].pos[1] += boxsize;
							adj_part[cnt].pos[2] -= boxsize;
							if(tmpidx>=0 && tmpidx<Ngrid) 
							{
								adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
								adj_part[cnt].group = tmpidx;
								dp->cuboid->count[tmpidx]++;
								if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
									printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
							}
							else
							{
								printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
							}
							cnt++;
							adj_part[cnt] = buffrecv[s][n];
						}
					}
					if(z==0)
					{
						tmpidx = getidx(x, y+Gside, z+Gside, Jn, Kn, minidx, npad);
						if (maxidx[2]==Gside-1 && (inbox(x,y+Gside,z+Gside,minidx,maxidx)) && dp->cuboid->tag[tmpidx]>=0)
						{
							adj_part[cnt].pos[1] += boxsize;
							adj_part[cnt].pos[2] += boxsize;
							if(tmpidx>=0 && tmpidx<Ngrid) 
							{
								adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
								adj_part[cnt].group = tmpidx;
								dp->cuboid->count[tmpidx]++;
								if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
									printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
							}
							else
							{
								printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
							}
							cnt++;
							adj_part[cnt] = buffrecv[s][n];
						}
					}
				}
			}
			if (z==Gside-1)
			{
				tmpidx = getidx(x, y, z-Gside, Jn, Kn, minidx, npad);
				if (minidx[2]==0 && (inbox(x,y,z-Gside,minidx,maxidx)) && dp->cuboid->tag[tmpidx]>=0)
				{
					adj_part[cnt].pos[2] -= boxsize;
					if(tmpidx>=0 && tmpidx<Ngrid) 
					{
						adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
						adj_part[cnt].group = tmpidx;
						dp->cuboid->count[tmpidx]++;
						if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
							printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
					}
					else
					{
						printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
					}
					cnt++;
					adj_part[cnt] = buffrecv[s][n];
				}
			}
			if(z==0)
			{
				tmpidx = getidx(x, y, z+Gside, Jn, Kn, minidx, npad);
				if (maxidx[2]==Gside-1 && (inbox(x,y,z+Gside,minidx,maxidx)) && dp->cuboid->tag[tmpidx]>=0)
				{
					adj_part[cnt].pos[2] += boxsize;
					if(tmpidx>=0 && tmpidx<Ngrid) 
					{
						adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
						adj_part[cnt].group = tmpidx;
						dp->cuboid->count[tmpidx]++;
						if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
							printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
					}
					else
					{
						printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
					}
					cnt++;
					adj_part[cnt] = buffrecv[s][n];
				}
			}
	  		if(inbox(x,y,z,minidx,maxidx) && tag[idx]>=0)
			{
				tmpidx = idx;
				if(tmpidx>=0 && tmpidx<Ngrid) 
				{
					adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
					adj_part[cnt].group = tmpidx;
					dp->cuboid->count[tmpidx]++;
					if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
						printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
				}
				else
				{
					printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
				}
				cnt++;
				adj_part[cnt] = buffrecv[s][n];
			}
//			if(minidx[0] == 0 && dp->cuboid->tag[idx-(Gside-1)*Jn*Kn] == TAG_FRONTIERS)
//				{
//					x-=Gside;
//					adj_part[cnt].pos[0] -= boxsize;
//				}
//			}
//			else if (x==0)
//			{
//				if(maxidx[0] == Gside-1 && dp->cuboid->tag[idx+(Gside-1)*Jn*Kn] == TAG_FRONTIERS)
//				{
//					x+=Gside;
//					adj_part[cnt].pos[0] += boxsize;
//				}
//			}
//			if (y==Gside-1)
//			{
//				if(minidx[1] == 0 && dp->cuboid->tag[idx-(Gside-1)*Kn] == TAG_FRONTIERS)
//				{
//					y-=Gside;
//					adj_part[cnt].pos[1] -= boxsize;
//				}
//			}
//			else if (y==0)
//			{
//				if(maxidx[1] == Gside-1 && dp->cuboid->tag[idx+(Gside-1)*Kn] == TAG_FRONTIERS)
//				{
//					y+=Gside;
//					adj_part[cnt].pos[1] += boxsize;
//				}
//			}
//			if (z==Gside-1)
//			{
//				if(minidx[2] == 0 && dp->cuboid->tag[idx-Gside+1] == TAG_FRONTIERS)
//				{
//					z-=Gside;
//					adj_part[cnt].pos[2] -= boxsize;
//				}
//			}
//			else if (z==0)
//			{
//				if(maxidx[2] == Gside-1 && dp->cuboid->tag[idx+Gside-1] == TAG_FRONTIERS)
//				{
//					z+=Gside;
//					adj_part[cnt].pos[2] += boxsize;
//				}
//			}
//			x = x - minidx[0] + npad;
//			y = y - minidx[1] + npad;
//			z = z - minidx[2] + npad;
			
// 			if (adj_part[cnt].pos[0]+1/nbox_boxsize >= boxsize)
//			{
//				if(minidx[0]/nbox_boxsize<=0.0)
//				{
//					adj_part[cnt].pos[0] -= boxsize;
//				}
//			}
//			else if (adj_part[cnt].pos[0]-1/nbox_boxsize <= 0.0)
//			{
//				if((maxidx[0]+1.000001)/nbox_boxsize>=boxsize)
//				{
//					adj_part[cnt].pos[0] += boxsize;
//				}
//			}
//			if (adj_part[cnt].pos[1]+1/nbox_boxsize >= boxsize)
//			{
//				if(minidx[1]/nbox_boxsize<=0.0)
//				{
//					adj_part[cnt].pos[1] -= boxsize;
//				}
//			}
//			else if (adj_part[cnt].pos[1]-1/nbox_boxsize <= 0.0)
//			{
//				if((maxidx[1]+1.000001)/nbox_boxsize>=boxsize)
//				{
//					adj_part[cnt].pos[1] += boxsize;
//				}
//			}
//			if (adj_part[cnt].pos[2]+1/nbox_boxsize >= boxsize)
//			{
//				if(minidx[2]/nbox_boxsize<=0.0)
//				{
//					adj_part[cnt].pos[2] -= boxsize;
//				}
//			}
//			else if (adj_part[cnt].pos[2]-1/nbox_boxsize <= 0.0)
//			{
//				if((maxidx[2]+1.000001)/nbox_boxsize>=boxsize)
//				{
//					adj_part[cnt].pos[2] += boxsize;
//				}
//			}
//
//			if(tmpidx>=0 && tmpidx<Ngrid) 
//			{
//				adj_part[cnt].tag = dp->cuboid->tag[tmpidx];
//				adj_part[cnt].group = tmpidx;
//				dp->cuboid->count[tmpidx]++;
//				if(dp->cuboid->tag[tmpidx] == TAG_OUTERCORE || dp->cuboid->tag[tmpidx] == TAG_FRONTIERS)
//					printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
//			}
//			else
//			{
//				printf("[%d]  box=%f, idx=%d, x=%f, y=%f, z=%f.\n", dp->rank, 1/nbox_boxsize, idx, adj_part[cnt].pos[0], adj_part[cnt].pos[1], adj_part[cnt].pos[2]);
//			}

//	        cnt++;
        }
    }
//	assert(cnt==recv_tot_part);
	part = (Body*)realloc((void*)dp->Part,sizeof(Body)*(npart+cnt));
    adj_part = part + npart;
	dp->Part = part;
	dp->NumPart += cnt;
	printf("[%d] Numpart = %d.\n", dp->rank, dp->NumPart);

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

