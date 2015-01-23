#include "global.h"
#include "domain.h"
#include <assert.h>

Domain* create_domain(Constants *constants, Status *status, System *system)
{
    Domain* domain;
    domain = (Domain*)xmalloc(sizeof(Domain),6001);

    Body *part = system->part;
    Int  npart = system->num_part;

	int rank, Nproc;
    double pos_min[3];
    double pos_max[3];
    double boxsize = constants->BOX_SIZE;
    int    cm_nside = (1<<(constants->MESH_BITS));
    int    peano_bit= constants->PEANO_BITS;
    int    mesh_bit = constants->MESH_BITS;
    double r_split = constants->SPLIT_SCALE;
    double mass    = constants->PART_MASS;
    double grav_const = constants->GRAV_CONST;

    double invdelta = (cm_nside)/boxsize;
    double delta = boxsize/cm_nside;
    domain->cm_bin = (double*)xmalloc(sizeof(double)*(cm_nside+1), 6002);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &Nproc);

	domain->rank = rank;
	domain->nproc = Nproc;

    int s, d, i, mesh_length;
    for (d=0; d<3; d++) {
        pos_min[d] = boxsize;
        pos_max[d] = 0.0;
    }
    for (i=0; i<cm_nside; i++)
        domain->cm_bin[i] = (double)((float)(i*delta));
    domain->cm_bin[cm_nside] = boxsize;

    float pos[3];
    Int n;
    for (n=0; n<npart; n++)
	{
        for (d=0; d<3; d++)
		{
            pos[d] = (part[n].pos[d]);

            if (pos[d]>=pos_max[d])
                pos_max[d] = pos[d];
            if (pos[d]<=pos_min[d])
                pos_min[d] = pos[d];
        }
    }
    int mesh_min[3];
    int mesh_max[3];
    int mesh_nside[3];

    for (d=0; d<3; d++)
	{
        mesh_min[d] = (int)(pos_min[d]*invdelta);
        mesh_max[d] = (int)(pos_max[d]*invdelta);
        if (pos_min[d] < domain->cm_bin[mesh_min[d]])
            mesh_min[d] -= 1;
        if (pos_max[d] < domain->cm_bin[mesh_max[d]])
            mesh_max[d] -= 1;

        if (mesh_min[d] == (cm_nside))
            mesh_min[d] = (cm_nside)-1;

        if (mesh_max[d] == (cm_nside))
            mesh_max[d] = (cm_nside)-1;

        assert(mesh_max[d] < cm_nside);

        mesh_nside[d] = 1+(mesh_max[d]-mesh_min[d]+1)+1;
    }

    mesh_length = mesh_nside[0]*mesh_nside[1]*mesh_nside[2];
    LOG_INFO(LOG_VB1, "[%d] mesh_min = %d %d %d\n", rank, mesh_min[0], mesh_min[1], mesh_min[2]);
    LOG_INFO(LOG_VB1, "[%d] mesh_max = %d %d %d\n", rank, mesh_max[0], mesh_max[1], mesh_max[2]);
    LOG_INFO(LOG_VB1, "[%d] mesh_nside = %d %d %d\n", rank, mesh_nside[0], mesh_nside[1], mesh_nside[2]);

    domain->count = (Int*)xmalloc(sizeof(Int)*mesh_length,6003);
    domain->tag = (int*)xmalloc(sizeof(int)*mesh_length,6003);
    domain->mroot = (Int*)xmalloc(sizeof(Int)*mesh_length,6003);
    domain->pidx = (Int*)xmalloc(sizeof(Int)*mesh_length,6003);
	domain->tree = NULL;

    for (i=0; i<mesh_length; i++) {
        domain->count[i] = 0;
        domain->tag[i] = 0;
        domain->mroot[i] = -1;
        domain->pidx[i] = -1;
    }

    domain->lower = status->lower_rank;
    domain->upper = status->upper_rank;

    domain->boxsize = boxsize;
    domain->split = r_split;
    domain->mass  = mass;
    domain->grav_const =  grav_const;

    domain->minimum[0] = mesh_min[0];
    domain->minimum[1] = mesh_min[1];
    domain->minimum[2] = mesh_min[2];
    domain->maximum[0] = mesh_max[0];
    domain->maximum[1] = mesh_max[1];
    domain->maximum[2] = mesh_max[2];
    domain->nSide[0] = mesh_nside[0];
    domain->nSide[1] = mesh_nside[1];
    domain->nSide[2] = mesh_nside[2];

    domain->nPad = 1;
    domain->peano_bit = peano_bit;
    domain->mesh_bit  =  mesh_bit;
    domain->cm_nside  =  cm_nside;
    domain->part  = system->part;
    domain->npart = system->num_part;

    return domain;
}


void mark_cuboid(Domain *domain)
{
    int rank, Nproc;
    int i, j, k;
    int x, y, z, idx, cnt;
    int In, Jn, Kn;

	rank = domain->rank;
	Nproc = domain->nproc;

    int minidx[3];
    int maxidx[3];
    int nside[3];
    minidx[0] = domain->minimum[0];
    minidx[1] = domain->minimum[1];
    minidx[2] = domain->minimum[2];
    maxidx[0] = domain->maximum[0];
    maxidx[1] = domain->maximum[1];
    maxidx[2] = domain->maximum[2];
    nside[0] = domain->nSide[0];
    nside[1] = domain->nSide[1];
    nside[2] = domain->nSide[2];
//    Peano lower = domain->lower[rank];
//    Peano upper = domain->upper[rank];
    int npad = domain->nPad;
    int nbits= domain->mesh_bit;
    int cm_nside = domain->cm_nside;

    Int *count  = domain->count;
    int *tag = domain->tag;

    In = nside[0];
    Jn = nside[1];
    Kn = nside[2];

    Peano cubekey;
    Peano *cdkey = (Peano*)xmalloc(sizeof(Peano)*In*Jn*Kn, 6004);
    int *label = (int*)xmalloc(sizeof(int)*In*Jn*Kn, 6005);
    int lower_mesh, upper_mesh;
    int shift_bit = 3*(domain->peano_bit - domain->mesh_bit);
    Peano lower, upper;
    lower = domain->lower[rank]>>shift_bit;
    upper = domain->upper[rank]>>shift_bit;

    upper_mesh = lower_mesh = -1;
    for ( i=0; i<In; i++ )
	{
        for ( j=0; j<Jn; j++ )
		{
            for ( k=0; k<Kn; k++ )
			{
                x = i - npad + minidx[0];
                y = j - npad + minidx[1];
                z = k - npad + minidx[2];
                idx = ( i*Jn + j )*Kn + k;
                //         count[idx] = 0;
                tag[idx] = TAG_IDLEFIELD;

                if ( i<npad || i>(In-npad-1) ||
                        j<npad || j>(Jn-npad-1) ||
                        k<npad || k>(Kn-npad-1)   )
                {
                    label[idx] = 0;
                    if ( i<npad || i>(In-npad-1) )
                        x = (x+cm_nside)%cm_nside;
                    if ( j<npad || j>(Jn-npad-1) )
                        y = (y+cm_nside)%cm_nside;
                    if ( k<npad || k>(Kn-npad-1) )
                        z = (z+cm_nside)%cm_nside;
                    cdkey[idx] = encoding_peano(x, y, z, nbits);
                }
                else
				{
                    cubekey = encoding_peano(x, y, z, nbits);
                    if ( cubekey >= lower && cubekey < upper )
                        label[idx] = 1;
                    else
                        label[idx] = 0;

                    cdkey[idx] = cubekey;
                }
//    printf(" lower = %d upper = %d \n", lower, upper);
                if (lower == cdkey[idx])
                    lower_mesh = idx;
                if ((upper-1) == cdkey[idx])
                    upper_mesh = idx;
            }
        }
    }
    LOG_INFO(LOG_VB1, "[%d] lower_mesh = %d upper_mesh = %d lower = %d upper = %d [%d, %d]\n", rank, lower_mesh, upper_mesh, lower, upper, cdkey[lower_mesh], cdkey[upper_mesh]);

    assert( 1==npad );

    for ( i=0; i<In; i++ )
	{
        for ( j=0; j<Jn; j++ )
		{
            for ( k=0; k<Kn; k++ )
			{
                idx = ( i*Jn + j )*Kn + k;

                cnt = 0;
                for (x=i-1; x<=i+1; x++)
				{
					if (x<0 || x>=In) continue;
                    for (y=j-1; y<=j+1; y++)
					{
						if (y<0 || y>=Jn) continue;
                        for (z=k-1; z<=k+1; z++)
						{
							if (z>=0 && z<Kn)
	                            cnt += label[(x*Jn + y )*Kn + z];
                        }
                    }
                }
                switch(cnt)
				{
					case 27:
						tag[idx] =  TAG_OUTERCORE;
						break;
					case  0:
						tag[idx] =  TAG_IDLEFIELD;
						break;
					default:
						if ( 1==label[idx] )
							tag[idx] = TAG_FRONTIERS;
						else if ( 0==label[idx] )
							tag[idx] = TAG_ADJOINING;
						break;
                }
            }
        }
    }

/////////////////////////////////////
    int nSplit, nLevel, splitLength;
    nSplit = Nproc;
    for (nLevel = 0; nSplit > (1<<nLevel); nLevel++);
    splitLength = 1<<nLevel;
    Peano *split = (Peano*) malloc(sizeof(Peano) * splitLength );

//   printf("nSplit = %d nLevel = %d\n", nSplit, nLevel);
    for (i=0; i<splitLength; i++) {
        if (i<nSplit)
            split[i] = ( domain->lower[i] >> shift_bit );
        else
            split[i] = ( domain->upper[nSplit-1] >> shift_bit);
//        printf(" split[%d] = %lu \n",i, split[i]);
    }

    int dom, l;
    int level = 1<< (nLevel);

    level >>= 1;
    for ( i=0; i<In; i++)
	{
        for ( j=0; j<Jn; j++)
		{
            for ( k=0; k<Kn; k++)
			{
                idx = (i*Jn+j)*Kn + k ;
                if ( TAG_ADJOINING == tag[idx] )
				{
                    cubekey = cdkey[idx];

                    /* determine which domin this cube belongs to */
                    l = dom = level ;
                    while ( l >>= 1 )
					{
                        if (cubekey < split[dom])
						{
                            dom -= l;
                        }
                        else
						{
                            dom += l;
                        }
                    }
                    if (cubekey < split[dom])
					{
                        dom -= 1;
                    }
                    /* determine which domin this cube belongs to */
                    tag[idx] = dom;
                }
            }
        }
    }


    MPI_Barrier(TREE_COMM_WORLD);
//    printf("mesh_bit=%d peano_bit=%d\n",domain->mesh_bit,domain->peano_bit);

    if ( -1 != upper_mesh ) {
        cubekey = cdkey[upper_mesh];

        /* determine which domin this cube belongs to */
        l = dom = level ;
        while ( l >>= 1 ) {
            if (cubekey < split[dom]) {
                dom -= l;
            }
            else {
                dom += l;
            }
        }
        if (cubekey < split[dom]) {
            dom -= 1;
        }
        /* determine which domin this cube belongs to */
//        printf("[%d] upper dom = %d\n", rank, dom);
        //  tag[idx] = dom;
    }
    if ( -1 != lower_mesh ) {
        cubekey = cdkey[lower_mesh];

        /* determine which domin this cube belongs to */
        l = dom = level ;
        while ( l >>= 1 ) {
            if (cubekey < split[dom]) {
                dom -= l;
            }
            else {
                dom += l;
            }
        }
        if (cubekey < split[dom]) {
            dom -= 1;
        }
        /* determine which domin this cube belongs to */
//        printf("[%d] lower dom = %d\n", rank, dom);
    }
// if rank != 0 && Nproc-1 { if (-1 != lower_mesh) send lower to rank-1, and if (-1 != upper_mesh) send upper to rank+1 }

    Body *part = domain->part;
//    int I, J, K;
    Int n;
    Int npart = domain->npart;
    double invdelta = (double)cm_nside/(domain->boxsize);
    double pos[3];
    double *bin  = domain->cm_bin;


    for (n=0; n<npart; n++)
	{
        pos[0] = (double)part[n].pos[0];
        pos[1] = (double)part[n].pos[1];
        pos[2] = (double)part[n].pos[2];

        x = (int)(pos[0]*invdelta);
        if (pos[0] < bin[x] )
            x -= 1;

        y = (int)(pos[1]*invdelta);
        if (pos[1] < bin[y] )
            y -= 1;

        z = (int)(pos[2]*invdelta);
        if (pos[2] < bin[z] )
            z -= 1;

        x += 1 - minidx[0];
        y += 1 - minidx[1];
        z += 1 - minidx[2];

        idx = (x*Jn + y) * Kn + z;
        part[n].tag = tag[idx];
        part[n].group = idx;
        count[idx]++;
    }

    free(split);
    free(cdkey);
    free(label);
}



void broadcast_frontiers(Domain *domain, System *sys)
{
    int minidx[3];
    int maxidx[3];
    int nside[3];
    int i, j, k, In, Jn, Kn;
    int x, y, z, dest, s, d, idx;
    Int n;

	int rank = domain->rank;
	int Nproc = domain->nproc;
    int npad = domain->nPad;
    int nbits= domain->mesh_bit;
    int cm_nside = domain->cm_nside;

    Int *count  = domain->count;
    int *tag = domain->tag;

    Int npart = domain->npart;
    Body *part = domain->part;

    minidx[0] = domain->minimum[0];
    minidx[1] = domain->minimum[1];
    minidx[2] = domain->minimum[2];
    maxidx[0] = domain->maximum[0];
    maxidx[1] = domain->maximum[1];
    maxidx[2] = domain->maximum[2];
    nside[0] = domain->nSide[0];
    nside[1] = domain->nSide[1];
    nside[2] = domain->nSide[2];


    int *numlock, *numpart;
    int *numsend, *numrecv;
    Body **buffsend, **buffrecv;
//    Body *adj_part;

    numsend = (int*) malloc(sizeof(int)*Nproc);
    numrecv = (int*) malloc(sizeof(int)*Nproc);
    numlock = (int*) malloc(sizeof(int)*Nproc);
    numpart = (int*) malloc(sizeof(int)*Nproc);

    buffsend = (Body**) malloc(sizeof(Body*)*Nproc);
    buffrecv = (Body**) malloc(sizeof(Body*)*Nproc);

    for (s=0; s<Nproc; s++)
	{
        numpart[s] = 0;
        numsend[s] = 0;
        numrecv[s] = 0;
        numlock[s] = 1;
        buffsend[s] = NULL;
        buffrecv[s] = NULL;
    }

    In = nside[0];
    Jn = nside[1];
    Kn = nside[2];
    int count_frontiers = 0;
    for (i=0; i<In; i++)
	{
        for (j=0; j<Jn; j++)
		{
            for (k=0; k<Kn; k++)
			{
                idx = (i*Jn+j)*Kn+k;
                if ( tag[idx] == TAG_FRONTIERS )
                    count_frontiers += count[idx];
            }
        }
    }
//    printf(" [%d] frontiers = %d\n", rank, count_frontiers);


    int cnt, active[27];
    for (i=0; i<In; i++)
        for (j=0; j<Jn; j++)
            for (k=0; k<Kn; k++)
			{
                idx = (i*Jn + j )*Kn + k;
                if (TAG_FRONTIERS == tag[idx])
				{
                    cnt = 0;
                    for (x=i-1; x<=i+1; x++)
                        for (y=j-1; y<=j+1; y++)
                            for (z=k-1; z<=k+1; z++)
							{
                                dest = tag[(x*Jn+y)*Kn+z];
                                if ( dest  >= 0 && numlock[dest] )
								{
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
    /*
        for (s=0; s<Nproc; s++) {
            printf(" send %10d part(s) from[%d]to[%d]\n", numsend[s], rank, s);
        }
    */
    MPI_Alltoall(numsend, 1, MPI_INT, numrecv, 1, MPI_INT, TREE_COMM_WORLD);

    for (s=0; s<Nproc; s++)
	{
        DBG_INFO(0, "[%d]: Recv %10d part(s) from [%d]\n", rank, numrecv[s], s);
    }

    int recv_tot_part, send_tot_part;

    send_tot_part = recv_tot_part = 0;
    for (s=0; s<Nproc; s++)
	{
        recv_tot_part += numrecv[s];
        send_tot_part += numsend[s];
    }
    LOG_INFO(LOG_VB1, " [%d] tot_send = %d  tot_recv = %d domain_part = %d\n", rank, send_tot_part, recv_tot_part, domain->npart);

    for (s=0; s<Nproc; s++)
	{
        if ( 0<numsend[s] )
            buffsend[s] = (Body*)malloc(sizeof(Body) * numsend[s] );
        if ( 0<numrecv[s] )
            buffrecv[s] = (Body*)malloc(sizeof(Body) * numrecv[s] );
    }

    int Ngrid = In*Jn*Kn;
    double boxsize = domain->boxsize;
    double nbox_boxsize = (1<<nbits)/boxsize;
    double invdelta = (double)cm_nside/(domain->boxsize);
    double *bin  = domain->cm_bin;
    double pos[3];

//    printf(" npart = %d napd = %d\n", npart, npad);
    /*

        for (n=0; n<npart; n++) {
            pos[0] = (double)part[n].pos[0];
            pos[1] = (double)part[n].pos[1];
            pos[2] = (double)part[n].pos[2];

            i = (int)(pos[0]*nbox_boxsize);
            if (pos[0] < bin[x] )
                i -= 1;

            j = (int)(pos[1]*nbox_boxsize);
            if (pos[1] < bin[y] )
                j -= 1;

            k = (int)(pos[2]*nbox_boxsize);
            if (pos[2] < bin[z] )
                k -= 1;

            i += npad - minidx[0];
            j += npad - minidx[1];
            k += npad - minidx[2];

    */
    int nb, cntpart, start, index_length;
    int flag;
    for (n=0; n<npart; n++)
	{
        if (TAG_FRONTIERS == part[n].tag)
		{
            pos[0] = (double)part[n].pos[0];
            pos[1] = (double)part[n].pos[1];
            pos[2] = (double)part[n].pos[2];

            i = (int)(pos[0]*nbox_boxsize);
            if (pos[0] < bin[i] )
                i -= 1;

            j = (int)(pos[1]*nbox_boxsize);
            if (pos[1] < bin[j] )
                j -= 1;

            k = (int)(pos[2]*nbox_boxsize);
            if (pos[2] < bin[k] )
                k -= 1;

            flag = 0; 
			/* flag to indicate the dimensions of the particle at the border.
			 * 1-6 bits means: Xmin Xmax Ymin Ymax Zmin Zmax.
			 * upper bits: the counter of the borders. */
            if (0==i)
			{
				flag |= 0x1;
				flag += 0x40;
			}
			if (cm_nside-1==i)
			{
				flag |= 0x2;
				flag += 0x40;
			}
            if (0==j)
			{
				flag |= 0x4;
				flag += 0x40;
			}
			if (cm_nside-1==j)
			{
				flag |= 0x8;
				flag += 0x40;
			}
			if (0==k)
			{
				flag |= 0x10;
				flag += 0x40;
			}
            if (cm_nside-1==k)
			{
				flag |= 0x20;
				flag += 0x40;
			}

            i += npad - minidx[0];
            j += npad - minidx[1];
            k += npad - minidx[2];

//            i = (int)(part[n].pos[0]*nbox_boxsize) - minidx[0] + npad;
//            j = (int)(part[n].pos[1]*nbox_boxsize) - minidx[1] + npad;
//            k = (int)(part[n].pos[2]*nbox_boxsize) - minidx[2] + npad;

            cnt =0 ;
            for (x=i-1; x<=i+1; x++)
                for (y=j-1; y<=j+1; y++)
                    for (z=k-1; z<=k+1; z++)
					{
                        dest = tag[(x*Jn+y)*Kn+z];
                        if ( dest >= 0 && numlock[dest] )
						{
                            nb = numpart[dest];

                            buffsend[dest][nb] = part[n];
                            buffsend[dest][nb].tag = flag;

                            numpart[dest]++;
                            numlock[dest] = 0;
                            active[cnt++] = dest;
                        }
                    }

            while ( 0< cnt-- )
                numlock[ active[cnt] ] = 1;
        }
    }


    for (s=0; s<Nproc; s++) {
        assert(numlock[s]);
        assert(numpart[s] == numsend[s]);
    }

    int data_btye = sizeof(Body);
    int recv, send;
    MPI_Status *status;
    MPI_Request *req_send, *req_recv;

    req_send = (MPI_Request*)malloc(sizeof(MPI_Request)*Nproc);
    req_recv = (MPI_Request*)malloc(sizeof(MPI_Request)*Nproc);
    status = (MPI_Status*)malloc(sizeof(MPI_Status)*Nproc);

    MPI_Datatype mpi_buffdata_type;
    MPI_Type_contiguous(sizeof(Body), MPI_CHAR, &mpi_buffdata_type );
    MPI_Type_commit(&mpi_buffdata_type);
    MPI_Barrier(TREE_COMM_WORLD);
    for (s=0; s<Nproc; s++) {
        send = (rank+s)%Nproc;
        recv = (rank-s+Nproc)%Nproc;
//        printf("%d:  %d ->[%d] :%d\n",s,recv, rank, numrecv[recv]);
        MPI_Sendrecv(buffsend[send], numsend[send], mpi_buffdata_type,
                     send, 0, buffrecv[recv], numrecv[recv], mpi_buffdata_type,
                     recv, 0, TREE_COMM_WORLD, &status[s]);
    }

//    int body_size = (int)sizeof(Body);
//    MPI_Barrier(TREE_COMM_WORLD);
//    for (s=0; s<Nproc; s++) {
//        send = (rank+s)%Nproc;
//        recv = (rank-s+Nproc)%Nproc;
//        printf("%d:  %d ->[%d] :%d\n",s,recv, rank, numrecv[recv]);
//        MPI_Sendrecv(buffsend[send], numsend[send]*body_size, MPI_CHAR,
//                    send, 0, buffrecv[recv], numrecv[recv]*body_size, MPI_CHAR,
//                    recv, 0, TREE_COMM_WORLD, &status[s]);
//    }

//    assert( cm_nside > 1 ); //  if ( 1 == cm_nside ) try 27-if scopes please;

	int rgflag = 0, gmask = 0, dcount = 0, pcor = 1;
	int xl, yl, zl;
	xl = yl = zl = 1;
	/* rgflag: exchange max and min flag of each dimention.
	 * gmask: the mask to hide the dimension with 1 mesh width.
	 * dcount: the counter of the dimension with 1 mesh width.
	 * pcor: coefficient of a particle should be copied. */

	if(minidx[0]==0) /* X dim in this grid touch the Box lower border. */
	{
		rgflag |= 0x2;
	}
	if(minidx[1]==0) /* Y dim in this grid touch the Box lower border. */
	{
		rgflag |= 0x8;
	}
	if(minidx[2]==0) /* Z dim in this grid touch the Box lower border. */
	{
		rgflag |= 0x20;
	}
	if(maxidx[0]==cm_nside-1) /* X dim in this grid touch the Box upper border. */
	{
		rgflag |= 0x1;
	}
	if(maxidx[1]==cm_nside-1) /* Y dim in this grid touch the Box upper border. */
	{
		rgflag |= 0x4;
	}
	if(maxidx[2]==cm_nside-1) /* Z dim in this grid touch the Box upper border. */
	{
		rgflag |= 0x10;
	}
	if (minidx[0] == maxidx[0])
	{
		gmask |= 0x3;
		xl = 3;
		dcount++;
		pcor *= 3;
	}
	if (minidx[1] == maxidx[1])
	{
		gmask |= 0xC;
		yl = 3;
		dcount++;
		pcor *= 3;
	}
	if (minidx[2] == maxidx[2])
	{
		gmask |= 0x30;
		zl = 3;
		dcount++;
		pcor *= 3;
	}
	gmask ^= 0x3F;
	rgflag &= gmask;

    Int np_adj_max = 0L;
    for (s=0; s<Nproc; s++)
	{
        for (nb=0; nb<numrecv[s]; nb++)
		{
			int ptag = buffrecv[s][nb].tag;
			np_adj_max += (1 + ((ptag & rgflag) ? ((1<<((ptag>>6) - dcount*2)) - 1) : 0)) * pcor;
        }
    }
    domain->part = xrealloc((void*)part, sizeof(Body)*(npart + np_adj_max), 6010);
	part = domain->part;
	sys->part = part;

//    int *tag = domain->tag;
    Int np_adj = 0;
    int np;
    Body *ghost = part + npart;

    Body tmp_body;
    int xg, yg, zg, xn, yn, zn;
    double dpos[3] = {0, boxsize, -boxsize};
	int dmesh[3] = {0, cm_nside, -cm_nside};
	int gminx, gminy, gminz, gmaxx, gmaxy, gmaxz;

	gminx = minidx[0] - 1;
	gminy = minidx[1] - 1;
	gminz = minidx[2] - 1;
	gmaxx = maxidx[0] + 1;
	gmaxy = maxidx[1] + 1;
	gmaxz = maxidx[2] + 1;

	if (dcount == 0) /* dcount == 0 */
	{
		for (s=0; s<Nproc; s++)
		{
			for (nb=0; nb<numrecv[s]; nb++)
			{
				flag = buffrecv[s][nb].tag;
				pos[0] = buffrecv[s][nb].pos[0];
				pos[1] = buffrecv[s][nb].pos[1];
				pos[2] = buffrecv[s][nb].pos[2];

				xg = (int)(pos[0]*nbox_boxsize);
				if (pos[0] < bin[xg] )
					xg -= 1;

				yg = (int)(pos[1]*nbox_boxsize);
				if (pos[1] < bin[yg] )
					yg -= 1;

				zg = (int)(pos[2]*nbox_boxsize);
				if (pos[2] < bin[zg] )
					zg -= 1;

				x = xg + 1 - minidx[0];
				y = yg + 1 - minidx[1];
				z = zg + 1 - minidx[2];

				np = 0;
				if ( !(flag & rgflag) ) /* !(flag & rgflag) */
				{
					idx = (x*Jn + y) * Kn + z;
					if (tag[idx]>=0)
					{
						count[idx] ++;
						ghost[np_adj+np] = buffrecv[s][nb];
						ghost[np_adj+np].tag = tag[idx];
						ghost[np_adj+np].group = idx;
						np++;
					}
				}
				else /* flag & rgflag */
				{
					int xflag = flag & rgflag & 0x3;
					int yflag = (flag & rgflag & 0xC) >> 2;
					int zflag = (flag & rgflag & 0x30) >> 4;

					xn = x + dmesh[xflag];
					yn = y + dmesh[yflag];
					zn = z + dmesh[zflag];

					idx = (x*Jn + y) * Kn + z;
					if (tag[idx]>=0 && xg>=gminx && xg<=gmaxx && yg>=gminy && yg<=gmaxy && zg>=gminz && zg<=gmaxz)
					{
						count[idx] ++;
						ghost[np_adj+np] = buffrecv[s][nb];
						ghost[np_adj+np].tag = tag[idx];
						ghost[np_adj+np].group = idx;
						np++;
					}
					if (xflag) /* xflag */
					{
						if (yg>=gminy && yg<=gmaxy && zg>=gminz && zg<=gmaxz)
						{
							idx = (xn*Jn + y) * Kn + z;
							if(tag[idx]>=0)
							{
								count[idx] ++;
								ghost[np_adj+np] = buffrecv[s][nb];
								ghost[np_adj+np].pos[0] += dpos[xflag];
								ghost[np_adj+np].tag = tag[idx];
								ghost[np_adj+np].group = idx;
								np++;
							}
						}
						if (yflag) /* xflag && yflag */
						{
							if (zg>=gminz && zg<=gmaxz)
							{
								idx = (xn*Jn + yn) * Kn + z;
								if(tag[idx]>=0)
								{
									count[idx] ++;
									ghost[np_adj+np] = buffrecv[s][nb];
									ghost[np_adj+np].pos[0] += dpos[xflag];
									ghost[np_adj+np].pos[1] += dpos[yflag];
									ghost[np_adj+np].tag = tag[idx];
									ghost[np_adj+np].group = idx;
									np++;
								}
							}
							if (zflag) /* xflag && yflag && zflag */
							{
								idx = (xn*Jn + yn) * Kn + zn;
								if(tag[idx]>=0)
								{
									count[idx] ++;
									ghost[np_adj+np] = buffrecv[s][nb];
									ghost[np_adj+np].pos[0] += dpos[xflag];
									ghost[np_adj+np].pos[1] += dpos[yflag];
									ghost[np_adj+np].pos[2] += dpos[zflag];
									ghost[np_adj+np].tag = tag[idx];
									ghost[np_adj+np].group = idx;
									np++;
								}
							} /* xflag && yflag && zflag */
						} /* xflag && yflag */
						if (zflag) /* xflag && zflag */
						{
							if (yg>=gminy && yg<=gmaxy)
							{
								idx = (xn*Jn + y) * Kn + zn;
								if(tag[idx]>=0)
								{
									count[idx] ++;
									ghost[np_adj+np] = buffrecv[s][nb];
									ghost[np_adj+np].pos[0] += dpos[xflag];
									ghost[np_adj+np].pos[2] += dpos[zflag];
									ghost[np_adj+np].tag = tag[idx];
									ghost[np_adj+np].group = idx;
									np++;
								}
							}
						} /* xflag && zflag */
					} /* xflag */
					if (yflag) /* yflag */
					{
						if (xg>=gminx && xg<=gmaxx && zg>=gminz && zg<=gmaxz)
						{
							idx = (x*Jn + yn) * Kn + z;
							if(tag[idx]>=0)
							{
								count[idx] ++;
								ghost[np_adj+np] = buffrecv[s][nb];
								ghost[np_adj+np].pos[1] += dpos[yflag];
								ghost[np_adj+np].tag = tag[idx];
								ghost[np_adj+np].group = idx;
								np++;
							}
						}
						if (zflag) /* yflag && zflag */
						{
							if (xg>=gminx && xg<=gmaxx)
							{
								idx = (x*Jn + yn) * Kn + zn;
								if(tag[idx]>=0)
								{
									count[idx] ++;
									ghost[np_adj+np] = buffrecv[s][nb];
									ghost[np_adj+np].pos[1] += dpos[yflag];
									ghost[np_adj+np].pos[2] += dpos[zflag];
									ghost[np_adj+np].tag = tag[idx];
									ghost[np_adj+np].group = idx;
									np++;
								}
							}
						} /* yflag && zflag */
					} /* yflag */
					if (zflag) /* zflag */
					{
						if (xg>=gminx && xg<=gmaxx && yg>=gminy && yg<=gmaxy)
						{
							idx = (x*Jn + y) * Kn + zn;
							if(tag[idx]>=0)
							{
								count[idx] ++;
								ghost[np_adj+np] = buffrecv[s][nb];
								ghost[np_adj+np].pos[2] += dpos[zflag];
								ghost[np_adj+np].tag = tag[idx];
								ghost[np_adj+np].group = idx;
								np++;
							}
						}
					} /* zflag */
				} /* flag & rgflag */
				np_adj += np;
			}
        }
    }
	else /* dcount > 0 */
	{
		for (s=0; s<Nproc; s++)
		{
			for (nb=0; nb<numrecv[s]; nb++)
			{
				flag = buffrecv[s][nb].tag;
				pos[0] = buffrecv[s][nb].pos[0];
				pos[1] = buffrecv[s][nb].pos[1];
				pos[2] = buffrecv[s][nb].pos[2];

				xg = (int)(pos[0]*nbox_boxsize);
				if (pos[0] < bin[xg] )
					xg -= 1;

				yg = (int)(pos[1]*nbox_boxsize);
				if (pos[1] < bin[yg] )
					yg -= 1;

				zg = (int)(pos[2]*nbox_boxsize);
				if (pos[2] < bin[zg] )
					zg -= 1;

				x = xg + 1 - minidx[0];
				y = yg + 1 - minidx[1];
				z = zg + 1 - minidx[2];

				np = 0;

				int xflag = flag & rgflag & 0x3;
				int yflag = (flag & rgflag & 0xC) >> 2;
				int zflag = (flag & rgflag & 0x30) >> 4;

				if (xg>=gminx && xg<=gmaxx && yg>=gminy && yg<=gmaxy && zg>=gminz && zg<=gmaxz)
				{
					for (i=0; i<xl; i++)
					{
						xn = x + dmesh[i];
						for (j=0; j<yl; j++)
						{
							yn = y + dmesh[j];
							for (k=0; k<zl; k++)
							{
								zn = z + dmesh[k];
								idx = (xn*Jn + yn) * Kn + zn;

								if(tag[idx]>=0)
								{
									count[idx] ++;
									ghost[np_adj+np] = buffrecv[s][nb];
									ghost[np_adj+np].pos[0] += dpos[i]; 
									ghost[np_adj+np].pos[1] += dpos[j]; 
									ghost[np_adj+np].pos[2] += dpos[k]; 
									ghost[np_adj+np].tag = tag[idx];
									ghost[np_adj+np].group = idx;
									np++;
								}
							}
						}
					}
				}

				if (xflag) /* xflag */
				{
					i = xflag;
					xn = x + dmesh[i];
					if (yg>=gminy && yg<=gmaxy && zg>=gminz && zg<=gmaxz)
					{
						for (j=0; j<yl; j++)
						{
							yn = y + dmesh[j];
							for (k=0; k<zl; k++)
							{
								zn = z + dmesh[k];
								idx = (xn*Jn + yn) * Kn + zn;

								if(tag[idx]>=0)
								{
									count[idx] ++;
									ghost[np_adj+np] = buffrecv[s][nb];
									ghost[np_adj+np].pos[0] += dpos[i]; 
									ghost[np_adj+np].pos[1] += dpos[j]; 
									ghost[np_adj+np].pos[2] += dpos[k]; 
									ghost[np_adj+np].tag = tag[idx];
									ghost[np_adj+np].group = idx;
									np++;
								}
							}
						}
					}
					if (yflag) /* xflag && yflag */
					{
						j = yflag;
						yn = y + dmesh[j];
						if (zg>=gminz && zg<=gmaxz)
						{
							for (k=0; k<zl; k++)
							{
								zn = z + dmesh[k];
								idx = (xn*Jn + yn) * Kn + zn;

								if(tag[idx]>=0)
								{
									count[idx] ++;
									ghost[np_adj+np] = buffrecv[s][nb];
									ghost[np_adj+np].pos[0] += dpos[i]; 
									ghost[np_adj+np].pos[1] += dpos[j]; 
									ghost[np_adj+np].pos[2] += dpos[k]; 
									ghost[np_adj+np].tag = tag[idx];
									ghost[np_adj+np].group = idx;
									np++;
								}
							}
						}
					} /* xflag && yflag */
					if (zflag) /* xflag && zflag */
					{
						k = zflag;
						zn = z + dmesh[k];
						if (yg>=gminy && yg<=gmaxy)
						{
							for (j=0; j<yl; j++)
							{
								yn = y + dmesh[j];
								idx = (xn*Jn + yn) * Kn + zn;

								if(tag[idx]>=0)
								{
									count[idx] ++;
									ghost[np_adj+np] = buffrecv[s][nb];
									ghost[np_adj+np].pos[0] += dpos[i]; 
									ghost[np_adj+np].pos[1] += dpos[j]; 
									ghost[np_adj+np].pos[2] += dpos[k]; 
									ghost[np_adj+np].tag = tag[idx];
									ghost[np_adj+np].group = idx;
									np++;
								}
							}
						}
					} /* xflag && zflag */
				} /* xflag */
				if (yflag) /* yflag */
				{
					j = yflag;
					yn = y + dmesh[j];
					if (xg>=gminx && xg<=gmaxx && zg>=gminz && zg<=gmaxz)
					{
						for (i=0; i<xl; i++)
						{
							xn = x + dmesh[i];
							for (k=0; k<zl; k++)
							{
								zn = z + dmesh[k];
								idx = (xn*Jn + yn) * Kn + zn;

								if(tag[idx]>=0)
								{
									count[idx] ++;
									ghost[np_adj+np] = buffrecv[s][nb];
									ghost[np_adj+np].pos[0] += dpos[i]; 
									ghost[np_adj+np].pos[1] += dpos[j]; 
									ghost[np_adj+np].pos[2] += dpos[k]; 
									ghost[np_adj+np].tag = tag[idx];
									ghost[np_adj+np].group = idx;
									np++;
								}
							}
						}
					}
					if (zflag) /* yflag && zflag */
					{
						k = zflag;
						zn = z + dmesh[k];
						if (xg>=gminx && xg<=gmaxx)
						{
							for (i=0; i<xl; i++)
							{
								xn = x + dmesh[i];
								idx = (xn*Jn + yn) * Kn + zn;

								if(tag[idx]>=0)
								{
									count[idx] ++;
									ghost[np_adj+np] = buffrecv[s][nb];
									ghost[np_adj+np].pos[0] += dpos[i]; 
									ghost[np_adj+np].pos[1] += dpos[j]; 
									ghost[np_adj+np].pos[2] += dpos[k]; 
									ghost[np_adj+np].tag = tag[idx];
									ghost[np_adj+np].group = idx;
									np++;
								}
							}
						}
					} /* yflag && zflag */
				}
				if (zflag) /* zflag */
				{
					k = zflag;
					zn = z + dmesh[k];
					if (xg>=gminx && xg<=gmaxx && yg>=gminy && yg<=gmaxy)
					{
						for (i=0; i<xl; i++)
						{
							xn = x + dmesh[i];
							for (j=0; j<yl; j++)
							{
								yn = y + dmesh[j];
								idx = (xn*Jn + yn) * Kn + zn;

								if(tag[idx]>=0)
								{
									count[idx] ++;
									ghost[np_adj+np] = buffrecv[s][nb];
									ghost[np_adj+np].pos[0] += dpos[i]; 
									ghost[np_adj+np].pos[1] += dpos[j]; 
									ghost[np_adj+np].pos[2] += dpos[k]; 
									ghost[np_adj+np].tag = tag[idx];
									ghost[np_adj+np].group = idx;
									np++;
								}
							}
						}
					}
				} /* zflag */
				np_adj += np;
            }
        }
    }
    LOG_INFO(LOG_VB1, "[%d] npart_adj = %d [max=%d]\n", rank, (int)np_adj, (int)np_adj_max);
    domain->npart_adj = np_adj;

	/* Indicate the first particle index of each mesh in the contiguous storage scheme. */
	npart = 0;
	np_adj = 0;
	for (i=0; i<Ngrid; i++)
	{
		if (tag[i] == TAG_OUTERCORE || tag[i] == TAG_FRONTIERS)
		{
			domain->pidx[i] = npart;
			npart += count[i];
		}
		else if (tag[i] >= 0)
		{
			domain->pidx[i] = domain->npart + np_adj;
			np_adj += count[i];
		}
		if (count[i]==0)
		{
			domain->pidx[i] = -1;
		}
//		DBG_INFO(DBG_MEM, "[%d] count[%d] = %d, pidx[%d] = %d.\n", rank, i, count[i], i, domain->pidx[i]);
	}

    LOG_INFO(LOG_VB1, "[%d] npart = %d, npart_adj = %d [max=%d]\n", rank, npart, (int)np_adj, (int)np_adj_max);

	assert(npart == domain->npart);
	assert(np_adj == domain->npart_adj);

    if (NULL != buffsend) {
        for (s=0; s<Nproc; s++) 
            if (NULL != buffsend[s])
                free(buffsend[s]);
    }
    free(buffsend);
    if (NULL != buffrecv) {
        for (s=0; s<Nproc; s++) 
            if (NULL != buffrecv[s])
                free(buffrecv[s]);
    }
    free(buffrecv);



    /*
       if (0==rank) {
           for (i=0; i<In; i++) {
               for (j=0; j<Jn; j++) {
                   for (k=0; k<Kn; k++){
                       idx = (i*Jn+j)*Kn+k ;
                       if (tag[idx] >= 0)
                           printf("%4d ", tag[idx] );
                       else
                           printf("%4d ", tag[idx] );
                   }
                   printf("\n");
               }
               printf("\n");
           }
       }
       */


    free(req_send);
    free(req_recv);
    free(status);

    MPI_Barrier(TREE_COMM_WORLD);
    free(numpart);
    free(numlock);
    free(numrecv);
    free(numsend);
}

void free_domain(Domain *domain)
{
    if (NULL!=domain) {
        if ( NULL != domain->count )
            free(domain->count);
        if ( NULL != domain->tag )
            free(domain->tag);
        if ( NULL != domain->cm_bin )
            free(domain->cm_bin);
        if ( NULL != domain->mroot )
            free(domain->mroot);
		if ( NULL != domain->tree )
			free(domain->tree);
		if ( NULL != domain->pidx )
			free(domain->pidx);
        free(domain);
    }
}

