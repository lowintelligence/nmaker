#include "partmesh.h"
#include <mpi.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <fftw.h>
#include <rfftw_mpi.h>
#include <assert.h>
#include <math.h>

static int *pm_local_start, *pm_local_xside;
static fftw_real *data, *work;
static fftw_complex *cdata;
static rfftwnd_mpi_plan forward_plan,  inverse_plan;
static int pm_nside;
static real pos_min[3];
static real pos_max[3];
static int mesh_min[3];
static int mesh_max[3];
static int mesh_nside[3];
static real *pm_real_field;
static real **recv_block;





void setup_partmesh_environment(Domain *dp, GlobalParam *gp)
{
    int rank, size;
    int num_pm_socket = gp->NumSocket;
    pm_nside = 6*(1<<gp->NumBits);
    int local_nx, local_x_start, local_ny_after_transpose,
        local_y_start_after_transpose, total_local_size;

    pm_local_start = (int*)malloc(sizeof(int)*num_pm_socket);
    pm_local_xside = (int*)malloc(sizeof(int)*num_pm_socket);

    MPI_Comm_rank(MESH_COMM_WORLD, &rank);
    MPI_Comm_size(MESH_COMM_WORLD, &size);

    forward_plan = rfftw3d_mpi_create_plan(MESH_COMM_WORLD,
                                           pm_nside, pm_nside, pm_nside,
                                           FFTW_REAL_TO_COMPLEX,
                                           FFTW_ESTIMATE);
     inverse_plan = rfftw3d_mpi_create_plan(MESH_COMM_WORLD,
                                            pm_nside, pm_nside, pm_nside,
                                            FFTW_COMPLEX_TO_REAL,
                                            FFTW_ESTIMATE);
    rfftwnd_mpi_local_sizes(forward_plan, &local_nx, &local_x_start,
                            &local_ny_after_transpose,
                            &local_y_start_after_transpose,
                            &total_local_size);

    MPI_Allgather(&local_x_start, 1, MPI_INT, pm_local_start,1, MPI_INT, MESH_COMM_WORLD);
    MPI_Allgather(&local_nx, 1, MPI_INT, pm_local_xside,1, MPI_INT, MESH_COMM_WORLD);

    data = (fftw_real*) malloc(sizeof(fftw_real) * total_local_size);
    work = (fftw_real*) malloc(sizeof(fftw_real) * total_local_size);
    cdata = (fftw_complex*) data;
    printf("[MESH%3d] local_x_start = %4d local_nx = %4d pm_nside = %5d\n", rank, local_x_start, local_nx, pm_nside);
}

typedef struct {
    Domain *dp;
    GlobalParam *gp;
} TaskCommParam;

typedef struct {
    int length;
    int index;
    //    real *buff;   // pointer to the buff
    int start[3]; // absoluted displacement
    int nside[3];


} BuffInfo;

void* convolution_gravity(void *param)
{
    int rank;
    int s;

    TaskCommParam* tcp = (TaskCommParam*)param;
    Domain *dp = tcp->dp;
    GlobalParam *gp = tcp->gp;
    real boxsize = gp->BoxSize;
    real mass;
    real delta = boxsize/pm_nside;
    real invdelta = 1.0/delta;
    Body *part = dp->Part;
    int num_pm_socket = gp->NumSocket;

    MPI_Comm_rank(MESH_COMM_WORLD, &rank);
    int n, d, npart;
    for (d=0; d<3; d++) {
        pos_min[d] = boxsize;
        pos_max[d] = 0.0;
    }
    npart = dp->NumPart;
    for (n=0; n<npart; n++) {
        for (d=0; d<3; d++) {
            if (part[n].pos[d]>pos_max[d])
                pos_max[d] = part[n].pos[d];
            if (part[n].pos[d]<pos_min[d])
                pos_min[d] = part[n].pos[d];
        }
    }
    for (d=0; d<3; d++) {
        mesh_min[d] = (int)(pos_min[d]*invdelta);
        mesh_max[d] = (int)(pos_max[d]*invdelta);
        mesh_nside[d] = 1+(mesh_max[d]-mesh_min[d]+1)+1;
    }
//printf(" %d - -%d %d %d -- %d\n", rank, mesh_nside[0], mesh_nside[1], mesh_nside[2],pm_nside);
    int field_length = mesh_nside[0]*mesh_nside[1]*mesh_nside[2];
    pm_real_field = (real*)malloc(sizeof(real)*field_length );
    
 assert(pm_real_field != NULL);
    /* cic */
    for (n=0; n<field_length; n++)
        pm_real_field[n] = 0.0;

//	printf("[%d] %d %d\n", rank, pm_local_start[rank], pm_local_xside[rank]);
/*
    if (0==rank) {
        for (s=0; s<num_pm_socket; s++)
            printf(" %d  c= [ %d , %d]\n", s, pm_local_start[s], pm_local_start[s] + pm_local_xside[s] - 1);
    }
*/
//MPI_Barrier(MPI_COMM_WORLD);
//exit(0);



    int I,J,K, In, Jn, Kn, NX, NY, NZ;
    real x, y, z, wx, wy, wz, wxn, wyn, wzn;
    NX = mesh_nside[0];
    NY = mesh_nside[1];
    NZ = mesh_nside[2];
real px, py, pz;
real pmx, pmy, pmz;

pmx = mesh_min[0]*delta;
pmy = mesh_min[1]*delta;
pmz = mesh_min[2]*delta;


    for (n=0; n<npart; n++) {
	px = part[n].pos[0] - pmx;
	py = part[n].pos[1] - pmy;
	pz = part[n].pos[2] - pmz;

	I = (int)(px*invdelta) + 1;
	J = (int)(py*invdelta) + 1;
	K = (int)(pz*invdelta) + 1;
//        I = (int)(part[n].pos[0]*invdelta) - mesh_min[0] + 1;
//        J = (int)(part[n].pos[1]*invdelta) - mesh_min[1] + 1;
//        K = (int)(part[n].pos[2]*invdelta) - mesh_min[2] + 1;
        x = (I + 0.5)*delta;
        y = (J + 0.5)*delta;
        z = (K + 0.5)*delta;

//        wxn = (part[n].pos[0] - x)*invdelta;
//        wyn = (part[n].pos[1] - y)*invdelta;
//        wzn = (part[n].pos[2] - z)*invdelta;

        wxn = (px - x)*invdelta;
        wyn = (py - y)*invdelta;
        wzn = (pz - z)*invdelta;

        if (wxn<0.0) {
            In = I - 1;
            wxn = -wxn;
        }
        else
            In = I + 1;

        wx = 1.0 - wxn;

        if (wyn<0.0) {
            Jn = I - 1;
            wyn = -wyn;
        }
        else
            Jn = J + 1;

        wy = 1.0 - wyn;

        if (wzn<0.0) {
            Kn = K - 1;
            wzn = -wzn;
        }
        else
            Kn = K + 1;

        wz = 1.0 - wzn;

        pm_real_field[(I *NY + J )*NZ + K ] += wx *wy *wz ;
        pm_real_field[(In*NY + J )*NZ + K ] += wxn*wy *wz ;
        pm_real_field[(I *NY + Jn)*NZ + K ] += wx *wyn*wz ;
        pm_real_field[(I *NY + J )*NZ + Kn] += wx *wy *wzn;
        pm_real_field[(In*NY + Jn)*NZ + K ] += wxn*wyn*wz ;
        pm_real_field[(In*NY + J )*NZ + Kn] += wxn*wy *wzn;
        pm_real_field[(I *NY + Jn)*NZ + Kn] += wx *wyn*wzn;
        pm_real_field[(In*NY + Jn)*NZ + Kn] += wxn*wyn*wzn;

 //       assert((I *NY + J )*NZ + K < field_length);
  //      assert((In*NY + Jn)*NZ + Kn< field_length);
    }

    mass = 1.0;
    real density = mass ;
    for (n=0; n<field_length; n++)
        pm_real_field[n] *= density;


//MPI_Barrier(MPI_COMM_WORLD);
//exit(0);

    BuffInfo *sendinfo; // for tree domain;
    BuffInfo *recvinfo; // for pm sockets;

    sendinfo = (BuffInfo*)malloc(sizeof(BuffInfo)*num_pm_socket);
    recvinfo = (BuffInfo*)malloc(sizeof(BuffInfo)*num_pm_socket);

    for (s=0; s<num_pm_socket; s++) {
        sendinfo[s].length = 0;
        sendinfo[s].start[0] = 0;
        sendinfo[s].start[1] = mesh_min[1];
        sendinfo[s].start[2] = mesh_min[2];
        sendinfo[s].nside[0] = 0;
        sendinfo[s].nside[1] = mesh_nside[1];
        sendinfo[s].nside[2] = mesh_nside[2];
        sendinfo[s].index = -1;

        recvinfo[s].length = 0;
        recvinfo[s].start[0] = 0;
        recvinfo[s].start[1] = 0;
        recvinfo[s].start[2] = 0;
        recvinfo[s].nside[0] = 0;
        recvinfo[s].nside[1] = 0;
        recvinfo[s].nside[2] = 0;
        recvinfo[s].index = -1;
    }

    int ab_xstart, ab_xend, s_xstart, s_xend, cut, width;
    ab_xstart = mesh_min[0] - 1;
    ab_xend   = mesh_max[0] + 1;
    cut = 0;
//MPI_Barrier(MPI_COMM_WORLD);
//exit(0);

//    printf(" [%d] ab_xstart =%d ab_xend = %d mesh_nside = %d\n", rank, ab_xstart, ab_xend, mesh_nside[0] );
    if (ab_xstart == -1) {
        s=num_pm_socket - 1;
        width = 1;
        sendinfo[s].start[0] = pm_nside -1 ;
        sendinfo[s].nside[0] = width;
        sendinfo[s].index = 0;
        sendinfo[s].length = width*mesh_nside[1]*mesh_nside[2];
        cut += width;
  //      printf(" [%d] s=%d, cut=%d ab_xstart =%d \n", rank,s, cut,ab_xstart );
    }
    for (s=0; s<num_pm_socket; s++) {
        s_xstart = pm_local_start[s];
        s_xend   = pm_local_start[s] + pm_local_xside[s] - 1;

        ab_xstart = cut + mesh_min[0] - 1;

        if (cut == mesh_nside[0] )
            break;
        if ( ab_xstart > s_xend)
            continue;
        if ( s_xstart<= ab_xstart && ab_xstart <= s_xend) {
            width = (s_xend > ab_xend ? ab_xend : s_xend ) - ab_xstart + 1;
            sendinfo[s].start[0] = ab_xstart;
            sendinfo[s].nside[0] = width;
            sendinfo[s].index = cut*mesh_nside[1]*mesh_nside[2];
            sendinfo[s].length = width*mesh_nside[1]*mesh_nside[2];
            cut += width;
            continue;
        }

    }
    if (ab_xend == pm_nside) {
        s=0;
        width = 1;
        sendinfo[s].start[0] = 0 ;
        sendinfo[s].nside[0] = width;
        sendinfo[s].index = cut*mesh_nside[1]*mesh_nside[2];
        sendinfo[s].length = width*mesh_nside[1]*mesh_nside[2];
        cut += width;
//        printf(" [%d] s=%d, cut=%d ab_xstart =%d \n", rank,s, cut,ab_xstart );
    }

    assert(cut == mesh_nside[0] );
    ///////////////////

//MPI_Barrier(MPI_COMM_WORLD);
//exit(0);

    BuffInfo template;
    MPI_Datatype mpi_buffinfo_type;
    MPI_Datatype type[4] = { MPI_INT, MPI_INT, MPI_INT, MPI_INT };
    int blocklen[4] = { 1, 1, 3, 3 };
    MPI_Aint disp[4], start_address, address;

    MPI_Get_address(&template, &start_address);

    MPI_Get_address(&template.length, &address);
    disp[0] = address - start_address;

    MPI_Get_address(&template.index, &address);
    disp[1] = address - start_address;

    MPI_Get_address(&template.start, &address);
    disp[2] = address - start_address;

    MPI_Get_address(&template.nside, &address);
    disp[3] = address - start_address;

    MPI_Type_create_struct(4, blocklen, disp, type, &mpi_buffinfo_type);
    MPI_Type_commit(&mpi_buffinfo_type);
//////////////////////////

//    for (s=0; s<num_pm_socket; s++) {
//        printf(" send [%d] s =%d length = %d\n", rank, s, sendinfo[s].length);
//    }
    
    MPI_Alltoall(sendinfo, 1, mpi_buffinfo_type, recvinfo, 1, mpi_buffinfo_type, MESH_COMM_WORLD);
    
 //   for (s=0; s<num_pm_socket; s++) {
 //       printf(" recv [%d] s =%d length = %d\n", rank, s, recvinfo[s].length);
 //   }

    recv_block = (real**)malloc(sizeof(real*)*num_pm_socket);
    for (s=0; s<num_pm_socket; s++) {
        if (recvinfo[s].length > 0)
            recv_block[s] = (real*)malloc(sizeof(real)*recvinfo[s].length);
        else
            recv_block[s] = NULL;
    }
    MPI_Status *status;
    status = (MPI_Status*)malloc(sizeof(MPI_Status)*num_pm_socket);
    int send, recv;
    for (s=0; s<num_pm_socket; s++) {
        send = (rank+s)%num_pm_socket;
        recv = (rank-s +num_pm_socket)%num_pm_socket;

        MPI_Sendrecv(&(pm_real_field[sendinfo[send].index]), 
                     sendinfo[send].length, MPI_INT, send, 1, 
                     recv_block[recv], recvinfo[recv].length, MPI_INT, recv, 1, 
                     MESH_COMM_WORLD, &status[s]);
    }    
    int LENGTH,DISP, X, Y, Z, XST, YST, ZST, Xab, Yab, Zab;
    real test =0.0;
    for (s=0; s<num_pm_socket; s++) {
        
        if (recvinfo[s].length > 0) {
            NX = recvinfo[s].nside[0];
            NY = recvinfo[s].nside[1];
            NZ = recvinfo[s].nside[2];
            LENGTH = recvinfo[s].length;
            DISP = pm_local_start[rank];
            XST = recvinfo[s].start[0];
            YST = recvinfo[s].start[1];
            ZST = recvinfo[s].start[2];
      //      printf(" [%d] << s =%d << x = %d, %d,%d %d %d \n",rank,  s, XST,DISP, NX, NY, NZ);
            for (X=0; X<NX; X++) {
                I = X + XST - DISP;
                for (Y=0; Y<NY; Y++) {
                    J = (Y + YST - 1 + NY)%NY;
                    for (Z=0; Z<NZ; Z++) {
                        K = (Z + ZST - 1 + NZ)%NZ;
                        
                        data[(I *pm_nside + J )*pm_nside + K ] 
                                  += recv_block[s][(X*NY+Y)*NZ+Z] ;
                                     
                    }
                }
                
            }
        }
    }


    int local_nx, local_x_start,local_ny_after_transpose,local_y_start_after_transpose,total_local_size;
    rfftwnd_mpi_local_sizes(forward_plan, &local_nx, &local_x_start,
                            &local_ny_after_transpose,
                            &local_y_start_after_transpose,
                            &total_local_size);


    rfftwnd_mpi(forward_plan, 1, data, work, FFTW_TRANSPOSED_ORDER);

    cdata = (fftw_complex*) data;
//////////////////////////////////////////
    real Kx, Ky, Kz, K2;
    real greenfunc, trun;
    
    trun = 2*acos(-1.0)*(1.2*delta)/boxsize;
    trun = (trun*trun);
    int  kx, ky, kz;
    for (Y = 0; Y < local_ny_after_transpose; ++Y) {
        ky = Y + local_y_start_after_transpose;
        if (ky > pm_nside/2)
            ky = pm_nside - ky;
        Ky = (real) ky;
        for (X = 0; X < pm_nside; ++X) {
            kx = X;
            if (kx>pm_nside/2)
                kx = pm_nside - kx;
            Kx = (real) kx;
            for (Z = 0; Z < (pm_nside/2+1); ++Z) {
                kz = Z;
                Kx = (real)kz;
                K2 = Kx*Kx + Ky*Ky + Kz*Kz;

                 greenfunc = -exp(-K2 * trun) / K2;
                cdata[(Y*pm_nside + X) * (pm_nside/2+1) + Z].im *= 2.0;
                
            }
        }
    }
/////////////////////////////////////////////
    rfftwnd_mpi(inverse_plan, 1, data, work, FFTW_TRANSPOSED_ORDER);

   for (s=0; s<num_pm_socket; s++) {
        
        if (recvinfo[s].length > 0) {
            NX = recvinfo[s].nside[0];
            NY = recvinfo[s].nside[1];
            NZ = recvinfo[s].nside[2];
            LENGTH = recvinfo[s].length;
            DISP = pm_local_start[rank];
            XST = recvinfo[s].start[0];
            YST = recvinfo[s].start[1];
            ZST = recvinfo[s].start[2];
            
            for (X=0; X<NX; X++) {
                I = X + XST - DISP;
                for (Y=0; Y<NY; Y++) {
                    J = (Y + YST - 1 + NY)%NY;
                    for (Z=0; Z<NZ; Z++) {
                        K = (Z + ZST - 1 + NZ)%NZ;
                        
                        // save back the data into recv_block
                        recv_block[s][(X*NY+Y)*NZ+Z] = data[(I *pm_nside + J )*pm_nside + K ];             
                    }
                }
                
            }
        }
    }    

    for (s=0; s<num_pm_socket; s++) {
        send = (rank+s)%num_pm_socket;
        recv = (rank-s +num_pm_socket)%num_pm_socket;
        MPI_Sendrecv(recv_block[recv], recvinfo[recv].length, MPI_INT, recv, 1, 
                     &(pm_real_field[sendinfo[send].index]), 
                     sendinfo[send].length, MPI_INT, send, 1, 
                     MESH_COMM_WORLD, &status[s]);
   
    }  
    /* pm_real_field is the potential we need */

    free(sendinfo);
    free(recvinfo);
    for (s=0; s<num_pm_socket; s++) 
        if (recv_block[s] != NULL)
            free(recv_block[s]);   
    free(recv_block);
    MPI_Type_free(&mpi_buffinfo_type);
        
    printf(" [%d] PM done \n", rank);
    
    free(work);
    free(data);

}



