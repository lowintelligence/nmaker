#include "partmesh.h"
#include "global.h"
#include <math.h>
#include <assert.h>

#ifndef FFTW3_LIB
    #include <fftw.h>
    #include <rfftw_mpi.h>
#else
    #include <fftw3.h>
    #include <fftw3-mpi.h>
#endif /* FFTW3_LIB */


void* partmesh_acc(void *param)
{
    printf("  partmesh start\n\n");
    /*
            pthread_t pm_thread;
            pthread_create(&pm_thread, NULL, convolution_gravity, param);
            pthread_join(pm_thread, 0);

            pthread_create(&pm_thread, NULL, pm_acceleration, param);
            pthread_join(pm_thread, 0);
    */
    convolution_gravity(param);
    pm_acceleration(param);


    printf("  partmesh end\n\n");
}

/*  create struct Partmesh compatiable with pthread calling,
    containing a mesh to determine the asignment of particles  */
Partmesh *create_partmesh(Constants *constants, Status *status, System *system)
{
    Partmesh *pm = xmalloc(sizeof(Partmesh), 4001);

    pm->part = system->part;
    pm->npart= system->num_part;
    pm->boxsize = constants->BOX_SIZE;
    pm->r_split = constants->SPLIT_SCALE;
    pm->pm_nside= constants->PM_NUM_SIDE;
    pm->mass    = constants->PART_MASS;
    pm->grav_const = constants->GRAV_CONST;

    Body *part = pm->part;
    double pos_min[3];
    double pos_max[3];
    double boxsize = pm->boxsize;
    int pm_nside = pm->pm_nside;
    pm->mesh_bin  = (double*)xmalloc(sizeof(double)*(pm_nside+1), 4002);
    pm->a_time = 1.0/(1.0+status->redshift);

    int mesh_min[3];
    int mesh_max[3];
    int mesh_nside[3];

    double invdelta = (pm_nside)/boxsize;
    double delta = boxsize/pm_nside;
    int s, d, i;
    for (d=0; d<3; d++) {
        pos_min[d] = boxsize;
        pos_max[d] = 0.0;
    }
    for (i=0; i<pm_nside; i++)
        pm->mesh_bin[i] = (double)((double)(i*delta));
    pm->mesh_bin[pm_nside] = boxsize;
    /*  generate array of mesh_bin to avoid asigning particle
        into an overrated grid index, induced by round-off error   */


    Int n, npart;
    double pos[3];
    npart = pm->npart;
    for (n=0; n<npart; n++) {
        for (d=0; d<3; d++) {
            pos[d] = (part[n].pos[d]);

            if (pos[d]>=pos_max[d])
                pos_max[d] = pos[d];
            if (pos[d]<=pos_min[d])
                pos_min[d] = pos[d];
        }
    }

    for (d=0; d<3; d++) {
        mesh_min[d] = (int)(pos_min[d]*invdelta);
        mesh_max[d] = (int)(pos_max[d]*invdelta);
        if (pos_min[d] < pm->mesh_bin[mesh_min[d]])
            mesh_min[d] -= 1;
        if (pos_max[d] < pm->mesh_bin[mesh_max[d]])
            mesh_max[d] -= 1;

        if (mesh_min[d] == (pm->pm_nside))
            mesh_min[d] = (pm->pm_nside)-1;

        if (mesh_max[d] == (pm->pm_nside))
            mesh_max[d] = (pm->pm_nside) - 1;

        assert(mesh_max[d] < pm->pm_nside);

        mesh_nside[d] = 1+(mesh_max[d]-mesh_min[d]+1)+1;
    }

    pm->potential = (Real*)xmalloc(sizeof(Real)*(mesh_nside[0]+4)*(mesh_nside[1]+4)*(mesh_nside[2]+4), 4003);

    pm->mesh_min[0] = mesh_min[0];
    pm->mesh_min[1] = mesh_min[1];
    pm->mesh_min[2] = mesh_min[2];
    pm->mesh_nside[0] = mesh_nside[0];
    pm->mesh_nside[1] = mesh_nside[1];
    pm->mesh_nside[2] = mesh_nside[2];
    pm->mesh_max[0] = mesh_max[0];
    pm->mesh_max[1] = mesh_max[1];
    pm->mesh_max[2] = mesh_max[2];

    return pm;
}

void free_partmesh(Partmesh *pm)
{
    if ( NULL != pm) {
        if ( NULL != (pm->potential) )
            free(pm->potential);  /* 4003 */

        if ( NULL != (pm->mesh_bin) )
            free(pm->mesh_bin);  /* 4002 */
        free(pm);  /* 4001 */
    }
}


typedef struct {
    int length;
    int index;
    int start[3]; /* global displacement */
    int nside[3];
} BuffInfo;


static int snap_index = 0;
static int snap_count = 0;

void* assign_particles(void *param)
{
    Partmesh *pm = (Partmesh*)param;

    int pm_nside = pm->pm_nside;
    int pm_nside_pad  = (pm_nside/2 + 1 )*2;
    int size_data;
    int local_nx, local_x_start, local_ny_after_transpose,
        local_y_start_after_transpose, total_local_size;
    long local_n0, local_n0_start, local_n1, local_n1_start;
    int *pm_local_start, *pm_local_xside;
    double *bin = pm->mesh_bin;

    double boxsize = (double)(pm->boxsize);
    double split   = (double)(pm->r_split);
    double mass  = (double)(pm->mass);
    double delta = boxsize/pm_nside;
    double invdelta = 1.0/delta;
    double gravconst = pm->grav_const;
    Body *part = pm->part;
    const Int npart = pm->npart;

    double pos_min[3];
    double pos_max[3];
    int mesh_min[3];
    int mesh_max[3];
    int mesh_nside[3];

    int s, d, i;
    Int n;

    mesh_min[0] = pm->mesh_min[0];
    mesh_min[1] = pm->mesh_min[1];
    mesh_min[2] = pm->mesh_min[2];

    mesh_max[0] = pm->mesh_max[0];
    mesh_max[1] = pm->mesh_max[1];
    mesh_max[2] = pm->mesh_max[2];

    mesh_nside[0] = pm->mesh_nside[0];
    mesh_nside[1] = pm->mesh_nside[1];
    mesh_nside[2] = pm->mesh_nside[2];

    int field_length = mesh_nside[0]*mesh_nside[1]*mesh_nside[2];
    pm->local_field = (double*)xmalloc(sizeof(double)*field_length,4010);

    for (i=0; i<field_length; i++)
        pm->local_field[i] = 0.0;

    int I,J,K, In, Jn, Kn, NX, NY, NZ;
    double px, py, pz, pmx, pmy, pmz;
    double x, y, z, wx, wy, wz, wxn, wyn, wzn;
    double pos[3];
    NX = mesh_nside[0];
    NY = mesh_nside[1];
    NZ = mesh_nside[2];

    pmx = mesh_min[0]*delta;
    pmy = mesh_min[1]*delta;
    pmz = mesh_min[2]*delta;

    /* asign particle into mesh by CIC */
    for (n=0; n<npart; n++) {
        pos[0] = part[n].pos[0];
        pos[1] = part[n].pos[1];
        pos[2] = part[n].pos[2];

        px = pos[0] - pmx;
        py = pos[1] - pmy;
        pz = pos[2] - pmz;

        I = (int)(pos[0]*invdelta);
        if (pos[0] < bin[I])
            I -= 1;

        J = (int)(pos[1]*invdelta);
        if (pos[1] < bin[J] )
            J -= 1;

        K = (int)(pos[2]*invdelta);
        if (pos[2] < bin[K])
            K -= 1;

        I += 1 - mesh_min[0];
        J += 1 - mesh_min[1];
        K += 1 - mesh_min[2];

        x = (I - 0.5)*delta;
        y = (J - 0.5)*delta;
        z = (K - 0.5)*delta;

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
            Jn = J - 1;
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

        pm->local_field[(I *NY+J )*NZ+K ] += part[n].mass * wx *wy *wz ;
        pm->local_field[(In*NY+J )*NZ+K ] += part[n].mass * wxn*wy *wz ;
        pm->local_field[(I *NY+J )*NZ+Kn] += part[n].mass * wx *wy *wzn;
        pm->local_field[(In*NY+Jn)*NZ+K ] += part[n].mass * wxn*wyn*wz ;
        pm->local_field[(In*NY+J )*NZ+Kn] += part[n].mass * wxn*wy *wzn;
        pm->local_field[(I *NY+Jn)*NZ+Kn] += part[n].mass * wx *wyn*wzn;
        pm->local_field[(I *NY+Jn)*NZ+K ] += part[n].mass * wx *wyn*wz ;
        pm->local_field[(In*NY+Jn)*NZ+Kn] += part[n].mass * wxn*wyn*wzn;
    }
    double vol = delta*delta*delta;
    double density = 1.0/vol;

    for (n=0; n<field_length; n++)
        pm->local_field[n] *= density;
    /* asign particle into mesh by CIC */
}

void* convolution_gravity(void *param)
{
    Partmesh *pm = (Partmesh*)param;

    int rank, Nproc;
    int pm_nside = pm->pm_nside;
    int pm_nside_pad  = (pm_nside/2 + 1 )*2;
    int size_data;
    int local_nx, local_x_start, local_ny_after_transpose,
        local_y_start_after_transpose, total_local_size;
    long local_n0, local_n0_start, local_n1, local_n1_start;
    int *pm_local_start, *pm_local_xside;
    double *bin = pm->mesh_bin;

    MPI_Comm_rank(PM_COMM_WORLD, &rank);
    MPI_Comm_size(PM_COMM_WORLD, &Nproc);

    pm_local_start = (int*)xmalloc(sizeof(int)*Nproc, 4004);
    pm_local_xside = (int*)xmalloc(sizeof(int)*Nproc, 4005);

#ifndef FFTW3_LIB
    fftw_real *data, *work;
    fftw_complex *cdata;
    rfftwnd_mpi_plan forward_plan,  inverse_plan;


    forward_plan = rfftw3d_mpi_create_plan(PM_COMM_WORLD,
                                           pm_nside, pm_nside, pm_nside,
                                           FFTW_REAL_TO_COMPLEX,
                                           FFTW_ESTIMATE);

    inverse_plan = rfftw3d_mpi_create_plan(PM_COMM_WORLD,
                                           pm_nside, pm_nside, pm_nside,
                                           FFTW_COMPLEX_TO_REAL,
                                           FFTW_ESTIMATE);

    rfftwnd_mpi_local_sizes(forward_plan, &local_nx, &local_x_start,
                            &local_ny_after_transpose,
                            &local_y_start_after_transpose,
                            &total_local_size);

    data = (fftw_real*) xmalloc(sizeof(fftw_real) * total_local_size, 4006);
    work = (fftw_real*) xmalloc(sizeof(fftw_real) * total_local_size, 4007);

    cdata = (fftw_complex*)data;

    size_data = sizeof(fftw_real);
#else
    double *data;
    fftw_complex *cdata;
    fftw_plan forward_plan,  inverse_plan;

//    fftw_mpi_local_size_3d_transposed((ptrdiff_t)pm_nside, (ptrdiff_t)pm_nside, (ptrdiff_t)pm_nside,
    fftw_mpi_local_size_3d((ptrdiff_t)pm_nside, (ptrdiff_t)pm_nside, (ptrdiff_t)pm_nside,
                           PM_COMM_WORLD,
                           (ptrdiff_t*)&local_n0, (ptrdiff_t*)&local_n0_start);
//                                      (ptrdiff_t*)&local_n1, (ptrdiff_t*)&local_n1_start);

    DBG_INFO(DBG_TMP, "%ld, %ld, PMWORLD=%ld, WORLD=%ld.\n", forward_plan, inverse_plan, PM_COMM_WORLD, MPI_COMM_WORLD);
    local_nx = (int)local_n0;
    local_x_start = (int)local_n0_start;
    local_ny_after_transpose = (int)local_n0;
    local_y_start_after_transpose = (int)local_n0_start;

    total_local_size = local_nx * pm_nside * 2 *(pm_nside/2+1);

    DBG_INFO(DBG_TMP, "%ld, %ld, %ld, %ld, %d, %ld.\n", local_n0, local_n0_start, local_n1, local_n1_start, local_nx, total_local_size);
    data = (double*) fftw_malloc(sizeof(double) * total_local_size);  /* 4008 */
    if ( !data ) {
        printf(" ERROR! : fail to alloc memory for data 4008 \n");
        system_exit(ERR_MEMORY);
    }
    cdata=(fftw_complex*)fftw_malloc(sizeof(double)*total_local_size);/* 4009 */
    /*    cdata = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)
                * local_nx * pm_nside *(pm_nside/2+1)); */

    if ( !cdata ) {
        printf(" ERROR! : fail to alloc memory for data 4009 \n");
        system_exit(ERR_MEMORY);
    }

    forward_plan = fftw_mpi_plan_dft_r2c_3d(pm_nside, pm_nside, pm_nside,
                                            data, cdata,
                                            PM_COMM_WORLD, FFTW_ESTIMATE);

    inverse_plan = fftw_mpi_plan_dft_c2r_3d(pm_nside, pm_nside, pm_nside,
                                            cdata, data,
                                            PM_COMM_WORLD, FFTW_ESTIMATE);

    DBG_INFO(DBG_TMP, "%ld, %ld.\n", forward_plan, inverse_plan);

    size_data = sizeof(double);
#endif /* FFTW3_LIB */

    MPI_Allgather(&local_x_start, 1, MPI_INT, pm_local_start,1, MPI_INT, PM_COMM_WORLD);
    MPI_Allgather(&local_nx, 1, MPI_INT, pm_local_xside,1, MPI_INT, PM_COMM_WORLD);

    double boxsize = (double)(pm->boxsize);
    double split   = (double)(pm->r_split);
    double mass  = (double)(pm->mass);
    double delta = boxsize/pm_nside;
    double invdelta = 1.0/delta;
    double gravconst = pm->grav_const;
    Body *part = pm->part;
    const Int npart = pm->npart;

    double pos_min[3];
    double pos_max[3];
    int mesh_min[3];
    int mesh_max[3];
    int mesh_nside[3];

    int s, d, i;
    Int n;

    mesh_min[0] = pm->mesh_min[0];
    mesh_min[1] = pm->mesh_min[1];
    mesh_min[2] = pm->mesh_min[2];

    mesh_max[0] = pm->mesh_max[0];
    mesh_max[1] = pm->mesh_max[1];
    mesh_max[2] = pm->mesh_max[2];

    mesh_nside[0] = pm->mesh_nside[0];
    mesh_nside[1] = pm->mesh_nside[1];
    mesh_nside[2] = pm->mesh_nside[2];

    int field_length = mesh_nside[0]*mesh_nside[1]*mesh_nside[2];

    double *pm_local_field = pm->local_field;
    int I,J,K;

    double *slab_send, *slab_recv;
    int slab_size, slab_size_max;
    int Ig, Jg, Kg;
    int conj;
    int pm_comm_mask;
    int pm_comm_loop, loop;
    int index, first_index, last_index;
    MPI_Status pm_status;
    slab_size_max = 0;
    for (s=0; s<Nproc; s++) {
        if (slab_size_max <=pm_local_xside[s]) {
            slab_size_max = pm_local_xside[s];
        }
    }
    slab_size_max *= (pm_nside_pad*pm_nside);
    slab_size = pm_nside * pm_nside_pad * pm_local_xside[rank];
    slab_send = (double*)xmalloc(sizeof(double)*slab_size_max, 4011);
    slab_recv = (double*)xmalloc(sizeof(double)*slab_size_max, 4012);

    pm_comm_loop = 1;
    while (pm_comm_loop < Nproc)
        pm_comm_loop<<=1;

    MPI_Barrier(PM_COMM_WORLD);
    int slab_size_rank, slab_size_conj;

    /* send distributed density to pm slab, according to pm rank */
    for (i=0; i<slab_size; i++)
        data[i] = 0.0;
    for (i=0; i<slab_size_max; i++) {
        slab_send[i] = 0.0;
    }
    first_index = pm_local_start[rank] * pm_nside_pad * pm_nside;
    for (I=0; I<mesh_nside[0]; I++) {
        Ig = ( I + mesh_min[0] - 1 + pm_nside)%pm_nside;
        for (J=0; J<mesh_nside[1]; J++) {
            Jg = ( J + mesh_min[1] - 1 + pm_nside)%pm_nside;
            for (K=0; K<mesh_nside[2]; K++) {
                Kg = ( K + mesh_min[2] - 1 + pm_nside)%pm_nside;
                index = (Ig * pm_nside + Jg ) * pm_nside_pad + Kg;
                index -= first_index;
                if (0<=index && index<slab_size)
                    slab_send[index] += pm_local_field[(I*mesh_nside[1]+J)*mesh_nside[2]+K];
            }
        }
    }

    for (i=0; i<slab_size; i++)
        data[i] += slab_send[i];

//	{
//		int nx, ny, nz;
//		double value, data_max = 0.0;
//		double data_min = data[3*pm_nside*pm_nside_pad];
//        for (nx = 3, ny = 0; ny<pm_nside; ny++) {
//            for (nz=0; nz<pm_nside; nz++) {
//                if( data[(nx*pm_nside + ny)*pm_nside_pad + nz] > data_max)
//                    data_max = data[(nx*pm_nside + ny)*pm_nside_pad + nz];
//                if( data[(nx*pm_nside + ny)*pm_nside_pad + nz] < data_min)
//                    data_min = data[(nx*pm_nside + ny)*pm_nside_pad + nz];
//            }
//        }
//		printf("[Den0 %d] %g %g\n", rank, data_max, data_min);
//	}

    for (loop=1, pm_comm_mask=1; loop<pm_comm_loop; loop++, pm_comm_mask++)
    {
        conj = rank^pm_comm_mask;
        if (conj < Nproc) {

            slab_size_rank = pm_nside * pm_nside_pad * pm_local_xside[rank];
            slab_size_conj = pm_nside * pm_nside_pad * pm_local_xside[conj];


            for (i=0; i<slab_size_max; i++)
                slab_send[i] = slab_recv[i] = 0.0;
            first_index = pm_local_start[conj] * pm_nside_pad * pm_nside;

            for (I=0; I<mesh_nside[0]; I++) {
                Ig = ( I + mesh_min[0] - 1 + pm_nside)%pm_nside;
                for (J=0; J<mesh_nside[1]; J++) {
                    Jg = ( J + mesh_min[1] - 1 + pm_nside)%pm_nside;
                    for (K=0; K<mesh_nside[2]; K++) {
                        Kg = ( K + mesh_min[2] - 1 + pm_nside)%pm_nside;
                        index = (Ig * pm_nside + Jg ) * pm_nside_pad + Kg;
                        index -= first_index;
                        if (0<=index && index<slab_size_conj)
                            slab_send[index] += pm_local_field[(I*mesh_nside[1]+J)*mesh_nside[2]+K];
                    }
                }
            }


            if (conj < rank)
                MPI_Send(slab_send,  slab_size_max*size_data, MPI_CHAR, conj, 1, PM_COMM_WORLD);
            else
                MPI_Recv(slab_recv,  slab_size_max*size_data, MPI_CHAR, conj, 1, PM_COMM_WORLD, &pm_status);

            if (conj < rank)
                MPI_Recv(slab_recv,  slab_size_max*size_data, MPI_CHAR, conj, 2, PM_COMM_WORLD, &pm_status);
            else
                MPI_Send(slab_send,  slab_size_max*size_data, MPI_CHAR, conj, 2, PM_COMM_WORLD);

            for (i=0; i<slab_size_rank; i++)
                data[i] += slab_recv[i];

        }
        MPI_Barrier(PM_COMM_WORLD);
    }
    /* send distributed density to pm slab, according to pm rank */

    free(pm_local_field);  /* 4010 */


#ifdef MAP_PPM
    char fname[50];
    FILE *fd;
    int nx, nz, ny;

    if (0==rank && 0 == snap_count % 1 ) {
        printf("   output  snap %d %d\n", snap_index, snap_count);
        sprintf(fname, "density_%d.ppm", snap_index);
        FILE *fp = fopen(fname, "w");
        (void) fprintf(fp, "P6\n%d %d\n255\n", pm_nside, pm_nside);

        double a_t = pm->a_time;
        double value, data_max = 0.0;
        for (nx = 3, ny = 0; ny<pm_nside; ny++) {
            for (nz=0; nz<pm_nside; nz++) {
                if( data[(nx*pm_nside + ny)*pm_nside_pad + nz] >= data_max)
                    data_max = data[(nx*pm_nside + ny)*pm_nside_pad + nz];
            }
        }
        //   data_max += 1.0;

        //    int color_max = (int)(sqrt(a_t)*256);
        int color_max = (int)(256);
        int color_value =0 ;
        for (ny = 0; ny<pm_nside; ny++) {
            for (nz=0; nz<pm_nside; nz++) {
                static unsigned char color[3];
                double weighted_value;
                weighted_value = value = 0.0;
                nx = 3;
                for (nx=0; nx<10; nx++)
                    weighted_value += data[(nx*pm_nside + ny)*pm_nside_pad + nz];
                weighted_value /= 20.0;
                value += weighted_value;

                for (nx=10; nx<20; nx++)
                    weighted_value += data[(nx*pm_nside + ny)*pm_nside_pad + nz];
                weighted_value /= 40.0;
                value += weighted_value;

                for (nx=20; nx<40; nx++)
                    weighted_value += data[(nx*pm_nside + ny)*pm_nside_pad + nz];
                weighted_value /= 100.0;
                value += weighted_value;

                value /= (1+0.5+0.2);
                //    value = (int)((double)color_max*log(value)/log(data_max));
                //    value += 1.0;
                color_value = (int)((double)color_max* ( (log(value) - log(data_max*0.001) )/ (log(data_max) - log(data_max*0.001) ) ) );
                if (color_value < 0)
                    color_value = 0;

                if (color_value > 255)
                    color_value = 255;

                //     printf("v[%e]=%f m[%e]=%f color = %d\n", value, log10(value), data_max, log10(data_max), color_value ) ;

                color[0] = 0;// (unsigned char)( color_value);
                color[1] = 0;// (unsigned char)( color_value);
                color[2] = (unsigned char)( color_value);

                (void) fwrite(color, 1, 3, fp);
            }
        }
        fclose(fp);
        snap_index ++;
    }
    snap_count ++;
#endif

    /*  convolution */
    MPI_Barrier(PM_COMM_WORLD);
//	{
//		int nx, ny, nz;
//		double value, data_max = 0.0;
//		double data_min = data[3*pm_nside*pm_nside_pad];
//        for (nx = 3, ny = 0; ny<pm_nside; ny++) {
//            for (nz=0; nz<pm_nside; nz++) {
//                if( data[(nx*pm_nside + ny)*pm_nside_pad + nz] > data_max)
//                    data_max = data[(nx*pm_nside + ny)*pm_nside_pad + nz];
//                if( data[(nx*pm_nside + ny)*pm_nside_pad + nz] < data_min)
//                    data_min = data[(nx*pm_nside + ny)*pm_nside_pad + nz];
//            }
//        }
//		printf("[Den1 %d] %g %g\n", rank, data_max, data_min);
//	}


#ifndef FFTW3_LIB
    rfftwnd_mpi(forward_plan, 1, data, work, FFTW_TRANSPOSED_ORDER);
    cdata = (fftw_complex*) data;
#else
    fftw_mpi_execute_dft_r2c(forward_plan, data, cdata);
#endif /* FFTW3_LIB */

    int pm_nside_h_pad = (pm_nside/2+1);
    int pm_nside_h = pm_nside/2;
    double Kx, Ky, Kz, K2;
    double greenfunc, scale2, m4piG;
    double sincx, sincy, sincz, gridgf, pi_nside;
    scale2 = 2.0 * M_PI * split/boxsize;
    scale2 = (scale2*scale2);

//    scale2 = 0.0;

    double pref = - gravconst * (delta*delta)/(M_PI * pm_nside);

    int  X, Y, Z, kx, ky, kz;

    pi_nside = (double) M_PI / pm_nside;

#ifndef FFTW3_LIB
    index = 0;
    double Ky2, Kxy2;
    int idx1, idx2;
    for (Y = 0; Y < local_ny_after_transpose; ++Y) {
        ky = (Y + local_y_start_after_transpose);
        if (ky > pm_nside_h)
            ky =  ky - pm_nside;
        Ky = (double) ky;
        Ky2 = Ky*Ky;

        sincy = ky * pi_nside;
        sincy = sin(sincy)/sincy;
        if (ky == 0)
            sincy = 1.0;

        for (Z = 0; Z <pm_nside_h_pad; ++Z) {
            Kz = (double)Z;
            K2 = Kz*Kz + Ky2;

            sincz = Z * pi_nside;
            sincz = sin(sincz)/sincz;
            if (Z == 0)
                sincz = 1.0;
            
            gridgf = sincy *sincz;
            gridgf = 1.0/(gridgf*gridgf*gridgf*gridgf);

            greenfunc = gridgf * pref * EXP(-K2 * scale2)/K2;

            index = (Y*pm_nside ) * pm_nside_h_pad + Z;

            cdata[index].re *= greenfunc;
            cdata[index].im *= greenfunc;
        }

        for (X = 1; X < pm_nside_h; ++X) {
            Kx = (double)X;
            Kxy2 = Kx*Kx + Ky2;

            sincx = Kx * pi_nside;
            sincx = sin(sincx)/sincx;

//#pragma ivdep
            for (Z = 0; Z <pm_nside_h_pad; ++Z) {
                Kz = (double)Z;
                K2 = Kz*Kz + Kxy2;

                sincz = Kz * pi_nside;
                sincz = sin(sincz)/sincz;
                if (Z == 0)
                    sincz = 1.0;

                gridgf = sincx * sincy * sincz;
                gridgf = 1.0/(gridgf*gridgf*gridgf*gridgf);

                greenfunc = gridgf * pref * EXP(-K2 * scale2)/K2;

                idx1 = (Y*pm_nside + X) * pm_nside_h_pad + Z;

                cdata[idx1].re *= greenfunc;
                cdata[idx1].im *= greenfunc;

                idx2 = (Y*pm_nside + pm_nside - X) * pm_nside_h_pad + Z;

                cdata[idx2].re *= greenfunc;
                cdata[idx2].im *= greenfunc;
            }
        }

        Kx = (double)pm_nside_h;
        Kxy2 = Kx*Kx + Ky2;
        sincx = Kx * pi_nside;
        sincx = sin(sincx)/sincx;

        for (Z = 0; Z <pm_nside_h_pad; ++Z) {
            Kz = (double)Z;
            K2 = Kz*Kz + Kxy2;

            sincz = Kz * pi_nside;
            sincz = sin(sincz)/sincz;
            if (Z == 0)
                sincz = 1.0;

            gridgf = sincx * sincy * sincz;
            gridgf = 1.0/(gridgf*gridgf*gridgf*gridgf);

            greenfunc = gridgf * pref * EXP(-K2 * scale2)/K2;

            index = (Y*pm_nside + X) * pm_nside_h_pad + Z;

            cdata[index].re *= greenfunc;
            cdata[index].im *= greenfunc;
        }
    }
    if ( 0== local_y_start_after_transpose) {
        cdata[0].re = 0.0;
        cdata[0].im = 0.0;
    }
//	{
//		int nx, ny, nz;
//		double value, data_max = 0.0;
//		double data_min = data[3*pm_nside*pm_nside_pad];
//        for (nx = 3, ny = 0; ny<pm_nside; ny++) {
//            for (nz=0; nz<pm_nside; nz++) {
//                if( data[(nx*pm_nside + ny)*pm_nside_pad + nz] > data_max)
//                    data_max = data[(nx*pm_nside + ny)*pm_nside_pad + nz];
//                if( data[(nx*pm_nside + ny)*pm_nside_pad + nz] < data_min)
//                    data_min = data[(nx*pm_nside + ny)*pm_nside_pad + nz];
//            }
//        }
//		printf("[Den3 %d] %g %g\n", rank, data_max, data_min);
//	}

    rfftwnd_mpi(inverse_plan, 1, data, work, FFTW_TRANSPOSED_ORDER);

//	{
//		int nx, ny, nz;
//		double value, data_max = 0.0;
//		double data_min = data[3*pm_nside*pm_nside_pad];
//        for (nx = 3, ny = 0; ny<pm_nside; ny++) {
//            for (nz=0; nz<pm_nside; nz++) {
//                if( data[(nx*pm_nside + ny)*pm_nside_pad + nz] > data_max)
//                    data_max = data[(nx*pm_nside + ny)*pm_nside_pad + nz];
//                if( data[(nx*pm_nside + ny)*pm_nside_pad + nz] < data_min)
//                    data_min = data[(nx*pm_nside + ny)*pm_nside_pad + nz];
//            }
//        }
//		printf("[Den4 %d] %g %g\n", rank, data_max, data_min);
//	}

    rfftwnd_mpi_destroy_plan(forward_plan);
    rfftwnd_mpi_destroy_plan(inverse_plan);

    free(work);  /* 4007 */

#else
    index = 0;
    double Ky2, Kx2;
    int idx1, idx2;
    for (Y = 0; Y < local_ny_after_transpose; ++Y) {
        ky = (Y + local_y_start_after_transpose);
        if (ky > pm_nside_h)
            ky =  ky - pm_nside;
        Ky = (double) ky;
        Ky2 = Ky*Ky;

        for (Z = 0; Z <pm_nside_h_pad; ++Z) {
            Kz = (double)Z;
            K2 = Kz*Kz + Ky2;
            greenfunc = pref * EXP(-K2 * scale2)/K2;

            index = (Y*pm_nside ) * pm_nside_h_pad + Z;

            cdata[index][0] *= greenfunc;
            cdata[index][1] *= greenfunc;
        }

        for (X = 1; X < pm_nside_h; ++X) {
            Kx = (double)X;
            Kx2 = Kx*Kx;
#pragma ivdep
            for (Z = 0; Z <pm_nside_h_pad; ++Z) {
                Kz = (double)Z;
                K2 = Kz*Kz + Ky2 + Kx2;

                greenfunc = pref * EXP(-K2 * scale2)/K2;

                idx1 = (Y*pm_nside + X) * pm_nside_h_pad + Z;

                cdata[idx1][0] *= greenfunc;
                cdata[idx1][1] *= greenfunc;

                idx2 = (Y*pm_nside + pm_nside - X) * pm_nside_h_pad + Z;

                cdata[idx2][0] *= greenfunc;
                cdata[idx2][1] *= greenfunc;
            }
        }

        Kx = (double)pm_nside_h;
        Kx2 = Kx*Kx;
        for (Z = 0; Z <pm_nside_h_pad; ++Z) {
            Kz = (double)Z;
            K2 = Kz*Kz + Ky2 + Kx2;
            greenfunc = pref * EXP(-K2 * scale2)/K2;

            index = (Y*pm_nside + X) * pm_nside_h_pad + Z;

            cdata[index][0] *= greenfunc;
            cdata[index][1] *= greenfunc;
        }
    }
    if ( 0== local_y_start_after_transpose) {
        cdata[0][0] = 0.0;
        cdata[0][1] = 0.0;
    }

    fftw_mpi_execute_dft_c2r(inverse_plan, cdata, data);

    fftw_destroy_plan(forward_plan);
    fftw_destroy_plan(inverse_plan);
#endif /* FFTW3_LIB */

    /*  convolution */


    double *pot_pad2;
    int NX_pad, NY_pad, NZ_pad, pot_pad2_size;
    NX_pad = (mesh_nside[0]+4);
    NY_pad = (mesh_nside[1]+4);
    NZ_pad = (mesh_nside[2]+4);
    pot_pad2_size = NX_pad * NY_pad * NZ_pad;

    pot_pad2 = (double*)xmalloc(sizeof(double)*pot_pad2_size, 4013);

    /* pull potential from pm_rank to domain_rank */
    for (i=0; i<pot_pad2_size; i++)
        pot_pad2[i] = 0.0;

    slab_size = pm_nside * pm_nside_pad * pm_local_xside[rank];

    for (i=0; i<slab_size_max; i++)
        slab_send[i] = 0.0;

    for (i=0; i<slab_size; i++)
        slab_send[i] = data[i];

    MPI_Barrier(PM_COMM_WORLD);

    first_index = pm_local_start[rank] * pm_nside_pad * pm_nside;
    for (I=0; I<NX_pad; I++) {
        Ig = ( I + mesh_min[0] - 3 + pm_nside)%pm_nside;
        for (J=0; J<NY_pad; J++) {
            Jg = ( J + mesh_min[1] - 3 + pm_nside)%pm_nside;
            for (K=0; K<NZ_pad; K++) {
                Kg = ( K + mesh_min[2] - 3 + pm_nside)%pm_nside;
                index = (Ig * pm_nside + Jg ) * pm_nside_pad + Kg;
                index -= first_index;
                if (0<=index && index<slab_size)
                    pot_pad2[(I*NY_pad+J)*NZ_pad+K] += slab_send[index];
            }
        }
    }

    for (loop=1, pm_comm_mask=1; loop<pm_comm_loop; loop++, pm_comm_mask++)
    {
        conj = rank^pm_comm_mask;

        if (conj < Nproc) {

            slab_size_rank = pm_nside * pm_nside_pad * pm_local_xside[rank];
            slab_size_conj = pm_nside * pm_nside_pad * pm_local_xside[conj];


            if (conj < rank)
                MPI_Send(slab_send,  slab_size_max*size_data, MPI_CHAR, conj, 1, PM_COMM_WORLD);
            else
                MPI_Recv(slab_recv,  slab_size_max*size_data, MPI_CHAR, conj, 1, PM_COMM_WORLD, &pm_status);

            if (conj < rank)
                MPI_Recv(slab_recv,  slab_size_max*size_data, MPI_CHAR, conj, 2, PM_COMM_WORLD, &pm_status);
            else
                MPI_Send(slab_send,  slab_size_max*size_data, MPI_CHAR, conj, 2, PM_COMM_WORLD);


            first_index = pm_local_start[conj] * pm_nside_pad * pm_nside;
            slab_size = pm_nside * pm_nside_pad * pm_local_xside[conj];
            for (I=0; I<NX_pad; I++) {
                Ig = ( I + mesh_min[0] - 3 + pm_nside)%pm_nside;
                for (J=0; J<NY_pad; J++) {
                    Jg = ( J + mesh_min[1] - 3 + pm_nside)%pm_nside;
                    for (K=0; K<NZ_pad; K++) {
                        Kg = ( K + mesh_min[2] - 3 + pm_nside)%pm_nside;
                        index  = (Ig * pm_nside + Jg ) * pm_nside_pad + Kg;
                        index -= first_index;
                        if (0<=index && index<slab_size)
                            pot_pad2[(I*NY_pad+J)*NZ_pad+K] += slab_recv[index];
                    }
                }
            }
        }
        MPI_Barrier(PM_COMM_WORLD);
    }
    /* pull potential from pm_rank to domain_rank */

    for (i=0; i<pot_pad2_size; i++)
        pm->potential[i] = (Real)pot_pad2[i];

    free(pot_pad2);   /* 4013 */

    free(slab_recv);  /* 4012 */
    free(slab_send);  /* 4011 */

#ifndef FFTW3_LIB
    free(data); /* 4006 */
#else
    fftw_free(cdata); /* 4009 */
    fftw_free(data);  /* 4008 */
#endif /* FFTW3_LIB */


    free(pm_local_xside); /* 4005 */
    free(pm_local_start); /* 4004 */

}


void* pm_acceleration(void* param)
{
    Partmesh *pm = (Partmesh*)param;
    Int n;
    double pos_min[3], boxsize, invdelta, delta;
    int pm_nside;
    int mesh_min[3];
    int mesh_max[3];
    int mesh_nside[3];
    int mesh_size;
    double *bin;
    Real *pot;
    Real *dpdx, *dpdy, *dpdz;
    Body *part;

    pm_nside = pm->pm_nside;
    boxsize = pm->boxsize;
    invdelta = (double)pm_nside/boxsize;
    delta = (double)boxsize/pm_nside;
    bin = pm->mesh_bin;

    Int npart = pm->npart;

    mesh_min[0] = pm->mesh_min[0];
    mesh_min[1] = pm->mesh_min[1];
    mesh_min[2] = pm->mesh_min[2];

    mesh_max[0] = pm->mesh_max[0];
    mesh_max[1] = pm->mesh_max[1];
    mesh_max[2] = pm->mesh_max[2];

    mesh_nside[0] = pm->mesh_nside[0];
    mesh_nside[1] = pm->mesh_nside[1];
    mesh_nside[2] = pm->mesh_nside[2];

    mesh_size = mesh_nside[0]*mesh_nside[1]*mesh_nside[2];
    dpdx = (Real*)xmalloc( sizeof(Real)*mesh_size, 4014 );
    dpdy = (Real*)xmalloc( sizeof(Real)*mesh_size, 4015 );
    dpdz = (Real*)xmalloc( sizeof(Real)*mesh_size, 4016 );

    pot = pm->potential;
    part= pm->part;

    MPI_Barrier(PM_COMM_WORLD);
    int rank;
    MPI_Comm_rank(PM_COMM_WORLD, &rank);

    int x0, x_2, x_1, x1, x2;
    int y0, y_2, y_1, y1, y2;
    int z0, z_2, z_1, z1, z2;
    int NX, NY, NZ, NX_PAD, NY_PAD, NZ_PAD, idx;

    NX = mesh_nside[0] ;
    NY = mesh_nside[1] ;
    NZ = mesh_nside[2] ;
    NX_PAD = NX + 4;
    NY_PAD = NY + 4;
    NZ_PAD = NZ + 4;
	
//	printf("[Pot %d] %f\n", rank, pot[0]);
    for (x0=2; x0<NX+2; x0++) {
        x_2 = x0-2;
        x_1 = x0-1;
        x1  = x0+1;
        x2  = x0+2;
        for (y0=2; y0<NY+2; y0++) {
            y_2 = y0-2;
            y_1 = y0-1;
            y1  = y0+1;
            y2  = y0+2;
            for (z0=2; z0<NZ+2; z0++) {
                z_2 = z0-2;
                z_1 = z0-1;
                z1  = z0+1;
                z2  = z0+2;

                idx = ((x0-2)*NY + y0-2)*NZ + z0-2;
                dpdx[idx] = -invdelta*((pot[(x1*NY_PAD + y0)*NZ_PAD + z0]
                                        - pot[(x_1*NY_PAD + y0)*NZ_PAD + z0])/1.5
                                       -(pot[(x2*NY_PAD + y0)*NZ_PAD + z0]
                                         - pot[(x_2*NY_PAD + y0)*NZ_PAD + z0])/12.0);
                dpdy[idx] = -invdelta*((pot[(x0*NY_PAD + y1)*NZ_PAD + z0]
                                        - pot[(x0*NY_PAD + y_1)*NZ_PAD + z0])/1.5
                                       -(pot[(x0*NY_PAD + y2)*NZ_PAD + z0]
                                         - pot[(x0*NY_PAD + y_2)*NZ_PAD + z0])/12.0);
                dpdz[idx] = -invdelta*((pot[(x0*NY_PAD + y0)*NZ_PAD + z1]
                                        - pot[(x0*NY_PAD + y0)*NZ_PAD + z_1])/1.5
                                       -(pot[(x0*NY_PAD + y0)*NZ_PAD + z2]
                                         - pot[(x0*NY_PAD + y0)*NZ_PAD + z_2])/12.0);

            }
        }
    }

    double pmx = mesh_min[0]*delta;
    double pmy = mesh_min[1]*delta;
    double pmz = mesh_min[2]*delta;
    double px, py, pz;
    double pos[3];
    int I, J, K, In, Jn, Kn;
    Real x, y, z, wx ,wy, wz, wxn, wyn, wzn;
    Real accx, accy, accz;

    Real acc_max, vel_max, vel_avg;
    acc_max = vel_max = vel_avg = 0.0;

    for (n=0; n<npart; n++) {
        pos[0] = (double)part[n].pos[0];
        pos[1] = (double)part[n].pos[1];
        pos[2] = (double)part[n].pos[2];

        assert( 0.0<=part[n].pos[0] && part[n].pos[0]<boxsize);
        assert( 0.0<=part[n].pos[1] && part[n].pos[1]<boxsize);
        assert( 0.0<=part[n].pos[2] && part[n].pos[2]<boxsize);

        px = pos[0] - pmx;
        py = pos[1] - pmy;
        pz = pos[2] - pmz;

        I = (int)(pos[0]*invdelta);
        if (pos[0] < bin[I] )
            I -= 1;

        J = (int)(pos[1]*invdelta);
        if (pos[1] < bin[J] )
            J -= 1;

        K = (int)(pos[2]*invdelta);
        if (pos[2] < bin[K] )
            K -= 1;

        I += 1 - mesh_min[0];
        J += 1 - mesh_min[1];
        K += 1 - mesh_min[2];

        x = (I-0.5)*delta;
        y = (J-0.5)*delta;
        z = (K-0.5)*delta;

        wxn = (px-x)*invdelta;
        wyn = (py-y)*invdelta;
        wzn = (pz-z)*invdelta;

        if (wxn<0.0) {
            In = I - 1;
            wxn = -wxn;
        }
        else
            In = I + 1;
        wx = 1.0 - wxn;

        if (wyn<0.0) {
            Jn = J - 1;
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

        accx= dpdx[(I *NY+J )*NZ+K ]*wx *wy *wz
            + dpdx[(In*NY+J )*NZ+K ]*wxn*wy *wz
            + dpdx[(I *NY+J )*NZ+Kn]*wx *wy *wzn
            + dpdx[(In*NY+Jn)*NZ+K ]*wxn*wyn*wz
            + dpdx[(In*NY+J )*NZ+Kn]*wxn*wy *wzn
            + dpdx[(I *NY+Jn)*NZ+Kn]*wx *wyn*wzn
            + dpdx[(I *NY+Jn)*NZ+K ]*wx *wyn*wz
            + dpdx[(In*NY+Jn)*NZ+Kn]*wxn*wyn*wzn;

        accy= dpdy[(I *NY+J )*NZ+K ]*wx *wy *wz
            + dpdy[(In*NY+J )*NZ+K ]*wxn*wy *wz
            + dpdy[(I *NY+J )*NZ+Kn]*wx *wy *wzn
            + dpdy[(In*NY+Jn)*NZ+K ]*wxn*wyn*wz
            + dpdy[(In*NY+J )*NZ+Kn]*wxn*wy *wzn
            + dpdy[(I *NY+Jn)*NZ+Kn]*wx *wyn*wzn
            + dpdy[(I *NY+Jn)*NZ+K ]*wx *wyn*wz
            + dpdy[(In*NY+Jn)*NZ+Kn]*wxn*wyn*wzn;

        accz= dpdz[(I *NY+J )*NZ+K ]*wx *wy *wz
            + dpdz[(In*NY+J )*NZ+K ]*wxn*wy *wz
            + dpdz[(I *NY+J )*NZ+Kn]*wx *wy *wzn
            + dpdz[(In*NY+Jn)*NZ+K ]*wxn*wyn*wz
            + dpdz[(In*NY+J )*NZ+Kn]*wxn*wy *wzn
            + dpdz[(I *NY+Jn)*NZ+Kn]*wx *wyn*wzn
            + dpdz[(I *NY+Jn)*NZ+K ]*wx *wyn*wz
            + dpdz[(In*NY+Jn)*NZ+Kn]*wxn*wyn*wzn;

        part[n].acc_pm[0] = accx;
        part[n].acc_pm[1] = accy;
        part[n].acc_pm[2] = accz;

    }

    free(dpdx); /* 4014 */
    free(dpdy); /* 4015 */
    free(dpdz); /* 4016 */

}

