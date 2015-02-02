#include "global.h"

#include "initial.h"
#include "partition.h"
#include "domain.h"
#include "partmesh.h"
#include "subtree.h"
#include "traversal.h"
#include "stepping.h"

#include <pthread.h>

#ifdef NMK_DEBUG
int dbg_level = DBG_TIME | DBG_TMP;
#endif
int log_level = LOG_ALL;

void warmUpMIC(const int micid, const int nppermic)
{
#ifdef __INTEL_OFFLOAD
	size_t buf1cnt = (size_t)1048576*1024*2/sizeof(CalcBody);
	size_t buf2cnt = (size_t)1048576*1024*1/sizeof(Node);
	double *part = (double *)xmemalign(buf1cnt*sizeof(CalcBody), 129);
	double *tree = (double *)xmemalign(buf2cnt*sizeof(Node), 130);
#pragma offload target(mic:micid) inout (part[0:buf1cnt]:alloc_if(1) free_if(1)) inout (tree[0:buf2cnt]:alloc_if(1) free_if(1)) signal(tree)
    {
        ;
    }
#pragma offload_wait target(mic:micid) wait(tree)
    free(part);
    free(tree);
#endif /* __INTEL_OFFLOAD */
} /* warmUpMIC() */

#ifdef NMK_REFINE_STEP
int main(int argc, char* argv[])
{
    Constants constants;
    Status status;
    System sys;
    Partmesh *pm_param;
    Domain  *domain;
    MPI_Init(&argc, &argv);


    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 2) {
        printf("miss the parameter files\n");
    }

    initialize_system(argv[1], &constants, &status, &sys);

    MPI_Barrier(MPI_COMM_WORLD);

#ifdef FFTW3_LIB
    fftw_mpi_init();
#endif /* FFTW3_LIB */

    int p, n, m, d, Nstep;
    Real time, t0, t1, dt;
    Stepping *step;

    Nstep = 120;
    step = create_stepping_fixed_loga(&constants, Nstep, 0.02, 0.4);

#ifdef TEST_STEP
    if (0 == rank)
        for (n=0; n<Nstep; n++) {
            CLOG("%d %f %f %f\n", n, step->kick1[n], step->kick2[n], step->draft[n]);
        }

#endif
    double start_pm, end_pm, dt_pm;


#ifdef __INTEL_OFFLOAD
//    warmUpMIC(rank%constants.NP_PER_NODE%constants.NMIC_PER_NODE, constants.NP_PER_NODE/constants.NMIC_PER_NODE);
    if (rank == 0)
        DBG_INFOL(DBG_MEMORY, "MIC memory warm-up finished.\n");
#endif

    MPI_Barrier(MPI_COMM_WORLD);

    double start_time, tstamp, end_time, start_step;
    start_time = dtime();

    int REFSTEP = 4;
    for (n=0; n<Nstep; n++) {
        MPI_Barrier(MPI_COMM_WORLD);

        if (0==rank) {
            CLOG("\n - - - - - - - - - - step %5d - - - - - - - - - - \n", n);
            CLOG(" z = %.5f\n", status.redshift);
        }

        if (status.redshift <= 0.0)
            break;

		if ( 0==n ) {
			partition_system(&constants, &status, &sys);
			domain = create_domain(&constants, &status, &sys);

			mark_cuboid(domain) ;
			broadcast_frontiers(domain, &sys) ;
			printf(" * redistribute particles! - \n");

			pm_param = create_partmesh(&constants, &status, &sys);
			pthread_t pm_thread1, pm_thread2, pm_thread3;

			pthread_create(&pm_thread1, NULL, assign_particles, pm_param);
			pthread_join(pm_thread1, 0);

			pthread_create(&pm_thread2, NULL, convolution_gravity, pm_param);


			build_subtrees(domain, &constants);
			printf( " * build subtree \n" );

			pthread_join(pm_thread2, 0);

			printf(" * PM gravity -  \n");

			pthread_create(&pm_thread3, NULL, pm_acceleration, pm_param);
			pthread_join(pm_thread3, 0);
			free_partmesh(pm_param);
		}
		kick1_system(step, &sys);
		printf(" * kick system 1\n");


		double rk1, rk2, rdr;
		rk1 = step->kick1[n]/REFSTEP;
		rk2 = step->kick2[n]/REFSTEP;
		rdr = step->draft[n]/REFSTEP;
		pthread_t pm_thrd;


		for (m=0; m<REFSTEP; m++)
		{
			printf(" --- n = %d, m = %d ---\n", n, m);

			if (0 == n && 0 == m) {
				tree_traversal(domain, &constants);
				printf(" * traversal - \n");
			}
			kick1_pp(rk1, &sys);
			printf(" * kick sub 1\n");
			draft_pp(rdr, &sys);
			printf(" * draft sub \n");

			if (REFSTEP - 1 == m) {
				/* periodic boundaries */
				double box = step->boxsize;
				for (p=0; p<sys.num_part; p++) {
					for (d=0; d<3; d++) {
						double pos = (double) sys.part[p].pos[d];

						if ( pos >= box )
						{
							pos -= box;
						}
						if ( pos < 0 )
						{
							pos += box;
					}
						sys.part[p].pos[d] = (Real) pos;

						if (sys.part[p].pos[d] >= box ) {
							printf(" bigger than bigger  %f box = %f %f\n", sys.part[p].pos[d], box, sys.part[p].vel[d] );
							system_exit(0);
						}
						if (sys.part[p].pos[d] < 0) {
							printf(" smaller than smaller  %f box = %f %f\n", sys.part[p].pos[d], box, sys.part[p].vel[d] );
							system_exit(0);
						}
					}
				}
				/* periodic boundaries */

				free_domain(domain);


				partition_system(&constants, &status, &sys);
				domain = create_domain(&constants, &status, &sys);

				mark_cuboid(domain) ;
				broadcast_frontiers(domain, &sys) ;
				printf(" * redistribute particles!\n");

				pm_param = create_partmesh(&constants, &status, &sys);
				pthread_t pm_thread1, pm_thread2, pm_thread3;

				pthread_create(&pm_thread1, NULL, assign_particles, pm_param);
				pthread_join(pm_thread1, 0);

				build_subtrees(domain, &constants);
				printf( " * build subtree \n" );

				pthread_create(&pm_thread2, NULL, convolution_gravity, pm_param);
				printf(" * traversal  \n");
				tree_traversal(domain, &constants);

				pthread_join(pm_thread2, 0);
				printf(" * PM gravity  \n");

				pthread_create(&pm_thread3, NULL, pm_acceleration, pm_param);
				pthread_join(pm_thread3, 0);
				free_partmesh(pm_param);
			}
			else {
				printf(" * traversal  \n");
				tree_traversal(domain, &constants);
			}

			kick2_pp(rk2, &sys);
			printf(" * kick sub 2 \n");
		}

		kick2_system(step, &sys);
		printf(" * kick system 2\n");
		update_stepping(step);
		status.redshift = 1.0/(step->time)-1.0;
	}

    free_domain(domain);
    /* main loop */
    free_stepping(step);

    write_snapshot_node(sys.part, sys.num_part, "output_snap", &constants, &status);

    finalize_system(&constants, &status, &sys);
	
    if (0==rank) {
        end_time = dtime();
        CLOG( " duration = %lf [sec]\n", tstamp-start_time);
    }

#ifdef FFTW3_LIB
    fftw_mpi_cleanup();
#endif /* FFTW3_LIB */

    MPI_Finalize();

    return 0;
} /* main() */

#else // NMK_REFINE_STEP

int main(int argc, char* argv[])
{
    Constants constants;
    Status status;
    System sys;
    Partmesh *pm_param;
    Domain  *domain;
    MPI_Init(&argc, &argv);


    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 2) {
        printf("miss the parameter files\n");
    }

    initialize_system(argv[1], &constants, &status, &sys);

    MPI_Barrier(MPI_COMM_WORLD);

#ifdef FFTW3_LIB
    fftw_mpi_init();
#endif /* FFTW3_LIB */

    int n, Nstep;
    Real time, t0, t1, dt;
	Stepping *step;

    Nstep = 120;
	step = create_stepping_fixed_loga(&constants, Nstep, 0.02, 0.4);

// The first PM
    pm_param = create_partmesh(&constants, &status, &sys);
	assign_particles(pm_param);
    convolution_gravity(pm_param);

    pm_acceleration(pm_param);
	free_partmesh(pm_param);

#ifdef TEST_STEP
    if (0 == rank)
    for (n=0; n<Nstep; n++) {
        CLOG("%d %f %f %f\n", n, step->kick1[n], step->kick2[n], step->draft[n]);
    }

#endif
    double start_pm, end_pm, dt_pm;


#ifdef __INTEL_OFFLOAD
//	warmUpMIC(rank%constants.NP_PER_NODE%constants.NMIC_PER_NODE, constants.NP_PER_NODE/constants.NMIC_PER_NODE);
	if (rank == 0)
		DBG_INFOL(DBG_MEMORY, "MIC memory warm-up finished.\n");
#endif

	MPI_Barrier(MPI_COMM_WORLD);

    double start_time, tstamp, end_time, start_step;
    start_time = dtime();

    for (n=0; n<Nstep; n++) {
        MPI_Barrier(MPI_COMM_WORLD);

        if (0==rank) {
            CLOG("\n - - - - - - - - - - step %5d - - - - - - - - - - \n", n);
            CLOG(" z = %.5f\n", status.redshift);
        }

        if (status.redshift <= 0.0)
            break;

        start_step = dtime();
		tstamp = start_step;

        partition_system(&constants, &status, &sys);

		DBG_INFO(DBG_TIME, "[Step %d, Pt A] Time of distributing init: %.4f sec\n", n, dtime() - tstamp);
		tstamp = dtime();

        domain = create_domain(&constants, &status, &sys);

		DBG_INFO(DBG_TIME, "[Step %d, Pt B] Time of domain decomposition: %.4f sec\n", n, dtime() - tstamp);
		tstamp = dtime();

        mark_cuboid(domain) ;
        broadcast_frontiers(domain, &sys) ;

		DBG_INFO(DBG_TIME, "[Step %d, Pt C] Time of subdomain creating: %.4f sec\n", n, dtime() - tstamp);
		tstamp = dtime();

        kick1_system(step, &sys);

		DBG_INFO(DBG_TIME, "[Step %d, Pt D] Time of force calculation: %.4f sec\n", n, dtime() - tstamp);
		tstamp = dtime();

#ifdef TEST_POT
        int i, point_mass_rank = 0;
        double hbox = 0.5*(step->boxsize);
        double pos_test[3];
        for (i=0; i<sys.num_part; i++) {
            sys.part[i].mass = 0.0;
        }
        if (point_mass_rank == rank) {
            sys.part[0].mass = 1.0;
            pos_test[0] = sys.part[0].pos[0];
            pos_test[1] = sys.part[0].pos[1];
            pos_test[2] = sys.part[0].pos[2];
        }

        MPI_Bcast(pos_test, 3*sizeof(double), MPI_CHAR, point_mass_rank, MPI_COMM_WORLD);
#endif /* TEST_POT */

		pm_param = create_partmesh(&constants, &status, &sys);
		pthread_t pm_thread1, pm_thread2, pm_thread3;

		pthread_create(&pm_thread1, NULL, assign_particles, pm_param);
		pthread_join(pm_thread1, 0);

		pthread_create(&pm_thread2, NULL, convolution_gravity, pm_param);
		pthread_join(pm_thread2, 0);

		pthread_create(&pm_thread3, NULL, pm_acceleration, pm_param);
		pthread_join(pm_thread3, 0);
		free_partmesh(pm_param);

		DBG_INFO(DBG_LOGIC, "[%d] PM done !\n", rank);

#ifdef TEST_POT
        if (0==n) {

            int sendcount;

            double *acc_test = (double*)malloc(sizeof(double)*sys.num_part);
            double *dis_test = (double*)malloc(sizeof(double)*sys.num_part);
            int  *recvcount  = (int*)malloc(sizeof(int)*size);
            int  *recvdispl  = (int*)malloc(sizeof(int)*size);

            recvdispl[0] = 0;
            for (i=1; i<size; i++) {
                recvdispl[i] = recvdispl[i-1] + status.num_part_rank[i-1];
            }

            for (i=0; i<size; i++) {
                recvcount[i] = status.num_part_rank[i]*sizeof(double);
                recvdispl[i] *= sizeof(double);
            }

            for (i=0; i<sys.num_part; i++) {
                double dx, dy, dz, ax, ay, az;
                dx = sys.part[i].pos[0] - pos_test[0];
                dy = sys.part[i].pos[1] - pos_test[1];
                dz = sys.part[i].pos[2] - pos_test[2];

                if (dx < 0.0)
                    dx = -dx;
                if (dx > hbox)
                    dx = (step->boxsize) - dx;

                if (dy < 0.0)
                    dy = -dy;
                if (dy > hbox)
                    dy = (step->boxsize) - dy;

                if (dz < 0.0)
                    dz = -dz;
                if (dz > hbox)
                    dz = (step->boxsize) - dz;

                ax = sys.part[i].acc[0];
                ay = sys.part[i].acc[1];
                az = sys.part[i].acc[2];

                dis_test[i] = sqrt(dx*dx + dy*dy + dz*dz);
                acc_test[i] = sqrt(ax*ax + ay*ay + az*az);
            }

            double *acc_total =(double*)malloc(sizeof(double)*constants.TOTAL_NUM_PART);
            double *dis_total =(double*)malloc(sizeof(double)*constants.TOTAL_NUM_PART);

            sendcount = sys.num_part*sizeof(double);

            MPI_Barrier(MPI_COMM_WORLD);

            MPI_Allgatherv(acc_test, sendcount, MPI_CHAR, acc_total, recvcount, recvdispl, MPI_CHAR, MPI_COMM_WORLD);
            MPI_Allgatherv(dis_test, sendcount, MPI_CHAR, dis_total, recvcount, recvdispl, MPI_CHAR, MPI_COMM_WORLD);

            if ( point_mass_rank == rank ) {
                FILE *fd_pot = fopen("dr_gravity.txt", "w");

                for (i=0; i<constants.TOTAL_NUM_PART; i++) {
                    fprintf(fd_pot, "%e %e\n", dis_total[i], acc_total[i]/(constants.GRAV_CONST) );
                }
                fclose(fd_pot);
            }

            free(recvcount);
            free(recvdispl);
            free(acc_total);
            free(dis_total);
            free(dis_test);
            free(acc_test);
        }
#endif /* TEST_POT */

		DBG_INFO(DBG_TIME, "[Step %d, Pt E] Time of PM calculation: %.4f sec\n", n, dtime() - tstamp);
		tstamp = dtime();

		build_subtrees(domain, &constants);
		DBG_INFO(DBG_TIME, "[Step %d, Pt F] Time of tree building: %.4f sec\n", n, dtime() - tstamp);
		tstamp = dtime();

		tree_traversal(domain, &constants);
		DBG_INFO(DBG_TIME, "[Step %d, Pt G] Time of force calculation: %.4f sec\n", n, dtime() - tstamp);
		tstamp = dtime();

        draft_system(step, &sys);
        kick2_system(step, &sys);

        free_domain(domain);

        update_stepping(step);
        status.redshift = 1.0/(step->time) - 1.0;

        if (0==rank) {
            tstamp =  dtime();
            CLOG( " step[%d] use %lf sec\n",n, tstamp-start_step );
            CLOG( " duration = %lf [sec]\n", tstamp-start_time);
        }

    } /* main loop */
    free_stepping(step);

    write_snapshot_node(sys.part, sys.num_part, "output_snap", &constants, &status);

    finalize_system(&constants, &status, &sys);
    if (0==rank) {
        end_time = dtime();
        CLOG( " duration = %lf [sec]\n", tstamp-start_time);
    }

#ifdef FFTW3_LIB
    fftw_mpi_cleanup();
#endif /* FFTW3_LIB */

    MPI_Finalize();

    return 0;
} /* main() */


#endif /* NMK_REFINE_STEP */

