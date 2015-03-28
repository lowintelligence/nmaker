#include "global.h"
#include "initial.h"
#include "iosnap.h"
#include <assert.h>
#include <string.h>
#include <math.h>

void readin_parameter_file(InputParam* param_ptr, char parameterfile[])
{
    char name[256];
    char value[256];
    FILE *fd;
    if ( !( fd = fopen(parameterfile, "r") ) ) {
        printf(" Cannot open parameterfile `%s'\n", parameterfile);
        system_exit(0);
    }

    while (2 == fscanf(fd, "%s %s", name, value) ) {
        if ( 0 == strcmp(name, "OMEGAM") ) {
            param_ptr->OMEGA_M0 = (Real) atof(value);
        }
        else if ( 0 == strcmp(name, "OMEGAX") ) {
            param_ptr->OMEGA_X0 = (Real) atof(value);
        }
        else if ( 0 == strcmp(name, "HUBBLE") ) {
            param_ptr->HUBBLE_0 = (Real) atof(value);
        }
        else if ( 0 == strcmp(name, "REDSHIFT") ) {
            param_ptr->REDSHIFT = (Real) atof(value);
        }
        else if ( 0 == strcmp(name, "OPENANGL") ) {
            param_ptr->OPEN_ANGLE = atof(value);
        }
        else if ( 0 == strcmp(name, "BOXSIZE") ) {
            param_ptr->BOX_SIZE  = atof(value);
        }
        else if ( 0 == strcmp(name, "SOFTING") ) {
            param_ptr->SOFTEN_SCALE = atof(value);
        }
        else if ( 0 == strcmp(name, "MESHPOW") ) {
            param_ptr->MESH_BITS  = atoi(value);
        }
        else if ( 0 == strcmp(name, "PEANOBIT") ) {
            param_ptr->PEANO_BITS = atoi(value);
        }
        else if ( 0 == strcmp(name, "NUMPART") ) {
            param_ptr->TOTAL_NUM_PART = atol(value);
        }
        else if ( 0 == strcmp(name, "ICPATH") ) {
            sprintf(param_ptr->IC_PATH_NAME, "%s", value);
        }
        else if ( 0 == strcmp(name, "NTHCPU") ) {
            param_ptr->NTHREAD_CPU = atoi(value);
        }
        else if ( 0 == strcmp(name, "NTHMIC") ) {
            param_ptr->NTHREAD_MIC = atoi(value);
        }
        else if ( 0 == strcmp(name, "NTHTREE") ) {
            param_ptr->NTHREAD_TREE = atoi(value);
        }
        else if ( 0 == strcmp(name, "NPNODE") ) {
            param_ptr->NP_PER_NODE = atoi(value);
        }
        else if ( 0 == strcmp(name, "MICNUM") ) {
            param_ptr->NMIC_PER_NODE = atoi(value);
        }
        else if ( 0 == strcmp(name, "NTEAMCPU") ) {
            param_ptr->NTEAM_CPU = atoi(value);
        }
        else if ( 0 == strcmp(name, "NTEAMMIC") ) {
            param_ptr->NTEAM_MIC = atoi(value);
        }
        else if ( 0 == strcmp(name, "RNDPART") ) {
            param_ptr->RND_PART = atoi(value);
        }
        else if ( 0 == strcmp(name, "CPURATIO") ) {
            param_ptr->CPU_RATIO = atof(value);
        }
        else if ( 0 == strcmp(name, "DYNLOAD") ) {
            param_ptr->DYNAMIC_LOAD = atoi(value);
        }
        else if ( 0 == strcmp(name, "TREEMAX") ) {
            param_ptr->MAX_TREE_LEVEL = atoi(value);
        }
        else if ( 0 == strcmp(name, "TREEMIN") ) {
            param_ptr->MIN_TREE_LEVEL = atoi(value);
        }
        else if ( 0 == strcmp(name, "NCRITI") ) {
            param_ptr->MAX_PACKAGE_SIZE = atoi(value);
        }
        else if ( 0 == strcmp(name, "PPSIZE") ) {
            param_ptr->MIC_PP_THRESHOLD = atoi(value);
        }
        else {
            printf("check the parameters %s, %s\n", name, value);
            exit(0);
        }
    }
}

void load_input_parameter_file(char filename[], InputParam* param_ptr)
{
    param_ptr->BOX_SIZE = 100000.0;
    param_ptr->OMEGA_M0 = 0.25;
    param_ptr->OMEGA_X0 = 0.75;
    param_ptr->HUBBLE_0 = 0.70;
    param_ptr->MESH_BITS   = 4; // 6 * ( 2^6 = 64 ) = 384
    param_ptr->PEANO_BITS  = 4;
    param_ptr->REDSHIFT    = 49.0;
    param_ptr->TOTAL_NUM_PART = 128*128*128;
    param_ptr->OPEN_ANGLE     =  0.3;
    param_ptr->SOFTEN_SCALE   = 10.0;
    sprintf( (param_ptr->IC_PATH_NAME), "./ic_128");
}

void init_random_particles(Constants *constants, System *sys, Long n_start, Long n_count, int myid)
{
	Int n;
	int m;

#ifdef NMK_VERIFY
    srand48(0);
	double b=0;
	for (n=0; n<n_count*3*myid; n++)
	{
		b+=drand48();		
	}
	printf("%f\n", b);
#else
    srand48(8888*myid);
#endif
	double velocity = constants->BOX_SIZE/((Real)(1<<(constants->MESH_BITS+1)));
    for (n=0; n<n_count; n++)
    {
		for (m=0; m<DIM; m++)
		{
#ifdef NMK_VERIFY
			sys->part[n].pos[m] = (Real) drand48()*24998.0+87501.0;
			sys->part[n].vel[m] = (Real) 0.0;
#else
			sys->part[n].pos[m] = (Real) drand48()*constants->BOX_SIZE;
			sys->part[n].vel[m] = (Real) drand48()*velocity;
#endif
			sys->part[n].acc[m]	= (Real) 0.0;
		}
		sys->part[n].id = n_start + n + 1;
		sys->part[n].mass = (Real) constants->PART_MASS;
    }
#if ((defined NMK_VERIFY) && (defined NMK_GEN_VER_FILE))
	if (myid == 0)
	{
		FILE *fp = fopen("pt.ori", "w+");
		for (n=0; n<n_count; n++)
		{
			fprintf(fp, "%d  %f  %f  %f\n", sys->part[n].id, sys->part[n].pos[0], sys->part[n].pos[1], sys->part[n].pos[2]);
		}
		fclose(fp);
	}
#endif
}

void initialize_system(char fnamepara[], Constants *constparam_ptr, Status* status_ptr, System * sys_ptr)
{
    int r, n, m;
    const int root_proc = 0;
    int rank_proc, num_proc;
    int rank_mesh, num_mesh, rank_tree, num_tree;
    InputParam *input_p;

    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);

    MPI_Comm_dup( MPI_COMM_WORLD, &PM_COMM_WORLD );
    MPI_Comm_size(PM_COMM_WORLD, &num_mesh);
    MPI_Comm_rank(PM_COMM_WORLD, &rank_mesh);

    MPI_Comm_dup( MPI_COMM_WORLD, &TREE_COMM_WORLD );
    MPI_Comm_size(TREE_COMM_WORLD, &num_tree);
    MPI_Comm_rank(TREE_COMM_WORLD, &rank_tree);

    input_p = create_input_parameters();

    if ( root_proc == rank_proc) {
//        load_input_parameter_file(fnamepara, input_p);
        readin_parameter_file(input_p, fnamepara);

    }
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Bcast( input_p, (int)sizeof(InputParam), MPI_CHAR, root_proc, MPI_COMM_WORLD);

    constparam_ptr->BOX_SIZE = input_p->BOX_SIZE;
    constparam_ptr->OMEGA_M0 = input_p->OMEGA_M0;
    constparam_ptr->OMEGA_X0 = input_p->OMEGA_X0;
    constparam_ptr->HUBBLE_0 = input_p->HUBBLE_0;
    constparam_ptr->NTHREAD_CPU  = input_p->NTHREAD_CPU;
    constparam_ptr->NTHREAD_MIC  = input_p->NTHREAD_MIC;
    constparam_ptr->NTHREAD_TREE  = input_p->NTHREAD_TREE;
    constparam_ptr->NP_PER_NODE  = input_p->NP_PER_NODE;
    constparam_ptr->NMIC_PER_NODE  = input_p->NMIC_PER_NODE;
    constparam_ptr->NTEAM_CPU  = input_p->NTEAM_CPU;
    constparam_ptr->NTEAM_MIC  = input_p->NTEAM_MIC;
    constparam_ptr->CPU_RATIO  = input_p->CPU_RATIO;
    constparam_ptr->DYNAMIC_LOAD  = input_p->DYNAMIC_LOAD;
    constparam_ptr->MAX_TREE_LEVEL  = input_p->MAX_TREE_LEVEL;
    constparam_ptr->MIN_TREE_LEVEL  = input_p->MIN_TREE_LEVEL;
    constparam_ptr->MAX_PACKAGE_SIZE  = input_p->MAX_PACKAGE_SIZE;
    constparam_ptr->MIC_PP_THRESHOLD  = input_p->MIC_PP_THRESHOLD;
	constparam_ptr->PEANO_BITS  = input_p->PEANO_BITS;
    constparam_ptr->MESH_BITS   = input_p->MESH_BITS;
    constparam_ptr->TOTAL_NUM_PART = input_p->TOTAL_NUM_PART;
    constparam_ptr->OPEN_ANGLE     = input_p->OPEN_ANGLE;
    constparam_ptr->SOFTEN_SCALE   = input_p->SOFTEN_SCALE;

//    printf(" mesh_level = %d\n", constparam_ptr->MESH_BITS);

    Real massp;
    massp = (3.0*(input_p->OMEGA_M0)*0.01)/1080886.296;
    massp *= pow(input_p->BOX_SIZE,3.0)/(input_p->TOTAL_NUM_PART);
    /* 10^10 solarmass/h */

    constparam_ptr->PART_MASS = (Real)massp;
    constparam_ptr->NUM_PROCESS = num_proc;
    constparam_ptr->PM_NUM_SIDE = 6*(1<<(constparam_ptr->MESH_BITS));
//    constparam_ptr->PM_NUM_SIDE = 4*(1<<(constparam_ptr->MESH_BITS));


    Real width_grid = (Real)(input_p->BOX_SIZE)/(constparam_ptr->PM_NUM_SIDE);

    constparam_ptr->SPLIT_SCALE  = 1.25*width_grid;
    constparam_ptr->CUTOFF_SCALE = 6.0*width_grid;
    constparam_ptr->GRAV_CONST = 43007.1;
	constparam_ptr->EPS2 = (constparam_ptr->SOFTEN_SCALE*constparam_ptr->SOFTEN_SCALE);

#ifdef NMK_PP_TAB
    Real x, dx, rs, G, eps2;
    Real sqrt_pi = sqrt(M_PI);
    constparam_ptr->delta=(Real)(constparam_ptr->CUTOFF_SCALE)/((double)TABLEN);
    dx = constparam_ptr->delta;
    rs = constparam_ptr->SPLIT_SCALE;
	G  = constparam_ptr->GRAV_CONST;
	eps2 = constparam_ptr->EPS2;
	
	LOG_INFO(LOG_VB1, " TAB = %d G = %f rs = %f width_grid = %f dx = %f eps2 = %f\n", TABLEN, G, rs, width_grid, dx, eps2 );

    for (n=0; n<TABLEN; n++) {
        x = n*dx;
        constparam_ptr->value[n] = G*(1.0-erf(x/2/rs)+(x/sqrt_pi/rs)*exp(-x*x/4/rs/rs))/pow( (x*x+rs*rs),1.5);
    }
    for (n=0; n<TABLEN-1; n++) {
        constparam_ptr->slope[n]=(constparam_ptr->value[n+1] - constparam_ptr->value[n])/dx;
    }
#endif 

	

//    Long start_read = (Long)(((float)( rank_proc )/(num_proc))*(constparam_ptr->TOTAL_NUM_PART));
//    Long npart_read = (Long)(((float)(rank_proc+1)/(num_proc))*(constparam_ptr->TOTAL_NUM_PART));
//    npart_read -= start_read;
	Long allpart = constparam_ptr->TOTAL_NUM_PART;
    Long start_read = allpart/num_proc*rank_proc+(allpart%num_proc+rank_proc)/num_proc*(rank_proc+allpart%num_proc-num_proc);
    Long npart_read = (allpart+rank_proc)/num_proc;

    sys_ptr->num_part = (Int)npart_read;
    sys_ptr->part = (Body*)xmalloc( npart_read*sizeof(Body), 1001 );

	if (input_p->RND_PART == 0)
	{
	    Snapshot *ic = create_snapshot(input_p->IC_PATH_NAME);
	    read_snapshot(sys_ptr->part, start_read, npart_read, ic);
	    free_snapshot(ic);
	}
	else
	{
		init_random_particles(constparam_ptr, sys_ptr, start_read, npart_read, rank_proc);
	}

    for (n=0; n<npart_read; n++)
	{
        sys_ptr->part[n].nstep = 1;
		for (m=0; m<DIM; m++)
		{
			if(sys_ptr->part[n].pos[m]>=constparam_ptr->BOX_SIZE)
			{
				sys_ptr->part[n].pos[m] = (Real) 0.0;
			}
		}
    }
//   write_snapshot_node(sys_ptr->part, npart_read, "snap1", constparam_ptr, status_ptr);

    status_ptr->redshift = input_p->REDSHIFT;
    status_ptr->timestep = 0.0;
    status_ptr->redshift_next_snap = 0.0;

    status_ptr->num_part_rank = (Int*)xmalloc(num_proc*sizeof(Int), 1003);
    status_ptr->lower_rank  = (Peano*)xmalloc(num_proc*sizeof(Peano), 1004);
    status_ptr->upper_rank  = (Peano*)xmalloc(num_proc*sizeof(Peano), 1005);
    status_ptr->weight_rank = (double*)xmalloc(num_proc*sizeof(Peano),1006);
    status_ptr->frac_sample = 1.0/0.0001;
//status_ptr->frac_sample = 1.0/0.001;

    Peano lower, upper, max;
    max = 1L<<(3*constparam_ptr->PEANO_BITS);
    lower = (Peano)(rank_proc*(max/num_proc) );
    upper = (Peano)( (rank_proc+1)*(max/num_proc) );

    double weight = 1.0;// + 0.1*rank_proc;

    MPI_Allgather(&(sys_ptr->num_part), sizeof(Int), MPI_CHAR,
                  status_ptr->num_part_rank, sizeof(Int), MPI_CHAR, MPI_COMM_WORLD);

    MPI_Allgather(&lower, sizeof(Peano), MPI_CHAR, status_ptr->lower_rank,
                  sizeof(Peano), MPI_CHAR, MPI_COMM_WORLD);

    MPI_Allgather(&upper, sizeof(Peano), MPI_CHAR, status_ptr->upper_rank,
                  sizeof(Peano), MPI_CHAR, MPI_COMM_WORLD);

    MPI_Allgather(&weight, 1, MPI_DOUBLE, status_ptr->weight_rank,
                  1, MPI_DOUBLE, MPI_COMM_WORLD);

//    printf("%u %lu %lu %d\n", sys_ptr->num_part, lower, upper, constparam_ptr->PEANO_BITS);

    MPI_Barrier(MPI_COMM_WORLD);

    free_input_parameters(input_p);


//   setup_partmesh_environment(constparam_ptr, status_ptr);

}

void finalize_system(Constants *constparam_ptr, Status* status_ptr, System * sys_ptr)
{
    if (sys_ptr) {
        if (sys_ptr->part)
            free(sys_ptr->part);
    }
    if (status_ptr) {
        if (status_ptr->num_part_rank)
            free(status_ptr->num_part_rank);
        if (status_ptr->lower_rank)
            free(status_ptr->lower_rank);
        if (status_ptr->upper_rank)
            free(status_ptr->upper_rank);
        if (status_ptr->weight_rank)
            free(status_ptr->weight_rank);
    }
}

InputParam* create_input_parameters()
{
    InputParam* param_ptr;
    param_ptr = (InputParam*)xmalloc( sizeof(InputParam), 1002);

    param_ptr->BOX_SIZE = 0.0;
    param_ptr->OMEGA_M0 = 0.0;
    param_ptr->OMEGA_X0 = 0.0;
    param_ptr->HUBBLE_0 = 0.0;
    param_ptr->REDSHIFT = 0.0;
    param_ptr->PEANO_BITS  = 0;
    param_ptr->MESH_BITS   = 0;
    param_ptr->TOTAL_NUM_PART = 0LL;
    param_ptr->OPEN_ANGLE     = 0.0;
    param_ptr->SOFTEN_SCALE   = 0.0;
    param_ptr->NTHREAD_CPU  = 1;
    param_ptr->NTHREAD_MIC  = 1;
    param_ptr->NTHREAD_TREE  = 1;
    param_ptr->NP_PER_NODE  = 1;
    param_ptr->NMIC_PER_NODE  = 1;
    param_ptr->NTEAM_CPU = 1;
    param_ptr->NTEAM_MIC = 1;
    param_ptr->CPU_RATIO  = 1.0;
    param_ptr->DYNAMIC_LOAD  = 0;
    param_ptr->MAX_TREE_LEVEL  = 5;
    param_ptr->MIN_TREE_LEVEL  = 1;
    param_ptr->MAX_PACKAGE_SIZE  = 128;
    param_ptr->MIC_PP_THRESHOLD  = 2048;
	param_ptr->RND_PART = 1;

    return param_ptr;
}

void free_input_parameters(InputParam* param_ptr)
{
    if (NULL != param_ptr)
        free(param_ptr);
    else {
        printf(" error: free_input_parameters \n");
        system_exit(0);
    }
}



