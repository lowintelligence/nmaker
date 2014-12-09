#include "driver.h"

double dtime();

MPI_Comm MESH_COMM_WORLD;
MPI_Comm TREE_COMM_WORLD;

typedef struct {
    Domain *dp;
    GlobalParam *gp;
} TaskCommParam;

void* thread_broad_front(void* param) {
    TaskCommParam* tcp = (TaskCommParam*)param;
    broadcast_frontiers(tcp->dp, tcp->gp) ;
}

void driver(void) {
    int myid, numprocs;
    GlobalParam allparam;
    Domain dom;

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#ifdef __INTEL_OFFLOAD
	int micid=myid%2;
	Body *part = malloc(16777216*2/numprocs*sizeof(Body));
	Node *tree = malloc(8388608*2/numprocs*sizeof(Node));
#pragma offload target(mic:micid) inout (part[0:16777216*2/numprocs]:alloc_if(1) free_if(1)) inout (tree[0:8388608*2/numprocs]:alloc_if(1) free_if(1)) signal(part)
	{
		;
	}
#pragma offload_wait target(mic:micid) wait(part)
	free(part);
	free(tree);
#endif
	MPI_Barrier(MPI_COMM_WORLD);

    setup_parameters(&allparam);
    print_parameters(&allparam);

    double tstart, tstart0, tstop, ttime;
	tstart0 = dtime();
	tstart = dtime();

	//    load_particle_into_domain(&dom, myid, allparam.NumSocket);
	init_particle_into_domain(&dom, &allparam, myid, allparam.NumSocket);
	tstop = dtime();
	ttime = tstop - tstart;
	printf("time A0 init:%lf\n",ttime);

	int i;
	for (i=0; i<100; i++)
	{
	setup_partmesh_environment(&dom, &allparam);
	tstop = dtime();
	ttime = tstop - tstart;
	printf("time A1 setup partmesh:%lf\n",ttime);

	if(i>0) tstart = dtime();
	decompose_domain(&dom, &allparam);
	tstop = dtime();
	ttime = tstop - tstart;
	printf("time A2 decompose:%lf\n",ttime);

	construct_subcuboid(&dom, &allparam);
	tstop = dtime();
	ttime = tstop - tstart;
	printf("time A3 subcuboid:%lf\n",ttime);

	long oldpartnum = dom.NumPart;

	//    broadcast_frontiers(&dom, &allparam) ;
	pthread_t taskcomm;
	pthread_t taskmesh;
	TaskCommParam taskparam;

	taskparam.dp = &dom;
	taskparam.gp = &allparam;

	pthread_create(&taskmesh, NULL, convolution_gravity, &taskparam);

	pthread_join(taskmesh, NULL);
	tstop = dtime();
	ttime = tstop - tstart;
	printf("time A6 PM:%lf\n",ttime);
	pthread_create(&taskcomm, NULL, thread_broad_front, &taskparam);   

	pthread_join(taskcomm, NULL);
	tstop = dtime();
	ttime = tstop - tstart;
	printf("time A4 broad front:%lf\n",ttime);

	build_subtree_on_subcuboid(&dom, &allparam, allparam.NumThreadPerSocket);
	tstop = dtime();
	ttime = tstop - tstart;
	printf("time A5 subtree:%lf\n",ttime);
	dtt_traversal(&dom, &allparam);
	tstop = dtime();
	ttime = tstop - tstart;
	printf("time A7 traverse:%lf\n",ttime);

	clean_particles(&dom, &allparam, oldpartnum);
	draft_kick_draft_fixed(dom.Part,dom.NumPart,get_stepping_increment(0));
	tstop = dtime();
	ttime = tstop - tstart;
	printf("Step %d, time A8 stepping:%lf\n", i, ttime);

	free(dom.domtree->tree);
	}

	tstop = dtime();
	ttime = tstop - tstart0;
	printf("time total:%lf\n",ttime);

}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm_dup( MPI_COMM_WORLD, &MESH_COMM_WORLD );
    MPI_Comm_dup( MPI_COMM_WORLD, &TREE_COMM_WORLD );

    driver();

    MPI_Finalize();

    return 0;
}

