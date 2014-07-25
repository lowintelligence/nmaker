#include "driver.h"

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

	setup_parameters(&allparam);
	print_parameters(&allparam);

	double time_start,time_end;
	time_start=dtime();

	//    load_particle_into_domain(&dom, myid, allparam.NumSocket);
	init_particle_into_domain(&dom, &allparam, myid, allparam.NumSocket);
	time_end=dtime();
	printf("time A0 init:%lf\n",(double)(time_end-time_start));
	setup_partmesh_environment(&dom, &allparam);
	time_end=dtime();
	printf("time A1 setup partmesh:%lf\n",(double)(time_end-time_start));

	decompose_domain(&dom, &allparam);
	time_end=dtime();
	printf("time A2 decompose:%lf\n",(double)(time_end-time_start));

	construct_subcuboid(&dom, &allparam);
	time_end=dtime();
	printf("time A3 subcuboid:%lf\n",(double)(time_end-time_start));

	//    broadcast_frontiers(&dom, &allparam) ;
	pthread_t taskcomm;
	pthread_t taskmesh;
	TaskCommParam taskparam;

	taskparam.dp = &dom;
	taskparam.gp = &allparam;

	pthread_create(&taskcomm, NULL, thread_broad_front, &taskparam);   

	pthread_join(taskcomm, NULL);
	time_end=dtime();
	printf("time A4 broad front:%lf\n",(double)(time_end-time_start));


	build_subtree_on_subcuboid(&dom, &allparam, allparam.NumThreadPerSocket);
	time_end=dtime();
	printf("time A5 subtree:%lf\n",(double)(time_end-time_start));

	//	convolution_gravity(&taskparam);


	pthread_create(&taskmesh, NULL, convolution_gravity, &taskparam);
	//    
	pthread_join(taskmesh, NULL);
	time_end=dtime();
	printf("time A6 PM:%lf\n",(double)(time_end-time_start));

	dtt_traversal(&dom, &allparam);
	time_end=dtime();
	printf("time total:%lf\n",(double)(time_end-time_start));
	//    inner_traversal(&dom, &allparam, allparam.NumThreadPerSocket);


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

