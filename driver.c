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
printf("A0\n");    

//    load_particle_into_domain(&dom, myid, allparam.NumSocket);
    init_particle_into_domain(&dom, &allparam, myid, allparam.NumSocket);
printf("A\n");    
	setup_partmesh_environment(&dom, &allparam);
printf("A1\n");    

    decompose_domain(&dom, &allparam);
printf("A2\n");    

    construct_subcuboid(&dom, &allparam);
printf("A3\n");    

//    broadcast_frontiers(&dom, &allparam) ;
    pthread_t taskcomm;
    pthread_t taskmesh;
    TaskCommParam taskparam;

    taskparam.dp = &dom;
    taskparam.gp = &allparam;

    pthread_create(&taskcomm, NULL, thread_broad_front, &taskparam);   

    pthread_join(taskcomm, NULL);
printf("A4\n");    
    

    build_subtree_on_subcuboid(&dom, &allparam, allparam.NumThreadPerSocket);

//	convolution_gravity(&taskparam);

    
    pthread_create(&taskmesh, NULL, convolution_gravity, &taskparam);
    
    pthread_join(taskmesh, NULL);
    
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

