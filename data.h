/* 
 * NBODY SIMULATION
 * PM + Tree + PP
 * qwang@nao.cas.cn
 */


#ifndef __INTERFACE_H
#define __INTERFACE_H

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <mpi.h>
#include <pthread.h>

#define TAG_IDLEFIELD -4
#define TAG_ADJOINING  0 // will be set to positive int as domain rank
#define TAG_FRONTIERS -1
#define TAG_OUTERCORE -2
#define TAG_INNERCORE -3

#define PROG_INITIAL_ITERATION 10
#define PROG_EXCHANGE_FRONTIER 30
#define PROG_BUILD_EXOTIC_TREE 40
#define PROG_PARTICLEMESH_DONE 50

#define DIM 3
typedef float real;   // persicion switch
typedef uint32_t count;
typedef real vect3d[DIM];

#include "fillcurve.h"

typedef struct {
    long int ID; // constant particle ID
    code key; // peano-hilbert key, which could be overlayed.
	int mortonkey; //Meng!
    int tag;
    int group; // comes from cuboid
    /*    real phi; // gravitional potential, which is optional. */
    vect3d pos; // position of particle
//    vect3d vel; // velocity of particle
//    vect3d acc; // acceration of particle
    real  mass;
} Body;  // Body of N-Body, You knew it.

/* Old traditional node structure.
typedef struct {
    int nPart;
    int sub[8]; // index  of sub-oct-tree
    real width;  // width of this node
    real mass;   // mass of this node
//    vect3d center;  // geometry center
    vect3d masscenter; // mass center
} Node;
*/

// Cao! New structure for morton key tree.
// We now have continous memory local trees!
typedef struct {
	int nPart;
	long int firstpart; // index of particle.
	long int firstchild; // index of the first child according to big node array.
	int childnum;
	int level;
	int mortonkey;
	real width;
	real mass;
	vect3d masscenter;
} Node;


typedef struct {
    int NumTree;
    int *root_tree; // as long as subcuboid, record the root index
    int *root_cell;
    Node *tree;
	int *partidxes; // Cao! to indicate the particle position.
	int *numparticles; // Cao! to store the particle numbers in meshes.
} DomainTree;



extern MPI_Comm MESH_COMM_WORLD;
extern MPI_Comm TREE_COMM_WORLD;


#endif /* __INTERFACE_H */

