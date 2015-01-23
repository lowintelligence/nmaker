#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include "global.h"

#define TAG_IDLEFIELD -4
#define TAG_ADJOINING  0 // will be set to positive int as domain rank
#define TAG_FRONTIERS -1
#define TAG_OUTERCORE -2
#define TAG_INNERCORE -3

typedef struct {
	int rank;
	int nproc;
    int minimum[3]; // step 1
    int maximum[3]; // step 1
    int nSide[3];   // step 1
    int nPad;   // 4/6
    int cm_nside;
    int peano_bit;
    int mesh_bit;
    Peano  *lower;
    Peano  *upper;
    double boxsize;
    double split;
    double mass;
    double grav_const;
    double *cm_bin;
    int  *tag;
    Int  *count;
    Int  npart;     // No. of particles
    Int  npart_adj; // No. of ghost particles
	Int  nnode;
    Body *part;
	Node *tree;
    Int  *mroot;
	Int  *pidx;
} Domain;

Domain* create_domain(Constants *constants, Status *status, System *system);
void mark_cuboid(Domain *domain);
void broadcast_frontiers(Domain *domain, System *sys);
void free_domain(Domain *domain);

#endif /* _DOMAIN_H_ */
