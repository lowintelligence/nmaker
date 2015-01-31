#ifndef _PARTMESH_H_
#define _PARTMESH_H_

#include "global.h"

typedef struct {
    Body *part;
    Int  npart;
    Real boxsize;
    Real r_split;
    int pm_nside;
    Real *potential;
    int  mesh_nside[3];
    int  mesh_max[3];
    int  mesh_min[3];
    double *mesh_bin;
    Real mass;
    Real grav_const;
    Real a_time;
    double *local_field;
} Partmesh;

Partmesh* create_partmesh(Constants *constants, Status *status, System *system);
void* partmesh_acc(void *param);
void free_partmesh(Partmesh *pm);

void* assign_particles(void *param);
void* convolution_gravity(void *param);
void* pm_acceleration(void *param);

#endif /* PARTMESH_H */
