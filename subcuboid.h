#ifndef SUBCUBOID_H
#define SUBCUBOID_H

#include "data.h"
#include "parameter.h"
#include "domain.h"

/* need to thread paralled construct */
void construct_subcuboid(Domain *dp, GlobalParam *gp);

void clean_particles(Domain *dp, GlobalParam *gp, long oldpartnum);

void free_subcuboid(SubCuboid* sp) ;

#endif
