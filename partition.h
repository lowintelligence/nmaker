#ifndef PARTITION_H
#define PARTITION_H

#include "domain.h"
#include "data.h"
#include <mpi.h>
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "parameter.h"

void decompose_domain(Domain* dp, GlobalParam *gp);


#endif /* PARTITION_H */

