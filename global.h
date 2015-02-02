#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include "offload.h"
#include "proto.h"
#include "nmk_utils.h"

#ifdef NMK_PP_TAB
#define TABLEN 128
#endif

typedef struct {
    /* phyical parameters */
    Real BOX_SIZE;
    Real OMEGA_M0;
    Real OMEGA_X0;
    Real HUBBLE_0;

    /* system parameters*/
    int NUM_PROCESS;
    int NTHREAD_CPU;
    int NTHREAD_MIC;
    int NTHREAD_TREE;
	int NP_PER_NODE;
	int NMIC_PER_NODE;
	int NTEAM_CPU;
    int NTEAM_MIC;
	double CPU_RATIO;
	int DYNAMIC_LOAD;

    /* runtime parameters */
    int  PEANO_BITS;
    int  MESH_BITS;
	int  MAX_TREE_LEVEL;
	int  MIN_TREE_LEVEL;
	int  MAX_PACKAGE_SIZE;
	int  MIC_PP_THRESHOLD;
    int  PM_NUM_SIDE;
#ifdef NMK_PP_TAB
    Real delta;
    Real value[TABLEN];
    Real slope[TABLEN];
#endif
    Long TOTAL_NUM_PART;
    Real PART_MASS;
    Real OPEN_ANGLE;
    Real SPLIT_SCALE;
    Real CUTOFF_SCALE;
    Real SOFTEN_SCALE;
    Real GRAV_CONST;
	Real EPS2;
} Constants; // for the whole system, never change

typedef struct {
    Real   redshift;
    Real   timestep;
    Real   redshift_next_snap;
    Int    *num_part_rank;
    Peano  *lower_rank;
    Peano  *upper_rank;

    double frac_sample;
    double *weight_rank;
} Status; // for the whole system, changed by step

typedef struct {
    Body *part;
    Int  num_part;
} System;

#include <mpi.h>
extern MPI_Comm PM_COMM_WORLD;
extern MPI_Comm TREE_COMM_WORLD;

#include <stdio.h>
#include <stdlib.h>

#define ERR_MEMORY 128
#define ALIGNCNT 64
__OffloadFunc_Macro__
void system_exit(int err_num);

__OffloadFunc_Macro__
void* xmalloc(size_t memo_size, int err_num);

__OffloadFunc_Macro__
void* xrealloc(void *mem, size_t memo_size, int err_num);

__OffloadFunc_Macro__
void* xmemalign(size_t memo_size, int err_num);

Peano encoding_peano(int x, int y, int z, int nBits);

#endif /* _GLOBAL_H_ */

