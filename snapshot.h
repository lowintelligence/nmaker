#ifndef _SNAPSHOT_H
#define _SNAPSHOT_H

#include "data.h"

typedef struct
{
    int npart[6];
    double mass[6];
    double time;
    double redshift;
    int flag_sfr;
    int flag_feedback;
    int npartTotal[6];
    int flag_cooling;
    int num_files;
    double BoxSize;
    double Omega0;
    double OmegaLambda;
    double HubbleParam;
    char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8];	/* fills to 256 Bytes */
} GadgetHeader;

Body* load_gadget2_snapshot(char *fname, int files, long int *npart, double *box);

void read_Particle(Body *pptr, char filename[], int n_start, int n_count);

#endif

