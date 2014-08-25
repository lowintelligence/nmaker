#include "stepping.h"


void kick_nbody(Body part[], int npart, real dKick) 
{
    int n, d;
    for (n=0; n<npart; n++) {
        for (d=0; d<DIM; d++)
            part[n].vel[d] += part[n].acc[d]*dKick;
    }
}


void draft_nbody(Body part[], int npart, real dDraft) 
{
    int n, d;
    for (n=0; n<npart; n++) {
        for (d=0; d<DIM; d++)
            part[n].pos[d] += part[n].vel[d]*dDraft;
    }
}

void draft_kick_draft_fixed(Body part[], int npart, real dt) 
{
    draft_nbody(part, npart, 0.5*dt); 

    kick_nbody(part, npart, dt);

    draft_nbody(part, npart, 0.5*dt);
}

real get_stepping_increment(int iStep) {
    return 0.001;
}

