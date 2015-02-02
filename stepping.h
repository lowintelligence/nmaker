#ifndef _STEPPING_H_
#define _STEPPING_H_

#include "global.h"

typedef struct {
    int step;
    int nsteps;
    double  time;
    double  time_start;
    double  time_end;
    double  dstep;
    double *kick1;
    double *kick2;
    double *draft;
    double  boxsize;
} Stepping;

void kick1_system(Stepping *st, System *sys);
void kick2_system(Stepping *st, System *sys);
void draft_system(Stepping *st, System *sys);
void kick1_pp(double k1, System *sys);
void kick2_pp(double k2, System *sys);
void draft_pp(double dr, System *sys);
void update_stepping(Stepping *st);

Stepping* create_stepping_fixed_loga(Constants *constants, int nsteps, double a_init, double a_end);
void free_stepping(Stepping *st);

#endif // _STEPPING_H_ 
