#include "global.h"
#include "stepping.h"
#include <math.h>

double a_flat_lcdm_t(double time);
double t_flat_lcdm_a(double a);

static double OmegaM, OmegaX;
static int Nblock;

void kick1_system(Stepping *st, System *sys)
{ 
    Int n, npart;
    npart = sys->num_part;
    Body *part = sys->part;
    int d;
    double kick1 = st->kick1[st->step];
    for (n=0; n<npart; n++) {
        for (d=0; d<3; d++) {
            part[n].vel[d] += kick1 * part[n].acc[d];
        }
    }
}

void kick2_system(Stepping *st, System *sys)
{ 
    Int n, npart;
    npart = sys->num_part;
    Body *part = sys->part;
    int d;
    double kick2 = st->kick2[st->step];
    for (n=0; n<npart; n++) {
        for (d=0; d<3; d++) {
            part[n].vel[d] += kick2 * part[n].acc[d];
        }
    }
}

void draft_system(Stepping *st, System *sys)
{
    int d;
    Int n, npart;
    npart = sys->num_part;
    Body *part = sys->part;
    double box = st->boxsize;
    double draft = st->draft[st->step];
    for (n=0; n<npart; n++) {
        for (d=0; d<3; d++) {
			double pos = (double) part[n].pos[d];
            pos += draft * (double) part[n].vel[d];

            if ( pos >= box )
			{
                pos -= box;
			}
            if ( pos < 0 )
			{
                pos += box;
			}
		
			part[n].pos[d] = (Real) pos;

            if ( part[n].pos[d] >= (Real) st->boxsize || part[n].pos[d] < 0)
                part[n].pos[d] = (Real) 0.0;
				
            if (part[n].pos[d] >= (Real) st->boxsize ) {
                printf(" bigger than bigger  %f box = %f %f\n", part[n].pos[d], box, part[n].vel[d] );
                system_exit(0);
            }
            if (part[n].pos[d] < 0) {
                printf(" smaller than smaller  %f box = %f %f\n", part[n].pos[d], box, part[n].vel[d] );
                system_exit(0);
            }
        }
    }
}

void update_stepping(Stepping *st)
{
    double ai = st->time_start;
    double af = st->time_end;
    st->step += 1;
    st->time = exp(log(ai)+(st->step)*(st->dstep));
}

double kick_loga(double loga_i, double loga_f) {
    int n;
    double kick_time = 0.0;
    double dloga = (loga_f - loga_i)/Nblock;
    double a_f = exp(loga_f);
    double a_i = exp(loga_i);
    double z1 = 1.0/(a_i);
    double h = 0.1*sqrt(OmegaM*z1*z1*z1 + OmegaX);
    kick_time = dloga*z1/h;
    for (n=1; n<Nblock; n++) {
        z1 = 1.0/(exp(loga_i+dloga*n));
        h = 0.1*sqrt(OmegaM*z1*z1*z1 + OmegaX);
        kick_time += 2.0*(1+n%2)*dloga*z1/h;
    }
    z1 = 1.0/(a_f);
    h = 0.1*sqrt(OmegaM*z1*z1*z1 + OmegaX);
    kick_time += dloga*z1/h;
    kick_time /= (3.0);
    return kick_time;
}

double draft_loga(double loga_i, double loga_f)
{
    int n;
    double kick_time = 0.0;
    double dloga = (loga_f - loga_i)/Nblock;
    double a_f = exp(loga_f);
    double a_i = exp(loga_i);
    double z1 = 1.0/(a_i);
    double h = 0.1*sqrt(OmegaM*z1*z1*z1 + OmegaX);
    kick_time = dloga*z1*z1/h;
    for (n=1; n<Nblock; n++) {
        z1 = 1.0/(exp(loga_i+dloga*n));
        h = 0.1*sqrt(OmegaM*z1*z1*z1 + OmegaX);
        kick_time += 2.0*(1+n%2)*dloga*z1*z1/h;
    }
    z1 = 1.0/(a_f);
    h = 0.1*sqrt(OmegaM*z1*z1*z1 + OmegaX);
    kick_time += dloga*z1*z1/h;
    kick_time /= (3.0);
    return kick_time;
}

Stepping* create_stepping_fixed_loga(Constants *constants, int nsteps, double a_init, double a_end)
{
    Stepping* st = (Stepping*)malloc(sizeof(Stepping));
    st->kick1 = (double*)malloc(sizeof(double)*nsteps);
    st->kick2 = (double*)malloc(sizeof(double)*nsteps);
    st->draft = (double*)malloc(sizeof(double)*nsteps);
    st->step = 0;
    st->nsteps = nsteps;
    st->time = a_init;
    st->time_start = a_init;
    st->time_end = a_end;

    Nblock = 2000;
    OmegaM = constants->OMEGA_M0;
    OmegaX = constants->OMEGA_X0;

    double loga_i, loga_f;
    loga_i = log(a_init);
    loga_f = log(a_end);
    double dloga = (loga_f - loga_i)/(nsteps);

    st->dstep = dloga;
    st->boxsize = constants->BOX_SIZE;

    int n;
    for (n=0; n<nsteps; n++) {
        st->kick1[n] = kick_loga(loga_i+n*dloga, loga_i+(n+0.5)*dloga )  ;
        st->kick2[n] = kick_loga(loga_i+(n+0.5)*dloga, loga_i+(n+1)*dloga )  ;
        st->draft[n] = draft_loga(loga_i+n*dloga, loga_i+(n+1)*dloga);
    }
    return st;
}

void free_stepping(Stepping *st)
{
    if ( NULL != st) {
        if (NULL != st->kick1)
            free(st->kick1);
        if (NULL != st->kick2)
            free(st->kick2);
        if (NULL != st->draft)
            free(st->draft);
        free(st);
    }
}

double a_flat_lcdm_t(double time)
{
    double t_star = 3.0*sqrt(OmegaX)/20.0;
    double kernel = sinh(t_star * time);
    double a = pow(kernel*kernel*OmegaM/OmegaX, 0.33333333f);
    return a;
}

double t_flat_lcdm_a(double a)
{
    double t_star = 3.0*sqrt(OmegaX)/20.0;
    double a3 = a*a*a;
    double f = OmegaX/OmegaM;
    double time = log( sqrt(f*a3) + sqrt(1.0+f*a3) )/t_star;
    return time;
}


