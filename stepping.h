#ifndef _STEPPING_H
#define _STEPPING_H


#include "data.h"
#include <stdlib.h>
#include <stdio.h>


void draft_kick_draft_fixed(Body part[], int npart, real dt) ;
real get_stepping_increment(int iStep);

#endif 
