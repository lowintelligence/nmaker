#ifndef INITIAL_H
#define INITIAL_H

#include "global.h"

typedef struct
{
    /* phyical parameters */
    double BOX_SIZE;
    double OMEGA_M0;
    double OMEGA_X0;
    double HUBBLE_0;
    double REDSHIFT;
    int    PEANO_BITS;

    /* system parameters*/
    int NTHREAD_CPU;
    int NTHREAD_MIC;
	int NTHREAD_TREE;
	int NP_PER_NODE;
	int NMIC_PER_NODE;
	int NTEAM_CPU;
    int NTEAM_MIC;
	double CPU_RATIO;
	int DYNAMIC_LOAD;
	int RND_PART;

    /* runtime parameters */
	int MAX_TREE_LEVEL;
	int MIN_TREE_LEVEL;
	int MAX_PACKAGE_SIZE;
	int MIC_PP_THRESHOLD;
    int  MESH_BITS;
    Long TOTAL_NUM_PART;
    double OPEN_ANGLE;
    double SOFTEN_SCALE;    
    Int    IC_NUM_FILE;
    char   IC_PATH_NAME[256];

} InputParam;

InputParam* create_input_parameters();
void free_input_parameters(InputParam* param_ptr);
void load_input_parameter_file(char filename[], InputParam* param_ptr);

void initialize_system(char fnamepara[], Constants *constparam_ptr, Status* status_ptr, System * sys_ptr) ;
void finalize_system(Constants *constparam_ptr, Status* status_ptr, System * sys_ptr);

#endif /* INITIAL_H */

