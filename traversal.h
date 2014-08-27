#ifndef _TRAVERSAL_H_
#define _TRAVERSAL_H_

#define _GNU_SOURCE
#include <sched.h>
#include <unistd.h>
#include <pthread.h>
#include <malloc.h>

#include "data.h"
#include "dtime.h"
#include "offload.h"
#include "domain.h"
#include "subtree.h"
#include "parameter.h"
#include "subcuboid.h"
#include "ppkernel.h"
#include "queue.h"

#define SQROOT3 1.73205080705688773

// Qiao's functions.
//int accepted_cell_branes(int one, int ns, Node *tree, Body *part);
//void naive_walk_cell(int add, int level, int body, Node *tree, Body *part);
//void inner_traversal(Domain *dp, GlobalParam *gp, int nThread);

// Cao!'s functions.
__OffloadFunc_Macro__
int accepted_cell_to_cell(int TA, int TB, double theta);

__OffloadFunc_Macro__
int dtt_process_cell(int TA, int TB, double theta, TWQ* PQ_treewalk, Block* pth);

void dtt_traversal(Domain *dp, GlobalParam *gp);

#endif /* TRAVERSAL_H */
