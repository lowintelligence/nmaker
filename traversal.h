#ifndef TRAVERSAL_H
#define TRAVERSAL_H

#include "data.h"
#include "domain.h"
#include "subtree.h"
#include "parameter.h"
#include "subcuboid.h"
#include "ppkernel.h"

#define SQROOT3 1.73205080705688773

// Qiao's functions.
int accepted_cell_branes(int one, int ns, Node *tree, Body *part);
void naive_walk_cell(int add, int level, int body, Node *tree, Body *part);
void inner_traversal(Domain *dp, GlobalParam *gp, int nThread);

// Cao!'s functions.
int accepted_cell_to_cell(int TA, int TB, double theta);
void dtt_process_cell(int TA, int TB, double theta);

#endif /* TRAVERSAL_H */
