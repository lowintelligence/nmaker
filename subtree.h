#ifndef SUBTREE_H
#define SUBTREE_H

#include "domain.h"
#include "data.h"
#include "parameter.h"

#define MAX_MORTON_LEVEL 5
#define MIN_TREE_LEVEL 1
#define MAX_PACKAGE_SIZE 2048

void build_subtree_on_subcuboid(Domain *dp, GlobalParam *gp, int nThread);



#endif /* subtree_h */
