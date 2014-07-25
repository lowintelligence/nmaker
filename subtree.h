#ifndef SUBTREE_H
#define SUBTREE_H

#include "domain.h"
#include "data.h"
#include "parameter.h"

#define MAX_MORTON_LEVEL 8
#define MIN_TREE_LEVEL 1
#define MAX_PACKAGE_SIZE 4096

void build_subtree_on_subcuboid(Domain *dp, GlobalParam *gp, int nThread);



#endif /* subtree_h */
