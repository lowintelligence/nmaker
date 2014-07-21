#include "traversal.h"
#include <time.h>
#include <math.h>

static Body *part;
static Node *tree; // thread number
static int FIRST_NODE;
double open_angle;

int accepted_cell_branes(int one, int ns, Node *tree, Body *part)
{
    int d, n, id;
    double dx[3], dr;
    if (one >= FIRST_NODE) {
        dr = 0.0;
        for (d=0; d<DIM; d++) {
            dx[d] = tree[ns].masscenter[d] - tree[one].masscenter[d];
            dr += dx[d]*dx[d];
        }
        dr = sqrt(dr) - tree[one].width ;
        if ( (tree[ns].width) < open_angle*dr )
            return 1;
        else
            return 0;
    }
    else if (one > -1)  {
        dr = 0.0;
        for (d=0; d<DIM; d++) {
            dx[d] = tree[ns].masscenter[d] - part[one].pos[d];
            dr += dx[d]*dx[d];
        }
        dr = sqrt(dr);
        if ( (tree[ns].width) < open_angle*dr )
            return 1;
        else
            return 0;
    }
    printf("error\n");
}

void naive_walk_cell(int add, int level, int body, Node *tree, Body *part)
{
    int id, n, d;
    double dx[3], dr3, dr2;
    if (add == body)
        return;

    if (add<FIRST_NODE) {
        dr2 = 0.0;
        for (d=0; d<DIM; d++) {
            dx[d] = part[add].pos[d] - part[body].pos[d];
            dr2 += dx[d]*dx[d];
        }
        dr3 = dr2 * sqrt(dr2);
        for (d=0; d<DIM; d++) {
            part[body].acc[d] += dx[d]/dr3;
        }
    }
    else
        for (n=0; n<8; n++) {
            id = tree[add].sub[n];
            if (id>-1) {
                if ( id >= FIRST_NODE ) {
                    if (accepted_cell_branes(body, id, tree, part)) {
                        dr2 = 0.0;
                        for (d=0; d<DIM; d++) {
                            dx[d] = tree[id].masscenter[d]-part[body].pos[d];
                            dr2 += dx[d]*dx[d];
                        }
                        dr3 = dr2 * sqrt(dr2);
                        for (d=0; d<DIM; d++) {
                            part[body].acc[d] += (tree[id].nPart) * dx[d]/dr3;
                        }
                    }
                    else
                        naive_walk_cell(id, level+1, body, tree, part);
                }
                else
                    naive_walk_cell(id, level+1, body, tree, part);
            }
        }
}


/* task pool */

void inner_traversal(Domain *dp, GlobalParam *gp, int nThread) {
    int i,j,k, ic, jc, kc;
    int In, Jn, Kn;
    int n, level;
    int group, index;
    int first_node;
    int npart;

    Node *tree;
    Body* part;
    DomainTree *dtp = dp->domtree;

    npart = dp->NumPart;
    part = dp->Part;
    tree = dp->domtree->tree;
    FIRST_NODE = dp->NumPart;

    open_angle = 0.3;
    level = gp->NumBits;
    int cnt = 0;
    In = dp->cuboid->nSide[0];
    Jn = dp->cuboid->nSide[1];
    Kn = dp->cuboid->nSide[2];

    for (n=0; n<npart; n++) {
        if (TAG_OUTERCORE == part[n].tag) {
            index = part[n].group;
            first_node = dtp->root_cell[index];
            naive_walk_cell(first_node, level, n, tree, part);
            cnt ++;
        }
    }
    printf("cnt = %d\n", cnt);

}


////////////////////////////////////////////////////////////////////////////////////////////////

/*
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define DIM 3
typedef double VECTOR[DIM];

typedef struct {
    VECTOR pos;
    VECTOR vel;
    VECTOR acc;
    double phi;
    long int ph_key;
    long int id;
} BODY;

typedef struct _node {
    int    nPart;
    int    sub[8];
    double width;
    VECTOR masscenter;
} NODE;


NODE *tree; // id of nodes
BODY *part;
static int NODELENGTH;
static int FIRST_NODE; // nPart
static int LAST_NODE;

void insert_particle_into_tree(VECTOR center,double size,int iNode,int iPart)
{
    int d, index, old, num;
    double qsize = 0.25*size;

    for (d=0, index=0; d<DIM; d++)
        if (part[iPart].pos[d] > center[d])
            index |= 1<<d ;

    num = tree[iNode].nPart ;
    tree[iNode].nPart += 1;
    tree[iNode].width = size;
    for (d=0; d<DIM; d++) {
        tree[iNode].masscenter[d] *= num;
        tree[iNode].masscenter[d] += part[iPart].pos[d];
        tree[iNode].masscenter[d] /= tree[iNode].nPart;
    }
    // this is empty
    if ( tree[iNode].sub[index] < 0) {
        tree[iNode].sub[index] = iPart;
        return;
    }
    for (d=0; d<DIM; d++)
        center[d] += (part[iPart].pos[d] > center[d]) ? qsize : -qsize ;
    // this is a particle
    if ( tree[iNode].sub[index] < FIRST_NODE ) {
        old = tree[iNode].sub[index]; // old particle;
        tree[iNode].sub[index] = ( ++LAST_NODE ) ; // new node
        insert_particle_into_tree(center, size/2, LAST_NODE, old);
    }
    // this is a node
    insert_particle_into_tree(center, size/2, tree[iNode].sub[index], iPart);
}

void construct_tree(int nPart, double size)
{
    int d, n;
    VECTOR center;
    FIRST_NODE = nPart; // nPart
    LAST_NODE  = nPart;
    NODELENGTH = nPart;

    tree = (NODE*)malloc( (int) sizeof(NODE) * NODELENGTH );
    for (n=0; n<NODELENGTH; n++) {
        tree[n].nPart = 0;
        for (d=0; d<(1<<DIM); d++)
            tree[n].sub[d] = -1;
    }
    tree -= FIRST_NODE;

    for (n=0; n<nPart; n++) {
        for (d=0; d<DIM; d++)
            center[d] = size/2;
        insert_particle_into_tree(center, size, FIRST_NODE, n);
    }
}

int *opencell;
int *accept;
double open_angle = 0.3;
double box ;

// only for particle-node
int accepted_cell_gadget(int one, int ns)
{
    double alpha = 0.001;
    int d, n, id;
    double dx[3], dr2, ext2;

        ext2 = (tree[ns].width) * (tree[ns].width) ;
        dr2 = 0.0;
        for (d=0; d<DIM; d++) {
            dx[d] = tree[ns].masscenter[d] - part[one].pos[d];
            dr2 += dx[d]*dx[d];
        }

        if ( (tree[ns].nPart)*ext2/dr2/dr2 <= alpha*(part[one].phi)) {
            return 1;
        }
        else
            return 0;

}

int accepted_cell_branes(int one, int ns)
{
    int d, n, id;
    double dx[3], dr;
    if (one >= FIRST_NODE) {
        dr = 0.0;
        for (d=0; d<DIM; d++) {
            dx[d] = tree[ns].masscenter[d] - tree[one].masscenter[d];
            dr += dx[d]*dx[d];
        }
        dr = sqrt(dr) - tree[one].width ;
        if ( (tree[ns].width/dr) < open_angle )
            return 1;
        else
            return 0;
    }
    else if (one > -1)  {
        dr = 0.0;
        for (d=0; d<DIM; d++) {
            dx[d] = tree[ns].masscenter[d] - part[one].pos[d];
            dr += dx[d]*dx[d];
        }
        dr = sqrt(dr);
        if ( (tree[ns].width/dr) < open_angle )
            return 1;
        else
            return 0;
    }
    printf("error\n");
}

int accepted_cell(int one, int ns) {
   // return accepted_cell_gadget(one, ns);
      return accepted_cell_branes(one, ns);
}

void gravity_sumover_accepted(int body,  int id_acc)
{
    int d, n, id;
    double dx[3], dr3, dr2;
    for (d=0; d<DIM; d++)
        part[body].acc[d] = 0.0;

//   printf("\n");
    for (n=0; n<id_acc; n++) {
        id = accept[n];
        //     printf("%d ", id);
        if (id != body)
            if (id >= FIRST_NODE ) {
                dr2 = 0.0;
                for (d=0; d<DIM; d++) {
                    dx[d] =  tree[id].masscenter[d] - part[body].pos[d];
                    dr2 += dx[d]*dx[d];
                }
                dr3 = dr2 * sqrt(dr2);
                for (d=0; d<DIM; d++) {
                    part[body].acc[d] += (tree[id].nPart) * dx[d]/dr3;
                }
            }
            else if (id > -1 ) {
                dr2 = 0.0;
                for (d=0; d<DIM; d++) {
                    dx[d] =  part[id].pos[d] - part[body].pos[d];
                    dr2 += dx[d]*dx[d];
                }
                dr3 = dr2 * sqrt(dr2);
                for (d=0; d<DIM; d++) {
                    part[body].acc[d] += dx[d]/dr3;
                }
            }
    }
}

void walk_cell(int one, int level, int os, int oe, int apt)
{
    int id, n, ns, ns2, cell;
    int new_os, new_oe, new_mid_oe, new_apt, new_mid_apt;

    if (one < 0)
        printf("ont a node or body\n");

    if ( os == oe && one < FIRST_NODE) {
        gravity_sumover_accepted(one, apt);
    }

    if ( os < oe && one < FIRST_NODE) {
        new_os = new_oe = oe+1;
        new_apt = apt;
        for (ns=os; ns<oe; ns++) {
            cell = opencell[ns] ;
            if ( accepted_cell(one, cell) ) {
                accept[new_apt] = cell;
                new_apt++;
            }
            else {
                for (n=0; n<8; n++) {
                    if (tree[cell].sub[n] >= FIRST_NODE) {
                        opencell[new_oe] = tree[cell].sub[n];
                        new_oe ++ ;
                    }
                    else if (tree[cell].sub[n] > -1 ) {
                        accept[new_apt] = tree[cell].sub[n];
                        new_apt++;
                    }
                }
            }
        }
        walk_cell(one, level, new_os, new_oe, new_apt);
    }

    if ( os == oe && one >= FIRST_NODE) {
        for (n=0; n<8; n++)
            if ( (id = tree[one].sub[n]) > -1 )
            {
                new_os = new_oe = oe+1;
                new_apt = apt;
                for (ns=n+1; ns<n+8; ns++) {
                    if (tree[one].sub[ns%8] >= FIRST_NODE) {
                        opencell[new_oe] = tree[one].sub[ns%8];
                        new_oe ++ ;
                    }
                    else if (tree[one].sub[ns%8] > -1) {
                        accept[new_apt] = tree[one].sub[ns%8];
                        new_apt++;
                    }
                }
                walk_cell(id, level+1, new_os, new_oe, new_apt);
            }
    }

    if ( os < oe && one >= FIRST_NODE) {
        for (n=0; n<8; n++)
            if ( (id = tree[one].sub[n]) > -1 ) {
                new_os = new_oe = oe+1;
                new_apt = apt;
                for (ns=n+1; ns<n+8; ns++) {
                    if (tree[one].sub[ns%8] >= FIRST_NODE) {
                        opencell[new_oe] = tree[one].sub[ns%8];
                        new_oe ++ ;
                    }
                    else if (tree[one].sub[ns%8] > -1) {
                        accept[new_apt] = tree[one].sub[ns%8];
                        new_apt++;
                    }
                }
                for (ns=os; ns<oe; ns++) {
                    cell = opencell[ns] ;
                    if ( accepted_cell(id, cell) ) {
                        accept[new_apt] = cell;
                        new_apt++;
                    }
                    else {
                        for (ns2=0; ns2<8; ns2++) {
                            if (tree[cell].sub[ns2] >= FIRST_NODE) {
                                opencell[new_oe] = tree[cell].sub[ns2];
                                new_oe ++ ;
                            }
                            else if (tree[cell].sub[ns2] > -1 ) {
                                accept[new_apt] = tree[cell].sub[ns2];
                                new_apt++;
                            }
                        }
                    }
                }
                walk_cell(id, level+1, new_os, new_oe, new_apt);
            }
    }
}

void walk_tree_for_gravity(int nPart)
{
    clock_t begin, end;
    double time_spent;
    begin = clock();
    opencell = (int*)malloc((int)sizeof(int)*2.5*nPart);
    accept   = (int*)malloc((int)sizeof(int)*2.5*nPart);
    opencell[0] = FIRST_NODE;
    walk_cell(FIRST_NODE, 0, 0, 0, 0);
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf(" recursive walk [%d] time = %lf[sec]\n", nPart, time_spent);
}

void naive_walk_cell(int add, int level, int body)
{
    int id, n, d;
    double dx[3], dr3, dr2;
    if (add == body)
        return;

    if (add<FIRST_NODE) {
        dr2 = 0.0;
        for (d=0; d<DIM; d++) {
            dx[d] = part[add].pos[d] - part[body].pos[d];
            dr2 += dx[d]*dx[d];
        }
        dr3 = dr2 * sqrt(dr2);
        for (d=0; d<DIM; d++) {
            part[body].acc[d] += dx[d]/dr3;
        }
    }
    else
        for (n=0; n<8; n++) {
            id = tree[add].sub[n];
            if (id>-1) {
                if ( id >= FIRST_NODE ) {
                    if (accepted_cell(body, id)) {
                        dr2 = 0.0;
                        for (d=0; d<DIM; d++) {
                            dx[d] = tree[id].masscenter[d]-part[body].pos[d];
                            dr2 += dx[d]*dx[d];
                        }
                        dr3 = dr2 * sqrt(dr2);
                        for (d=0; d<DIM; d++) {
                            part[body].acc[d] += (tree[id].nPart) * dx[d]/dr3;
                        }
                    }
                    else
                        naive_walk_cell(id, level+1, body);
                }
                else
                    naive_walk_cell(id, level+1, body);
            }
        }
}

void naive_walk_tree_for_gravity(int nPart)
{
    int d, n;
    clock_t begin, end;
    double time_spent;
    begin = clock();
    for (n=0; n<nPart; n++) {
        for (d=0; d<DIM; d++)
            part[n].acc[d] = 0.0;

        naive_walk_cell(FIRST_NODE, 0, n);
    }
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf(" naive_walk [%d] time = %lf[sec]\n", nPart, time_spent);
}

void walk_tree(int add, int level)
{
    int id, n;
    printf("[%d]%d included %d width = %f\n", level, add, tree[add].nPart, tree[add].width);
    for (n=0; n<8; n++) {
        id = tree[add].sub[n];

        if ( id > FIRST_NODE )  // nodes
        {
            walk_tree(id, level+1) ;
        }
        else if ( id >= 0) // particles
        {
            printf(" [%d.%d]-> %d : %f %f %f \n", add, n, id, part[id].pos[0], part[id].pos[1], part[id].pos[2]);
        }
    }
}

float ran3(long *idum);
void direct_addition(BODY *bptr, int nPart);

int main(void)
{
    int d, n, nNode;
    int nPart;
    int level;
    long int seed;
    double center[3], box;
    clock_t begin, end;
    double time_spent;
    FILE *fd;

    nPart = 20000000;
    box = 1.0;
    seed = 4658714;
    part = (BODY*)malloc( (int) sizeof (BODY)*nPart );
    for (n=0; n<nPart; n++) {
        for (d=0; d<DIM; d++) {
            part[n].pos[d] = box*ran3(&seed);
        }
    }

    for (n=0; n<nPart; n++) {
        part[n].phi = 0.0;
        for (d=0; d<DIM; d++) {
            part[n].phi += part[n].acc[d]* part[n].acc[d];
        }
        part[n].phi = sqrt(part[n].phi);
    }

    begin = clock();
    construct_tree(nPart, box) ;
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("construct tree [%d] time = %lf[sec]\n", nPart, time_spent);

    walk_tree_for_gravity(nPart);
    fd = fopen("nbody_recursive_treewalk.txt","w");
    fprintf(fd, "%d %lf\n", nPart, box);
    for (n=0; n<nPart; n++) {
        fprintf(fd, "%lf %lf %lf %e %e %e\n", (part[n].pos[0]), (part[n].pos[1]), (part[n].pos[2]), (part[n].acc[0]), (part[n].acc[1]), (part[n].acc[2]));
    }
    fclose(fd);

    naive_walk_tree_for_gravity(nPart) ;

    fd = fopen("nbody_iterative_treewalk.txt","w");
    fprintf(fd, "%d %lf\n", nPart, box);
    for (n=0; n<nPart; n++) {
        fprintf(fd, "%lf %lf %lf %e %e %e\n",(part[n].pos[0]), (part[n].pos[1]), (part[n].pos[2]), (part[n].acc[0]), (part[n].acc[1]), (part[n].acc[2]));
    }
    fclose(fd);

//   free(part);
    return 0;
}
*/

// Cao! Use this critertian for cell/cell acceptance judging. 
int accepted_cell_to_cell(int TA, int TB, double theta/* Here theta is the open_angle  */)
{
	double delta, dr;
	int d;

//# Strict substraction of R for correct acceptance:
//	double Rmax = SQROOT3 * (pTA->bmax + pTB->bmax); 
//#	Alternative subtraction of R:
	double Rmax = tree[TA]->width + tree[TB]->width;

	dr = 0.0;

#pragma unroll
	for (d=0; d<DIM; d++)
	{
		delta = tree[TB]->masscenter[d] - tree[TA]->masscenter[d];
		dr += delta*delta;
	}

	dr = SQRT(dr) - Rmax;
	if ( Rmax < theta * dr )
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

// Cao! Walk the tree using Dual Tree Traversel algorithm.
void dtt_process_cell(int TA, int TB, double theta)
{
	assert ((tree[TA]->level >=0) && (tree[TA]->level <= MAX_MORTON_LEVEL));
	assert ((tree[TB]->level >=0) && (tree[TB]->level <= MAX_MORTON_LEVEL));

	if(TA == TB)
	{
		if (tree[TA]->childnum == 0)
		{
			/* Cao! Here we send the tree node, and the Queue      */
			/*      processing function should pack the particles  */
			/*      and do some further optimizatoins.             */
			/*      PQ_ppnode: node pp packs, !0 means accepted    */
			/*      PQ_particle: save the packed pps.              */
			EnqueueP_PPnode(PQ_ppnode, TA, TB, 0);
		}
		else
		{
			int i;
			for (i=0; i<tree[TA]->childnum; i++)
			{
				EnqueueP_Cell(PQ_treewalk, tree[TA]->firstchild+i, TB, theta);
			}
		}
	}
	else
	{
		if (accepted_cell_to_cell(TA, TB, theta))
		{
			/* Cao! See comments above. */
			EnqueueP_PPnode(PQ_ppnode, TA, TB, 1);
		}
		else
		{
			if ( tree[TA]->chlidnum==0 && tree[TB]->childnum==0 )
			{
				/* Cao! See comments above. */
				EnqueueP_PPnode(PQ_ppnode, TA, TB, 0);
			}
			else if ( tree[TB]->childnum==0 || (tree[TA]->childnum>0 && tree[TA]->width>=tree[TB]->width) )
			{
				int i;
				for (i=0; i<tree[TA]->childnum; i++)
				{
					EnqueueP_Cell(PQ_treewalk, tree[TA]->firstchild+i, TB, theta);
				}
			}
			else
			{
				int i;
				for (i=0; i<tree[TB]->childnum; i++)
				{
					EnqueueP_Cell(PQ_treewalk, TA, tree[TB]->firstchild+i, theta);
				}
			}
		}
	}
}

