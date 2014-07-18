#include "subtree.h"
#include <assert.h>

static Body *part;
static Node *tree;


static double frac_part_node = 0.7;
// new space to build tree as a uint to sent to the other domain or node
// build top-subtree
// the tree must be smoothed during tree-traveral
// divide part int NumThread to




typedef struct {

} ExternTree;


typedef struct {
    int **buff; // need to count the number of particle to the other domain.
} ExchangeFrontier;


void free_domaintree(DomainTree *dtp)
{
    if ( NULL != dtp) {
        if ( NULL != dtp->root_tree)
            free(dtp->root_tree);
        if ( NULL != dtp->root_cell)
            free(dtp->root_cell);

        free(dtp);
    }
}


void insert_particle_into_tree(vect3d center, double size, int iNode,int iPart, int first, int *plast)
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
    if ( tree[iNode].sub[index] < first ) {
        old = tree[iNode].sub[index]; // old particle;
        (*plast) ++;
        tree[iNode].sub[index] = (*plast) ; // new node
        insert_particle_into_tree(center, size/2, (*plast), old, first, plast);
    }
    // this is a node
    insert_particle_into_tree(center, size/2, tree[iNode].sub[index], iPart, first, plast);
}


typedef struct {
    int rank;
    int lower;
    int upper;
    int npart;
    int first;
    int length;
    int totpart;
    int level;
    int min[3];
    double box;
//    Body *part;
//    Tree *tree;
    int  *root_cell;
} TaskBuildSubTree;



void* mt_build_subtree(void *param) {
    int n, x, y, z, idx;
    int lower, upper, rank, domainpart;
    int first, last, length;
    double  size;
    double nbox_boxsize;
    int *root;
    Node *cells;
//   Body *part;
    TaskBuildSubTree *taskp = (TaskBuildSubTree*)param;
    vect3d center;

    lower = taskp->lower;
    upper = taskp->upper;
    rank  = taskp->rank;
//   cells = sub->tree->cell;
    first = last = taskp->first;
    length= taskp->length;
    domainpart = taskp->totpart;
//    part =  sub->part;
    size =  taskp->box/(1<<(taskp->level) );
    nbox_boxsize = 1.0/size;
    root = taskp->root_cell;

//    printf(" [%d] npart = %d, first = %d\n", rank, domainpart, first);

    for (n=0; n<domainpart; n++) {
        if ( last-first>=length ) {
            printf("error\n");
            exit(0);
        }
        idx = part[n].group;
        if ( idx>=lower && idx<upper ) {
            // insert_a_particle
            x = (int) ( part[n].pos[0] * nbox_boxsize );
            y = (int) ( part[n].pos[1] * nbox_boxsize );
            z = (int) ( part[n].pos[2] * nbox_boxsize );

            center[0] = (0.5+x)*size;
            center[1] = (0.5+y)*size;
            center[2] = (0.5+z)*size;
            if ( -1 == root[idx] )
                root[idx] = last;

            insert_particle_into_tree(center, size, root[idx], n, first, &last);

        }
    }
    printf(" (%d) first = %d , last = %d, frac=%.3f\n",rank, first, last, (double)(last-first)/length);

}


void build_subtree_on_subcuboid(Domain *dp, GlobalParam *gp, int nThread) {
    int i, j, k, In, Jn, Kn, n, d, m, q, Ngrid, accum_thread;
    int *accum, *lower, *upper, *npart;
    int node_length;
    SubCuboid* sub = dp->cuboid;
    In = sub->nSide[0];
    Jn = sub->nSide[1];
    Kn = sub->nSide[2];
    Ngrid = In*Jn*Kn;
    assert(nThread > 0);
    accum_thread = (dp->NumPart) / nThread;

    accum = (int*)malloc( sizeof(int)*In*Jn*Kn);

    accum[0] =  sub->count[0];
    for (n=1; n<Ngrid; n++) {
        accum[n] = accum[n-1] + sub->count[n];
    }
    lower = (int*)malloc(sizeof(int)*nThread);
    upper = (int*)malloc(sizeof(int)*nThread);
    npart = (int*)malloc(sizeof(int)*nThread);

    lower[0] = 0;
    upper[nThread-1]=Ngrid;

    for (m=1, n=0; n<Ngrid; n++) {
        if ( accum[n] > m*accum_thread ) {
            lower[m] = upper[m-1] = n;
            m++;
        }
    }

    npart[0] = accum[upper[0]-1];
    for (m=1; m<nThread; m++) {
        npart[m] = accum[upper[m]-1] - accum[lower[m]-1];
    }
    free(accum);

    for (m=0; m<nThread; m++) {
        printf("%d:lower=%d upper=%d npart=%d\n",m,lower[m],upper[m],npart[m]);
    }

    /* construct domaintree */
    int index_length = In * Jn * Kn;
    int first_node = dp->NumPart;
    dp->domtree = (DomainTree*)malloc(sizeof(DomainTree));
    DomainTree *dtp = dp->domtree;
    dtp->NumTree = nThread;

    dtp->root_tree = (int*)malloc(sizeof(int)*index_length );
    dtp->root_cell = (int*)malloc(sizeof(int)*index_length );

    /* allocate space for subtree */
    node_length = (int) ( (dp->NumPart) * frac_part_node );

    dtp->tree = (Node*)malloc(sizeof(Node)*node_length);

    for (n=0; n<node_length; n++) {
        dtp->tree[n].nPart = 0;
        for (d=0; d<8; d++)
            dtp->tree[n].sub[d] = -1;
    }
    dtp->tree -= first_node ;
    tree = (dtp->tree );

    /*
        first_node = (dp->NumPart);
        for (q=0; q<nThread; q++) {
            dtp->tree[q].first_node = first_node ;
            dtp->tree[q].last_node  = first_node ;
            node_length = (int) ( npart[q] * frac_part_node );
            dtp->tree[q].node_length = node_length;

            dtp->tree[q].cell = (dtp->gtree - first_node);
       //     printf();
            first_node += node_length;
        }
    //  dtp->gtree -= dp->NumPart;
        */

    for (i=0, m=0; i<index_length; i++) {
        dtp->root_tree[i] = m;
        dtp->root_cell[i] = -1;
        if (i==upper[m])
            m++;
    }

    /* build subtree */
    pthread_t *task;
    TaskBuildSubTree *buildsub;

    task = (pthread_t*)malloc(sizeof(pthread_t)*nThread);
    buildsub = (TaskBuildSubTree*)malloc(sizeof(TaskBuildSubTree)*nThread);

    /* set global variable for particles */
    part = dp->Part;

    first_node = (dp->NumPart);
    for (q=0; q<nThread; q++) {
        buildsub[q].rank  = q;
        buildsub[q].lower = lower[q];
        buildsub[q].upper = upper[q];
        buildsub[q].npart = npart[q];
        buildsub[q].first = first_node ;
        buildsub[q].length= (int) ( npart[q] * frac_part_node );
        buildsub[q].totpart = dp->NumPart;
//        buildsub[q].tree  = &(dtp->tree[q]);
        buildsub[q].root_cell = dtp->root_cell;
//       buildsub[q].part  = dp->Part;
        buildsub[q].level = gp->NumBits;
        buildsub[q].box   = gp->BoxSize;
        buildsub[q].min[0]= sub->minimum[0];
        buildsub[q].min[1]= sub->minimum[1];
        buildsub[q].min[2]= sub->minimum[2];

        first_node +=  buildsub[q].length ;
    }

    for (q=0; q<nThread; q++) {
        pthread_create(&task[q], NULL, mt_build_subtree, &buildsub[q]);
    }

    for (q=0; q<nThread; q++) {
        pthread_join(task[q], NULL);
    }
    free(lower);
    free(upper);
    free(npart);

    free(task);
    free(buildsub);

}

void packing_group(DomainTree *dtp) {

}


// Cao! Multithreaded build morton based subtree kernel.
void mt_build_subtree_morton(void *param) {
    int m, n;
    int rank;
    int first, last, length;
    int *root;
	double size;

	TaskBuildSubTree *taskp = (TaskBuildSubTree*)param;

// Now first and last should be the mesh index.
// Here we should pass from the caller.
	rank  = taskp->rank;
	first = taskp->first;
	last = taskp->last;
	root = taskp->root_cell;
	parts = taskp->partidxes; // This should be added in data.h to indicate the particle position.
	numpart = taskp->numparticles; // Number of particles in meshes, should be added in data.h.
	size = taskp->box/(1<<(taskp->level));
	
	for (m=first; m<last; m++)
	{
		int ipart = parts[m];
		int mpart = numpart[m];

		float sumcenterx, sumcentery, sumcenterz, invnpart;

		// Calculate the root node.
		sumcenterx = sumcentery = sumcenterz = 0;

		Node* pnode = &tree[root[m]];

		pnode->level = 0;
		pnode->firstpart = ipart;
		pnode->childnum = 0;
		pnode->firstchild = -1;
		pnode->mortonkey = 0;
		pnode->width = size;
		pnode->mass = 0;

		if(mpart == 0) continue;

		for (n=0; n<mpart; n++)
		{
			pnode->nPart++;
			pnode->mass += part[ipart+n].mass;
			sumcenterx += part[ipart+n].pos[0];
			sumcentery += part[ipart+n].pos[1];
			sumcenterz += part[ipart+n].pos[2];
		}
		invnpart = 1.0f / (float)pnode->nPart;
		pnode->masscenter[0] = sumcenterx * invnpart;
		pnode->masscenter[1] = sumcentery * invnpart;
		pnode->masscenter[2] = sumcenterz * invnpart;
		pnode++;

		int level = 1;
		int endmask = 0;
		Node* prnode = &tree[root[m]];

        printf(" [%d] building the mesh(%d), particles in mesh = %d.\n", rank, mpart);

		while (!endmask || level>MAX_MORTON_LEVEL/* Define the number please */)
		{
			endmask = 1;
			levelshift = (MAX_MORTON_LEVEL-level)*3;
			
			n = 0;
			// Find and put the first particle into the tree.
			if(level > MIN_TREE_LEVEL/* Define */)
			{
				// If the package is small enough, 
				// We move to the next up-level node.
				while(n<mprt && prnode->nPart <= MIN_PACKAGE_SIZE/* Define  */)
				{
					n += prnode->nPart;
					prnode++;
				}
			}
			if (n<mprt)
			{
				assert (ipart+n == prnode->firstpart);
				// We have a certain new node of this level.
				// So we cannot stop scanning at this level.
				endmask = 0;
				curcode = part[ipart+n].mortonkey >> levelshift;

				prnode->childnum++;
				prnode->firstchild = pnode - root[m];
				pnode->level = level;
				pnode->mortonkey = part[ipart+n].mortonkey & ((~0) << levelshift);
				pnode->nPart = 1;
				pnode->firstpart = ipart+n;
				pnode->firstchild = -1;
				pnode->childnum = 0;
				pnode->width = prnode->width * 0.5f; 
				pnode->mass = 0;

				// Scan all the node to build this level tree stucture.
				while (1)
				{
					sumcenterx = sumcentery = sumcenterz = 0;
					while (n<mpart && ((part[ipart+n].mortonkey >> levershift) == curcode))
					{
						pnode->nPart++;
						pnode->mass += part[ipart+n].mass;
						sumcenterx += part[ipart+n].pos[0];
						sumcentery += part[ipart+n].pos[1];
						sumcenterz += part[ipart+n].pos[2];
						n++;
					}
					invnpart = 1.0f / (float)pnode->nPart;
					pnode->masscenter[0] = sumcenterx * invnpart;
					pnode->masscenter[1] = sumcentery * invnpart;
					pnode->masscenter[2] = sumcenterz * invnpart;
					pnode++;

					if ( n>=mpart ) 
					{
						// Here n>=mpart means we have scanned the last particle of this lever.
						// So we should jump to a deeper level scan.
						// prnode++ also changes the up-level.
						prnode++;
						break;
					}

					curcode = part[ipart+n].mortonkey >> levelshift;

					if ( (curcode^(prnode->mortonkey>>levelshift)) >> 3 )
					{
						// We finished all children of one node,
						// then shift to the next up-level node.
						prnode++;
						if(level > MIN_TREE_LEVEL)
						{
							// If the package is small enough, 
							// We move to the next up-level node.
							while(n<mpart && prnode->pNode <= MIN_PACKAGE_SIZE)
							{
								n += prnode->nPart;
								prnode++;
							}
						}
						// If n>=mpart, we need to jump to a deeper level scan.
						// Otherwise, continue in this level.
						if (n>=mpart) break;
						curcode = part[ipart+n].mortonkey >> levelshift;
					}	
						
					prnode->firstchild = pnode - root[m];
					prnode->childnum++;
					pnode->level = level;
					pnode->mortonkey = curcode << levelshift;
					pnode->nPart = 1;
					pnode->firstpart = ipart+n;
					pnode->firstchild = -1;
					pnode->childnum = 0;
					pnode->width = prnode->width * 0.5f;
					pnode->mass = 0;
				}
			}
			level++;
		}
	}
}
