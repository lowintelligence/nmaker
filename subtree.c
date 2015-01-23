#include "global.h"
#include "domain.h"
#include "subtree.h"
#include <assert.h>

static double factor_part_node = 0.3;
// new space to build tree as a uint to sent to the other domain or node
// build top-subtree
// the tree must be smoothed during tree-traveral
// divide part int NumThread to

void check_subtree(Domain *dp, Constants *constants)
{
	int In, Jn, Kn, Ngrid;
	Int i, j, k;
	long int incorrect;
	Int n1, n2;
	Body *part = dp->part;
	Node *tree = dp->tree;

	In = dp->nSide[0];
	Jn = dp->nSide[1];
	Kn = dp->nSide[2];

	Ngrid = In*Jn*Kn;
	incorrect = 0;
	n1 = 0;
	n2 = 0;
	for (i=0; i<Ngrid; i++)
	{
		int rooti = dp->mroot[i];
		int n;

		// Check empty mesh.
		if (dp->count[i]==0)
		{
			if(dp->pidx[i]!=-1)
			{
				incorrect += 0x1;
			}

			if(tree[rooti].nPart!=0 || tree[rooti].nChildren!=0 || tree[rooti].firstChild!=-1 || tree[rooti].firstPart!=-1)
			{
				incorrect += 0x100000000L;
			}

			continue;
		}

		// Check first particle position.
		if ((dp->tag[i]>=0 && dp->pidx[i]!=dp->npart+n2) || (dp->tag[i]<0 && dp->pidx[i]!=n1) || dp->pidx[i]!=tree[rooti].firstPart)
		{
			incorrect += 0x100L;
			DBG_INFO(DBG_MEMORY, "i=%ld, partindex=%ld, numpart=%ld, treefirstPart=%ld.\n", i, dp->pidx[i],  dp->count[i], tree[rooti].firstPart);
		}
		n = dp->tag[i]<0 ? n1 : dp->npart+n2;
		// Check the mesh group id.
		if (part[n].group!=i)
		{
			incorrect += 0x10000;
		}

		// Check the particle number.
		Int ncount = 0;
		Int npart = dp->npart+dp->npart_adj;
		if (dp->tag[i]<0)
		{
			for(; n1<npart && part[n1].group==i; n1++)
			{
				ncount ++;
			}
		}
		else
		{
			for(; n2<dp->npart_adj && part[dp->npart+n2].group==i; n2++)
			{
				ncount ++;
			}
		}

		if (ncount != dp->count[i] || tree[rooti].nPart != dp->count[i])
		{
			incorrect += 0x1000000L;
		}

		// Check subtree.
		// Broad travel.
		int js=tree[rooti].firstChild;
		if(rooti!=js-1)
		{
			incorrect += 0x10000000000L;
//			printf("i=%ld, numparts=%ld, mroot=%ld, rpart=%ld, firstChild=%ld, childpart=%ld.\n", i, dp->count[i],  dp->mroot[i], tree[dp->mroot[i]].firstPart, tree[rooti].firstChild, tree[tree[rooti].firstChild].firstPart);
		}
		int je=tree[rooti].firstChild+tree[rooti].nChildren;
		int level=1;
		int j;
		int skipnp=0;
		int firstnl=js; // First tree node of next level.
		while (firstnl!=-1 && level<=constants->MAX_TREE_LEVEL)
		{
			int numnl = 0; // Next level tree nodes.
			int numl = skipnp; // This level particles in the tree.
			firstnl = -1; 
			if(tree[js-1].level!=level-1)
			{
				incorrect += 0x10000000000L;
			}
			for (j=js; j<je; j++)
			{
				if(tree[j].level!=level)
				{
					incorrect += 0x10000000000L;
//					printf("i=%ld, j=%ld, je=%ld, js=%ld, root=%ld, rnumchild=%d, rnump=%d, rp=%d, ilevel=%ld, jlevel=%ld, jnpart=%d, jp=%d, meshnp=%d, meshp=%d.\n", i, j, je, js, rooti, tree[rooti].nChildren, tree[rooti].nPart, tree[rooti].firstPart, level, tree[j].level, tree[j].nPart, tree[j].firstPart, dp->count[i], dp->pidx[i]);
				}
				if(tree[j].nPart==0 || (tree[j].nChildren==0 && tree[j].nPart>constants->MAX_PACKAGE_SIZE))
				{
					incorrect += 0x1000000L;
				}
				numl+=tree[j].nPart;
				if(tree[j].firstPart<dp->pidx[i] || tree[j].firstPart+tree[j].nPart>dp->pidx[i]+dp->count[i]) // Particle pointer error.
				{
//					printf("i=%ld, j=%ld, je=%ld, js=%ld, ilevel=%d, jlevel=%d, jnpart=%d, jp=%d, meshnp=%d, meshp=%d.\n", i, j, je, js, level, tree[j].level, tree[j].nPart, tree[j].firstPart, dp->count[i], dp->pidx[i]);
//					sleep(1);
					incorrect += 0x100;
				}
				if(tree[j].nChildren>0)
				{
					numnl+=tree[j].nChildren;
					if(firstnl==-1) firstnl = tree[j].firstChild;
				}
				else
				{
					skipnp+=tree[j].nPart;
				}
			}
			if(numl!=dp->count[i])
			{
				incorrect += 0x1000000L;
			}
			if(firstnl!=-1)
			{
				if(j!=firstnl)
				{
					incorrect += 0x10000000000L;
				}
				js = firstnl;
				je = firstnl+numnl;
			}
			level++;
		}
		if(firstnl!=-1 || (i<Ngrid-1 && j>(dp->mroot[i+1])))
		{
			incorrect += 0x1000000000000L;
		}
	}
	
	DBG_INFO(DBG_TMP, "[%d], Check subtree, %ld/%ld meshes, %ld/%d particles, incorrect code = %lx.\n", dp->rank, i, Ngrid, n1+n2, dp->npart+dp->npart_adj, incorrect);
}

Morton mortonEncode(Real x, Real y, Real z, int Bits, double meshsize){
	Morton xi, yi, zi, X, Y, Z, code;
	int i;
	double size, nmesh_meshsize, nbox_boxsize;
	nbox_boxsize = 1.0/meshsize;
    size = meshsize/(1<<Bits);//each grid cover how much kpc
    nmesh_meshsize = 1.0/size;//1kpc need how many grids
	
	xi = (Morton) ( x * nbox_boxsize );
    yi = (Morton) ( y * nbox_boxsize );
    zi = (Morton) ( z * nbox_boxsize );

	x = x - meshsize*xi;
	y = y - meshsize*yi;
	z = z - meshsize*zi;

	X = (Morton) ( x * nmesh_meshsize );
    Y = (Morton) ( y * nmesh_meshsize );
    Z = (Morton) ( z * nmesh_meshsize );
	
	code = 0;
	for(i=0;i<Bits;i++)
		code|=(((X&(1<<i))<<2*i)|((Y&(1<<i))<<2*i+1)|((Z&(1<<i))<<2*i+2));
	return code;
}

void exchange(Body* part, Int i, Int j)
{
	Body temp;
	temp    = part[i];
	part[i] = part[j];
	part[j] = temp;	
}

void bsort2p(Body *part, Domain *dp)
{

	Int i, j;
	Body temp,temp1;
	int meshid;
	Int nPart = dp->npart+dp->npart_adj;
	int *count = dp->count;
	int *pidx = dp->pidx;
	int Ngrid = dp->nSide[0]*dp->nSide[1]*dp->nSide[2];

	Int *ppart=(Int*)malloc(sizeof(Int)*Ngrid);

	for(i=0;i<Ngrid;i++)
	{
		ppart[i] = pidx[i];
	}

	for(i=0;i<nPart;)
	{
		meshid = part[i].group;
		
		if(i>=pidx[meshid] && i<pidx[meshid]+count[meshid])
		{
			i++;
//		if(i%10000 == 0)
//			DBG_INFO(DBG_LOGIC, "[%d]: bucked, i=%d.\n", dp->rank, i);
		}
		else
		{
			j = ppart[meshid];//current pointer for each mesh	
			temp = part[j];
			part[j] = part[i];
			ppart[meshid]++;	
			while(1){
				meshid = temp.group;
				if(i>=pidx[meshid] && i<pidx[meshid]+count[meshid])
				{
					part[i] = temp;
					i++;
//		if(i%100 == 0)
//			DBG_INFO(DBG_LOGIC, "[%d]: bucked, i=%d.\n", dp->rank, i);
					break;
				}
				j = ppart[meshid];
				temp1 = part[j];
				ppart[meshid]++;
				part[j] = temp;
				temp = temp1;
			}
		}
	}

	free(ppart);		
}

void mortonSort(Body* part, Int numparts, int Bits)
{
	Int i, m;
	Morton j, k;
	Body temp, temp1;
	Morton ngrids = 1 << (3*Bits);
	Int *pparti = (Int*)xmalloc(sizeof(Int)*ngrids, 7007);
	Int *ppartj = (Int*)xmalloc(sizeof(Int)*ngrids, 7007);
	Int *numpartis = (Int*)xmalloc(sizeof(Int)*ngrids, 7007);

	for(j=0;j<ngrids;j++)
	{
		pparti[j] = -1;
		ppartj[j] = -1;
		numpartis[j] = 0;
	}

	for(i=0;i<numparts;i++)
	{
		numpartis[part[i].mkey]++;	
	}

	m=0;
	for(j=0;j<ngrids;j++)
	{
		if(numpartis[j]>0)
		{
			pparti[j] = m;
			ppartj[j] = m;
			m += numpartis[j];
		}
	}

	for(i=0;i<numparts;)
	{
		Morton mtid = part[i].mkey;
		m = ppartj[mtid];
		if(i>=m && i<m+numpartis[mtid])
		{
			i++;
			continue;
		}
		else
		{
			j = pparti[mtid];//current pointer for each mesh	
			temp = part[j];
			part[j] = part[i];
			pparti[mtid]++;	
			while(1)
			{
				mtid = temp.mkey;
				if(i>=ppartj[mtid] && i<ppartj[mtid]+numpartis[mtid])
				{
					part[i] = temp;
					i++;
					break;
				}
				j = pparti[mtid];
				temp1 = part[j];
				pparti[mtid]++;
				part[j] = temp;
				temp = temp1;
			}
		}
	}

	free(pparti);	
	free(ppartj);	
	free(numpartis);	
}

int CompKey(const void *p1, const void *p2)
{
	return ((Body*)p1)->mkey - ((Body*)p2)->mkey;
}

typedef struct
{
    int rank;
	int tid;
    int lower;
    int upper;
    int level;
    Int npart;
    Int totpart;
    double box;
	int maxlevel;
	int minlevel;
	int maxsize;
    Int *mroot;
	Int *pidx;
	Int *count;
	Body *part;
	Node *tree;
} TaskBuildSubTree;

// Multithreaded build morton based subtree kernel.
void* build_tree_morton(void *param)
{
    int m, n;
    int rank, tid;
    int first, last, length;
	int maxlevel, minlevel, maxsize;
    Int *root, *pidx, *count;
	double size;
	unsigned int levelshift;
	Morton curcode;

	TaskBuildSubTree *taskp = (TaskBuildSubTree*)param;
	Body *part = taskp->part;
	Node *tree = taskp->tree;

// Now first and last should be the mesh index.
// Here we should pass from the caller.
	rank = taskp->rank;
	tid = taskp->tid;
	first = taskp->lower;
	last = taskp->upper;
	root = taskp->mroot;  //???
	pidx = taskp->pidx; 
	count = taskp->count; 
	size = taskp->box/(1<<(taskp->level));
	maxlevel = taskp->maxlevel;
	minlevel = taskp->minlevel;
	maxsize = taskp->maxsize;
	double ttree = 0.0;
	double tstamp;
		
	for (m=first; m<last; m++)
	{
		int ipart = pidx[m];
		int mpart = count[m];
		Node *pnode, *prnode;

		for(n=0;n<mpart;n++)
			part[ipart+n].mkey = mortonEncode(part[ipart+n].pos[0], part[ipart+n].pos[1], part[ipart+n].pos[2], maxlevel, size);

		tstamp=dtime();
		//  if(maxlevel>5)
		if(mpart<(1<<(maxlevel*2+4)))  
		{
			qsort(part+ipart, mpart, sizeof(Body), CompKey);
			ttree += dtime() - tstamp;
//			printf("[%d] qsort %d of %d particles finished, time=%.4f.\n", tid, m, mpart, time_end-time_start);
		}
		else
		{
			mortonSort(part+ipart, mpart, maxlevel);
			ttree += dtime() - tstamp;
//			printf("[%d] bsort %d of %d particles finished, time=%.4f.\n", tid, m, mpart, time_end-time_start);
		}

		double sumcenterx, sumcentery, sumcenterz, invnpart;

		// Calculate the root node.
		sumcenterx = sumcentery = sumcenterz = 0;

		pnode = tree+root[m];
//		if (m==211) printf ("[%d]Meshloop: m=%d, tree=%lx, pnode=%lx, mpart=%d, first=%d, last=%d, ipart=%d, root=%d.\n", tid, m, tree, pnode, mpart, first, last, ipart, root[m]);
		
		pnode->nPart = 0;
		pnode->level = 0;
		pnode->firstPart = ipart;
		pnode->nChildren = 0;
		pnode->firstChild = -1;
		pnode->mortonKey = 0;
		pnode->width = size;
		pnode->mass = 0;

		if(mpart == 0) continue;

		for (n=0; n<mpart; n++)
		{ // Here we need to confirm particles have the same mass.
			pnode->mass += part[ipart+n].mass;
			sumcenterx += part[ipart+n].pos[0];
			sumcentery += part[ipart+n].pos[1];
			sumcenterz += part[ipart+n].pos[2];
		} //comp the massCenter
		pnode->nPart = mpart;
		invnpart = 1 / (Real) pnode->nPart;
		pnode->massCenter[0] = sumcenterx * invnpart;
		pnode->massCenter[1] = sumcentery * invnpart;
		pnode->massCenter[2] = sumcenterz * invnpart;
		Real drmax = (Real) 0.000026;
		for (n=0; n<pnode->nPart; n++)
		{
			Real dx, dy, dz, dr;
			dx = part[ipart+n].pos[0] - pnode->massCenter[0];
			dy = part[ipart+n].pos[1] - pnode->massCenter[1];
			dz = part[ipart+n].pos[2] - pnode->massCenter[2];
			dr = dx*dx + dy*dy + dz*dz;
			drmax = dr>drmax ? dr : drmax;
		}
		pnode->rmax = SQRT(drmax);
		pnode++;

		int level = 1;
		int endmask = 0;
		prnode = tree+root[m];

//        printf(" [%d] building the mesh(%d), particles in mesh = %d.\n", tid,m,mpart);

		while (!endmask && level<maxlevel/* Define the number please */)
		{
			int lpmask = mpart;
			endmask = 1;
			levelshift = (maxlevel-level)*DIM; //for computing the key on this level
//			if (m==211)	printf ("[%d]Newlevel, m=%d, prnode=%lx, root=%d, l=%d, 1p=%d, nP=%d, n=%d, ip=%d, mp=%d.\n", tid, m, prnode, root[m], prnode->level, prnode->firstPart, prnode->nPart, n, ipart, mpart);
			
			n = prnode->firstPart - ipart;
			// Find and put the first particle into the tree.
			if(level > minlevel/* Define */)
			{
				// If the package is small enough, 
				// We move to the next up-level node.
				while(n<mpart && prnode->nPart <= maxsize/* Define  */)
				{
//					if (m==211)	printf ("[%d]Bscan, m=%d, prnode=%lx, root=%d, l=%d, level=%d, 1p=%d, nP=%d, n=%d, ip=%d, mp=%d.\n", tid, m, prnode, root[m], prnode->level, level, prnode->firstPart, prnode->nPart, n, ipart, mpart);
					n = prnode->firstPart - ipart + prnode->nPart;
					prnode++;
				}
//				if (m==211)	printf ("[%d]Afterscan, m=%d, prnode=%lx, root=%d, l=%d, level=%d, 1p=%d, nP=%d, n=%d, ip=%d, mp=%d.\n", tid, m, prnode, root[m], prnode->level, level, prnode->firstPart, prnode->nPart, n, ipart, mpart);
			}

			if (n<mpart)
			{
//				if (ipart+n != prnode->firstPart) printf ("[%d]Err, m=%d, prnode=%lx, root=%d, l=%d, 1p=%d, nP=%d, n=%d, ip=%d, mp=%d.\n", tid, m, prnode, root[m], prnode->level, prnode->firstPart, prnode->nPart, n, ipart, mpart);
				n = prnode->firstPart - ipart;
				assert (ipart+n == prnode->firstPart);
				
				// We have a certain new node of this level.
				// So we cannot stop scanning at this level.
				endmask = 0;//if need a new tree node;
				curcode = part[ipart+n].mkey >> levelshift;

				prnode->nChildren++;
				prnode->firstChild = pnode - tree;//Cao!
				pnode->level = level;
				pnode->mortonKey = part[ipart+n].mkey & ((~0) << levelshift);
				pnode->nPart = 0;
				pnode->firstPart = ipart+n;
				pnode->firstChild = -1;
				pnode->nChildren = 0;
				pnode->width = prnode->width * 0.5f; 
				pnode->mass = 0;

				// Scan all the node to build this level tree stucture.
				while (1)
				{
					Int i, oldn = n;
					sumcenterx = sumcentery = sumcenterz = 0;
					while (n<mpart && ((part[ipart+n].mkey >> levelshift) == curcode))
					{ // Here we need to confirm particles have the same mass.
						pnode->nPart++;
						pnode->mass += part[ipart+n].mass;
						sumcenterx += part[ipart+n].pos[0];
						sumcentery += part[ipart+n].pos[1];
						sumcenterz += part[ipart+n].pos[2];
						n++;
					}
//					if (m == 211) printf ("[%d]Child, m=%d, pnode=%lx, l=%d, 1p=%d, nP=%d, n=%d, ip=%d, mp=%d.\n", tid, m, prnode, pnode->level, pnode->firstPart, pnode->nPart, n, ipart, mpart);
//					if (m==211)	printf ("[%d]ChildParent, m=%d, prnode=%lx, root=%d, l=%d, level=%d, 1p=%d, nP=%d, n=%d, ip=%d, mp=%d.\n", tid, m, prnode, root[m], prnode->level, level, prnode->firstPart, prnode->nPart, n, ipart, mpart);
					lpmask = n;
					invnpart = 1 / (Real) pnode->nPart;
					pnode->massCenter[0] = sumcenterx * invnpart;
					pnode->massCenter[1] = sumcentery * invnpart;
					pnode->massCenter[2] = sumcenterz * invnpart;
					Real drmax = (Real) 0.000026;
					for (i=0; i<pnode->nPart; i++)
					{
						Real dx, dy, dz, dr;
						dx = part[ipart+oldn+i].pos[0] - pnode->massCenter[0];
						dy = part[ipart+oldn+i].pos[1] - pnode->massCenter[1];
						dz = part[ipart+oldn+i].pos[2] - pnode->massCenter[2];
						dr = dx*dx + dy*dy + dz*dz;
						drmax = dr>drmax ? dr : drmax;
					}
					pnode->rmax = SQRT(drmax);
					pnode++;

					if ( n>=mpart ) 
					{
						// Here n>=mpart means we have scanned the last particle of this lever.
						// So we should jump to a deeper level scan.
						// prnode++ also changes the up-level.
						prnode++;
						break;
					}

					curcode = part[ipart+n].mkey >> levelshift;

					if ( (curcode^(prnode->mortonKey>>levelshift)) >> DIM )
					{
						// We finished all children of one node,
						// then shift to the next up-level node.
						prnode++;
						n = prnode->firstPart - ipart;
						if(level > minlevel)
						{
							// If the package is small enough, 
							// We move to the next up-level node.
							while(n<mpart && prnode->nPart <= maxsize)
							{
								n = prnode->firstPart - ipart + prnode->nPart;
								prnode++;
//								if (m==211)	printf ("[%d]Blockjump, m=%d, prnode=%lx, root=%d, l=%d, level=%d, 1p=%d, nP=%d, n=%d, ip=%d, mp=%d.\n", tid, m, prnode, root[m], prnode->level, level, prnode->firstPart, prnode->nPart, n, ipart, mpart);
							}
						}
						// If n>=mpart, we need to jump to a deeper level scan.
						// Otherwise, continue in this level.
						if (n>=mpart) break;
						n = prnode->firstPart - ipart;
						curcode = part[ipart+n].mkey >> levelshift;
						prnode->firstChild = pnode - tree;
					}	
						
					prnode->nChildren++;
					pnode->level = level;
					pnode->mortonKey = curcode << levelshift;
					pnode->nPart = 0;
					pnode->firstPart = ipart+n;
					pnode->firstChild = -1;
					pnode->nChildren = 0;
					pnode->width = prnode->width * 0.5f;
					pnode->mass = 0;
				}
			}
			if (lpmask < mpart) mpart = lpmask;
			level++;
		}
//		printf("[%d] mesh %d finished, level=%d.\n", tid, m, level);
	}
	DBG_INFO(DBG_TIME, "[%d-%d] bsort or qsort used %.4f seconds\n", rank, tid, ttree);
}


void build_subtrees(Domain *dp, Constants *constants)
{
    int i, j, k, In, Jn, Kn, n, d, m, q, Ngrid;
	Int accum_thread;
    Int *accum, *lower, *upper, *npart;
    Int nNode;

//    SubCuboid* sub = dp->cuboid;

	Int nThread = constants->NTHREAD_TREE;

    In = dp->nSide[0];
    Jn = dp->nSide[1];
    Kn = dp->nSide[2];

    Ngrid = In*Jn*Kn;

	double tstamp=dtime();
    DBG_INFO(DBG_LOGIC, "[%d] Before bucket sort for mesh.\n", dp->rank);
	bsort2p(dp->part, dp);//bucketsort by the mesh index. And get count and pidx.
	DBG_INFO(DBG_TIME, "[%d] Time bucketsort: %.4f\n", dp->rank, dtime()-tstamp);

    assert(nThread > 0);

    /* allocate space for subtree */
    nNode = (Int) ((double)(dp->npart+dp->npart_adj) * factor_part_node) + Ngrid*((1<<(constants->MIN_TREE_LEVEL*3) + 1));

    dp->tree = (Node*)xmalloc(sizeof(Node)*nNode, 7000);
	dp->nnode = nNode;

    for (n=0; n<nNode; n++) {
        dp->tree[n].nPart = 0;
        dp->tree[n].nChildren = 0;
        dp->tree[n].firstChild = -1;
        dp->tree[n].firstPart = -1;
    }
	DBG_INFO(DBG_LOGIC, "[%d] Tree initialization.\n", dp->rank);

	/* Init the root node of each grid. */
	dp->mroot[0] = 0;
    for (i=1; i<Ngrid; i++)
	{ 
        dp->mroot[i] = dp->mroot[i-1]+(int)((double)dp->count[i-1]*factor_part_node)+(1<<(constants->MIN_TREE_LEVEL*3))+1;
	}

	/* Build local trees on meshes with local particles. */
    accum = (int*)xmalloc(sizeof(int)*In*Jn*Kn, 7000);

    lower = (int*)xmalloc(sizeof(int)*nThread, 7000);
    upper = (int*)xmalloc(sizeof(int)*nThread, 7000);
    npart = (int*)xmalloc(sizeof(int)*nThread, 7000);

    accum_thread = (dp->npart+dp->npart_adj) / nThread;

    accum[0] = 0;
    for (n=1; n<Ngrid; n++) {
		accum[n] = accum[n-1] + dp->count[n];
    }

    lower[0] = 0;

    for (m=1, n=0; n<Ngrid; n++) {
        if ( accum[n] > m*accum_thread ) {
            lower[m] = upper[m-1] = n;
            m++;
        }
    }
    upper[nThread-1]=Ngrid;

    npart[0] = accum[upper[0]-1];
    for (m=1; m<nThread; m++) {
        npart[m] = accum[upper[m]-1] - accum[lower[m]-1];
    }
    free(accum);

    for (m=0; m<nThread; m++) {
        DBG_INFO(DBG_DATA, "[%d-%d]:lower=%d upper=%d npart=%d\n",dp->rank,m,lower[m],upper[m],npart[m]);
    }

    /* build subtree */
    pthread_t *task;
    TaskBuildSubTree *buildsub;

    task = (pthread_t*)xmalloc(sizeof(pthread_t)*nThread, 7004);
    buildsub = (TaskBuildSubTree*)xmalloc(sizeof(TaskBuildSubTree)*nThread, 7004);
	

	/* Build multi-local-trees with multithread. */
    int level;
    for (q=0; q<nThread; q++)
	{
        buildsub[q].rank  = dp->rank;
        buildsub[q].tid   = q;
        buildsub[q].lower = lower[q];
        buildsub[q].upper = upper[q];
        buildsub[q].npart = npart[q];
        buildsub[q].totpart = dp->npart + dp->npart_adj;
        buildsub[q].mroot = dp->mroot; //Root nodes
        buildsub[q].level = constants->MESH_BITS;
        buildsub[q].box   = constants->BOX_SIZE;
		buildsub[q].maxlevel  = constants->MAX_TREE_LEVEL;
		buildsub[q].minlevel  = constants->MIN_TREE_LEVEL;
		buildsub[q].maxsize = constants->MAX_PACKAGE_SIZE;
		buildsub[q].pidx  = dp->pidx;
		buildsub[q].count = dp->count;
		buildsub[q].part = dp->part;
		buildsub[q].tree = dp->tree;
    }
	DBG_INFO(DBG_LOGIC, "[%d]: Before subtree pthread.\n", dp->rank);
    for (q=0; q<nThread; q++) {
        pthread_create(&task[q], NULL, build_tree_morton, &buildsub[q]);
    }

    for (q=0; q<nThread; q++) {
        pthread_join(task[q], NULL);
    }
	DBG_INFO(DBG_LOGIC, "[%d]: After subtree pthread.\n", dp->rank);
	free(lower);
    free(upper);
    free(npart);

    free(task);
    free(buildsub);

	check_subtree(dp, constants);
}
