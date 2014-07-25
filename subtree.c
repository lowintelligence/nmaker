#include "subtree.h"
#include <assert.h>

static Body *part;
static Node *tree;


static double frac_part_node = 0.5;
// new space to build tree as a uint to sent to the other domain or node
// build top-subtree
// the tree must be smoothed during tree-traveral
// divide part int NumThread to




typedef struct {

} ExternTree;


typedef struct {
    int **buff; // need to count the number of particle to the other domain.
} ExchangeFrontier;


int mortonEncode(float x,float y,float z,double meshsize){
	int xi,yi,zi,X,Y,Z,code;
	int i;
	double size,nmesh_meshsize,nbox_boxsize;
	int Bits=MAX_MORTON_LEVEL;
	nbox_boxsize=1.0/meshsize;
    size =  meshsize/(1<<Bits);//each grid cover how much kpc
    nmesh_meshsize = 1.0/size;//1kpc need how many grids
	
	xi = (int) ( x * nbox_boxsize );
    yi = (int) ( y * nbox_boxsize );
    zi = (int) ( z * nbox_boxsize );
	x=x-meshsize*xi;
	y=y-meshsize*yi;
	z=z-meshsize*zi;

	X = (int) ( x * nmesh_meshsize );
    Y = (int) ( y * nmesh_meshsize );
    Z = (int) ( z * nmesh_meshsize );
	
	code=0;
	for(i=0;i<Bits;i++)
		code|=(((X&(1<<i))<<2*i)|((Y&(1<<i))<<2*i+1)|((Z&(1<<i))<<2*i+2));
	return code;
}

void exchange(Body* part,int i,int j){
	Body temp=part[i];
	part[i]=part[j];
	part[j]=temp;	
}

void bucketsort(Body* part,int numparts,int* numparticles,int* partidxes,int Ngrid){
	int i,j,k,m;
	int* ppart=(int*)malloc(sizeof(int)*Ngrid);
	for(i=0;i<numparts;i++){
		numparticles[part[i].group]++;	
	}

	m=0;
	for(i=0;i<Ngrid;i++){
		ppart[i]=-1;
		partidxes[i]=-1;	
		if(numparticles[i]==0)
			continue;
		partidxes[i]=m;	
		ppart[i]=m;
		m+=numparticles[i];
	}
//	printf("540: partidx %d, ppart %d, numpart %d.\n", partidxes[540], ppart[540], numparticles[540]);

	for(i=0;i<numparts;i++){
		int meshid=part[i].group;
		m = partidxes[meshid];
		j=ppart[meshid];//current pointer for each mesh	
		if(j>=numparts)
		{
			printf("Part %d, ppart %d, group %d, meshid %d, meshnum %d, partidx %d.\n", i, ppart[meshid], part[i].group, meshid, numparticles[meshid], partidxes[meshid]);
		}
		if(meshid==540 && ppart[meshid]%100000==0) printf("i=%d, gi=%d, j=%d, gj=%d, ppart=%d\n ", i, part[i].group, j, part[j].group, ppart[meshid]);
		if(i>=m && i<m+numparticles[meshid])
			continue;
		else
		{
			exchange(part,i,j);
			ppart[meshid]++;	
		}
	}
	free(ppart);		
}

void mortonSort(Body* part,int numparts,int Bits){
	int i,j,k,m;
	int ngrids=1<<Bits;
	ngrids=ngrids*ngrids*ngrids;
	int* pparti=(int*)malloc(sizeof(int)*ngrids);
	int* ppartj=(int*)malloc(sizeof(int)*ngrids);
	int* numpartis=(int*)malloc(sizeof(int)*ngrids);
	for(i=0;i<ngrids;i++){
		pparti[i]=-1;
		ppartj[i]=-1;
		numpartis[i]=0;
	}
	for(i=0;i<numparts;i++){
		numpartis[part[i].mortonkey]++;	
	}

	m=0;
	for(i=0;i<ngrids;i++){
		if(numpartis[i]==0)
			continue;
		pparti[i]=m;
		ppartj[i]=m;
		m+=numpartis[i];
	}

	for(i=0;i<numparts;i++){
		int mtid=part[i].mortonkey;
		m=ppartj[mtid];
		j=pparti[mtid];//current pointer for each mesh	
		if(i>=m && i<m+numpartis[mtid])
			continue;
		else
		{
			exchange(part,i,j);
			pparti[mtid]++;	
		}
	}

	free(pparti);	
	free(numpartis);	
}

int CompKey(const void*p1,const void*p2){
	return ((Body*)p1)->mortonkey-((Body*)p2)->mortonkey;
}

void packing_group(Domain *dp) {
	int i;
	float j,k;
	Node* ptree=tree;
	int node_length=(dp->NumPart)*frac_part_node;
	for(i=dp->NumPart;i<node_length+dp->NumPart;i++)
		if(tree[i].nPart==2){
			//printf("center=%lf,%lfd,%lf,width=%f,npart=%d\n",tree[i].masscenter[0],tree[i].masscenter[1],tree[i].masscenter[2],tree[i].width,tree[i].nPart);
			k=100000/tree[i].width/16;
			if(k>j)
				j=k;
		}
			printf("width=%f\n",j);
	return;
}

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

/*  
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
*/

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

	int* partidxes;
	int* numparticles;

} TaskBuildSubTree;

// Cao! Multithreaded build morton based subtree kernel.
void* mt_build_subtree_morton(void *param) {
    int m, n;
    int rank;
    int first, last, length;
    int *root,*parts,*numpart;
	double size;
	unsigned int levelshift;
	int curcode;

	TaskBuildSubTree *taskp = (TaskBuildSubTree*)param;

// Now first and last should be the mesh index.
// Here we should pass from the caller.
	rank  = taskp->rank;
	first = taskp->lower;
	last = taskp->upper;
	root = taskp->root_cell;  //???
	parts = taskp->partidxes; 
	numpart = taskp->numparticles; 
	size = taskp->box/(1<<(taskp->level));
	double mengchen=0;
  clock_t time_start,time_end;
		
	for (m=first; m<last; m++)
	{
		int ipart = parts[m];
		int mpart = numpart[m];
		Node *pnode, *prnode;
	//!Meng
	    for(n=0;n<mpart;n++)
			part[ipart+n].mortonkey=mortonEncode(part[ipart+n].pos[0],part[ipart+n].pos[1],part[ipart+n].pos[2],size);
	
  time_start=clock();
  if(MAX_MORTON_LEVEL>5)
		qsort(part+ipart,mpart,sizeof(Body),CompKey);
  else
		mortonSort(part+ipart,mpart,MAX_MORTON_LEVEL);
  time_end=clock();
  mengchen+=(double)(time_end-time_start);

		float sumcenterx, sumcentery, sumcenterz, invnpart;

		// Calculate the root node.
		sumcenterx = sumcentery = sumcenterz = 0;

		pnode = tree+root[m];
//		if (m==211) printf ("[%d]Meshloop: m=%d, tree=%lx, pnode=%lx, mpart=%d, first=%d, last=%d, ipart=%d, root=%d.\n", rank, m, tree, pnode, mpart, first, last, ipart, root[m]);
		
		pnode->nPart = 0;
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
		} //comp the masscenter
		invnpart = 1.0f / (float)pnode->nPart;
		pnode->masscenter[0] = sumcenterx * invnpart;
		pnode->masscenter[1] = sumcentery * invnpart;
		pnode->masscenter[2] = sumcenterz * invnpart;
		pnode++;

		int level = 1;
		int endmask = 0;
		prnode = tree+root[m];

//        printf(" [%d] building the mesh(%d), particles in mesh = %d.\n", rank,m,mpart);

		while (!endmask && level<MAX_MORTON_LEVEL/* Define the number please */)
		{
			int lpmask = mpart;
			endmask = 1;
			levelshift = (MAX_MORTON_LEVEL-level)*DIM; //for computing the key on this level
//			if (m==211)	printf ("[%d]Newlevel, m=%d, prnode=%lx, root=%d, l=%d, 1p=%d, nP=%d, n=%d, ip=%d, mp=%d.\n", rank, m, prnode, root[m], prnode->level, prnode->firstpart, prnode->nPart, n, ipart, mpart);
			
			n = prnode->firstpart - ipart;
			// Find and put the first particle into the tree.
			if(level > MIN_TREE_LEVEL/* Define */)
			{
				// If the package is small enough, 
				// We move to the next up-level node.
				while(n<mpart && prnode->nPart <= MAX_PACKAGE_SIZE/* Define  */)
				{
//					if (m==211)	printf ("[%d]Bscan, m=%d, prnode=%lx, root=%d, l=%d, level=%d, 1p=%d, nP=%d, n=%d, ip=%d, mp=%d.\n", rank, m, prnode, root[m], prnode->level, level, prnode->firstpart, prnode->nPart, n, ipart, mpart);
					n = prnode->firstpart - ipart + prnode->nPart;
					prnode++;
				}
//				if (m==211)	printf ("[%d]Afterscan, m=%d, prnode=%lx, root=%d, l=%d, level=%d, 1p=%d, nP=%d, n=%d, ip=%d, mp=%d.\n", rank, m, prnode, root[m], prnode->level, level, prnode->firstpart, prnode->nPart, n, ipart, mpart);
			}

			if (n<mpart)
			{
//				if (ipart+n != prnode->firstpart) printf ("[%d]Err, m=%d, prnode=%lx, root=%d, l=%d, 1p=%d, nP=%d, n=%d, ip=%d, mp=%d.\n", rank, m, prnode, root[m], prnode->level, prnode->firstpart, prnode->nPart, n, ipart, mpart);
				n = prnode->firstpart - ipart;
				assert (ipart+n == prnode->firstpart);
				
				// We have a certain new node of this level.
				// So we cannot stop scanning at this level.
				endmask = 0;//if need a new tree node;
				curcode = part[ipart+n].mortonkey >> levelshift;

				prnode->childnum++;
				prnode->firstchild = pnode - tree;//Cao!
				pnode->level = level;
				pnode->mortonkey = part[ipart+n].mortonkey & ((~0) << levelshift);
				pnode->nPart = 0;
				pnode->firstpart = ipart+n;
				pnode->firstchild = -1;
				pnode->childnum = 0;
				pnode->width = prnode->width * 0.5f; 
				pnode->mass = 0;

				// Scan all the node to build this level tree stucture.
				while (1)
				{
					sumcenterx = sumcentery = sumcenterz = 0;
					while (n<mpart && ((part[ipart+n].mortonkey >> levelshift) == curcode))
					{
						pnode->nPart++;
						pnode->mass += part[ipart+n].mass;
						sumcenterx += part[ipart+n].pos[0];
						sumcentery += part[ipart+n].pos[1];
						sumcenterz += part[ipart+n].pos[2];
						n++;
					}
//					if (m == 211) printf ("[%d]Child, m=%d, pnode=%lx, l=%d, 1p=%d, nP=%d, n=%d, ip=%d, mp=%d.\n", rank, m, prnode, pnode->level, pnode->firstpart, pnode->nPart, n, ipart, mpart);
//					if (m==211)	printf ("[%d]ChildParent, m=%d, prnode=%lx, root=%d, l=%d, level=%d, 1p=%d, nP=%d, n=%d, ip=%d, mp=%d.\n", rank, m, prnode, root[m], prnode->level, level, prnode->firstpart, prnode->nPart, n, ipart, mpart);
					lpmask = n;
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

					if ( (curcode^(prnode->mortonkey>>levelshift)) >> DIM )
					{
						// We finished all children of one node,
						// then shift to the next up-level node.
						prnode++;
						n = prnode->firstpart - ipart;
						if(level > MIN_TREE_LEVEL)
						{
							// If the package is small enough, 
							// We move to the next up-level node.
							while(n<mpart && prnode->nPart <= MAX_PACKAGE_SIZE)
							{
								n = prnode->firstpart - ipart + prnode->nPart;
								prnode++;
//								if (m==211)	printf ("[%d]Blockjump, m=%d, prnode=%lx, root=%d, l=%d, level=%d, 1p=%d, nP=%d, n=%d, ip=%d, mp=%d.\n", rank, m, prnode, root[m], prnode->level, level, prnode->firstpart, prnode->nPart, n, ipart, mpart);
							}
						}
						// If n>=mpart, we need to jump to a deeper level scan.
						// Otherwise, continue in this level.
						if (n>=mpart) break;
						n = prnode->firstpart - ipart;
						curcode = part[ipart+n].mortonkey >> levelshift;
					}	
						
					prnode->firstchild = pnode - tree;
					prnode->childnum++;
					pnode->level = level;
					pnode->mortonkey = curcode << levelshift;
					pnode->nPart = 0;
					pnode->firstpart = ipart+n;
					pnode->firstchild = -1;
					pnode->childnum = 0;
					pnode->width = prnode->width * 0.5f;
					pnode->mass = 0;
				}
			}
			if (lpmask < mpart) mpart = lpmask;
			level++;
		}
		//printf("[%d] mesh %d finished, level=%d.\n", rank, m, level);
	}
	printf("mengchen qsort=%lf\n", mengchen/CLOCKS_PER_SEC);
}

/*  
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
/*	int code;
	code=mortonEncode(6249,6249,6249,size);
	printf("code=%o\n",code);
	return 0;*/

//    printf(" [%d] npart = %d, first = %d\n", rank, domainpart, first);
	
/*	int mybool=1;
		Body temp;
		temp.group=part[10].group;
		part[10].group=part[domainpart-1].group;
		part[domainpart-1].group=temp.group;##/
    for (n=0; n<domainpart; n++) {
        if ( last-first>=length ) {
            printf("error\n");
            exit(0);
        }
        idx = part[n].group;
        if ( idx>=lower && idx<upper ) {
/*            if(mybool==1){
				printf("shock!\n");
				mybool=0;
			}##/
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
/*		else
			mybool=1;##/
    }
    printf(" (%d) first = %d , last = %d, frac=%.3f\n",rank, first, last, (double)(last-first)/length);

}*/


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
    node_length = (int) ( (dp->NumPart) * frac_part_node + index_length);

    dtp->tree = (Node*)malloc(sizeof(Node)*node_length);

    for (n=0; n<node_length; n++) {
        dtp->tree[n].nPart = 0;
        dtp->tree[n].childnum = 0;
        dtp->tree[n].firstchild = -1;
        dtp->tree[n].firstpart = -1;
    }
//    dtp->tree -= first_node ;
    tree = (dtp->tree );
	printf("tree\n");

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
        dtp->root_tree[i] = m;//this grid is belong to which thread? 
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
	dtp->partidxes=(int*)malloc(Ngrid*sizeof(int));
	dtp->numparticles=(int*)malloc(Ngrid*sizeof(int));
	for(m=0;m<Ngrid;m++){
		dtp->partidxes[m]=-1;
		dtp->numparticles[m]=0;
	}
  
//  clock_t time_start,time_end;
//  time_start=clock();
//    printf("[] Before bucket sort for mesh.\n");
	bucketsort(part,dp->NumPart,dtp->numparticles,dtp->partidxes,Ngrid);//bucketsort by the mesh index. And get numparticles and partidxes.
  
//  time_end=clock();
//  printf("time bucketsort:%lf\n",(double)(time_end-time_start)/CLOCKS_PER_SEC);
//  sleep(10);	
	printf("tree1\n");

/* 	m=part[0].group;
	dtp->partidxes[m]=0;
	for(i=1;i<dp->NumPart;i++)
		if(part[i].group!=m){
			dtp->numparticles[m]=i-dtp->partidxes[m];
			m=part[i].group;
			dtp->partidxes[m]=i;
		}
	dtp->numparticles[m]=dp->NumPart-dtp->partidxes[m];*/
	
	dtp->root_cell[0] = 0;
    for (i=1; i<Ngrid; i++){ 
        dtp->root_cell[i] = dtp->root_cell[i-1]+dtp->numparticles[i-1]*frac_part_node+1;
	}

/*	printf("nSide=%d,%d,%d",dp->cuboid->nSide[0],dp->cuboid->nSide[1],dp->cuboid->nSide[2]);

	int mcc=0,myrank;
	char mc;
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	while(myrank==0){
		scanf("%c",&mc);
	printf("%d.................\n",part[mcc].group);
	mcc++;
	}*/



//    first_node = (dp->NumPart);
    first_node = 0;
    for (q=0; q<nThread; q++) {
        buildsub[q].rank  = q;
        buildsub[q].lower = lower[q];
        buildsub[q].upper = upper[q];
        buildsub[q].npart = npart[q];
        buildsub[q].first = first_node ;
        buildsub[q].length= (int) ( npart[q] * frac_part_node );//node_length
        buildsub[q].totpart = dp->NumPart;
//        buildsub[q].tree  = &(dtp->tree[q]);
        buildsub[q].root_cell = dtp->root_cell; //???
//       buildsub[q].part  = dp->Part;
        buildsub[q].level = gp->NumBits;
        buildsub[q].box   = gp->BoxSize;
        buildsub[q].min[0]= sub->minimum[0];
        buildsub[q].min[1]= sub->minimum[1];
        buildsub[q].min[2]= sub->minimum[2];
		buildsub[q].partidxes=dtp->partidxes;
		buildsub[q].numparticles=dtp->numparticles;

        first_node +=  buildsub[q].length ;
    }
	printf("Before subtree pthread\n");
    for (q=0; q<nThread; q++) {
        pthread_create(&task[q], NULL, mt_build_subtree_morton, &buildsub[q]);
    }

    for (q=0; q<nThread; q++) {
        pthread_join(task[q], NULL);
    }
	printf("After subtree pthread\n");
	free(lower);
    free(upper);
    free(npart);

    free(task);
    free(buildsub);

}
