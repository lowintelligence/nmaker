#include "subcuboid.h"
#include <assert.h>
#include <stdio.h>


SubCuboid* get_frame_subcuboid(Domain* dp, GlobalParam *gp, int npad, int nbits)
{
    int n, m, i, j, k;
    int x, y, z, idx;
    int minidx[3];
    int maxidx[3];
    int nside[3];
    code cubekey;
    code lower, upper;
    double boxsize = gp->BoxSize;
    double nbox_boxsize = (double)(1<<nbits)/boxsize;
    unsigned long int npart = dp->NumPart;
    real minreal[3], maxreal[3];
    Body *part = dp->Part;
    SubCuboid *sub;
//    printf(" nbits = %d\n", nbits);

    for (m=0; m<3; m++) {
        minreal[m] = boxsize;
        maxreal[m] = 0.0;
    }
    for (n=0; n<npart; n++) {
        for ( m=0; m<3; m++ ) {
            if (part[n].pos[m] < minreal[m])
                minreal[m] = part[n].pos[m];
            if (part[n].pos[m] > maxreal[m])
                maxreal[m] = part[n].pos[m];
        }
    }

    for (m=0; m<3; m++) {
        minidx[m] = (int)(minreal[m]*nbox_boxsize);
        maxidx[m] = (int)(maxreal[m]*nbox_boxsize);
        nside[m] = maxidx[m] - minidx[m] + 1 + 2*npad; /* 001 */

#ifndef SILENCE
        printf(" pos[%d]=[%f, %f],[%d<<%d (...%d...) %d>>%d]\n", m,
               minreal[m],maxreal[m],npad,minidx[m],nside[m],maxidx[m],npad);
#endif
    }

    sub = (SubCuboid*) malloc( sizeof(SubCuboid) );
    sub->nPadding = npad;
    sub->nBits = nbits;

    for (m=0; m<3; m++) {
        sub->minimum[m] = minidx[m];
        sub->maximum[m] = maxidx[m];
        sub->nSide[m]   =  nside[m];
    }
    lower = dp->thisdom->LowerHilbertKey;
    upper = dp->thisdom->UpperHilbertKey;

    printf(" [%d]sub->Key:%ld %ld\n", dp->rank, lower, upper);

    sub->tag   = (int*)malloc(sizeof(int)*nside[0]*nside[1]*nside[2]);
    sub->count = (int*)malloc(sizeof(int)*nside[0]*nside[1]*nside[2]);
    sub->mesh  = NULL;

    return sub;
}

/***************** identify the cuboid within boundary *******************/
void mark_cuboid(Domain *dp, SubCuboid *sub)
{
    int n, m, i, j, k;
    int x, y, z, X, Y, Z, idx, Idx, cnt;
    int In, Jn, Kn;
    int minidx[3];
    int maxidx[3];
    int nside[3];
    int npad = sub->nPadding;
    int nbits= sub->nBits;
    code cubekey, *cdkey;
    code lower, upper;
    int *label;

    for (m=0; m<3; m++) {
        minidx[m] = sub->minimum[m];
        maxidx[m] = sub->maximum[m];
        nside[m] = sub->nSide[m];
//    printf("%d %d %d %d\n", m, minidx[m], maxidx[m], nside[m]);
    }
//    printf("npad = %d , nbits = %d\n", npad, nbits);

    printf(" - min : %d %d %d\n", minidx[0], minidx[1], minidx[2]);
    printf(" - max : %d %d %d\n", maxidx[0], maxidx[1], maxidx[2]);
    printf(" nside : %d %d %d\n", nside[0], nside[1], nside[2]);
    In = nside[0];
    Jn = nside[1];
    Kn = nside[2];
    label = (int*)malloc(sizeof(int)*In*Jn*Kn);
    cdkey = (code*)malloc(sizeof(code)*In*Jn*Kn);

    assert( In>0 && Jn>0 && Kn>0 );

    lower = dp->thisdom->LowerHilbertKey;
    upper = dp->thisdom->UpperHilbertKey;

    for ( i=0; i<In; i++ ) {
        for ( j=0; j<Jn; j++ ) {
            for ( k=0; k<Kn; k++ ) {
                x = i - npad + minidx[0];
                y = j - npad + minidx[1];
                z = k - npad + minidx[2];
                idx = ( i*Jn + j )*Kn + k;
                sub->count[idx] = 0;
                sub->tag[idx] = TAG_IDLEFIELD;
                if ( i<npad || i>(In-npad-1)
                        || j<npad || j>(Jn-npad-1)
                        || k<npad || k>(Kn-npad-1) ) {
                    label[idx] = 0;
                    x = (x+In)%In;
                    y = (y+Jn)%Jn;
                    z = (z+Kn)%Kn;
                    cdkey[idx] = coding_3d_filling_curve(x, y, z, nbits);
                }
                else {
                    cubekey = coding_3d_filling_curve(x, y, z, nbits);
                    if ( cubekey>=lower && cubekey<upper )
                        label[idx] = 1 ;
                    else
                        label[idx] = 0 ;

                }
            }
        }
    }

    /*
    for ( i=npad; i<In-npad; i++ ) {
        for ( j=npad; j<Jn-npad; j++ ) {
            for ( k=npad; k<Kn-npad; k++ ) {
                Idx = ( i*Jn + j )*Kn + k;
                if (label[Idx]) {
                    cnt = 0;
                    for (x=i-1; x<=i+1; x++)
                        for (y=j-1; y<=j+1; y++)
                            for (z=k-1; z<=k+1; z++) {
                            cnt += label[(x*Jn + y )*Kn + z];
                        }
                    if ( 0 == cnt ) {
                        sub->tag[Idx] = -1;
                    } else
                    if ( 27 > cnt ) {
                        sub->tag[Idx] = 1; // boundary grid
                    }
                    else
                        sub->tag[Idx] = 2; // inner grid

                }
            }
        }
    }
    */

    assert( 1==npad );

    for ( i=0; i<In; i++ ) {
        for ( j=0; j<Jn; j++ ) {
            for ( k=0; k<Kn; k++ ) {
                Idx = ( i*Jn + j )*Kn + k;

                cnt = 0;
                for (x=i-1; x<=i+1; x++) {
                    X = x;
                    if (x<0)
                        X += In;
                    if (x>=In)
                        X -= In;
                    for (y=j-1; y<=j+1; y++) {
                        Y = y;
                        if ( y<0)
                            Y += Jn;
                        if (y>=Jn)
                            Y -= Jn;
                        for (z=k-1; z<=k+1; z++) {
                            Z = z;
                            if (z < 0)
                                Z += Kn;
                            if (z>=Kn)
                                Z -= Kn;
                            cnt += label[(X*Jn + Y )*Kn + Z];
                        }
                    }
                }

                switch(cnt) {
                case 27:
                    sub->tag[Idx] =  TAG_OUTERCORE;
                    break;
                case  0:
                    sub->tag[Idx] =  TAG_IDLEFIELD;
                    break;
                default:
                    if ( 1==label[Idx] )
                        sub->tag[Idx] = TAG_FRONTIERS;
                    else if ( 0==label[Idx] )
                        sub->tag[Idx] = TAG_ADJOINING;
                    break;
                }
                /*
                if ( 27 == cnt ) {
                    sub->tag[Idx] = 2; // inner grid
                //     printf("here\n");
                }
                if (  0 == cnt ) {
                    sub->tag[Idx] = -1; // too far
                }
                else if ( 27  > cnt ) {
                    if (label[Idx])
                        sub->tag[Idx] = 1; // boundary grid
                    else
                        sub->tag[Idx] = 0; //neigbor grid

                }
                */
            }
        }
    }


    int nSplit, nLevel, splitLength;
    nSplit = dp->NumDom;
    for (nLevel = 0; nSplit > (1<<nLevel); nLevel++);
    splitLength = 1<<nLevel;
    code *split = (code*) malloc( sizeof(code) * splitLength );

//   printf("nSplit = %d nLevel = %d\n", nSplit, nLevel);
    for (i=0; i<splitLength; i++) {
        if (i<nSplit)
            split[i] = dp->DomList[i].LowerHilbertKey;
        else
            split[i] =  dp->DomList[nSplit-1].UpperHilbertKey;
//        printf(" split[%d] = %d \n",i, split[i]);
    }

    int dom, l;
    int level = 1<< (nLevel);

    level >>= 1;
    for ( i=0; i<In; i++) {
        for ( j=0; j<Jn; j++) {
            for ( k=0; k<Kn; k++) {
                idx = (i*Jn+j)*Kn + k ;
                if ( TAG_ADJOINING == sub->tag[idx]  ) {

                    cubekey = cdkey[idx];

                    /* determine which domin this cube belongs to */
                    l = dom = level ;
                    while ( l >>= 1 ) {
                        if (cubekey < split[dom]) {
                            dom -= l;
                        }
                        else {
                            dom += l;
                        }
                    }
                    if (cubekey < split[dom]) {
                        dom -= 1;
                    }
                    /* determine which domin this cube belongs to */

                    sub->tag[idx] = dom;

                }

            }
        }
    }

    free(split);
    /*
        int myrank;
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        if (myrank == 0) {
            FILE *fd;
            //  assert(fd=fopen("test_boundary_label.dat", "w")) ;
            assert(fd=fopen("test_boundary_tag.dat", "w")) ;
            for ( i=0; i<In; i++) {
                for ( j=0; j<Jn; j++) {
                    for ( k=0; k<Kn; k++) {
                        //              fprintf(fd, "%d ", label[ (i*Jn+j)*Kn + k] );
                        fprintf(fd, "%2d ", sub->tag[ (i*Jn+j)*Kn + k ] );
                    }
                    fprintf(fd, "\n");
                }
                fprintf(fd, "\n");
            }
            fclose(fd);
        }

    */
    free(label);
    free(cdkey);
}

void mark_particle_by_cuboid(SubCuboid* sub, Domain* dp, double boxsize, int nbits)
{
    int n, npart = dp->NumPart;
    int x, y, z, idx;
    int npad = sub->nPadding;
    int real2int;
    int minidx[3];
    int In, Jn, Kn;
    double nbox_boxsize = (1<<nbits)/boxsize;
    Body *part = dp->Part;

    In = sub->nSide[0];
    Jn = sub->nSide[1];
    Kn = sub->nSide[2];

    minidx[0] = sub->minimum[0];
    minidx[1] = sub->minimum[1];
    minidx[2] = sub->minimum[2];

    for (n=0; n<npart; n++) {
        x = (int)(part[n].pos[0]*nbox_boxsize) - minidx[0] + npad;
        y = (int)(part[n].pos[1]*nbox_boxsize) - minidx[1] + npad;
        z = (int)(part[n].pos[2]*nbox_boxsize) - minidx[2] + npad;
        idx = (x*Jn + y) * Kn + z;
        part[n].tag = sub->tag[idx];
        part[n].group = idx;
        sub->count[idx]++;
    }
}

#include <stdio.h>
/* check the particle number */
void statistic_subcuboid(SubCuboid* sub, Domain *dp, int rank) {
#define NCEIL 5
    int i,j,k, n;
    int In, Jn, Kn, idx;
    int c, ceil[NCEIL] = { 8, 64, 512, 4096, 10000 };
    int accu[NCEIL] = {0,0,0,0,0};
    In = sub->nSide[0];
    Jn = sub->nSide[1];
    Kn = sub->nSide[2];
    FILE *fd;
    char fname[100];

    sprintf(fname,"clip%03d.dat", rank);

    for (i=0; i<In; i++) {
        for (j=0; j<Jn; j++) {
            for (k=0; k<Kn; k++) {
                idx = (i*Jn+j)*Kn+k;
                for (c=0; c<NCEIL; c++ )
                    if (sub->count[idx] > ceil[c])
                        accu[c]++;
            }
        }
    }

    int cnt_outercore, cnt_frontiers;
    Body *part = dp->Part;

    cnt_outercore = cnt_frontiers = 0;

    for (n=0; n<dp->NumPart; n++) {
        if ( TAG_FRONTIERS == part[n].tag )
            cnt_frontiers++;
        if ( TAG_OUTERCORE == part[n].tag )
            cnt_outercore++;
    }
    printf(" [%d]  frontiers=%d outercore=%d ratio = %f\n", rank, cnt_frontiers, cnt_outercore, (float)cnt_frontiers/cnt_outercore );

    for (c=0; c<NCEIL; c++ )
        printf("%8d: %8d\n", ceil[c], accu[c]);
    /*
    assert(fd = fopen(fname,"w") );
    i = In/2;
        for (j=0; j<Jn; j++) {
            for (k=0; k<Kn; k++) {
                idx = (i*Jn+j)*Kn+k;
             //   if (sub->count[idx] > 0)
                    fprintf(fd,"%d %d %d\n", j,k, sub->count[idx]);
            }
            fprintf(fd,"\n");
        }

    fclose(fd);
    */
    /*
    sprintf(fname,"domain%03d.part", rank);
    assert(fd = fopen(fname,"w") );

    for (n=0; n<dp->NumPart; n++) {
     //   assert( 0!= part[n].tag );
        if ( 1== part[n].tag )
            fprintf(fd,"%f %f %f\n", part[n].pos[0], part[n].pos[1], part[n].pos[2]);
    }

    fclose(fd);
    */

}
#include <math.h>

void statistic_packing(SubCuboid* sub, Domain *dp, int rank) {
    int  n, i, max, min, db, Nbin;
    int In, Jn, Kn, index_total;
    Nbin = 50000;
    max = 50000;
    min = 0;
    db = (max-min)/Nbin;
    int *numcube, *count;


    numcube = (int*) malloc(sizeof(int) * Nbin );
    for (n=0; n<Nbin; n++) {
        numcube[n] = 0;
    }
    In = sub->nSide[0];
    Jn = sub->nSide[1];
    Kn = sub->nSide[2];
    index_total = In * Jn * Kn;

    count = sub->count;
    for (i=0; i<index_total; i++) {
        for (n=0; n<Nbin; n++) {
            //    printf("%d %d\n", count[i], db*n+min);
            if ( count[i] > db*n+min ) {
                numcube[n] ++;
            }
            else
                break;
        }
    }
    char fname[100];
    sprintf(fname,"test_packing_%d", rank);
    FILE *fd = fopen(fname,"w");
    for (n=0; n<Nbin; n++) {
//        fprintf(fd, "%e %e\n", (double)(db*n+min)/(512.0*512*512), (double)numcube[n]/index_total );
        fprintf(fd, "%d %d %f\n",(db*n+min), numcube[n], 1.0/sqrt(numcube[n]) );
    }
    fclose(fd);
    free(numcube);
}


void construct_subcuboid(Domain *dp, GlobalParam *gp)
{
    SubCuboid* sub;

    sub = get_frame_subcuboid(dp, gp, 1, gp->NumBits);
    mark_cuboid(dp, sub);
    mark_particle_by_cuboid(sub, dp, gp->BoxSize, gp->NumBits);

//   if ( 0 == dp->rank)
//    statistic_subcuboid(sub, dp, dp->rank);
//    statistic_packing(sub, dp, dp->rank);
    dp->cuboid = sub;
}

void free_subcuboid(SubCuboid* sp)
{
    if ( NULL != sp->tag )
        free(sp->tag);
    if ( NULL != sp->count)
        free(sp->count);
    if ( NULL != sp->mesh )
        free(sp->mesh);
    free(sp);
}


/* sperate mesh into several threads */

/* 4 to 6 grid */
/* 4 to 3 grid */

/* boundary hash table : the first element is LowerHilbertKey */



