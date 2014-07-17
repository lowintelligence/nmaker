#ifndef DOMAIN_H
#define DOMAIN_H

#include "data.h"
#include "fillcurve.h"
#include "snapshot.h"

typedef struct {
    unsigned long int NumPart;
    code LowerHilbertKey;
    code UpperHilbertKey;
} DomainInfo;

typedef struct {
    int minimum[3]; // step 1
    int maximum[3]; // step 1
    int nSide[3];   // step 1
    int nPadding;   // 4/6
    int nBits;
//    code LowerKey;  // step 1
//    code UpperKey;  // step 1
//////////    int nbits;
    int  *tag;      // step 5 boundary? core? other domain?
    int  *count;    // step 2 --> first group to MIC -->tag the particle iGroup
//   code *cubekey;
//////////    code *key;      // step 4 --> combine minimum+(i,j,k)
    double *mesh;   // step 3 --> cic scheme
} SubCuboid;

//struct BufferExchangeFrontier;

typedef struct {
    int rank;
    long int NumPart;
    Body *Part;

    int NumDom;
    DomainInfo* DomList;
    DomainInfo* thisdom;
    SubCuboid  *cuboid;
    DomainTree *domtree;
    int progress;
    struct RecvFront *buff;
    // send-recv communication rank/buffer
} Domain;


void load_particle_into_domain(Domain* dp, int myid, int numprocs);

#endif /* DOMAIN_H */
