#include "global.h"
#include <malloc.h>

MPI_Comm PM_COMM_WORLD;
MPI_Comm TREE_COMM_WORLD;

void system_exit(int err_num)
{
#ifdef __INTEL_OFFLOAD
	exit(err_num);
#else
    MPI_Abort(MPI_COMM_WORLD, err_num);
#endif
}

void* xmalloc(size_t memo_size, int err_num)
{
    void *p;
    p =malloc(memo_size);
    if (!p) {
        LOG_INFOL(LOG_FATAL, " Allocate error [ %d ] ! \n", err_num );
        system_exit(0);
    }
    return p;
}

void* xrealloc(void *mem, size_t memo_size, int err_num)
{
    void *p;
    p =realloc(mem, memo_size);
    if (!p) {
        LOG_INFOL(LOG_FATAL, " Re-allocate error [ %d ] ! \n", err_num );
        system_exit(0);
    }
    return p;
}

void* xmemalign(size_t memo_size, int err_num)
{
    void *p;
    p = memalign(ALIGNCNT, memo_size);
    if (!p) {
        LOG_INFOL(LOG_FATAL, " Aligned allocate error [ %d ] ! \n", err_num );
        system_exit(0);
    }
    return p;
}

static int rotcube[8][8] =
{
	{0,7,6,1,2,5,4,3},
    {0,3,4,7,6,5,2,1},
    {0,3,4,7,6,5,2,1},
    {2,3,0,1,6,7,4,5},
    {2,3,0,1,6,7,4,5},
    {6,5,2,1,0,3,4,7},
    {6,5,2,1,0,3,4,7},
    {4,3,2,5,6,1,0,7}
};

Peano encoding_peano(int x, int y, int z, int nBits)
{
    Peano key = 0;
    int d, i3, rank, temp[8], dirc[8] = { 0,1,3,2,6,7,5,4 };
    while (nBits--) {
        i3 = ( ((x>>nBits)&1)*2 + ((y>>nBits)&1) )*2 + ((z>>nBits)&1);
        for (rank=0; dirc[rank] != i3; rank++);
        key += rank;
        if (nBits ==0)
            break;
        key <<= 3;
        for (d=0; d<8; d++)
            temp[d] = dirc[d];
        for (d=0; d<8; d++)
            dirc[rotcube[rank][d]] = temp[d];
    }
    return key;
}

