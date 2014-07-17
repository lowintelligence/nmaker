#include "fillcurve.h"

static int rotcube[8][8] =
{   {0,7,6,1,2,5,4,3},
    {0,3,4,7,6,5,2,1},
    {0,3,4,7,6,5,2,1},
    {2,3,0,1,6,7,4,5},
    {2,3,0,1,6,7,4,5},
    {6,5,2,1,0,3,4,7},
    {6,5,2,1,0,3,4,7},
    {4,3,2,5,6,1,0,7}
};

code coding_3d_filling_curve(int x, int y, int z, int nBits) {
    code key = 0LL;
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

