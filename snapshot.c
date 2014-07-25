#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "snapshot.h"

/*
void read_Header(Header *hptr, char filename[])
{
    int dummy;
    FILE* fsnap;

    if ( !(fsnap=fopen(filename,"r")) )
    {
        printf("[readsnap.c] cannot open file %s \n", filename);
        exit(0);
    }

    fread(&dummy, sizeof(dummy), 1, fsnap);
    fread(hptr, sizeof(Header), 1, fsnap);
    fread(&dummy, sizeof(dummy), 1, fsnap);
}
*/
void read_Particle(Body *pptr, char filename[], int n_start, int n_count)
{
    int n;
    int dummy;
    int disp, num_part;
    int id;

    FILE* fsnap;
    GadgetHeader head;
    if ( !(fsnap=fopen(filename,"r")) )
    {
        printf("[readsnap.c] cannot open file %s \n", filename);
        exit(0);
    }
//   printf("reading file `%s' \n", filename);
    fread(&dummy, sizeof(dummy), 1, fsnap);
    fread(&head, sizeof(GadgetHeader), 1, fsnap);
    fread(&dummy, sizeof(dummy), 1, fsnap);
    num_part = head.npart[1];

    disp  = (2*sizeof(dummy) + sizeof(GadgetHeader));
    disp += (sizeof(dummy) + n_start*3*sizeof(float));
    fseek(fsnap, disp,0);

    for (n=0; n<n_count; n++)
    {
        fread(&(pptr[n].pos[0]), sizeof(float), 3, fsnap);
    }

    disp  += (2*sizeof(dummy) + num_part*3*sizeof(float));
    fseek(fsnap, disp,0);
    for (n=0; n<n_count; n++)
    {
//        fread(&(pptr[n].vel[0]), sizeof(float), 3, fsnap);
    }

    disp   = (2*sizeof(dummy) + sizeof(GadgetHeader));
    disp  += (2*sizeof(dummy) + 2*num_part*3*sizeof(float) );
    disp  += (2*sizeof(dummy) + n_start*sizeof(int) );
    fseek(fsnap, disp,0);
    for (n=0; n<n_count; n++)
    {
        fread(&id, sizeof(int), 1, fsnap);
        pptr[n].ID = (long)id;
    }
}
/*
void read_block(Particle *pptr, char filename[], int iblock, int nblock)
{
    int n;
    int dummy;
    int disp, num_part;
    int n_start, n_end, n_count;
    Header head;
    FILE* fsnap;

    if ( !(fsnap=fopen(filename,"r")) )
    {
        printf("[readsnap.c] cannot open file %s \n", filename);
        exit(0);
    }

    fread(&dummy, sizeof(dummy), 1, fsnap);
    fread(&head, sizeof(Header), 1, fsnap);
    fread(&dummy, sizeof(dummy), 1, fsnap);
    num_part = head.npart[1];

    n_start = iblock * (int)(num_part/nblock);
    n_end = (iblock+1) * (int)(num_part/nblock);
    if ( (1+iblock) == nblock)
        n_end = num_part;
    n_count = n_end - n_start;

    disp  = (2*sizeof(dummy) + sizeof(Header));
    disp += (sizeof(dummy) + n_start*3*sizeof(float));
    fseek(fsnap, disp,0);

    for (n=0; n<n_count; n++)
    {
        fread(&(pptr[n].pos[0]), sizeof(float), 3, fsnap);
//      printf("%d  %f\n", n+n_start, pptr[n].pos[0]);
    }

    disp  += (2*sizeof(dummy) + num_part*3*sizeof(float));
    fseek(fsnap, disp,0);
    for (n=0; n<n_count; n++)
    {
        fread(&(pptr[n].vel[0]), sizeof(float), 3, fsnap);
//      printf("%d  %f\n", n+n_start, pptr[n].vel[0]);
    }
}
*/

Body* load_gadget2_snapshot(char *fname, int files, long int *npart, double *box)
{
    FILE *fd;
    char buf[200];
    int i, j, k, dummy, ntot_withmasses;
    int t, n, off, pc, pc_new, pc_sph;

    long int NumPart, Ngas;
    double Time, Redshift;
    GadgetHeader header1;
    Body *P;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

    for(i = 0, pc = 1; i < files; i++, pc = pc_new)
    {
        if(files > 1)
            sprintf(buf, "%s.%d", fname, i);
        else
            sprintf(buf, "%s", fname);

        if(!(fd = fopen(buf, "r")))
        {
            printf("can't open file `%s`\n", buf);
            exit(0);
        }
#ifndef SILENCE
        printf("reading `%s' ...\n", buf);
        fflush(stdout);
#endif
        fread(&dummy, sizeof(dummy), 1, fd);
        fread(&header1, sizeof(header1), 1, fd);
        fread(&dummy, sizeof(dummy), 1, fd);

        if(files == 1)
        {
            for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
                NumPart += header1.npart[k];
            Ngas = header1.npart[0];
        }
        else
        {
            for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
                NumPart += header1.npartTotal[k];
            Ngas = header1.npartTotal[0];
        }

        for(k = 0, ntot_withmasses = 0; k < 6; k++)
        {
            if(header1.mass[k] == 0)
                ntot_withmasses += header1.npart[k];
        }

        if(i == 0) {
            *npart = NumPart;
            *box = header1.BoxSize;
            P = (Body*) malloc(sizeof(Body)*NumPart);
            if (P == NULL)
                exit(0);
        }
#ifndef SILENCE
        if(i == 0) {
            printf("NumPart = %ld  Ngas = %ld\n", NumPart, Ngas);
            if (NumPart > 10000000)
                printf("allocate memory %.2f GB \n", sizeof(Body)*NumPart/(1024.*1024.*1024.0) );
            else if (NumPart > 10000)
                printf("allocate memory %.2f MB \n", sizeof(Body)*NumPart/(1024.*1024.0) );
            else
                printf("allocate memory %.3f KB \n", sizeof(Body)*NumPart/(1024.0) );
        }
#endif

        SKIP;
        for(k = 0, pc_new = pc; k < 6; k++)
        {
            for(n = 0; n < header1.npart[k]; n++)
            {
                fread(&P[pc_new].pos[0], sizeof(float), 3, fd);
                pc_new++;
            }
        }
        SKIP;

        SKIP;
        for(k = 0, pc_new = pc; k < 6; k++)
        {
            for(n = 0; n < header1.npart[k]; n++)
            {
//                fread(&P[pc_new].vel[0], sizeof(float), 3, fd);
                pc_new++;
            }
        }
        SKIP;


        SKIP;
        for(k = 0, pc_new = pc; k < 6; k++)
        {
            for(n = 0; n < header1.npart[k]; n++)
            {
                fread( &P[pc_new].ID, sizeof(int), 1, fd);

                //   printf("%ld\n", P[pc_new].ID );
                pc_new++;
            }
        }
        SKIP;


        if(ntot_withmasses > 0)
            SKIP;
        for(k = 0, pc_new = pc; k < 6; k++)
        {
            for(n = 0; n < header1.npart[k]; n++)
            {
                if(header1.mass[k] == 0)
                    fread(&P[pc_new].mass, sizeof(float), 1, fd);
                else
                    P[pc_new].mass = header1.mass[k];
                pc_new++;
            }
        }
        if(ntot_withmasses > 0)
            SKIP;


        if(header1.npart[0] > 0)
        {
            printf("error1 without gas now!\n");
            exit(0);
        }

        fclose(fd);
    }


    Time = header1.time;
    Redshift = header1.time;

    return P;
}



