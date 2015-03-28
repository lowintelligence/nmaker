#include "global.h"
#include "iosnap.h"

Snapshot* create_snapshot(char path_name[]) {
    Snapshot *snap;
    const int root_proc = 0;
    int rank_proc;
    int n, nfile, *npart;
    char content_name[260];

    if ( 256 < sizeof(path_name)) {
        printf("error: file name is too long!\n");
        system_exit(0);
    }

    sprintf(content_name, "%s.cnt", path_name);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);

    if (root_proc == rank_proc) {
        FILE *fin;
        fin = fopen(content_name, "r");
        if ( !fin ) {
            printf("error: cannot create snapshot! %s\n", content_name);
            system_exit(0);
        }
        fscanf(fin, "%d", &(nfile) );
        npart = xmalloc(nfile*sizeof(int) ,1011) ;

        for (n=0; n<nfile; n++) {
            fscanf(fin, "%d", &(npart[n]) );
        }
        fclose(fin);
    }
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Bcast( &nfile, 1, MPI_INT, root_proc, MPI_COMM_WORLD);

    snap = (Snapshot*)xmalloc( sizeof(Snapshot), 1012);
    snap->num_file = nfile;
    snap->num_part = (int*)xmalloc( nfile*sizeof(int) ,1013);
    sprintf((snap->path_name), "%s", path_name );

    if (root_proc == rank_proc) {
        MPI_Bcast( npart, nfile, MPI_INT, root_proc, MPI_COMM_WORLD);
        for (n=0; n<nfile; n++) {
            snap->num_part[n] = npart[n];
        }
        free(npart);
    }
    else {
        MPI_Bcast(&(snap->num_part[0]),nfile,MPI_INT,root_proc,MPI_COMM_WORLD);
    }

    //for (n=0; n<nfile; n++) printf(" [%d] num_part[%d/%d] = %d\n", rank_proc, n, snap->num_file, snap->num_part[n]);

    return snap;
}

void free_snapshot(Snapshot* snap_ptr) {
    if (NULL != snap_ptr) {
        free(snap_ptr->num_part);
        free(snap_ptr);
    }
    else {
        printf(" error: free_snap \n");
        system_exit(0);
    }
}

typedef struct io_header
{
    int      npart[6];
    double   mass[6];
    double   time;
    double   redshift;
    int      flag_sfr;
    int      flag_feedback;
    int      npartTotal[6];
    int      flag_cooling;
    int      num_files;
    double   BoxSize;
    double   Omega0;
    double   OmegaLambda;
    double   HubbleParam;
    char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];
    /* fills to 256 Bytes */
} Header;

void read_body_from_file(Body *sys, char filename[], Long n_start, Long n_body)
{
    int n;
    int dummy;
    long int disp, num_part;
    int data_size, id_size;

//    unsigned int id;
    int id;
    float data[3];

    id_size   = (int)sizeof(int);
    data_size = (int)sizeof(float);

    FILE* fsnap;
    Header head;
    if ( !(fsnap=fopen(filename,"r")) )
    {
        printf(" error: cannot open file %s \n", filename);
        exit(0);
    }
    printf("reading file `%s' \n", filename);
    fread(&dummy, sizeof(dummy), 1, fsnap);
    fread(&head, sizeof(Header), 1, fsnap);
    fread(&dummy, sizeof(dummy), 1, fsnap);
    num_part = head.npart[1];

    printf(" readin mass_p = %f\n", head.mass[1]);
    printf(" readin hubble = %f\n", head.HubbleParam);
    printf(" readin omega0 = %f\n", head.Omega0);
    printf(" readin boxsize = %f\n", head.BoxSize);
    printf(" readin numpart = %d\n", head.npartTotal[1]);
    printf(" readin redshift= %f\n", head.redshift);


    disp  = (3*sizeof(dummy) + sizeof(Header) + (n_start )*3*data_size);
    fseek(fsnap, disp, 0);

    for (n=0; n<n_body; n++)
    {
        fread( &(data[0]), data_size, 3, fsnap);
        sys[n].pos[0] = (Real)data[0];
        sys[n].pos[1] = (Real)data[1];
        sys[n].pos[2] = (Real)data[2];
    }
    disp  = (5*sizeof(dummy) + sizeof(Header) + (n_start + num_part)*3*data_size);
    fseek(fsnap, disp, 0);
    for (n=0; n<n_body; n++)
    {
        fread( &(data[0]), data_size, 3, fsnap);
        sys[n].vel[0] = (Real)data[0];
        sys[n].vel[1] = (Real)data[1];
        sys[n].vel[2] = (Real)data[2];
    }

    disp  = (7*sizeof(dummy) + sizeof(Header) + n_start*id_size + num_part*6*data_size);
    fseek(fsnap, disp, 0);

    for (n=0; n<n_body; n++)
    {
        fread(&id, id_size, 1, fsnap);
        sys[n].id = (Long)id;
        sys[n].mass = head.mass[1];
    }
    fclose(fsnap);

    double gdt2unit = pow( 1.0/(1.0+head.redshift), 1.5); //for gadget2 format
    for (n=0; n<n_body; n++) {
        sys[n].vel[0] *= gdt2unit;
        sys[n].vel[1] *= gdt2unit;
        sys[n].vel[2] *= gdt2unit;
    }

}

void read_snapshot(Body *sys, Long n_start, Long n_count, Snapshot *snap)
{
    int n;
    int rank_proc;
    char fname[260];
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);

    if (NULL == sys) {
        printf("error: NULL == Body* \n");
        system_exit(0);
    }

    Long n_end = n_start + n_count - 1;
    Long n_s, n_e, n_c;
    Long start_f, end_f;

    for (n=0, end_f=0; n<(snap->num_file); n++)
    {
        start_f = end_f;
        end_f = start_f + snap->num_part[n];
        if ( start_f > n_end )
            continue;
        if ( end_f < n_start )
            continue;

        if (n_start < start_f)
            n_s = start_f;
        else
            n_s = n_start;

        if (n_end > end_f)
            n_e = end_f;
        else
            n_e = n_end;

        n_c = n_e - n_s + 1;
        n_s -= start_f;

        sprintf(fname, "%s.%d", snap->path_name, n);
        read_body_from_file(sys, fname, n_s, n_c);
    }

}

void write_snapshot_node(Body *sys, Long n_count, char fname[], Constants *cst, Status *sta)
{
    int n;
    int rank_proc;
    const int root = 0;
    char filename[200];
    int rank, numproc;
    int npart = (int)n_count;

	CLOG("Write %d particle into file '%s'\n", npart, fname);
    Header head;
//////////////////////////////////////
    head.npart[0] = 0;
    head.npart[1] = npart;
    head.npart[2] = 0;
    head.npart[3] = 0;
    head.npart[4] = 0;
    head.npart[5] = 0;

    head.npartTotal[0] = 0;
    head.npartTotal[1] = (int)(cst->TOTAL_NUM_PART);
    head.npartTotal[2] = 0;
    head.npartTotal[3] = 0;
    head.npartTotal[4] = 0;
    head.npartTotal[5] = 0;

    head.mass[0] = 0;
    head.mass[1] = (double)(cst->PART_MASS);
    head.mass[2] = 0;
    head.mass[3] = 0;
    head.mass[4] = 0;
    head.mass[5] = 0;

    head.time = sta->redshift;
    head.redshift  = sta->redshift;

    head.flag_sfr  = 0;
    head.flag_feedback = 0;
    head.flag_cooling  = 0;

    head.BoxSize   = cst->BOX_SIZE;
    head.Omega0    = cst->OMEGA_M0;
    head.OmegaLambda = cst->OMEGA_X0;
    head.HubbleParam = cst->HUBBLE_0;
    head.num_files = cst->NUM_PROCESS;
///////////////////////////////////
    if (NULL == sys) {
        printf("error: NULL == Body* \n");
        system_exit(0);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numproc);

    FILE* fsnap;
    sprintf(filename, "%s.%d", fname, rank);
    if ( !(fsnap=fopen(filename,"w")) )
    {
        printf(" error: cannot open file %s \n", filename);
        exit(0);
    }
    int id, dummy = 0;
    float data[3];

    fwrite(&dummy, sizeof(int), 1, fsnap);
    fwrite(&head, sizeof(Header), 1, fsnap);
    fwrite(&dummy, sizeof(int), 1, fsnap);

    fwrite(&dummy, sizeof(int), 1, fsnap);
    for (n=0; n<npart; n++) {
        data[0] = (float)sys[n].pos[0];
        data[1] = (float)sys[n].pos[1];
        data[2] = (float)sys[n].pos[2];
        fwrite(data, sizeof(float), 3, fsnap);
    }
    fwrite(&dummy, sizeof(int), 1, fsnap);


    double gdt2unit = pow( 1.0/(1.0+head.redshift), 1.5); //for gadget2 format

    fwrite(&dummy, sizeof(int), 1, fsnap);
    for (n=0; n<npart; n++) {
        data[0] = (float)(sys[n].vel[0]/gdt2unit);
        data[1] = (float)(sys[n].vel[1]/gdt2unit);
        data[2] = (float)(sys[n].vel[2]/gdt2unit);
        fwrite(data, sizeof(float), 3, fsnap);
    }
    fwrite(&dummy, sizeof(int), 1, fsnap);

    fwrite(&dummy, sizeof(int), 1, fsnap);
    for (n=0; n<npart; n++) {
        id = (int)sys[n].id;
        fwrite(&id, sizeof(int), 1, fsnap);
    }
    fwrite(&dummy, sizeof(int), 1, fsnap);

    fclose(fsnap);
    Int *num_part;
    if ( root == rank)
    {
        num_part = (int*)xmalloc(numproc*sizeof(int), 1014);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Gather(&npart, sizeof(Int), MPI_CHAR, num_part, sizeof(Int), MPI_CHAR, root, MPI_COMM_WORLD);

    if ( root ==  rank)
    {
        FILE *fd;
        sprintf(filename, "%s.cnt", fname);
        fd = fopen(filename, "w");
        fprintf(fd, "%d\n", head.num_files);
        for (n=0; n<head.num_files; n++)
            fprintf(fd, "%d\n", num_part[n]);
        fclose(fd);
        free(num_part);
    }

}


