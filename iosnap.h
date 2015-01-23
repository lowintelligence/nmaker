#ifndef _IOSNAP_H_
#define _IOSNAP_H_

#include "global.h"

typedef struct
{
    char path_name[1024];
    int  num_file;
    int *num_part;
} Snapshot;

Snapshot* create_snapshot(char path_name[]);
void write_snapshot_node(Body *sys, Long n_count, char fname[], Constants *cst, Status *sta);

#endif /* _IOSNAP_H_ */
