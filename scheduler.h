#ifndef SCHEDULER_H
#define SCHEDULER_H

#include "data.h"
#include "domain.h"
#include "parameter.h"

typedef struct RecvFront {

} RecvBuffer;

typedef struct {
    real mass;
    real pos[3];
} BuffData;

typedef struct BufferFrontier {
    int nPart;
    BuffData *buff;
} BuffFront;

#endif /* SCHEDULER_H */
