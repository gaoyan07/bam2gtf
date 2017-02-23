#ifndef _SPLICE_GRAPH_H
#define _SPLICE_GRAPH_H
#include "gtf.h"

typedef struct SG {
    exon_t e;
    struct SG *child; int child_n, child_m;
} SGnode; // node of splicing-graph

#endif
