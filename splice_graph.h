#ifndef _SPLICE_GRAPH_H
#define _SPLICE_GRAPH_H
#include "gtf.h"

typedef struct {
    uint32_t node_id; // unique id in corresponding gene-locus
    exon_t e;
    uint32_t *next_id; int next_n, next_m;
    uint32_t *pre_id; int pre_n, pre_m;
    uint32_t *pre_domn; int pre_domn_n, pre_domn_m;
    uint32_t *post_domn; int post_domn_n, post_domn_m;
} SGnode; // node of splicing-graph

typedef struct {
    uint32_t don_id, acc_id;
    uint8_t is_rev;
    double cov;
} SGedge; // edge of splicing-graph, splice junction

typedef struct {
    SGnode v;        // virtual start and end node
    SGnode *node;    // sort by e.start and e.end
    int node_n, node_m;
    uint32_t root_id;        // id of root node
    SGedge *edge;    // sort by (acc.start-don.end)
    int edge_n, edge_m;
    int32_t start, end;
    // for calculate max-edge
    uint8_t **path_map; // size: node_n * node_n, but only 1/2 of the bit-map is valid
                        // path_map[i][j] == 1 : path from i to j
                        // 0: none path
                        // 1: original edge
                        // 2: indirect path
} SG;

typedef struct {
    uint32_t SG_id;     // generate from which SG
    uint32_t v_start, v_end; // virtual start and end node
    uint32_t *node_id; int node_n, node_m;
    uint32_t *edge_id; int edge_n, edge_m;
} SGasm;

typedef struct {
    SGasm **sg_asm;
    int sg_asm_n, sg_asm_m;
} SGasm_group;
 
typedef struct {
    SG **SG;
    int SG_n, SG_m;
} SG_group;

int sg(int argc, char *argv[]);

#endif
