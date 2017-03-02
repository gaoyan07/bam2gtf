#ifndef _SPLICE_GRAPH_H
#define _SPLICE_GRAPH_H
#include "gtf.h"

typedef struct {
    uint32_t node_id; // unique id in corresponding gene-locus
    exon_t e;
    uint32_t *next_id;  int next_n, next_m;
    uint32_t *pre_id;    int pre_n, pre_m;
    uint32_t *pre_domn;  int pre_domn_n, pre_domn_m;
    uint32_t *post_domn; int post_domn_n, post_domn_m;
} SGnode; // node of splicing-graph

typedef struct {
    uint32_t site_id;
    int32_t site;
    uint32_t *exon_id; int exon_n, exon_m;
    uint8_t type; // 0: 5'(donor), 1: 3'(acceptor)
} SGsite; // splice-site of splicing-graph

typedef struct {
    uint32_t don_site_id, acc_site_id;
    uint8_t is_rev;
    double cov;
} SGedge; // edge of splicing-graph, splice junction

typedef struct {
    SGnode v;  // virtual start and end node
    SGnode *node; int node_n, node_m; // sort by e.start and e.end 
    SGsite *site; int site_n, site_m; // sort by site
    SGedge *edge; int edge_n, edge_m; // sort by don_id and acc_id
    int32_t start, end;
    // for calculate max-edge
    uint8_t **path_map; // size: node_n * node_n, but only 1/2 of the bit-map is valid
                        // path_map[i][j] == 1 : path from i to j
                        // 0: none path
                        // 1: original edge
                        // 2: indirect path
} SG;

typedef struct {
    uint32_t SG_id;
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
