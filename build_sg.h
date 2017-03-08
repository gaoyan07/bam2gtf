#ifndef _SPLICE_GRAPH_H
#define _SPLICE_GRAPH_H
#include "gtf.h"
#include "utils.h"

#define _insert(v, p, n, m, type) { \
    int i, flag=0;                  \
    for (i = 0; i < n; ++i) {       \
        if (p[i] == v) {            \
            flag = 1;               \
            break;                  \
        }                           \
    }                               \
    if (flag == 0) {                \
        if (n == m) {               \
            _realloc(p, m, type)    \
        }                           \
        p[n++] = v;                 \
    }                               \
}

typedef struct {
    uint32_t node_id; // unique id in corresponding gene-locus
    exon_t e;
    uint32_t *next_id;  int32_t next_n, next_m;
    uint32_t *pre_id;    int32_t pre_n, pre_m;
    uint32_t *pre_domn;  int32_t pre_domn_n, pre_domn_m;
    uint32_t *post_domn; int32_t post_domn_n, post_domn_m;
} SGnode; // node of splicing-graph

typedef struct {
    uint32_t site_id;
    int32_t site;
    uint32_t *exon_id; int32_t exon_n, exon_m;
    uint8_t type; // 0: 5'(donor), 1: 3'(acceptor)
} SGsite; // splice-site of splicing-graph

typedef struct {
    uint32_t don_site_id, acc_site_id;
    uint8_t is_rev;
    int32_t cov;
    uint8_t motif, is_anno;
    int32_t uniq_c, multi_c, max_over;
} SGedge; // edge of splicing-graph, splice junction

typedef struct {
    //SGnode v;  // virtual start and end node
    // virtual_start: node[0]; virtual_end: node[node_n-1]
    SGnode *node; int32_t node_n, node_m; // sort by e.start and e.end 
    SGsite *don_site; int32_t don_site_n, don_site_m; // sort by site
    SGsite *acc_site; int32_t acc_site_n, acc_site_m; // sort by site
    SGedge *edge; int32_t edge_n, edge_m; // sort by don_id and acc_id
    int32_t tid; uint8_t is_rev;
    // boundaries of splice-sites
    int32_t start, end; 
    // for calculate max-edge XXX
    uint8_t **path_map; // size: node_n * node_n, but only 1/2 of the bit-map is valid
                        // path_map[i][j] == 1 : path from i to j
                        // 0: none path
                        // 1: original edge
                        // 2: indirect path
} SG;

typedef struct {
    uint32_t SG_id;
    uint32_t v_start, v_end; // virtual start and end node
    uint32_t *node_id; int32_t node_n, node_m;
    uint32_t *edge_id; int32_t edge_n, edge_m;
} SGasm;

typedef struct {
    SGasm **sg_asm;
    int sg_asm_n, sg_asm_m;
} SGasm_group;
 
typedef struct {
    SG **SG;
    int32_t SG_n, SG_m;
} SG_group;

SG *sg_init_node(SG *sg);
SG *sg_init_site(SG *sg);
SG_group *sg_init_group(int g_n);
void sg_free_group(SG_group *sg_g);

int sg_update_node(SG *sg, exon_t e);
int sg_update_site(SG *sg, int32_t site, uint8_t type);

int sg_bin_sch_node(SG sg, exon_t e, int *hit);
int sg_bin_sch_site(SGsite *site, int32_t site_n, int32_t s, int *hit);
int sg_bin_sch_edge(SG sg, uint32_t don_site_id, uint32_t acc_site_id, int *hit);

void cal_pre_domn(SG *sg);
void cal_post_domn(SG *sg);

SG_group *construct_SpliceGraph(FILE *gtf_fp, chr_name_t *cname);

int build_sg(int argc, char *argv[]);
void sg_dump(SG_group sg_g, const char *sg_name);
SG_group *sg_restore(const char *sg_name);

#endif
