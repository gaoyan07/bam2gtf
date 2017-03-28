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

#define _bin_insert(v, p, n, m, type) { \
    int flag=0,k_i=-1,left=0,right=n-1,mid;    \
    type mid_v, tmp_v;                 \
    if (right == -1) k_i = 0;   \
    else {                      \
        while (left <= right) { \
            mid = (left+right) >> 1;    \
            mid_v = p[mid];             \
            if (mid_v == v) {           \
                flag = 1; break;        \
            } else if (mid_v > v) {     \
                if (mid != 0) {         \
                    tmp_v = p[mid-1];   \
                }                       \
                if (mid == 0 || v > tmp_v)  \
                    k_i = mid;              \
            } else left = mid+1;        \
        }                               \
    }                                   \
    if (k_i == -1) k_i = n;         \
                                    \
    if (flag == 0) {                \
        if (n == m) {               \
            _realloc(p, m, type)    \
        }                           \
        if (k_i <= n-1)             \
            memmove(p+k_i+1, p+k_i, (n-k_i)*sizeof(type));  \
        (p)[k_i] = v;               \
        (n)++;                      \
    }                               \
}

typedef struct {
    exon_t up, se, down;
} SE_t;   // skipped exon

typedef struct {
    exon_t lon, shor, down;
} A5SS_t; // alternative 3' splice site

typedef struct {
    exon_t up, lon, shor;
} A3SS_t; // alternative 3' splice site

typedef struct {
    exon_t up, fir, sec, down;
} MXE_t; // mutually exclusive exon

typedef struct {
    exon_t up, down;
} RI_t;  // retained intron

typedef struct {
    SE_t *se; int32_t se_n, se_m;
    A5SS_t *a5ss; int32_t a5ss_n, a5ss_m;
    A3SS_t *a3ss; int32_t a3ss_n, a3ss_m;
    MXE_t *mxe; int32_t mxe_n, mxe_m;
    RI_t *ri; int32_t ri_n, ri_m;
} ASE_t;

typedef struct {
    uint32_t node_id; // unique id in corresponding gene-locus
    int32_t start, end; /* real exon */ exon_t node_e;    // node in splice-graph
    uint8_t is_init, is_termi;
    uint8_t is_asm; uint32_t uniq_c, multi_c;
    uint32_t *next_id;  int32_t next_n, next_m;
    uint32_t *pre_id;    int32_t pre_n, pre_m;
    uint32_t *pre_domn;  int32_t pre_domn_n, pre_domn_m;
    uint32_t *post_domn; int32_t post_domn_n, post_domn_m;
} SGnode; // node of splicing-graph

typedef struct {
    uint32_t site_id;
    int32_t site;
    uint32_t *exon_id; int32_t exon_n, exon_m;
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
} SG;

typedef struct {
    uint32_t SG_id;
    uint32_t v_start, v_end; // virtual start and end node
    uint32_t *node_id; int32_t node_n, node_m;
    uint32_t *edge_id; int32_t edge_n, edge_m;
    int32_t start, end;
} SGasm;

typedef struct {
    SGasm **sg_asm;
    int sg_asm_n, sg_asm_m;
} SGasm_group;
 
typedef struct {
    SG **SG;
    int32_t SG_n, SG_m;
    chr_name_t *cname;
} SG_group;

SG *sg_init_node(SG *sg);
SG *sg_init_site(SG *sg);
SG_group *sg_init_group(int g_n);
void sg_free_group(SG_group *sg_g);

int sg_update_node(SG *sg, exon_t e, int32_t start, int32_t end);
int sg_update_site(SG *sg, int32_t site, uint8_t type);

int sg_bin_sch_node(SG *sg, exon_t e, int *hit);
int sg_bin_sch_site(SGsite *site, int32_t site_n, int32_t s, int *hit);
int sg_bin_sch_edge(SG *sg, uint32_t don_site_id, uint32_t acc_site_id, int *hit);

void cal_pre_domn(SG *sg);
void cal_post_domn(SG *sg);

SG_group *construct_SpliceGraph(FILE *gtf_fp, chr_name_t *cname);

int build_sg(int argc, char *argv[]);

#endif
