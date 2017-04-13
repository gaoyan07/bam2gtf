#ifndef _BUILD_SG_H
#define _BUILD_SG_H
#include <stdlib.h>
#include <string.h>
#include "gtf.h"
#include "utils.h"

#define _insert(v, p, n, m, type) { \
    int _i, _flag=0;                  \
    for (_i = 0; _i < n; ++_i) {       \
        if (p[_i] == v) {            \
            _flag = 1;               \
            break;                  \
        }                           \
    }                               \
    if (_flag == 0) {                \
        if (n == m) {               \
            _realloc(p, m, type)    \
        }                           \
        p[n++] = v;                 \
    }                               \
}

#define _bin_insert(v, p, n, m, type) { \
    int _flag=0,_k_i=-1,_left=0,_right=n-1,_mid;    \
    type _mid_v, _tmp_v;                 \
    if (_right == -1) _k_i = 0;   \
    else {                      \
        while (_left <= _right) { \
            _mid = (_left+_right) >> 1;    \
            _mid_v = p[_mid];             \
            if (_mid_v == v) {           \
                _flag = 1; break;        \
            } else if (_mid_v > v) {     \
                if (_mid != 0) {         \
                    _tmp_v = p[_mid-1];   \
                }                       \
                if (_mid == 0 || v > _tmp_v) { \
                    _k_i = _mid;          \
                    break;              \
                }                       \
                else _right = _mid-1;     \
            } else _left = _mid+1;        \
        }                               \
    }                                   \
    if (_k_i == -1) _k_i = n;         \
                                    \
    if (_flag == 0) {                \
        if (n == m) {               \
            _realloc(p, m, type)    \
        }                           \
        if (_k_i <= n-1)             \
            memmove(p+_k_i+1, p+_k_i, (n-_k_i)*sizeof(type));  \
        (p)[_k_i] = v;               \
        (n)++;                      \
    }                               \
}

#define _bin_search(v, p, n, type, hit, i) { \
    int _left =0,_right=n-1,_mid; \
    type _mid_v;    \
    hit = 0;               \
    if (_right == -1) hit=0;   \
    else {  \
        while (_left <= _right) {   \
            _mid = (_left+_right) >> 1; \
            _mid_v = p[_mid];       \
            if (_mid_v == v) {  \
                i = _mid;   \
                hit = 1;   \
                break;      \
            } else if (_mid_v > v) {   \
                _right = _mid-1;    \
            } else {    \
                _left = _mid+1; \
            }   \
        }   \
    }   \
}

#define _node_len(n, i) ((n)[i].end-(n)[i].start+1)

typedef struct {
    int up, se, down;
    int asm_i, sg_i;
    int up_c, down_c, both_c, skip_c;
} SE_t;   // skipped exon

typedef struct {
    int lon, shor, down;
    int asm_i, sg_i;
    int lon_c, shor_c;
} A5SS_t; // alternative 3' splice site

typedef struct {
    int up, lon, shor;
    int asm_i, sg_i;
    int lon_c, shor_c;
} A3SS_t; // alternative 3' splice site

typedef struct {
    int up, fir, sec, down;
    int asm_i, sg_i;
    int fir_up_c, fir_down_c, fir_both_c, sec_up_c, sec_down_c, sec_both_c;
} MXE_t; // mutually exclusive exon

typedef struct {
    int up, down, in;
    int asm_i, sg_i;
    int sj_c;
} RI_t;  // retained intron

typedef struct {
    SE_t *se; int se_n, se_m;
    A5SS_t *a5ss; int a5ss_n, a5ss_m;
    A3SS_t *a3ss; int a3ss_n, a3ss_m;
    MXE_t *mxe; int mxe_n, mxe_m;
    RI_t *ri; int ri_n, ri_m;
} ASE_t;

typedef struct {
    int node_id, s_site_id, e_site_id; // unique id in corresponding gene-locus
    int start, end; /* real exon */ exon_t node_e;    // node in splice-graph
    uint8_t is_init, is_termi;
    uint8_t is_asm; int uniq_c, multi_c;
    int *next_id, next_n, next_m;
    int *pre_id, pre_n, pre_m;
    int *pre_domn, pre_domn_n, pre_domn_m;
    int *post_domn, post_domn_n, post_domn_m;
} SGnode; // node of splicing-graph

typedef struct {
    int site_id;
    int site;
    int *exon_id; int exon_n, exon_m;
} SGsite; // splice-site of splicing-graph

typedef struct {
    int don_site_id, acc_site_id;
    uint8_t is_rev;
    uint8_t motif, is_anno;
    int uniq_c, multi_c, max_over;
    anc_t *anc; // for each uniq-junction
} SGedge; // edge of splicing-graph, splice junction

typedef struct {
    //SGnode v;  // virtual start and end node
    // virtual_start: node[0]; virtual_end: node[node_n-1]
    SGnode *node; int node_n, node_m; // sort by e.start and e.end 
    SGsite *don_site; int don_site_n, don_site_m; // sort by site
    SGsite *acc_site; int acc_site_n, acc_site_m; // sort by site
    SGedge *edge; int edge_n, edge_m; // sort by don_id and acc_id
    int tid; uint8_t is_rev;
    // boundaries of splice-sites
    int start, end; 
} SG;

typedef struct {
    int SG_id;
    int v_start, v_end; // virtual start and end node
    int *node_id; int node_n, node_m;
    int *edge_id; int edge_n, edge_m;
    int start, end;
} SGasm;

typedef struct {
    SGasm **sg_asm;
    int sg_asm_n, sg_asm_m;
} SGasm_group;
 
typedef struct {
    SG **SG;
    int SG_n, SG_m;
    chr_name_t *cname;
} SG_group;

typedef struct {
    int sam_n, tol_rep_n, *rep_n;
    char **in_name;
    int no_novel_sj, no_novel_com, only_novel;
    int use_multi, read_type, anchor_len, intron_len;
    int merge_out;
} sg_para;

int comp_sj_sg(sj_t sj, SG sg);

#define sg_add_edge(ed, ei, ed_n, ed_m, _don_site_id, _acc_site_id, _is_rev, _is_anno) { \
    if (ed_n++ >= ed_m) _realloc(ed, ed_m, SGedge) \
    /* copy edge */ \
    if (ei <= ed_n-2) memmove(ed+ei+1, ed+ei, (ed_n-ei-1) * sizeof(SGedge)); \
    /* set edge */ \
    ed[ei].don_site_id = _don_site_id, ed[ei].acc_site_id = _acc_site_id,   \
    ed[ei].is_rev = _is_rev; ed[ei].is_anno = _is_anno;  \
    ed[ei].motif=0; ed[ei].uniq_c=0; ed[ei].multi_c=0; ed[ei].max_over=0; \
    ed[ei].anc=NULL; \
}

#define PAIR "paried"
#define SING "single"
#define PAIR_T 1
#define SING_T 0

int sg_par_input(sg_para *sgp, char *in);
sg_para *sg_init_para(void);
void sg_free_para(sg_para *sgp);

SG *sg_init_node(SG *sg);
SG *sg_init_site(SG *sg);
SG_group *sg_init_group(int g_n);
void sg_free_group(SG_group *sg_g);

int sg_update_node(SG *sg, exon_t e, int start, int end);
int sg_update_site(SG *sg, int site, uint8_t type);

int sg_bin_sch_node(SG *sg, exon_t e, int *hit);
int err_sg_bin_sch_node(const char *func, const int line, SG *sg, exon_t e, int *hit);
#define _err_sg_bin_sch_node(sg, e, hit) err_sg_bin_sch_node(__func__, __LINE__, sg, e, hit)
int sg_bin_sch_site(SGsite *site, int site_n, int s, int *hit);
int err_sg_bin_sch_site(const char *func, const int line, SGsite *site, int site_n, int s, int *hit);
#define _err_sg_bin_sch_site(site, site_n, s, hit) err_sg_bin_sch_site(__func__, __LINE__, site, site_n, s, hit)
int sg_bin_sch_edge(SG *sg, int don_site_id, int acc_site_id, int *hit);
int err_sg_bin_sch_edge(const char *func, const int line, SG *sg, int don_site_id, int acc_site_id, int *hit);
#define _err_sg_bin_sch_edge(sg, don_site_id, acc_site_id, hit) err_sg_bin_sch_edge(__func__, __LINE__, sg, don_site_id, acc_site_id, hit)

void cal_pre_domn(SG *sg);
void cal_post_domn(SG *sg);

SG_group *construct_SpliceGraph(FILE *gtf_fp, chr_name_t *cname);

int build_sg(int argc, char *argv[]);

#endif
