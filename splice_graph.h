#ifndef _SG_H
#define _SG_H
#include <stdlib.h>
#include <string.h>
#include "gtf.h"
#include "utils.h"

//XXX uint16_t for exon-number
// exon-number should not exceed pow(2,16) (65536)
#define gtf_exon_cnt_t uint16_t
#define gec_t gtf_exon_cnt_t

#define _node_len(n) ((n).end-(n).start+1)

typedef struct {
    uint32_t tid:16, intv_n:14, is_uniq:1, is_splice:1;
    int32_t start, end;
    int32_t *exon_end, *intr_end; // exon_end[intv_n]: exonic, intr_end[intv_n-1]: intronic
} ad_t;   // alignment details: start, end, intv_n, intv[]

typedef struct {
    int up, se, down;
    int asm_i, sg_i;
    int up_c, down_c, ud_both_c, skip_c;
    int body_c;
} SE_t;   // skipped exon

typedef struct {
    int lon, shor, down;
    int asm_i, sg_i;
    int shor_c, pj_c, lon_c, pl_both_c;
    int body_c;
} A5SS_t; // alternative 3' splice site

typedef struct {
    int up, lon, shor;
    int asm_i, sg_i;
    int lon_c, pj_c, lp_both_c, shor_c;
    int body_c;
} A3SS_t; // alternative 3' splice site

typedef struct {
    int up, fir, sec, down;
    int asm_i, sg_i;
    int fir_up_c, fir_down_c, fir_both_c, sec_up_c, sec_down_c, sec_both_c;
    int fir_body_c, sec_body_c;
} MXE_t; // mutually exclusive exon

typedef struct {
    int up, down, in;
    int asm_i, sg_i;
    int ej_c, pj1_c, pj2_c, pj_both_c;
    int body_c;
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
    int start, end; // real exon
    exon_t node_e;  // node in splice-graph
                    // NOW: no differences between (start,end) and node_n
    uint8_t is_init:1, is_termi:1, is_asm:1;
    int uniq_c, multi_c;
    gec_t *next_id, next_n; int next_m;
    gec_t *pre_id, pre_n; int pre_m;
    gec_t *pre_domn, pre_domn_n; int pre_domn_m;
    gec_t *post_domn, post_domn_n; int post_domn_m;
} SGnode; // node of splicing-graph

typedef struct {
    int site_id;
    int site;
    int *exon_id; int exon_n, exon_m; // in new definition, one site <-> one node
} SGsite; // splice-site of splicing-graph

typedef struct {
    int don_id, acc_id; // don/acc node id
    uint8_t is_rev:2, is_anno:2, motif:4;
    int uniq_c, multi_c, max_over;
} SGedge; // edge of splicing-graph: splice junction

typedef struct {
    // virtual_start: node[0]; virtual_end: node[node_n-1]
    SGnode *node; int node_n, node_m; // sort by e.start and e.end 
    SGsite *don_site; int don_site_n, don_site_m; // sort by site
    SGsite *acc_site; int acc_site_n, acc_site_m; // sort by site
    SGedge *edge; int edge_n, edge_m; // sort by don_id and acc_id
    int tid; uint8_t is_rev;
    int start, end; // boundaries of splice-sites
} SG;

typedef struct {
    int SG_id;
    int v_start, v_end; // virtual start and end node
    gec_t *node_id, node_n; int node_m;
    int *edge_id, edge_n; int edge_m;
    int start, end;
} SGasm;

typedef struct {
    int ASM_id, SG_id;
    gec_t v_start, v_end; // virtual start and end node
    gec_t *node_id;
    int *uniq_sj_c, *uniq_tot_c, *multi_sj_c, *multi_tot_c;
} SGiso; // each ASM has one SGiso

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
    int n_threads;
    int sam_n, tot_rep_n, *rep_n;
    char **in_name;
    uint8_t no_novel_sj:1, no_novel_com:1, only_novel:1, use_multi:1, read_type:1, merge_out:1, rm_edge:1;
    int intron_len, edge_wt;
    uint64_t iso_cnt_max; int asm_exon_max, iso_read_cnt_min;
    int anchor_len[5]; // [anno, non-canonical, GT/AG, GC/AG, AT/AC]
    int uniq_min[5];   // [anno, non-canonical, GT/AG, GC/AG, AT/AC]
    int all_min[5];    // [anno, non-canonical, GT/AG, GC/AG, AT/AC]
} sg_para;

int comp_sj_sg(sj_t sj, SG sg);

#define sg_add_edge(ed, ei, ed_n, ed_m, _don_id, _acc_id, _is_rev, _is_anno) { \
    if (ed_n++ >= ed_m) _realloc(ed, ed_m, SGedge) \
    /* copy edge */ \
    if (ei <= ed_n-2) memmove(ed+ei+1, ed+ei, (ed_n-ei-1) * sizeof(SGedge)); \
    /* set edge */ \
    ed[ei].don_id = _don_id, ed[ei].acc_id = _acc_id,   \
    ed[ei].is_rev = _is_rev; ed[ei].is_anno = _is_anno;  \
    ed[ei].motif=0; ed[ei].uniq_c=0; ed[ei].multi_c=0; ed[ei].max_over=0; \
}

#define PAIR "paried"
#define SING "single"
#define PAIR_T 1
#define SING_T 0

sg_para *sg_init_para(void);
void sg_free_para(sg_para *sgp);

SG *sg_init_node(SG *sg);
SG *sg_init_site(SG *sg);
SG_group *sg_init_group(int g_n);
void sg_free_group(SG_group *sg_g);

int sg_update_node(SG *sg, exon_t e, int start, int end);
int sg_update_site(SG *sg, int site, uint8_t type);

int sg_bin_sch_node(SG *sg, exon_t e, int *hit);
int err_sg_bin_sch_node(const char *func, const int line, SG *sg, exon_t e);
#define _err_sg_bin_sch_node(sg, e ) err_sg_bin_sch_node(__func__, __LINE__, sg, e)
int sg_bin_sch_site(SGsite *site, int site_n, int s, int *hit);
int err_sg_bin_sch_site(const char *func, const int line, SGsite *site, int site_n, int s);
#define _err_sg_bin_sch_site(site, site_n, s) err_sg_bin_sch_site(__func__, __LINE__, site, site_n, s)
int sg_bin_sch_edge(SG *sg, int don_id, int acc_id, int *hit);
int err_sg_bin_sch_edge(const char *func, const int line, SG *sg, int don_id, int acc_id);
#define _err_sg_bin_sch_edge(sg, don_id, acc_id) err_sg_bin_sch_edge(__func__, __LINE__, sg, don_id, acc_id)

void cal_pre_domn(SG *sg);
void cal_post_domn(SG *sg);

SG_group *construct_SpliceGraph(FILE *gtf_fp, chr_name_t *cname);

int build_sg(int argc, char *argv[]);

#endif
