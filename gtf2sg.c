#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "splice_graph.h"
#include "utils.h"
#include "gtf.h"

extern char PROG[20];
int sg_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s build-sg [option] <in.gtf>\n\n", PROG);
    err_printf("Options:\n\n");
    err_printf("         -f --prefix  [STR]    file name to store splice-graph. [in.gtf]\n");
    err_printf("\n");

    return 0;
}

/***************************
 *     alloc and free      *
 ***************************/
SG *sg_init_node(SG *sg)
{
    int i;
    for (i = 0; i < sg->node_n; ++i) {
        sg->node[i].node_id = i; sg->node[i].s_site_id = sg->node[i].e_site_id = -1;
        sg->node[i].is_asm = 0; sg->node[i].uniq_c = 0; sg->node[i].multi_c = 0; sg->node[i].is_init = 0; sg->node[i].is_termi = 0;
        sg->node[i].next_n = 0; sg->node[i].next_m = 1;
        sg->node[i].next_id = (gec_t*)_err_malloc(sizeof(gec_t));
        sg->node[i].pre_n = 0; sg->node[i].pre_m = 1;
        sg->node[i].pre_id = (gec_t*)_err_malloc(sizeof(gec_t));
        sg->node[i].pre_domn_n = 1; sg->node[i].pre_domn_m = 2;
        sg->node[i].pre_domn = (gec_t*)_err_malloc(2 * sizeof(gec_t)); sg->node[i].pre_domn[0] = i;
        sg->node[i].post_domn_n = 1; sg->node[i].post_domn_m = 2;
        sg->node[i].post_domn = (gec_t*)_err_malloc(2 * sizeof(gec_t)); sg->node[i].post_domn[0] = i;
    }
    return sg;
}

SG *sg_init_site(SG *sg)
{
    int i;
    for (i = 0; i < sg->don_site_n; ++i) {
        sg->don_site[i].site_id = i;
        sg->don_site[i].exon_n = 0; sg->don_site[i].exon_m = 1;
        sg->don_site[i].exon_id = (int*)_err_malloc(sizeof(int));
    }
    for (i = 0; i < sg->acc_site_n; ++i) {
        sg->acc_site[i].site_id = i;
        sg->acc_site[i].exon_n = 0; sg->acc_site[i].exon_m = 1;
        sg->acc_site[i].exon_id = (int*)_err_malloc(sizeof(int));
    }
    return sg;
}

SG *sg_init(void)
{
    SG *sg = (SG*)_err_malloc(sizeof(SG));
    sg->node_n = 0, sg->node_m = 2;
    sg->node = (SGnode*)_err_malloc(2 * sizeof(SGnode));
    sg->don_site_n = 0, sg->don_site_m = 2;
    sg->don_site = (SGsite*)_err_malloc(2 * sizeof(SGsite));
    sg->acc_site_n = 0, sg->acc_site_m = 2;
    sg->acc_site = (SGsite*)_err_malloc(2 * sizeof(SGsite));
    sg->edge_n = 0, sg->edge_m = 2;
    sg->edge = (SGedge*)_err_malloc(2 * sizeof(SGedge));

    sg->start = MAX_SITE, sg->end = 0;
    return sg;
}

SG_group *sg_init_group(int g_n)
{
    SG_group *sg_g = (SG_group*)_err_malloc(sizeof(SG_group));
    sg_g->SG_n = g_n, sg_g->SG_m = g_n;
    sg_g->SG = (SG**)_err_malloc(g_n * sizeof(SG*));
    sg_g->cname = chr_name_init();
    int i; for (i = 0; i < g_n; ++i) sg_g->SG[i] = sg_init();
    return sg_g;
}

SG_group *sg_realloc_group(SG_group *sg_g)
{
    sg_g->SG_m <<= 1;
    sg_g->SG = (SG**)_err_realloc(sg_g->SG, sg_g->SG_m * sizeof(SG*));
    int i; for (i = (sg_g->SG_m >> 1); i < sg_g->SG_m; ++i) sg_g->SG[i] = sg_init();
    return sg_g;
}

void sg_free_node(SG *sg) 
{ 
    int i;
    for (i = 0; i < sg->node_n; ++i) {
        free(sg->node[i].next_id); free(sg->node[i].pre_id); 
        free(sg->node[i].pre_domn); free(sg->node[i].post_domn); 
    }
    free(sg->node);
}

void sg_free_site(SG *sg)
{
    int i; 
    for (i = 0; i < sg->don_site_n; ++i) free(sg->don_site[i].exon_id);
    for (i = 0; i < sg->acc_site_n; ++i) free(sg->acc_site[i].exon_id);
    free(sg->don_site); free(sg->acc_site);
}

void sg_free(SG *sg)
{
    sg_free_node(sg); sg_free_site(sg); 
    free(sg->edge); free(sg);
}

void sg_free_group(SG_group *sg_g)
{
    int i; for (i = 0; i < sg_g->SG_m; i++) sg_free(sg_g->SG[i]);
    free(sg_g->SG); chr_name_free(sg_g->cname);
    free(sg_g);
}
/***************************/

/****************************************
 * construct splice graph from GTF file *
 ****************************************/
// binary search node
// compare start, then end
int sg_bin_sch_node(SG *sg, exon_t e, int *hit)
{
    *hit = 0;
    int start = e.start, end = e.end, mid_s, mid_e, tmp_s, tmp_e;
    int left = 0, right = sg->node_n-1, mid;
    if (right == -1) return 0;

    while (left <= right) {
        mid = ((left + right) >> 1);
        mid_s = sg->node[mid].node_e.start, mid_e = sg->node[mid].node_e.end;
        if (mid_s == start && mid_e == end) { *hit = 1; return mid; }
        else if (mid_s > start || (mid_s == start && mid_e > end)) { // [mid] is bigger than query
            if (mid != 0) {
                tmp_s = sg->node[mid-1].node_e.start, tmp_e = sg->node[mid-1].node_e.end;
            }
            if (mid == 0 || (start > tmp_s || (start == tmp_s && end > tmp_e))) {
                return mid;
            } else right = mid-1;
        } else left = mid + 1;
    }
    return sg->node_n;
}

int err_sg_bin_sch_node(const char *func, const int line, SG *sg, exon_t e)
{
    int hit;
    int id = sg_bin_sch_node(sg, e, &hit);
    char head[100]; sprintf(head, "%s:%d", func, line);
    if (hit == 0) _err_fatal_simple(head, "Can not hit node.\n");
    return id;
}

int sg_update_node(SG *sg, exon_t e, int start, int end)
{
    int hit = 0;
    int n_i = sg_bin_sch_node(sg, e, &hit);
    if (hit == 0) { // insert new node
        if (sg->node_n++ >= sg->node_m) _realloc(sg->node, sg->node_m, SGnode)
        // copy node
        if (n_i <= sg->node_n-2)
            memmove(sg->node+n_i+1, sg->node+n_i, (sg->node_n-n_i-1) * sizeof(SGnode));
        // set node
        sg->node[n_i].node_e = e;
        sg->node[n_i].start = start;
        sg->node[n_i].end = end;
        if (start != 0 && start < sg->start) sg->start = start;
        if (end != MAX_SITE && end > sg->end) sg->end = end;
    } else {
        if (start < sg->node[n_i].start) sg->node[n_i].start = start;
        if (end > sg->node[n_i].end) sg->node[n_i].end = end;
        if (start != 0 && start < sg->start) sg->start = start;
        if (end != MAX_SITE && end > sg->end) sg->end = end;
    }
    return 0;
}

int sg_bin_sch_site(SGsite *site, int site_n, int s, int *hit)
{
    *hit = 0;
    int mid_s, tmp_s;
    int left = 0, right = site_n-1, mid;
    if (right == -1) return 0;

    while (left <= right) {
        mid = ((left + right) >> 1);
        mid_s = site[mid].site;
        if (mid_s == s) { *hit = 1; return mid; }
        else if (mid_s > s) { // [mid] is bigger than query
            if (mid != 0) tmp_s = site[mid-1].site;
            if (mid == 0 || s > tmp_s ) return mid;
            else right = mid-1;
        } else left = mid + 1;
    }
    return site_n;
}

int err_sg_bin_sch_site(const char *func, const int line, SGsite *site, int site_n, int s)
{
    int hit;
    int id = sg_bin_sch_site(site, site_n, s, &hit);
    char head[100]; sprintf(head, "%s:%d", func, line);
    if (hit == 0) _err_fatal_simple(head, "Can not hit site.\n");
    return id;
}

int sg_update_site(SG *sg, int site, uint8_t type)
{
    int hit = 0;
    if (type == DON_SITE_F) {
        int s_i = sg_bin_sch_site(sg->don_site, sg->don_site_n, site, &hit);
        if (hit == 0) {
            //if (site < sg->start) sg->start = site;

            if (sg->don_site_n++ >= sg->don_site_m) _realloc(sg->don_site, sg->don_site_m, SGsite)
                // copy site
                if (s_i <= sg->don_site_n-2)
                    memmove(sg->don_site+s_i+1, sg->don_site+s_i, (sg->don_site_n-s_i-1) * sizeof(SGsite));
            // set site
            sg->don_site[s_i].site = site;
        }
    } else {
        int s_i = sg_bin_sch_site(sg->acc_site, sg->acc_site_n, site, &hit);
        if (hit == 0) {
            //if (site > sg->end) sg->end = site;

            if (sg->acc_site_n++ >= sg->acc_site_m) _realloc(sg->acc_site, sg->acc_site_m, SGsite)
                // copy site
                if (s_i <= sg->acc_site_n-2)
                    memmove(sg->acc_site+s_i+1, sg->acc_site+s_i, (sg->acc_site_n-s_i-1) * sizeof(SGsite));
            // set site
            sg->acc_site[s_i].site = site;
        }
    }
    return 0;
}

int sg_bin_sch_edge(SG *sg, int don_id, int acc_id, int *hit)
{
    *hit = 0;
    int mid_d, mid_a, tmp_d, tmp_a;
    int left = 0, right = sg->edge_n-1, mid;
    if (right == -1) return 0;

    while (left <= right) {
        mid = ((left + right) >> 1);
        mid_d = sg->edge[mid].don_id, mid_a = sg->edge[mid].acc_id;
        if (mid_d == don_id && mid_a == acc_id) { *hit = 1; return mid; }
        else if (mid_d > don_id || (mid_d == don_id && mid_a > acc_id)) { // [mid] is bigger than query
            if (mid != 0) {
                tmp_d = sg->edge[mid-1].don_id, tmp_a = sg->edge[mid-1].acc_id;
            }
            if (mid == 0 || (don_id > tmp_d || (don_id == tmp_d && acc_id > tmp_a))) {
                return mid;
            } else right = mid - 1;
        } else left = mid + 1;
    }
    return sg->edge_n;
}
int err_sg_bin_sch_edge(const char *func, const int line, SG *sg, int don_id, int acc_id)
{
    int hit;
    int id = sg_bin_sch_edge(sg, don_id, acc_id, &hit);
    char head[100]; sprintf(head, "%s:%d", func, line);
    if (hit == 0) _err_fatal_simple(head, "Can not hit edge.\n");
    return id;
}

// update edge of splicing-graph
int sg_update_edge(SG *sg, int don_id, int acc_id, int don_site_id, int acc_site_id, uint8_t is_rev)
{
    int hit = 0;
    int e_i = sg_bin_sch_edge(sg, don_id, acc_id, &hit);
    if (hit == 0) { // insert new edge
        uint8_t is_anno=1;
        sg_add_edge(sg->edge, e_i, (sg->edge_n), (sg->edge_m), don_id, acc_id, is_rev, is_anno)
    }
    // set next/pre
    _bin_insert(acc_id, sg->node[don_id].next_id, sg->node[don_id].next_n, sg->node[don_id].next_m, gec_t)
    _bin_insert(don_id, sg->node[acc_id].pre_id, sg->node[acc_id].pre_n, sg->node[acc_id].pre_m, gec_t)
    // set site
    _bin_insert(don_id, sg->don_site[don_site_id].exon_id, sg->don_site[don_site_id].exon_n, sg->don_site[don_site_id].exon_m, int)
    _bin_insert(acc_id, sg->acc_site[acc_site_id].exon_id, sg->acc_site[acc_site_id].exon_n, sg->acc_site[acc_site_id].exon_m, int)
    return 0;
}

// order: 1, pre; 2, post
void intersect_domn(gec_t **com, gec_t *new_domn, gec_t *com_n, gec_t new_n, int order)
{
    int i, j, domn_i=0;
    for (i=0, j=0; i<*com_n && j<new_n; ) {
        if ((*com)[i] == new_domn[j]) {
            (*com)[domn_i++] = (*com)[i];
            i++, j++;
        } else if ((*com)[i] < new_domn[j]) {
            if (order == 1) j++;
            else i++;
        } else {
            if (order == 1) i++;
            else j++;
        }
    }
    *com_n = domn_i;
}

// calculate predominate nodes with edge weight
// XXX require at least one rep has >0 weight
void cal_pre_domn(SG *sg, double **rep_W, uint8_t ***asm_matrix)
{
    int i, j, first;
    for (i = 0; i < sg->node_n; ++i) {
        if (sg->node[i].pre_n == 0) continue;

        gec_t com_n, com_m, *com;

        for (first = 0; first < sg->node[i].pre_n; ++first) {
            if (sg->node[i].pre_id[first] == 0 || ((i == sg->node_n-1 || rep_W[sg->node[i].pre_id[first]][i] > 0) && sg->node[sg->node[i].pre_id[first]].pre_domn_n > 1))
                goto pre_domn;
        }
        continue;

pre_domn:
        com_n = sg->node[sg->node[i].pre_id[first]].pre_domn_n, com_m = sg->node[sg->node[i].pre_id[first]].pre_domn_m;
        com = (gec_t*)_err_malloc(com_m * sizeof(gec_t));

        for (j = 0; j < com_n; ++j) com[j] = sg->node[sg->node[i].pre_id[first]].pre_domn[j];

        asm_matrix[0][sg->node[i].pre_id[first]][i] = 1;
        for (j = first+1; j < sg->node[i].pre_n; ++j) {
            if (sg->node[i].pre_id[j] == 0 || ((i == sg->node_n-1 || rep_W[sg->node[i].pre_id[j]][i] > 0) && sg->node[sg->node[i].pre_id[j]].pre_domn_n > 1))
                goto pre_domn_x;
            else continue;
pre_domn_x:
            asm_matrix[0][sg->node[i].pre_id[j]][i] = 1;
            intersect_domn(&com, sg->node[sg->node[i].pre_id[j]].pre_domn, &com_n, sg->node[sg->node[i].pre_id[j]].pre_domn_n, 1);
        }
        if (com_n+1 > sg->node[i].pre_domn_m) {
            sg->node[i].pre_domn_m = com_n+1; 
            sg->node[i].pre_domn = (gec_t*)_err_realloc(sg->node[i].pre_domn, (com_n+1) * sizeof(gec_t));
        }
        for (j = 0; j < com_n; ++j) sg->node[i].pre_domn[j+1] = com[j];
        sg->node[i].pre_domn_n = 1+com_n;
        free(com);
    }
}

void cal_post_domn(SG *sg, double **rep_W, uint8_t ***asm_matrix)
{
    int i, j, first;
    for (i = sg->node_n-1; i >= 0; --i) {
        if (sg->node[i].next_n == 0) continue;
            
        gec_t com_n, com_m, *com;
        first = 0;
        for (first = 0; first < sg->node[i].next_n; ++first) {
            if (sg->node[i].next_id[first] == sg->node_n-1 || ((i == 0 || rep_W[i][sg->node[i].next_id[first]] > 0) && sg->node[sg->node[i].next_id[first]].post_domn_n > 1))
                goto post_domn;
        }
        continue;

post_domn:
        com_n = sg->node[sg->node[i].next_id[first]].post_domn_n, com_m = sg->node[sg->node[i].next_id[first]].post_domn_m;
        com = (gec_t*)_err_malloc(com_m * sizeof(gec_t));
        for (j = 0; j < com_n; ++j) com[j] = sg->node[sg->node[i].next_id[first]].post_domn[j];

        asm_matrix[0][i][sg->node[i].next_id[first]] |= 2;
        for (j = first+1; j < sg->node[i].next_n; ++j) {
            if (sg->node[i].next_id[j] == sg->node_n-1 || ((i == 0 || rep_W[i][sg->node[i].next_id[j]] > 0) && sg->node[sg->node[i].next_id[j]].post_domn_n > 1))
                goto post_domn_x;
            else continue;
post_domn_x:
            asm_matrix[0][i][sg->node[i].next_id[j]] |= 2;
            intersect_domn(&com, sg->node[sg->node[i].next_id[j]].post_domn, &com_n, sg->node[sg->node[i].next_id[j]].post_domn_n, 2);
        }

        if (com_n+1 > sg->node[i].post_domn_m) {
            sg->node[i].post_domn_m = com_n+1; 
            sg->node[i].post_domn = (gec_t*)_err_realloc(sg->node[i].post_domn, (com_n+1) * sizeof(gec_t));
        }
        for (j = 0; j < com_n; ++j) sg->node[i].post_domn[j+1] = com[j];
        sg->node[i].post_domn_n = com_n+1;
        free(com);
    }
}

typedef struct {
    int start, end;
    int n, m, *s;
} coor_pair_t;

// 1st: merge all connected-exons
// 2nd: use node to split exon
// 3rd: use boundary to split node
// 4th: use splited node to update sg
int merge_exon(coor_pair_t **coor, int *coor_n, int *coor_m, int start, int end) {
    if (*coor_n == 0) {
        (*coor)[0].start = start, (*coor)[0].end = end;
        (*coor)[0].n = 0, (*coor)[0].m = 1; (*coor)[0].s = (int*)_err_malloc(sizeof(int));
        *coor_n = 1;
    } else {
        if (start <= (*coor)[*coor_n-1].end+1) {
            if (end > (*coor)[*coor_n-1].end) (*coor)[*coor_n-1].end = end;
        }
        else { // add new coor
           if (*coor_n == *coor_m) _realloc(*coor, *coor_m, coor_pair_t)
           (*coor)[*coor_n].start = start, (*coor)[*coor_n].end = end;
           (*coor)[*coor_n].n = 0; (*coor)[*coor_n].m = 1, (*coor)[*coor_n].s = (int*)_err_malloc(sizeof(int));
           (*coor_n)++;
        }
    }
    return *coor_n;
}

// split exon into smaller exons
void coor_push_s(coor_pair_t *coor, int start1, int start2) {
    _bin_insert(start1, coor->s, coor->n, coor->m, int)
    _bin_insert(start2, coor->s, coor->n, coor->m, int)
}

// if node is contained in coor, cut coor
void split_exon(coor_pair_t *coor, int coor_n, int start, int end, int *coor_i) {
    int i;
    for (i = *coor_i; i < coor_n; ++i) {
        if (coor[i].start <= start && coor[i].end >= end) {
            coor_push_s(coor+i, start, end+1);
            *coor_i = i;
            break;
        }
    }
}

void gen_split_node(SG *sg, coor_pair_t *coor, int coor_n) {
    sg->node_n = 0;
    sg_update_node(sg, (exon_t){sg->tid, sg->is_rev, 0, 0}, 0, 0);
    int i, j;
    for (i = 0; i < coor_n; ++i) {
        for (j = 0; j < coor[i].n-1; ++j) {
            sg_update_node(sg, (exon_t){sg->tid, sg->is_rev, coor[i].s[j], coor[i].s[j+1]-1}, coor[i].s[j], coor[i].s[j+1]-1); 
            sg_update_site(sg, coor[i].s[j]-1, ACC_SITE_F);
            sg_update_site(sg, coor[i].s[j+1], DON_SITE_F);
        }
    }
    sg_update_node(sg, (exon_t){sg->tid, sg->is_rev, MAX_SITE, MAX_SITE}, MAX_SITE, MAX_SITE);
}

void sg_merge_end(gene_t *gene) {
    int *end = (int*)_err_malloc(5 * sizeof(int)), n1 = 0, m1 = 5;
    int *start = (int*)_err_malloc(5 * sizeof(int)), n2 = 0, m2 = 5;
    int i, j, min, max;
    for (i = 0; i < gene->trans_n; ++i) {
        _bin_insert(gene->trans[i].exon[0].end, end, n1, m1, int)
        _bin_insert(gene->trans[i].exon[gene->trans[i].exon_n-1].start, start, n2, m2, int)
    }
    for (i = 0; i < n1; ++i) {
        min = MAX_SITE;
        for (j = 0; j < gene->trans_n; ++j) {
            if (gene->trans[j].exon[0].end != end[i]) continue;
            min = MIN_OF_TWO(min, gene->trans[j].exon[0].start);
        }
        for (j = 0; j < gene->trans_n; ++j) { // merge init-exon's start, which share same end
            if (gene->trans[j].exon[0].end != end[i]) continue;
            gene->trans[j].exon[0].start = min;
        }
    }
    for (i = 0; i < n2; ++i) {
        max = 0;
        for (j = 0; j < gene->trans_n; ++j) {
            if (gene->trans[j].exon[gene->trans[j].exon_n-1].start!= start[i]) continue;
            max = MAX_OF_TWO(max, gene->trans[j].exon[gene->trans[j].exon_n-1].end);
        }
        for (j = 0; j < gene->trans_n; ++j) { // merge term-exon' end, which share same start
            if (gene->trans[j].exon[gene->trans[j].exon_n-1].start!= start[i]) continue;
            gene->trans[j].exon[gene->trans[j].exon_n-1].end = max;
        }
    }
    free(start), free(end);
}

void gen_gene_split_exon(coor_pair_t *coor, int coor_n, gene_t *gene) {
    int i, j, k, p, q;
    for (i = 0; i < gene->trans_n; ++i) {
        for (j = 0; j < gene->trans[i].exon_n; ++j) {
            int start = gene->trans[i].exon[j].start, end = gene->trans[i].exon[j].end;
            for (k = 0; k < coor_n; ++k) {
                if (start > coor[k].end) continue;
                else if (end < coor[k].start) break;
                for (p = 0; p < coor[k].n; ++p) {
                    if (coor[k].s[p] > gene->trans[i].exon[j].start) {
                        gene->trans[i].exon[j].end = coor[k].s[p]-1;
                        exon_t e = gene->trans[i].exon[j];
                        for (q = p; q < coor[k].n-1 && coor[k].s[q] < end; ++q) {
                            e.start = coor[k].s[q], e.end = coor[k].s[q+1]-1;
                            if (gene->trans[i].exon_n == gene->trans[i].exon_m) _realloc(gene->trans[i].exon, gene->trans[i].exon_m, exon_t)
                            gene->trans[i].exon[gene->trans[i].exon_n++] = e;
                        }
                        goto next_exon;
                    }
                }
            }
next_exon:;
        }
        sort_exon(gene->trans+i);
    }
}

int gen_sg_node(SG *sg, gene_t *gene)
{
    int i, j; exon_t e;
    // merge init/term exon
    // XXX should be optional
    sg_merge_end(gene);
    // generate node
    for (i = 0; i < gene->trans_n; ++i) {
        for (j = 0; j < gene->trans[i].exon_n; ++j) {
            e = gene->trans[i].exon[j];
            sg_update_node(sg, e, e.start, e.end);
        }
    }
    int coor_n = 0, coor_i, coor_m = gene->trans[0].exon_n; coor_pair_t *coor = (coor_pair_t*)_err_malloc(coor_m * sizeof(coor_pair_t));
    // merege into exonic region
    for (i = 0; i < sg->node_n; ++i)
        merge_exon(&coor, &coor_n, &coor_m, sg->node[i].start, sg->node[i].end);
    // split gene-exon into splice-exon
    coor_i = 0;
    for (i = 0; i < sg->node_n; ++i)
        split_exon(coor, coor_n, sg->node[i].start, sg->node[i].end, &coor_i);
    // generate split-node
    gen_split_node(sg, coor, coor_n);

    // split exon in gene
    gen_gene_split_exon(coor, coor_n, gene);

    for (i = 0; i < coor_n; ++i) free(coor[i].s); free(coor);
    return 0;
}

// construct splice-graph for each gene
// new definition for exon, might be partial exon, not the real exon
void build_SpliceGraph_core(SG *sg, gene_t *gene)
{
    int i, j, don_id, acc_id, don_site_id=0, acc_site_id=0; exon_t e;
    
    sg->tid = gene->tid, sg->is_rev = gene->is_rev;
    // gen sg node (pseu-exon)
    gen_sg_node(sg, gene);

    // alloc for next_id/pre_id/pre_domn/post_domn
    sg_init_node(sg); sg_init_site(sg);

    // search node and generate edge 
    for (i = 0; i < gene->trans_n; ++i) {
        // keep 1 exon transcript ???
        // XXX if (gene->trans[i].exon_n == 1) continue;
        e = gene->trans[i].exon[0];

        acc_id = don_id = _err_sg_bin_sch_node(sg, e);
        acc_site_id = _err_sg_bin_sch_site(sg->acc_site, sg->acc_site_n, e.start-1);
        don_site_id = _err_sg_bin_sch_site(sg->don_site, sg->don_site_n, e.end+1);
        sg->node[don_id].s_site_id = acc_site_id;
        sg->node[don_id].e_site_id = don_site_id;

        // set next_id of v_start
        _bin_insert(don_id, sg->node[0].next_id, sg->node[0].next_n, sg->node[0].next_m, gec_t)
        _bin_insert(0, sg->node[don_id].pre_id, sg->node[don_id].pre_n, sg->node[don_id].pre_m, gec_t)

        // set edge
        int hit = 0;
        int e_i = sg_bin_sch_edge(sg, 0, acc_id, &hit);
        if (hit == 0) // insert new edge
            sg_add_edge(sg->edge, e_i, (sg->edge_n), (sg->edge_m), 0, acc_id, sg->is_rev, 1);
        // set site
        _bin_insert(acc_id, sg->acc_site[acc_site_id].exon_id, sg->acc_site[acc_site_id].exon_n, sg->acc_site[acc_site_id].exon_m, int)

        sg->node[don_id].is_init = 1;

        for (j = 1; j < gene->trans[i].exon_n; ++j) {
            e = gene->trans[i].exon[j];

            acc_id = _err_sg_bin_sch_node(sg, e);
            acc_site_id = _err_sg_bin_sch_site(sg->acc_site, sg->acc_site_n, e.start-1);
            sg->node[acc_id].s_site_id = acc_site_id;

            sg_update_edge(sg, don_id, acc_id, don_site_id, acc_site_id, gene->is_rev);

            don_site_id = _err_sg_bin_sch_site(sg->don_site, sg->don_site_n, e.end+1);
            sg->node[acc_id].e_site_id = don_site_id;
            don_id = acc_id;
        }

        // set pre_id of v_end
        _bin_insert(acc_id, sg->node[sg->node_n-1].pre_id, sg->node[sg->node_n-1].pre_n, sg->node[sg->node_n-1].pre_m, gec_t)
        _bin_insert((gec_t)sg->node_n-1, sg->node[acc_id].next_id, sg->node[acc_id].next_n, sg->node[acc_id].next_m, gec_t)

        // set edge
        hit = 0;
        e_i = sg_bin_sch_edge(sg, acc_id, sg->node_n-1, &hit);
        if (hit == 0) // insert new edge
            sg_add_edge(sg->edge, e_i, (sg->edge_n), (sg->edge_m), acc_id, sg->node_n-1, sg->is_rev, 1);
        // set site
        _bin_insert(don_id, sg->don_site[don_site_id].exon_id, sg->don_site[don_site_id].exon_n, sg->don_site[don_site_id].exon_m, int)


        sg->node[acc_id].is_termi = 1;
    }
    // after update!!! cal pre/post domn
    // cal_pre_domn(sg); cal_post_domn(sg); 
}

SG_group *construct_SpliceGraph(FILE *gtf, chr_name_t *cname)
{
    err_func_format_printf(__func__, "constructing splice-graph with GTF file ...\n", __func__);
    gene_group_t *gg = gene_group_init();
    // XXX optional: merge overlapping genes into 1 complex gene
    int g_n = read_gene_group(gtf, cname, gg);
    SG_group *sg_g = sg_init_group(g_n);
    int i;
    if (cname->chr_n > sg_g->cname->chr_m) {
        sg_g->cname->chr_name = (char**)_err_realloc(sg_g->cname->chr_name, cname->chr_n * sizeof(char*));
        for (i = sg_g->cname->chr_m; i < cname->chr_n; ++i) sg_g->cname->chr_name[i] = (char*)_err_malloc(100 * sizeof(char));
        sg_g->cname->chr_m = cname->chr_n;
    }
    sg_g->cname->chr_n = cname->chr_n;
    for (i = 0; i < sg_g->cname->chr_n; ++i) strcpy(sg_g->cname->chr_name[i], cname->chr_name[i]);

    for (i = 0; i < gg->gene_n; ++i) build_SpliceGraph_core(sg_g->SG[i], gg->g+i);

    gene_group_free(gg);
    err_func_format_printf(__func__, "constructing splice-graph with GTF file done!\n");
    return sg_g;
}

/****************************************/

#ifdef __BUILD_SG_MAIN__
const struct option sg_long_opt [] = {
    {"prefix", 1, NULL, 'f' },
    {0, 0, 0, 0}
};

int main(int argc, char *argv[])
{
    int c;
    char prefix[1024] = "";
    while ((c = getopt_long(argc, argv, "f:", sg_long_opt, NULL)) >= 0) {
        switch (c) {
            case 'f': strcpy(prefix, optarg); break;
            default: err_printf("Error: unknown option: %s.\n", optarg);
                     return sg_usage();
        }
    }
    if (argc - optind != 1) return sg_usage();
    FILE *gtf_fp = xopen(argv[optind], "r");
    if (strlen(prefix) == 0) strcpy(prefix, argv[optind]);

    chr_name_t *cname = chr_name_init();
    SG_group *sg_g = construct_SpliceGraph(gtf_fp, cname);

    sg_free_group(sg_g); err_fclose(gtf_fp); chr_name_free(cname); 
    return 0;
}
#endif
