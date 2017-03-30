#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "build_sg.h"
#include "utils.h"

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
        sg->node[i].next_id = (uint32_t*)_err_malloc(sizeof(uint32_t));
        sg->node[i].pre_n = 0; sg->node[i].pre_m = 1;
        sg->node[i].pre_id = (uint32_t*)_err_malloc(sizeof(uint32_t));
        sg->node[i].pre_domn_n = 1; sg->node[i].pre_domn_m = 2;
        sg->node[i].pre_domn = (uint32_t*)_err_malloc(2 * sizeof(uint32_t)); sg->node[i].pre_domn[0] = i;
        sg->node[i].post_domn_n = 1; sg->node[i].post_domn_m = 2;
        sg->node[i].post_domn = (uint32_t*)_err_malloc(2 * sizeof(uint32_t)); sg->node[i].post_domn[0] = i;
    }
    return sg;
}

SG *sg_init_site(SG *sg)
{
    int i;
    for (i = 0; i < sg->don_site_n; ++i) {
        sg->don_site[i].site_id = i;
        sg->don_site[i].exon_n = 0; sg->don_site[i].exon_m = 1;
        sg->don_site[i].exon_id = (uint32_t*)_err_malloc(sizeof(uint32_t));
    }
    for (i = 0; i < sg->acc_site_n; ++i) {
        sg->acc_site[i].site_id = i;
        sg->acc_site[i].exon_n = 0; sg->acc_site[i].exon_m = 1;
        sg->acc_site[i].exon_id = (uint32_t*)_err_malloc(sizeof(uint32_t));
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
    sg_free_node(sg); sg_free_site(sg); free(sg->edge);
    free(sg);
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
    int32_t start = e.start, end = e.end, mid_s, mid_e, tmp_s, tmp_e;
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
int err_sg_bin_sch_node(const char *func, const int line, SG *sg, exon_t e, int *hit)
{
    int id = sg_bin_sch_node(sg, e, hit);
    char head[100]; sprintf(head, "%s:%d", func, line);
    if (*hit == 0) _err_fatal_simple(head, "Can not hit node.\n");
    return id;
}

int sg_update_node(SG *sg, exon_t e, int32_t start, int32_t end)
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

int sg_bin_sch_site(SGsite *site, int site_n, int32_t s, int *hit)
{
    *hit = 0;
    int32_t mid_s, tmp_s;
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
int err_sg_bin_sch_site(const char *func, const int line, SGsite *site, int site_n, int32_t s, int *hit)
{
    int id = sg_bin_sch_site(site, site_n, s, hit);
    char head[100]; sprintf(head, "%s:%d", func, line);
    if (*hit == 0) _err_fatal_simple(head, "Can not hit site.\n");
    return id;
}

int sg_update_site(SG *sg, int32_t site, uint8_t type)
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

int sg_bin_sch_edge(SG *sg, uint32_t don_site_id, uint32_t acc_site_id, int *hit)
{
    *hit = 0;
    uint32_t mid_d, mid_a, tmp_d, tmp_a;
    int left = 0, right = sg->edge_n-1, mid;
    if (right == -1) return 0;

    while (left <= right) {
        mid = ((left + right) >> 1);
        mid_d = sg->edge[mid].don_site_id, mid_a = sg->edge[mid].acc_site_id;
        if (mid_d == don_site_id && mid_a == acc_site_id) { *hit = 1; return mid; }
        else if (mid_d > don_site_id || (mid_d == don_site_id && mid_a > acc_site_id)) { // [mid] is bigger than query
            if (mid != 0) {
                tmp_d = sg->edge[mid-1].don_site_id, tmp_a = sg->edge[mid-1].acc_site_id;
            }
            if (mid == 0 || (don_site_id > tmp_d || (don_site_id == tmp_d && acc_site_id > tmp_a))) {
                return mid;
            } else right = mid - 1;
        } else left = mid + 1;
    }
    return sg->edge_n;
}
int err_sg_bin_sch_edge(const char *func, const int line, SG *sg, uint32_t don_site_id, uint32_t acc_site_id, int *hit)
{
    int id = sg_bin_sch_edge(sg, don_site_id, acc_site_id, hit);
    char head[100]; sprintf(head, "%s:%d", func, line);
    if (*hit == 0) _err_fatal_simple(head, "Can not hit edge.\n");
    return id;
}
// update edge of splicing-graph
int sg_update_edge(SG *sg, uint32_t don_id, uint32_t acc_id, uint32_t don_site_id, uint32_t acc_site_id, uint8_t is_rev)
{
    int hit = 0;
    int e_i = sg_bin_sch_edge(sg, don_site_id, acc_site_id, &hit);
    if (hit == 0) { // insert new edge
        if (sg->edge_n++ >= sg->edge_m) _realloc(sg->edge, sg->edge_m, SGedge)
        // copy edge
        if (e_i <= sg->edge_n-2)
            memmove(sg->edge+e_i+1, sg->edge+e_i, (sg->edge_n-e_i-1) * sizeof(SGedge));
        // set edge
        sg->edge[e_i].don_site_id = don_site_id, sg->edge[e_i].acc_site_id = acc_site_id;
        sg->edge[e_i].is_rev = is_rev;
    }
    // set next/pre
    _insert(acc_id, sg->node[don_id].next_id, sg->node[don_id].next_n, sg->node[don_id].next_m, uint32_t)
    _insert(don_id, sg->node[acc_id].pre_id, sg->node[acc_id].pre_n, sg->node[acc_id].pre_m, uint32_t)
    // set site
    _insert(don_id, sg->don_site[don_site_id].exon_id, sg->don_site[don_site_id].exon_n, sg->don_site[don_site_id].exon_m, uint32_t)
    _insert(acc_id, sg->acc_site[acc_site_id].exon_id, sg->acc_site[acc_site_id].exon_n, sg->acc_site[acc_site_id].exon_m, uint32_t)
    return 0;
}

// order:
// 1: pre
// 2: post
void intersect_domn(uint32_t **com, uint32_t *new_domn, int *com_n, int new_n, int order)
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

void cal_pre_domn(SG *sg)
{
    int i, j;
    for (i = 0; i < sg->node_n; ++i) {
        if (sg->node[i].pre_n == 0) continue;

        int com_n = sg->node[sg->node[i].pre_id[0]].pre_domn_n, com_m = sg->node[sg->node[i].pre_id[0]].pre_domn_m;
        uint32_t *com = (uint32_t*)_err_malloc(com_m * sizeof(uint32_t));
        for (j = 0; j < com_n; ++j) com[j] = sg->node[sg->node[i].pre_id[0]].pre_domn[j];
        for (j = 1; j < sg->node[i].pre_n; ++j) {
            intersect_domn(&com, sg->node[sg->node[i].pre_id[j]].pre_domn, &com_n, sg->node[sg->node[i].pre_id[j]].pre_domn_n, 1);
        }
        if (com_n+1 > sg->node[i].pre_domn_m) {
            sg->node[i].pre_domn_m = com_n+1; 
            sg->node[i].pre_domn = (uint32_t*)_err_realloc(sg->node[i].pre_domn, (com_n+1) * sizeof(uint32_t));
        }
        for (j = 0; j < com_n; ++j) sg->node[i].pre_domn[j+1] = com[j];
        sg->node[i].pre_domn_n = 1+com_n;
        free(com);
    }
}

void cal_post_domn(SG *sg)
{
    int i, j;
    for (i = sg->node_n-1; i >= 0; --i) {
        if (sg->node[i].next_n == 0) continue;

        int com_n = sg->node[sg->node[i].next_id[0]].post_domn_n, com_m = sg->node[sg->node[i].next_id[0]].post_domn_m;
        uint32_t *com = (uint32_t*)_err_malloc(com_m * sizeof(uint32_t));
        for (j = 0; j < com_n; ++j) com[j] = sg->node[sg->node[i].next_id[0]].post_domn[j];
        for (j = 1; j < sg->node[i].next_n; ++j) 
            intersect_domn(&com, sg->node[sg->node[i].next_id[j]].post_domn, &com_n, sg->node[sg->node[i].next_id[j]].post_domn_n, 2);
        if (com_n+1 > sg->node[i].post_domn_m) {
            sg->node[i].post_domn_m = com_n+1; 
            sg->node[i].post_domn = (uint32_t*)_err_realloc(sg->node[i].post_domn, (com_n+1) * sizeof(uint32_t));
        }
        for (j = 0; j < com_n; ++j) sg->node[i].post_domn[j+1] = com[j];
        sg->node[i].post_domn_n = com_n+1;
        free(com);
    }
}

// construct splice-graph for each gene
void construct_SpliceGraph_core(SG *sg, gene_t gene)
{
    int i, j, hit; uint32_t don_id, acc_id, don_site_id, acc_site_id; exon_t e; int32_t start, end;
    
    sg->tid = gene.tid, sg->is_rev = gene.is_rev;
    // generate node
    // update v_start
    sg_update_node(sg, (exon_t){gene.tid, gene.is_rev, 0, 0}, 0, 0);
    for (i = 0; i < gene.trans_n; ++i) {
        for (j = 0; j < gene.trans[i].exon_n; ++j) {
            e = gene.trans[i].exon[j];
            start = e.start, end = e.end;
            if ((gene.is_rev == 0 && j == 0) || (gene.is_rev == 1 && j == gene.trans[i].exon_n-1)) e.start = 0;
            else sg_update_site(sg, e.start-1, ACC_SITE_F);
            if ((j == gene.trans[i].exon_n-1 && gene.is_rev == 0) || (j == 0 && gene.is_rev == 1)) e.end = MAX_SITE;
            else sg_update_site(sg, e.end+1, DON_SITE_F);
            sg_update_node(sg, e, start, end);
        }
    }
    // update v_end
    sg_update_node(sg, (exon_t){gene.tid, gene.is_rev, MAX_SITE, MAX_SITE}, MAX_SITE, MAX_SITE);

    // alloc for next_id/pre_id/pre_domn/post_domn
    sg_init_node(sg); sg_init_site(sg);

    // search node and generate edge 
    if (gene.is_rev == 0) { // forward strand
        for (i = 0; i < gene.trans_n; ++i) {
            if (gene.trans[i].exon_n == 1) continue;
            e = gene.trans[i].exon[0]; e.start = 0;

            don_id = _err_sg_bin_sch_node(sg, e, &hit);
            don_site_id = _err_sg_bin_sch_site(sg->don_site, sg->don_site_n, e.end+1, &hit);
            sg->node[don_id].s_site_id = -1;
            sg->node[don_id].e_site_id = don_site_id;

            // set next_id of v_start
            _insert(don_id, sg->node[0].next_id, sg->node[0].next_n, sg->node[0].next_m, uint32_t)
            _insert(0, sg->node[don_id].pre_id, sg->node[don_id].pre_n, sg->node[don_id].pre_m, uint32_t)

            sg->node[don_id].is_init = 1;

            for (j = 1; j < gene.trans[i].exon_n; ++j) {
                e = gene.trans[i].exon[j];
                if (j == gene.trans[i].exon_n-1) e.end = MAX_SITE;

                acc_id = _err_sg_bin_sch_node(sg, e, &hit);
                acc_site_id = _err_sg_bin_sch_site(sg->acc_site, sg->acc_site_n, e.start-1, &hit);
                sg->node[acc_id].s_site_id = acc_site_id;

                sg_update_edge(sg, don_id, acc_id, don_site_id, acc_site_id, gene.is_rev);

                if (j == gene.trans[i].exon_n-1) {
                    sg->node[acc_id].e_site_id = -1;
                    break;
                }
                don_site_id = _err_sg_bin_sch_site(sg->don_site, sg->don_site_n, e.end+1, &hit);
                sg->node[acc_id].e_site_id = don_site_id;
                don_id = acc_id;
            }
            // set pre_id of v_end
            _insert(acc_id, sg->node[sg->node_n-1].pre_id, sg->node[sg->node_n-1].pre_n, sg->node[sg->node_n-1].pre_m, uint32_t)
            _insert((uint32_t)sg->node_n-1, sg->node[acc_id].next_id, sg->node[acc_id].next_n, sg->node[acc_id].next_m, uint32_t)

            sg->node[acc_id].is_termi = 1;
        }
    } else { // reverse strand
        for (i = 0; i < gene.trans_n; ++i) {
            if (gene.trans[i].exon_n == 1) continue;
            e = gene.trans[i].exon[0]; e.end = MAX_SITE;

            acc_id = _err_sg_bin_sch_node(sg, e, &hit);
            acc_site_id = _err_sg_bin_sch_site(sg->acc_site, sg->acc_site_n, e.start-1, &hit);
            sg->node[acc_id].e_site_id = -1;
            sg->node[acc_id].s_site_id = acc_site_id;
            // set pre_id of v_end
            _insert(acc_id, sg->node[sg->node_n-1].pre_id, sg->node[sg->node_n-1].pre_n, sg->node[sg->node_n-1].pre_m, uint32_t)
            _insert((uint32_t)sg->node_n-1, sg->node[acc_id].next_id, sg->node[acc_id].next_n, sg->node[acc_id].next_m, uint32_t)

            sg->node[acc_id].is_termi = 1;

            for (j = 1; j < gene.trans[i].exon_n; ++j) {
                e = gene.trans[i].exon[j]; 
                if (j == gene.trans[i].exon_n-1) e.start = 0;

                don_id = _err_sg_bin_sch_node(sg, e, &hit);
                don_site_id = _err_sg_bin_sch_site(sg->don_site, sg->don_site_n, e.end+1, &hit);
                sg->node[don_id].e_site_id = don_site_id;

                sg_update_edge(sg, don_id, acc_id, don_site_id, acc_site_id, gene.is_rev);

                if (j == gene.trans[i].exon_n-1) { 
                    sg->node[don_id].s_site_id = -1;
                    break;
                }
                acc_site_id = _err_sg_bin_sch_site(sg->acc_site, sg->acc_site_n, e.start-1, &hit);
                sg->node[don_id].s_site_id = acc_site_id;
                acc_id = don_id;
            }
            // set next_id of v_start
            _insert(don_id, sg->node[0].next_id, sg->node[0].next_n, sg->node[0].next_m, uint32_t)
            _insert(0, sg->node[don_id].pre_id, sg->node[don_id].pre_n, sg->node[don_id].pre_m, uint32_t)

            sg->node[don_id].is_init = 1;
        }
    }
    // cal pre/post domn
    cal_pre_domn(sg); cal_post_domn(sg); 
}

SG_group *construct_SpliceGraph(FILE *gtf, chr_name_t *cname)
{
    print_format_time(stderr); err_printf("[%s] constructing splice-graph with GTF file ...\n", __func__);
    gene_group_t *gg = gene_group_init();
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

    for (i = 0; i < gg->gene_n; ++i) construct_SpliceGraph_core(sg_g->SG[i], gg->g[i]);

    gene_group_free(gg);
    print_format_time(stderr); err_printf("[%s] constructing splice-graph with GTF file done!\n", __func__);
    return sg_g;
}

/****************************************/
int trav_SpliceGraph()
{
    return 0;
}

const struct option sg_long_opt [] = {
    {"prefix", 1, NULL, 'f' },
    {0, 0, 0, 0}
};

int build_sg(int argc, char *argv[])
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
