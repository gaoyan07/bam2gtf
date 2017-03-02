#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "splice_graph.h"
#include "utils.h"

extern char PROG[20];
int sg_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s sg [option] <in.gtf>\n\n", PROG);
    err_printf("Options:\n\n");
    //err_printf("         -e --exon-min    [INT]    minimum length of internal exon. [%d]\n", INTER_EXON_MIN_LEN);
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
        sg->node[i].node_id = i;
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
    for (i = 0; i < sg->site_n; ++i) {
        sg->site[i].site_id = i;
        sg->site[i].exon_n = 0; sg->site[i].exon_m = 1;
        sg->site[i].exon_id = (uint32_t*)_err_malloc(sizeof(uint32_t));
    }
    return sg;
}

SG *sg_init(void)
{
    SG *sg = (SG*)_err_malloc(sizeof(SG));
    sg->v.next_n = 0; sg->v.next_m = 1;
    sg->v.next_id = (uint32_t*)_err_malloc(sizeof(uint32_t));
    sg->v.pre_n = 0; sg->v.pre_m = 1;
    sg->v.pre_id = (uint32_t*)_err_malloc(sizeof(uint32_t));
    sg->node_n = 0, sg->node_m = 2;
    sg->node = (SGnode*)_err_malloc(2 * sizeof(SGnode));
    sg->site_n = 0, sg->site_m = 2;
    sg->site = (SGsite*)_err_malloc(2 * sizeof(SGsite));
    sg->edge_n = 0, sg->edge_m = 2;
    sg->edge = (SGedge*)_err_malloc(2 * sizeof(SGedge));
    // path_map will be alloced when SG is done
    return sg;
}

SG_group *sg_init_group(int g_n)
{
    SG_group *sg_g = (SG_group*)_err_malloc(sizeof(SG_group));
    sg_g->SG_n = g_n, sg_g->SG_m = g_n;
    sg_g->SG = (SG**)_err_malloc(g_n * sizeof(SG*));
    int i; for (i = 0; i < g_n; ++i) sg_g->SG[i] = sg_init();
    return sg_g;
}

SGasm *sg_init_asm(uint32_t sg_id, uint32_t v_start, uint32_t v_end)
{
    SGasm *sg_asm = (SGasm*)_err_malloc(sizeof(SGasm));
    sg_asm->SG_id = sg_id;
    sg_asm->v_start = v_start; sg_asm->v_end = v_end;
    sg_asm->node_n = 0; sg_asm->node_m = 1;
    sg_asm->node_id = (uint32_t*)_err_malloc(sizeof(uint32_t));
    sg_asm->edge_n = 0; sg_asm->edge_m = 1;
    sg_asm->edge_id = (uint32_t*)_err_malloc(sizeof(uint32_t));
    return sg_asm;
}

SGasm_group *sg_init_asm_group(void)
{
    SGasm_group *asm_g = (SGasm_group*)_err_malloc(sizeof(SGasm_group));
    asm_g->sg_asm_n = 0; asm_g->sg_asm_m = 1;
    asm_g->sg_asm = (SGasm**)_err_malloc(sizeof(SGasm*));
    asm_g->sg_asm[0] = sg_init_asm(0, 0, 0);
    return asm_g;
}

SG_group *sg_realloc_group(SG_group *sg_g)
{
    sg_g->SG_m <<= 1;
    sg_g->SG = (SG**)_err_realloc(sg_g->SG, sg_g->SG_m * sizeof(SG*));
    int i; for (i = (sg_g->SG_m >> 1); i < sg_g->SG_m; ++i) sg_g->SG[i] = sg_init();
    return sg_g;
}

SGasm_group *sg_realloc_asm_group(SGasm_group *asm_g)
{
    asm_g->sg_asm_m <<= 1;
    asm_g->sg_asm = (SGasm**)_err_realloc(asm_g->sg_asm, asm_g->sg_asm_m * sizeof(SGasm*));
    int i;
    for (i = asm_g->sg_asm_m >> 1; i <= asm_g->sg_asm_m; ++i) asm_g->sg_asm[i] = sg_init_asm(0, 0, 0);
    return asm_g;
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
    int i; for (i = 0; i < sg->site_m; ++i) free(sg->site[i].exon_id);
}

void sg_free(SG *sg)
{
    free(sg->v.next_id); free(sg->v.pre_id);
    sg_free_node(sg); sg_free_site(sg); free(sg->edge);
    if (sg->node_n > 0) {
        int i; for (i = 0; i < sg->node_n; ++i) free(sg->path_map[i]);
        free(sg->path_map);
    }
    free(sg);
}

void sg_free_group(SG_group *sg_g)
{
    int i; for (i = 0; i < sg_g->SG_m; i++) sg_free(sg_g->SG[i]);
    free(sg_g->SG); free(sg_g);
}

void sg_free_asm(SGasm *sg_asm)
{
    free(sg_asm->node_id); free(sg_asm->edge_id);
    free(sg_asm);
}

void sg_free_asm_group(SGasm_group *asm_g)
{
    int i; for (i = 0; i < asm_g->sg_asm_m; ++i) sg_free_asm(asm_g->sg_asm[i]);
    free(asm_g->sg_asm);
    free(asm_g);
}
/***************************/

/****************************************
 * construct splice graph from GTF file *
 ****************************************/
// binary search node
// compare start, then end
int sg_bin_sch_node(SG sg, exon_t e, int *hit)
{
    *hit = 0;
    int32_t start = e.start, end = e.end, mid_s, mid_e, tmp_s, tmp_e;
    int left = 0, right = sg.node_n-1, mid;
    if (right == -1) return 0;

    while (left <= right) {
        mid = ((left + right) >> 1);
        mid_s = sg.node[mid].e.start, mid_e = sg.node[mid].e.end;
        if (mid_s == start && mid_e == end) { *hit = 1; return mid; }
        else if (mid_s > start || (mid_s == start && mid_e > end)) { // [mid] is bigger than query
            if (mid != 0) {
                tmp_s = sg.node[mid-1].e.start, tmp_e = sg.node[mid-1].e.end;
            }
            if (mid == 0 || (start > tmp_s || (start == tmp_s && end > tmp_e))) {
                return mid;
            } else right = mid-1;
        } else left = mid + 1;
    }
    return sg.node_n;
}

int sg_update_node(SG *sg, exon_t e)
{
    int hit = 0;
    int n_i = sg_bin_sch_node(*sg, e, &hit);
    if (hit == 0) { // insert new node
        if (sg->node_n++ >= sg->node_m) _realloc(sg->node, sg->node_m, SGnode)
        // copy node
        if (n_i <= sg->node_n-2)
            memmove(sg->node+n_i+1, sg->node+n_i, (sg->node_n-n_i-1) * sizeof(SGnode));
        // set node
        sg->node[n_i].e = e;
    }
    return 0;
}

int sg_bin_sch_site(SG sg, int32_t s, int *hit)
{
    *hit = 0;
    int32_t mid_s, tmp_s;
    int left = 0, right = sg.site_n-1, mid;
    if (right == -1) return 0;

    while (left <= right) {
        mid = ((left + right) >> 1);
        mid_s = sg.site[mid].site;
        if (mid_s == s) { *hit = 1; return mid; }
        else if (mid_s > s) { // [mid] is bigger than query
            if (mid != 0) tmp_s = sg.site[mid-1].site;
            if (mid == 0 || s > tmp_s ) return mid;
            else right = mid-1;
        } else left = mid + 1;
    }
    return sg.site_n;
}

int sg_update_site(SG *sg, int32_t site, uint8_t type)
{
    int hit = 0;
    int s_i = sg_bin_sch_site(*sg, site, &hit);
    if (hit == 0) {
        if (sg->site_n++ >= sg->site_m) _realloc(sg->site, sg->site_m, SGsite)
        // copy site
        if (s_i <= sg->site_n-2)
            memmove(sg->site+s_i+1, sg->site+s_i, (sg->site_n-s_i-1) * sizeof(SGsite));
        // set site
        sg->site[s_i].site = site;
        sg->site[s_i].type = type;
    }
    return 0;
}

int sg_bin_sch_edge(SG sg, uint32_t don_site_id, uint32_t acc_site_id, int *hit)
{
    *hit = 0;
    uint32_t mid_d, mid_a, tmp_d, tmp_a;
    int left = 0, right = sg.edge_n-1, mid;
    if (right == -1) return 0;

    while (left <= right) {
        mid = ((left + right) >> 1);
        mid_d = sg.edge[mid].don_site_id, mid_a = sg.edge[mid].acc_site_id;
        if (mid_d == don_site_id && mid_a == acc_site_id) { *hit = 1; return mid; }
        else if (mid_d > don_site_id || (mid_d == don_site_id && mid_a > acc_site_id)) { // [mid] is bigger than query
            if (mid != 0) {
                tmp_d = sg.edge[mid-1].don_site_id, tmp_a = sg.edge[mid-1].acc_site_id;
            }
            if (mid == 0 || (don_site_id > tmp_d || (don_site_id == tmp_d && acc_site_id > tmp_a))) {
                return mid;
            } else right = mid - 1;
        } else left = mid + 1;
    }
    return sg.edge_n;
}

// update edge of splicing-graph
int sg_update_edge(SG *sg, uint32_t don_id, uint32_t acc_id, uint32_t don_site_id, uint32_t acc_site_id, uint8_t is_rev)
{
    int hit = 0;
    int e_i = sg_bin_sch_edge(*sg, don_site_id, acc_site_id, &hit);
    if (hit == 0) { // insert new edge
        if (sg->edge_n++ >= sg->edge_m) _realloc(sg->edge, sg->edge_m, SGedge)
        // copy edge
        if (e_i <= sg->edge_n-2)
            memmove(sg->edge+e_i+1, sg->edge+e_i, (sg->edge_n-e_i-1) * sizeof(SGedge));
        // set edge
        sg->edge[e_i].don_site_id = don_site_id, sg->edge[e_i].acc_site_id = acc_site_id;
        sg->edge[e_i].is_rev = is_rev;
        // set child
        if (sg->node[don_id].next_n >= sg->node[don_id].next_m) _realloc(sg->node[don_id].next_id, sg->node[don_id].next_m, uint32_t)
        sg->node[don_id].next_id[sg->node[don_id].next_n++] = acc_id;
        if (sg->node[acc_id].pre_n >= sg->node[acc_id].pre_m) _realloc(sg->node[acc_id].pre_id, sg->node[acc_id].pre_m, uint32_t)
        sg->node[acc_id].pre_id[sg->node[acc_id].pre_n++] = don_id;
        // set site
        if (sg->site[don_site_id].exon_n >= sg->site[don_site_id].exon_m) _realloc(sg->site[don_site_id].exon_id, sg->site[don_site_id].exon_m, uint32_t)
        sg->site[don_site_id].exon_id[sg->site[don_site_id].exon_n++] = don_id;
        if (sg->site[acc_site_id].exon_n >= sg->site[acc_site_id].exon_m) _realloc(sg->site[acc_site_id].exon_id, sg->site[acc_site_id].exon_m, uint32_t)
        sg->site[acc_site_id].exon_id[sg->site[acc_site_id].exon_n++] = don_id;
    }
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

void sg_node_add_next(SGnode *node, uint32_t id)
{
    int i;
    for (i = 0; i < node->next_n; ++i) {
        if (node->next_id[i] == id) return;
    }
    if (node->next_n == node->next_m) _realloc(node->next_id, node->next_m, uint32_t)
    node->next_id[node->next_n++] = id;
}

void sg_node_add_pre(SGnode *node, uint32_t id)
{
    int i;
    for (i = 0; i < node->pre_n; ++i) {
        if (node->pre_id[i] == id) return;
    }
    if (node->pre_n == node->pre_m) _realloc(node->pre_id, node->pre_m, uint32_t)
    node->pre_id[node->pre_n++] = id;
}

// construct splice-graph for each gene
void construct_SpliceGraph_core(SG *sg, gene_t gene)
{
    int i, j, hit; uint32_t don_id, acc_id, don_site_id, acc_site_id; exon_t e;
    
    // generate node
    for (i = 0; i < gene.trans_n; ++i) {
        for (j = 0; j < gene.trans[i].exon_n; ++j) {
            e = gene.trans[i].exon[j];
            if ((j == 0 && gene.is_rev == 0) || (j == gene.trans[i].exon_n-1 && gene.is_rev == 1)) e.start = 0;
            if ((j == gene.trans[i].exon_n-1 && gene.is_rev == 0) || (j == 0 && gene.is_rev == 1)) e.end = MAX_END;
            sg_update_node(sg, e);
            if (e.start != 0) sg_update_site(sg, e.start-1, ACC_SITE_F);
            if (e.end != MAX_END) sg_update_site(sg, e.end+1, DON_SITE_F);
        }
    }
    // alloc for next_id/pre_id/pre_domn/post_domn
    sg_init_node(sg); sg_init_site(sg);
    // alloc path_map
    sg->path_map = (uint8_t**)_err_malloc(sg->node_n * sizeof(uint8_t*));
    for (i = 0; i < sg->node_n; ++i) sg->path_map[i] = (uint8_t*)_err_calloc(sg->node_n, sizeof(uint8_t));

    // search node and generate edge 
    for (i = 0; i < gene.trans_n; ++i) {
        e = gene.trans[i].exon[0]; 
        if (gene.is_rev == 0 || (gene.trans[i].exon_n-1==0 && gene.is_rev==1)) e.start = 0;
        if ((gene.trans[i].exon_n-1==0 && gene.is_rev == 0) || (gene.is_rev==1)) e.end = MAX_END;

        don_id = sg_bin_sch_node(*sg, e, &hit);
        if (hit == 0) err_fatal_simple("Can not hit node.(1)\n");
        don_site_id = sg_bin_sch_site(*sg, e.end+1, &hit);
        if (hit == 0) err_fatal_simple("Can not hit edge.(1)\n");

        // set next_id of s
        sg_node_add_next(&(sg->v), don_id);
        
        for (j = 1; j < gene.trans[i].exon_n; ++j) {
            e = gene.trans[i].exon[j]; 
            if ((j == 0 && gene.is_rev == 0) || (j == gene.trans[i].exon_n-1 && gene.is_rev==1)) e.start = 0;
            if ((j == gene.trans[i].exon_n-1 && gene.is_rev == 0) || (j == 0 && gene.is_rev==1)) e.end = MAX_END;
 
            acc_id = sg_bin_sch_node(*sg, e, &hit);
            if (hit == 0) err_fatal_simple("Can not hit node.(2)\n");
            acc_site_id = sg_bin_sch_site(*sg, e.start-1, &hit);
            if (hit == 0) err_fatal_simple("Can not hit site.(2)\n");
            sg_update_edge(sg, don_id, acc_id, don_site_id, acc_site_id, gene.is_rev);
            // path map
            sg->path_map[don_id][acc_id] = 1;

            don_id = acc_id;
            don_site_id = sg_bin_sch_site(*sg, e.end+1, &hit);
            // XXX if (hit == 0) err_fatal_simple("Can not hit site.(3)\n");
        }
        // set pre_id of e
        sg_node_add_pre(&(sg->v), don_id);
    }
    // update path_map
    int k;
    for (i = 0; i < sg->node_n; ++i) {
        for (j = i+2; j < sg->node_n; ++j) {
            if (sg->path_map[i][j] == 1) continue;
            for (k = i+1; k < j; ++k) {
                if (sg->path_map[i][k] > 0 && sg->path_map[k][j] > 0) 
                    sg->path_map[i][j] = 2;
            }
        }
    }
    // cal pre/post domn
    cal_pre_domn(sg); cal_post_domn(sg); 
}

SG_group *construct_SpliceGraph(FILE *gtf, chr_name_t *cname)
{
    gene_group_t *gg = gene_group_init();
    int g_n = read_gene_group(gtf, cname, gg);
    SG_group *sg_g = sg_init_group(g_n);

    int i; for (i = 0; i < gg->gene_n; ++i) construct_SpliceGraph_core(sg_g->SG[i], gg->g[i]);

    gene_group_free(gg);
    return sg_g;
}
/****************************************/

/***************************************************************
 * update splice graph with splice-junctions (short-read data) *
 ***************************************************************/
int predict_SpliceGraph_core(SG_group sg_g, sj_t *sj_group, int sj_n, SG_group *sr_sg_g, int no_novel_sj)
{
    int i, j, hit; uint32_t don_site_id, acc_site_id, edge_id; int don_sg_i=0, acc_sg_i=0, sg_i;
    for (i = 0; i < sj_n; ++i) {
        // for sorted sj_group
        while (don_sg_i < sg_g.SG_n && sj_group[i].don > sg_g.SG[don_sg_i]->end) don_sg_i++;
        while (acc_sg_i < sg_g.SG_n && sj_group[i].acc > sg_g.SG[don_sg_i]->end) don_sg_i++;

        if (don_sg_i == sg_g.SG_n || acc_sg_i == sg_g.SG_n) break;
        if (don_sg_i != acc_sg_i || sj_group[i].don < sg_g.SG[don_sg_i]->start || sj_group[i].acc < sg_g.SG[acc_sg_i]->start) continue; // splice-site span gene-loci OR novel splice-site

        sg_i = don_sg_i;
        // 1. (don_site, acc_site) => (don_site_id, acc_site_id)
        don_site_id = sg_bin_sch_site(*(sg_g.SG[sg_i]), sj_group[i].don, &hit); if (hit == 0) continue;
        acc_site_id = sg_bin_sch_site((*sg_g.SG[sg_i]), sj_group[i].acc, &hit); if (hit == 0) continue;
        // 2. search_edge(don_site_id, acc_site_id) => edge_id
        edge_id = sg_bin_sch_edge(*(sg_g.SG[sg_i]), don_site_id, acc_site_id, &hit); if (hit == 0 && no_novel_sj == 1) continue;
        
        // 3. update_site(don_site, acc_site)
        sg_update_site(sr_sg_g->SG[sg_i], sj_group[i].don, DON_SITE_F);
        sg_update_site(sr_sg_g->SG[sg_i], sj_group[i].acc, ACC_SITE_F);
        // 4. update_node(edge[edge_id].don_exon[], edge[edge_id].acc_exon[])
        for (j = 0; j < sg_g.SG[sg_i]->site[don_site_id].exon_n; ++j) {
            int e_i = sg_g.SG[sg_i]->site[don_site_id].exon_id[j];
            exon_t e = sg_g.SG[sg_i]->node[e_i].e;
            sg_update_node(sr_sg_g->SG[sg_i], e);
        }
        for (j = 0; j < sg_g.SG[sg_i]->site[acc_site_id].exon_n; ++j) {
            int e_i = sg_g.SG[sg_i]->site[acc_site_id].exon_id[j];
            exon_t e = sg_g.SG[sg_i]->node[e_i].e;
            sg_update_node(sr_sg_g->SG[sg_i], e);
        }
    }
    // alloc mem for node&site
    for (i = 0; i < sg_g.SG_n; ++i) {
        sg_init_node(sr_sg_g->SG[i]); sg_init_site(sr_sg_g->SG[i]);
    }
    // 5. updage_edge(edge[edge_id])
    return 0;
}

// sr_sg_g.SG_n == sg_g.SG_n
int predict_SpliceGraph(SG_group sg_g, FILE *sj_p, chr_name_t *cname, SG_group *sr_sg_g, int no_novel_sj)
{
    sj_t *sj_group = (sj_t*)_err_malloc(10000 * sizeof(sj_t));
    int sj_n, sj_m = 10000;
    sj_n = read_sj_group(sj_p, cname, &sj_group, sj_m);
    predict_SpliceGraph_core(sg_g, sj_group, sj_n, sr_sg_g, no_novel_sj);
    free(sj_group);
    return 0;
}
/***************************************************************/


/*****************************
 *       generate ASM        *
 *****************************/
// sort edge by interval XXX
/*int sg_max_edge(SG *sg)
{
    int i, j;
    for (i = 1; i < sg->edge_n; ++i) {
        for (j = 0; j < i; ++j) {
            if ([j] > [i])
                [i] is NOT max
        }
    }
    return 0;
}*/

void cal_cand_node(SG sg, uint32_t **entry, uint32_t **exit, int *entry_n, int *exit_n)
{
    int i, n1=0, n2=0;
    for (i = 0; i < sg.node_n; ++i) {
        if (sg.node[i].next_n > 1) n1++;
        if (sg.node[i].pre_n > 1) n2++;
    }
    *entry_n = n1, *exit_n = n2;
    *entry = (uint32_t*)_err_malloc(n1 * sizeof(uint32_t));
    *exit = (uint32_t*)_err_malloc(n2 * sizeof(uint32_t));

    n1 = 0, n2 = 0;
    for (i = 0; i < sg.node_n; ++i) {
        if (sg.node[i].next_n > 1) (*entry)[n1++] = sg.node[i].node_id;
        if (sg.node[i].pre_n > 1) (*exit)[n2++] = sg.node[i].node_id;;
    }
}

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
            m <<= 1;                \
            p = (type*)_err_realloc(p, m * sizeof(type));   \
        }                           \
        p[n++] = v;                 \
    }                               \
}

void sg_update_asm(SG sg, SGasm *sg_asm, uint32_t pre_id, uint32_t cur_id)
{
    int hit = 0;
    uint32_t edge_i = sg_bin_sch_edge(sg, pre_id, cur_id, &hit);
    if (hit == 0) err_fatal_simple("Can not hit edge.(3)\n");
    _insert(cur_id, sg_asm->node_id, sg_asm->node_n, sg_asm->node_m, uint32_t)
    _insert(edge_i, sg_asm->edge_id, sg_asm->edge_n, sg_asm->edge_m, uint32_t)
}

void sg_update_asm_edge(SG sg, SGasm *sg_asm, uint32_t pre_id, uint32_t cur_id)
{
    int hit = 0;
    uint32_t edge_i = sg_bin_sch_edge(sg, pre_id, cur_id, &hit);
    if (hit == 0) err_fatal_simple("Can not hit edge.(3)\n");
    _insert(edge_i, sg_asm->edge_id, sg_asm->edge_n, sg_asm->edge_m, uint32_t)
}

void sub_splice_graph(SG sg, SGasm *sg_asm, uint32_t cur_id, uint32_t e_id)
{
    //if (sg.node[cur_id].next_n == 0) err_fatal_simple("No next node!\n");

    if (cur_id == e_id) return;
    int i;
    for (i = 0; i < sg.node[cur_id].next_n; ++i) {
        if (sg.node[cur_id].next_id[i] == e_id) {
            sg_update_asm_edge(sg, sg_asm, cur_id, e_id);
            continue;
        } else {
            sg_update_asm(sg, sg_asm, cur_id, sg.node[cur_id].next_id[i]); 
            sub_splice_graph(sg, sg_asm, sg.node[cur_id].next_id[i], e_id);
        }
    }
}

int sg_asm_group_add(SGasm_group *asm_g, SGasm *sg_asm)
{
    if (asm_g->sg_asm_n == asm_g->sg_asm_m) sg_realloc_asm_group(asm_g);
    SGasm *a = asm_g->sg_asm[asm_g->sg_asm_n];

    a->SG_id = sg_asm->SG_id;
    a->v_start = sg_asm->v_start; a->v_end = sg_asm->v_end;
    a->node_n = sg_asm->node_n; a->edge_n = sg_asm->edge_n;
    a->node_id = (uint32_t*)_err_realloc(a->node_id, sg_asm->node_n * sizeof(uint32_t));
    a->edge_id = (uint32_t*)_err_realloc(a->edge_id, sg_asm->edge_n * sizeof(uint32_t));
    int i;
    for (i = 0; i < sg_asm->node_n; ++i) a->node_id[i] = sg_asm->node_id[i];
    for (i = 0; i < sg_asm->edge_n; ++i) a->edge_id[i] = sg_asm->edge_id[i];

    asm_g->sg_asm_n++;
    return asm_g->sg_asm_n;
}

void gen_ASM(SG sg, uint32_t sg_id, SGasm_group *asm_g)
{
    int entry_n, exit_n;
    uint32_t *entry, *exit;
    cal_cand_node(sg, &entry, &exit, &entry_n, &exit_n);
    if (entry_n == 0 || exit_n == 0) return;

    int i, j, hit;
    for (i = 0; i < entry_n; ++i) {
        hit = 0;
        for (j = 0; j < exit_n; ++j) {
            int post_domn_n = sg.node[entry[i]].post_domn_n;
            int pre_domn_n = sg.node[exit[j]].pre_domn_n;
            if (post_domn_n > 1 && pre_domn_n > 1 
            && sg.node[entry[i]].post_domn[1] == exit[j] 
            && sg.node[exit[j]].pre_domn[1] == entry[i]) {
                SGasm *sg_asm = sg_init_asm(sg_id, entry[i], exit[j]);
                sub_splice_graph(sg, sg_asm, entry[i], exit[j]);
                sg_asm_group_add(asm_g, sg_asm);
                //asm_sg = rm_max_edge(asm_sg);
                //get_ASM(asm_sg, asm_g);
                sg_free_asm(sg_asm);
                hit = 1; break;
            }
        }
        if (hit) continue;
    }
    free(entry); free(exit);
}
/*****************************/

int trav_SpliceGraph()
{
    return 0;
}

const struct option sg_long_opt [] = {
    {0, 0, 0, 0}
};

int sg(int argc, char *argv[])
{
    int c;
    while ((c = getopt_long(argc, argv, "", sg_long_opt, NULL)) >= 0) {
        switch (c) {
            default: err_printf("Error: unknown option: %s.\n", optarg);
                     return sg_usage();
        }
    }
    if (argc - optind != 1) return sg_usage();
    FILE *gtf = fopen(argv[optind], "r");
    chr_name_t *cname = chr_name_init();

    SG_group *sg_g = construct_SpliceGraph(gtf, cname);

    int i, j;
    for (i = 0; i < sg_g->SG[0]->v.next_n; ++i) printf("%d\t", sg_g->SG[0]->v.next_id[i]); printf("\n");
    for (i = 0; i < sg_g->SG[0]->v.pre_n; ++i) printf("%d\t", sg_g->SG[0]->v.pre_id[i]); printf("\n");
    for (i = 0; i < sg_g->SG[0]->node_n; ++i) printf("\t%d", i+1); printf("\n");

    for (i = 0; i < sg_g->SG[0]->node_n; ++i) {
        printf("%d\t", i+1);
        for (j = 0; j < i; ++j) {
            printf("%d\t", sg_g->SG[0]->path_map[j][i]);
        }
        printf("\n");
    }

    SGasm_group *asm_g = sg_init_asm_group();
    for (i = 0; i < sg_g->SG_n; ++i) gen_ASM(*(sg_g->SG[i]), i, asm_g);
    printf("ASM: %d\n", asm_g->sg_asm_n);

    sg_free_group(sg_g); sg_free_asm_group(asm_g);
    fclose(gtf); chr_name_free(cname);
    return 0;
}
