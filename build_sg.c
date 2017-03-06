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
    err_printf("Usage:   %s sg [option] <in.gtf> <in.sj>\n\n", PROG);
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

    sg->start = CHR_MAX_END, sg->end = 0;
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
    int i; for (i = 0; i < sg->site_n; ++i) free(sg->site[i].exon_id);
    free(sg->site);
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
        if (sg->is_rev == 0) {
            if (type == DON_SITE_F && site < sg->start) sg->start = site;
            else if (type == ACC_SITE_F && site > sg->end) sg->end = site;
        } else {
            if (type == DON_SITE_F && site > sg->end) sg->end = site;
            else if (type == ACC_SITE_F && site < sg->start) sg->start = site;
        }

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
        // set next/pre
        if (sg->node[don_id].next_n >= sg->node[don_id].next_m) _realloc(sg->node[don_id].next_id, sg->node[don_id].next_m, uint32_t)
        sg->node[don_id].next_id[sg->node[don_id].next_n++] = acc_id;
        if (sg->node[acc_id].pre_n >= sg->node[acc_id].pre_m) _realloc(sg->node[acc_id].pre_id, sg->node[acc_id].pre_m, uint32_t)
        sg->node[acc_id].pre_id[sg->node[acc_id].pre_n++] = don_id;
        // set site
        _insert(don_id, sg->site[don_site_id].exon_id, sg->site[don_site_id].exon_n, sg->site[don_site_id].exon_m, uint32_t)
        _insert(acc_id, sg->site[acc_site_id].exon_id, sg->site[acc_site_id].exon_n, sg->site[acc_site_id].exon_m, uint32_t)
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

// construct splice-graph for each gene
void construct_SpliceGraph_core(SG *sg, gene_t gene)
{
    int i, j, hit; uint32_t don_id, acc_id, don_site_id, acc_site_id; exon_t e;
    
    sg->tid = gene.tid, sg->is_rev = gene.is_rev;
    // generate node
    for (i = 0; i < gene.trans_n; ++i) {
        for (j = 0; j < gene.trans[i].exon_n; ++j) {
            e = gene.trans[i].exon[j];
            if ((j == 0 && gene.is_rev == 0) || (j == gene.trans[i].exon_n-1 && gene.is_rev == 1)) e.start = 0;
            if ((j == gene.trans[i].exon_n-1 && gene.is_rev == 0) || (j == 0 && gene.is_rev == 1)) e.end = CHR_MAX_END;
            sg_update_node(sg, e);
            if (e.start != 0) sg_update_site(sg, e.start-1, ACC_SITE_F);
            if (e.end != CHR_MAX_END) sg_update_site(sg, e.end+1, DON_SITE_F);
        }
    }

    // alloc for next_id/pre_id/pre_domn/post_domn
    sg_init_node(sg); sg_init_site(sg);
    // XXX alloc path_map
    sg->path_map = (uint8_t**)_err_malloc(sg->node_n * sizeof(uint8_t*));
    for (i = 0; i < sg->node_n; ++i) sg->path_map[i] = (uint8_t*)_err_calloc(sg->node_n, sizeof(uint8_t));

    // search node and generate edge 
    if (gene.is_rev == 0) {
        for (i = 0; i < gene.trans_n; ++i) {
            if (gene.trans[i].exon_n == 1) continue;
            e = gene.trans[i].exon[0];
            e.start = 0;

            don_id = sg_bin_sch_node(*sg, e, &hit); if (hit == 0) err_fatal_simple("Can not hit node.(1)\n");
            don_site_id = sg_bin_sch_site(*sg, e.end+1, &hit); if (hit == 0) err_fatal_simple("Can not hit site.(1)\n");

            // set next_id of s
            _insert(don_id, sg->v.next_id, sg->v.next_n, sg->v.next_m, uint32_t)

            for (j = 1; j < gene.trans[i].exon_n; ++j) {
                e = gene.trans[i].exon[j];
                if (j == gene.trans[i].exon_n-1) e.end = CHR_MAX_END;

                acc_id = sg_bin_sch_node(*sg, e, &hit); if (hit == 0) err_fatal_simple("Can not hit node.(2)\n");
                acc_site_id = sg_bin_sch_site(*sg, e.start-1, &hit); if (hit == 0) err_fatal_simple("Can not hit site.(2)\n");

                sg_update_edge(sg, don_id, acc_id, don_site_id, acc_site_id, gene.is_rev);
                // XXX path map
                sg->path_map[don_id][acc_id] = 1;

                don_id = acc_id;
                don_site_id = sg_bin_sch_site(*sg, e.end+1, &hit);
                // XXX if (hit == 0) err_fatal_simple("Can not hit site.(3)\n");
            }
            // set pre_id of e
            _insert(don_id, sg->v.pre_id, sg->v.pre_n, sg->v.pre_m, uint32_t)
        }
    } else {
        for (i = 0; i < gene.trans_n; ++i) {
            if (gene.trans[i].exon_n == 1) continue;
            e = gene.trans[i].exon[0]; 
            e.end = CHR_MAX_END;

            acc_id = sg_bin_sch_node(*sg, e, &hit); if (hit == 0) err_fatal_simple("Can not hit node.(3)\n");
            acc_site_id = sg_bin_sch_site(*sg, e.start-1, &hit); if (hit == 0) err_fatal_simple("Can not hit site.(3)\n");

            // set next_id of s
            _insert(acc_id, sg->v.pre_id, sg->v.pre_n, sg->v.pre_m, uint32_t)

            for (j = 1; j < gene.trans[i].exon_n; ++j) {
                e = gene.trans[i].exon[j]; 
                if (j == gene.trans[i].exon_n-1) e.start = 0;

                don_id = sg_bin_sch_node(*sg, e, &hit); if (hit == 0) err_fatal_simple("Can not hit node.(4)\n");
                don_site_id = sg_bin_sch_site(*sg, e.end+1, &hit); if (hit == 0) err_fatal_simple("Can not hit site.(4)\n");

                sg_update_edge(sg, don_id, acc_id, don_site_id, acc_site_id, gene.is_rev);
                // XXX path map
                sg->path_map[don_id][acc_id] = 1;

                acc_id = don_id;
                acc_site_id = sg_bin_sch_site(*sg, e.start-1, &hit);
                // XXX if (hit == 0) err_fatal_simple("Can not hit site.(3)\n");
            }
            // set pre_id of e
            _insert(acc_id, sg->v.next_id, sg->v.next_n, sg->v.next_m, uint32_t)
        }
    }
    // XXX update path_map
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
int trav_SpliceGraph()
{
    return 0;
}

const struct option sg_long_opt [] = {
    {0, 0, 0, 0}
};

int build_sg(int argc, char *argv[])
{
    int c;
    while ((c = getopt_long(argc, argv, "", sg_long_opt, NULL)) >= 0) {
        switch (c) {
            default: err_printf("Error: unknown option: %s.\n", optarg);
                     return sg_usage();
        }
    }
    if (argc - optind != 2) return sg_usage();
    FILE *gtf_fp = xopen(argv[optind], "r");
    chr_name_t *cname = chr_name_init();

    SG_group *sg_g = construct_SpliceGraph(gtf_fp, cname);

    /*FILE *sj_fp = xopen(argv[optind+1], "r");
    SG_group *sr_sg_g = sg_init_group(1);

    int i, j;
    for (i = 0; i < sr_sg_g->SG[0]->v.next_n; ++i) printf("%d\t", sr_sg_g->SG[0]->v.next_id[i]); printf("\n");
    for (i = 0; i < sr_sg_g->SG[0]->v.pre_n; ++i) printf("%d\t", sr_sg_g->SG[0]->v.pre_id[i]); printf("\n");
    for (i = 0; i < sr_sg_g->SG[0]->node_n; ++i) printf("\t%d", i+1); printf("\n");

    for (i = 0; i < sr_sg_g->SG[0]->node_n; ++i) {
        printf("%d\t", i+1);
        for (j = 0; j < i; ++j) {
            printf("%d\t", sr_sg_g->SG[0]->path_map[j][i]);
        }
        printf("\n");
    }*/

    sg_free_group(sg_g); err_fclose(gtf_fp); chr_name_free(cname); 
    //err_fclose(sj_fp); sg_free_group(sr_sg_g);
    return 0;
}
