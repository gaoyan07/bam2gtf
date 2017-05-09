#include <stdio.h>
#include <stdlib.h>
#include "gtf.h"
#include "build_sg.h"
#include "utils.h"

#define sg_update_edge_wei(ed, sj) { \
    ed.motif = sj->motif; \
    ed.uniq_c = sj->uniq_c; \
    ed.multi_c = sj->multi_c; \
    ed.max_over = sj->max_over; \
}

// update based on splice-junction
int update_SpliceGraph_edge(SG_group *sg_g, sj_t *sj_group, int sj_n, sg_para *sgp)
{
    err_func_format_printf(__func__, "updating splice-graph with splice-junctions ...\n");
    int no_novel_sj = sgp->no_novel_sj; //, no_novel_com = sgp->no_novel_com; // XXX no need for no_novel_com
    int hit, i, j;
    int sj_i = 0, last_sg_i = 0, sg_i;
    while (sj_i < sj_n && last_sg_i < sg_g->SG_n) {
        if (sg_g->SG[last_sg_i]->node_n <= 3 || sg_g->SG[last_sg_i]->start == MAX_SITE || sg_g->SG[last_sg_i]->end == 0) { ++last_sg_i; continue; }
        sj_t *sj = sj_group+sj_i;
        if (sgp->rm_edge && sj->uniq_c < sgp->edge_wt) { ++sj_i; continue; } // remove sj/edge
        int comp_res = comp_sj_sg(*sj, *(sg_g->SG[last_sg_i]));
        if (comp_res < 0) { ++sj_i; continue; }
        else if (comp_res > 0) { ++last_sg_i; continue; }
        for (sg_i = last_sg_i; sg_i < sg_g->SG_n; ++sg_i) {
            if (comp_sj_sg(*sj, *(sg_g->SG[sg_i])) < 0) break;
            SG *sg = sg_g->SG[sg_i]; SGnode *node = sg->node;
            SGsite *don_site = sg->don_site, *acc_site = sg->acc_site; int32_t acc_n = sg->acc_site_n, don_n = sg->don_site_n;
            // 0. search site/edge: (GTF_don_site_id, GTF_acc_site_id) => GTF_edge_id
            int GTF_don_site_id = sg_bin_sch_site(don_site, don_n, sj->don, &hit); if (hit == 0) continue;
            int GTF_acc_site_id = sg_bin_sch_site(acc_site, acc_n, sj->acc, &hit); if (hit == 0) continue;
            int GTF_edge_id = sg_bin_sch_edge(sg, GTF_don_site_id, GTF_acc_site_id, &hit);
            // 1. for valid sj
            if (hit == 1) { // update edge weight
                sg_update_edge_wei(sg->edge[GTF_edge_id], sj)
            } else if (hit == 0 && no_novel_sj == 0) { // add new edge
                int is_anno = 0, is_rev = sg->is_rev;
                // add edge between all candidate nodes
                for (i = 0; i < don_site[GTF_don_site_id].exon_n; ++i) {
                    int don_id = don_site[GTF_don_site_id].exon_id[i];
                    for (j = 0; j < acc_site[GTF_acc_site_id].exon_n; ++j) {
                        int acc_id = acc_site[GTF_acc_site_id].exon_id[j];
                        _bin_insert(acc_id, node[don_id].next_id, node[don_id].next_n, node[don_id].next_m, gec_t)
                        _bin_insert(don_id, node[acc_id].pre_id, node[acc_id].pre_n, node[acc_id].pre_m, gec_t)
                    }
                }

                sg_add_edge(sg->edge, GTF_edge_id, (sg->edge_n), (sg->edge_m), GTF_don_site_id, GTF_acc_site_id, is_rev, is_anno)
                sg_update_edge_wei(sg->edge[GTF_edge_id], sj)
            }
        }
        // one sj could be used for multi gene
        sj_i++;
    }
    // cal domn after remove 0-weight edge
    if (sgp->rm_edge == 0) {
        for (i = 0; i < sg_g->SG_n; ++i) {
            SG *sg = sg_g->SG[i];
            cal_pre_domn(sg), cal_post_domn(sg);
        }
    }
    err_func_format_printf(__func__, "updating splice-graph with splice-junctions done!\n");
    return 0;
}

void remove_SpliceGraph_edge(SG_group *sg_g, int edge_wt)
{
    int sg_i, edge_i, i, j;
    for (sg_i = 0; sg_i < sg_g->SG_n; ++sg_i) {
        SG *sg = sg_g->SG[sg_i];
        SGnode *node = sg->node; SGsite *don = sg->don_site, *acc = sg->acc_site; SGedge *edge = sg->edge;
        for (edge_i = 0; edge_i < sg->edge_n; ++edge_i) {
            if (edge[edge_i].uniq_c < edge_wt) {
                int don_id = edge[edge_i].don_site_id, acc_id = edge[edge_i].acc_site_id;
                int *don_exon_id = don[don_id].exon_id, don_exon_n = don[don_id].exon_n;
                int acc_site = acc[acc_id].site;
                for (i = 0; i < don_exon_n; ++i) {
                    int cur_id = don_exon_id[i];
                    for (j = 0; j < node[cur_id].next_n; ++j) {
                        int next_id = node[cur_id].next_id[j];
                        if (node[next_id].start == acc_site + 1) { // remove edge
                            int m, n, hit=0;
                            // remove id->next_id, next_n
                            node[cur_id].next_n--;
                            for (m = j; m < node[cur_id].next_n; ++m) node[cur_id].next_id[m] = node[cur_id].next_id[m+1];
                            // remove next_id->pre_id, pre_n
                            if (node[next_id].pre_n < 1)
                                err_printf("Can not hit pre_id. (%d %d)\n", cur_id, next_id);

                            for (m = 0; m < node[next_id].pre_n; ++m) {
                                if (node[next_id].pre_id[m] == cur_id) {
                                    hit = 1;
                                    node[next_id].pre_n--;
                                    for (n = m; n < node[next_id].pre_n; ++n) {
                                        node[next_id].pre_id[n] = node[next_id].pre_id[n+1];
                                    }
                                    break;
                                }
                            }
                            if (hit == 0)
                                err_printf("Can not hit pre_id. (%d %d)\n", cur_id, next_id);
                        }
                    }
                }
            }
        }

        cal_pre_domn(sg), cal_post_domn(sg);
    }
}

int update_SpliceGraph(SG_group *sg_g, sj_t *sj_group, int sj_n, sg_para *sgp)
{
    update_SpliceGraph_edge(sg_g, sj_group, sj_n, sgp);
    if (sgp->rm_edge) remove_SpliceGraph_edge(sg_g, sgp->edge_wt);
    return 0;
}
