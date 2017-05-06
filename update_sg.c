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
    print_format_time(stderr); err_printf("[%s] updating splice-graph with splice-junctions ...\n", __func__);
    int no_novel_sj = sgp->no_novel_sj; //, no_novel_com = sgp->no_novel_com; // XXX no need for no_novel_com
    int hit, i, j;
    int sj_i = 0, last_sg_i = 0, sg_i;
    while (sj_i < sj_n && last_sg_i < sg_g->SG_n) {
        if (sg_g->SG[last_sg_i]->node_n <= 3 || sg_g->SG[last_sg_i]->start == MAX_SITE || sg_g->SG[last_sg_i]->end == 0) { ++last_sg_i; continue; }
        sj_t *sj = sj_group+sj_i;
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
    for (i = 0; i < sg_g->SG_n; ++i) {
        SG *sg = sg_g->SG[i];
        cal_pre_domn(sg), cal_post_domn(sg);
    }
    print_format_time(stderr); err_printf("[%s] updating splice-graph with splice-junctions done!\n", __func__);
    return 0;
}

int update_SpliceGraph(SG_group *sg_g, sj_t *sj_group, int sj_n, sg_para *sgp)
{
    update_SpliceGraph_edge(sg_g, sj_group, sj_n, sgp);
    // remove_SpliceGraph_edge(sg_g);
    return 0;
}

