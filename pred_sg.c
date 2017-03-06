#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gtf.h"
#include "build_sg.h"

/****************************************************************
 * predict splice graph with splice-junctions (short-read data) *
 * based on reference-splice-graph                              *
 ****************************************************************/
int sg_update_edge_pred(SG *sg, sj_t sj, uint32_t don_site_id, uint32_t acc_site_id)
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
        if (sj.strand == 0) sg->edge[e_i].is_rev = 2;
        else sg->edge[e_i].is_rev = sj.strand - 1; // 0:+, 1:-, 2:undefined
        sg->edge[e_i].motif = sj.motif, sg->edge[e_i].is_anno = sj.is_anno;
        sg->edge[e_i].uniq_c = sj.uniq_c, sg->edge[e_i].multi_c = sj.multi_c, sg->edge[e_i].max_over = sj.max_over;
    }
    return 0;
}

int comp_sj_sg(sj_t sj, SG sg)
{
    if (sj.tid < sg.tid) return -1;
    else if (sj.tid > sg.tid) return 1;
    else { // tid ==
        if (sj.acc < sg.start) return -1;
        else if (sj.don > sg.end) return 1;
        else return 0;
    }
}

int predict_SpliceGraph_core(SG_group sg_g, sj_t *sj_group, int sj_n, SG_group *sr_sg_g, int no_novel_sj)
{
    int i, j, k, hit; uint32_t don_site_id, acc_site_id;
    uint32_t GTF_don_site_id, GTF_acc_site_id, GTF_edge_id;
    uint8_t **node_map = (uint8_t**)_err_malloc(sg_g.SG_n * sizeof(uint8_t*));
    for (i = 0; i < sg_g.SG_n; ++i) {
        sr_sg_g->SG[i]->tid = sg_g.SG[i]->tid;
        sr_sg_g->SG[i]->is_rev = sg_g.SG[i]->is_rev;
        node_map[i] = (uint8_t*)_err_calloc(sg_g.SG[i]->node_n, sizeof(uint8_t));
        for (j = 0; j < sg_g.SG[i]->node_n; ++j) {
            if (sg_g.SG[i]->node[j].e.start == 0) node_map[i][j] |= 2;
            if (sg_g.SG[i]->node[j].e.end == CHR_MAX_END) node_map[i][j] |= 1;
        }
    }
    int last_sg_i=0, sg_i;
    
    //NEW
    i = 0, last_sg_i = 0;
    while (i < sj_n && last_sg_i < sg_g.SG_n) {
        if (sg_g.SG[last_sg_i]->node_n <= 2) { last_sg_i++; continue; }
        int comp_res = comp_sj_sg(sj_group[i], *(sg_g.SG[last_sg_i]));
        if (comp_res < 0) { i++; continue; }
        else if (comp_res > 0) { last_sg_i++; continue; }
        else { 
            for (sg_i = last_sg_i; sg_i < sg_g.SG_n; ++sg_i) {
                if (comp_sj_sg(sj_group[i], *(sg_g.SG[sg_i])) < 0) break;
                SG *sg = sg_g.SG[sg_i], *sr_sg = sr_sg_g->SG[sg_i];
                // 0. search site/edge: (GTF_don_site_id, GTF_acc_site_id) => GTF_edge_id
                GTF_don_site_id = sg_bin_sch_site(*sg, sj_group[i].don, &hit); if (hit == 0) continue;
                GTF_acc_site_id = sg_bin_sch_site(*sg, sj_group[i].acc, &hit); if (hit == 0) continue;
                GTF_edge_id = sg_bin_sch_edge(*sg, GTF_don_site_id, GTF_acc_site_id, &hit); if (hit == 0 && no_novel_sj == 1) continue;
                // 1. update site(don_site, acc_site)
                sg_update_site(sr_sg, sj_group[i].don, DON_SITE_F);
                sg_update_site(sr_sg, sj_group[i].acc, ACC_SITE_F);
                // 2. update node
                //    2.1 map[exon] = 1
                for (j = 0; j < sg->site[GTF_don_site_id].exon_n; ++j) {
                    uint32_t GTF_don_id = sg->site[GTF_don_site_id].exon_id[j];
                    node_map[sg_i][GTF_don_id] |= 1;
                    if (node_map[sg_i][GTF_don_id] == 3) { // update node
                        exon_t e = sg->node[GTF_don_id].e;
                        sg_update_node(sr_sg, e);
                        node_map[sg_i][GTF_don_id] |= 4;
                    }
                }
                for (j = 0; j < sg->site[GTF_acc_site_id].exon_n; ++j) {
                    uint32_t GTF_acc_id = sg->site[GTF_acc_site_id].exon_id[j];
                    node_map[sg_i][GTF_acc_id] |= 2;
                    if (node_map[sg_i][GTF_acc_id] == 3) { // update node
                        exon_t e = sg->node[GTF_acc_id].e;
                        sg_update_node(sr_sg, e);
                        node_map[sg_i][GTF_acc_id] |= 4;
                    }
                }
                break;
            }
            i++;
        }
    }
    // NEW-END

    // OLD
    /*for (i = 0; i < sj_n; ++i) {
        // for sorted sj_group
        while (last_sg_i < sg_g.SG_n && ((sg_g.SG[last_sg_i]->tid < sj_group[i].tid) || (sg_g.SG[last_sg_i]->tid == sj_group[i].tid && sg_g.SG[last_sg_i]->end < sj_group[i].don))) {
            last_sg_i++;
        }
        if (last_sg_i >= sg_g.SG_n) 
            break;

        sg_i = last_sg_i;
        while (sg_i < sg_g.SG_n && sj_group[i].tid == sg_g.SG[sg_i]->tid && sj_group[i].don >= sg_g.SG[sg_i]->start && sj_group[i].acc <= sg_g.SG[sg_i]->end) {
            SG *sg = sg_g.SG[sg_i], *sr_sg = sr_sg_g->SG[sg_i];
            // 0. search site/edge: (GTF_don_site_id, GTF_acc_site_id) => GTF_edge_id
            GTF_don_site_id = sg_bin_sch_site(*sg, sj_group[i].don, &hit); if (hit == 0) goto NEXT1;
            GTF_acc_site_id = sg_bin_sch_site(*sg, sj_group[i].acc, &hit); if (hit == 0) goto NEXT1;
            GTF_edge_id = sg_bin_sch_edge(*sg, GTF_don_site_id, GTF_acc_site_id, &hit); if (hit == 0 && no_novel_sj == 1) goto NEXT1;
            // 1. update site(don_site, acc_site)
            sg_update_site(sr_sg, sj_group[i].don, DON_SITE_F);
            sg_update_site(sr_sg, sj_group[i].acc, ACC_SITE_F);
            // 2. update node
            //    2.1 map[exon] = 1
            for (j = 0; j < sg->site[GTF_don_site_id].exon_n; ++j) {
                uint32_t GTF_don_id = sg->site[GTF_don_site_id].exon_id[j];
                node_map[sg_i][GTF_don_id] |= 1;
                if (node_map[sg_i][GTF_don_id] == 3) { // update node
                    exon_t e = sg->node[GTF_don_id].e;
                    sg_update_node(sr_sg, e);
                    node_map[sg_i][GTF_don_id] |= 4;
                }
            }
            for (j = 0; j < sg->site[GTF_acc_site_id].exon_n; ++j) {
                uint32_t GTF_acc_id = sg->site[GTF_acc_site_id].exon_id[j];
                node_map[sg_i][GTF_acc_id] |= 2;
                if (node_map[sg_i][GTF_acc_id] == 3) { // update node
                    exon_t e = sg->node[GTF_acc_id].e;
                    sg_update_node(sr_sg, e);
                    node_map[sg_i][GTF_acc_id] |= 4;
                }
            }
            break;
NEXT1: sg_i++;
        }
    }*/
    // OLD-END
    // 3. alloc mem for node&site
    for (i = 0; i < sr_sg_g->SG_n; ++i) {
        sg_init_node(sr_sg_g->SG[i]); sg_init_site(sr_sg_g->SG[i]);
        if (sr_sg_g->SG[i]->node_n == 0) continue;
        sr_sg_g->SG[i]->path_map = (uint8_t**)_err_malloc(sr_sg_g->SG[i]->node_n * sizeof(uint8_t*));
        for (j = 0; j < sr_sg_g->SG[i]->node_n; ++j) sr_sg_g->SG[i]->path_map[j] = (uint8_t*)_err_calloc(sr_sg_g->SG[i]->node_n, sizeof(uint8_t));
    }
    // 4. updage_edge(edge[edge_id])
    // NEW
    i = 0, last_sg_i = 0;
    while (i < sj_n && last_sg_i < sr_sg_g->SG_n) {
        if (sr_sg_g->SG[last_sg_i]->node_n == 0) { ++last_sg_i; continue; }
        int comp_res = comp_sj_sg(sj_group[i], *(sr_sg_g->SG[last_sg_i]));
        if (comp_res < 0) { ++i; continue; }
        else if (comp_res > 0) { ++last_sg_i; continue; }
        else {
            for (sg_i = last_sg_i; sg_i < sg_g.SG_n; ++sg_i) {
                SG *sg = sg_g.SG[sg_i], *sr_sg = sr_sg_g->SG[sg_i];
                if (sr_sg->node_n == 0 || sr_sg->site_n == 0) continue;
                if (comp_sj_sg(sj_group[i], *sr_sg) < 0) break;
                // 4.0. (don_site, acc_site) => (don_site_id, acc_site_id)
                don_site_id = sg_bin_sch_site(*sr_sg, sj_group[i].don, &hit); if (hit == 0) continue;
                acc_site_id = sg_bin_sch_site(*sr_sg, sj_group[i].acc, &hit); if (hit == 0) continue;
                GTF_don_site_id = sg_bin_sch_site(*sg, sj_group[i].don, &hit); if (hit == 0) continue;
                GTF_acc_site_id = sg_bin_sch_site(*sg, sj_group[i].acc, &hit); if (hit == 0) continue;
                GTF_edge_id = sg_bin_sch_edge(*sg, GTF_don_site_id, GTF_acc_site_id, &hit); if (hit == 0 && no_novel_sj == 1) continue;
                // XXX 4.1. update edge(sj_group[i])
                sg_update_edge_pred(sr_sg, sj_group[i], don_site_id, acc_site_id);
                // 4.2. update node()
                uint32_t GTF_don_site_id = sg_bin_sch_site(*sg, sj_group[i].don, &hit);
                for (j = 0; j < sg->site[GTF_don_site_id].exon_n; ++j) {
                    uint32_t GTF_don_id = sg->site[GTF_don_site_id].exon_id[j];
                    if (node_map[sg_i][GTF_don_id] == 7) {
                        exon_t e = sg->node[GTF_don_id].e;
                        uint32_t don_id = sg_bin_sch_node(*sr_sg, e, &hit);
                        _insert(don_id, sr_sg->site[don_site_id].exon_id, sr_sg->site[don_site_id].exon_n, sr_sg->site[don_site_id].exon_m, uint32_t)
                        if (e.start == 0) _insert(don_id, sr_sg->v.next_id, sr_sg->v.next_n, sr_sg->v.next_m, uint32_t)
                    }
                }
                uint32_t GTF_acc_site_id = sg_bin_sch_site(*sg, sj_group[i].acc, &hit);
                for (j = 0; j < sg->site[GTF_acc_site_id].exon_n; ++j) {
                    uint32_t GTF_acc_id = sg->site[GTF_acc_site_id].exon_id[j];
                    if (node_map[sg_i][GTF_acc_id] == 7) {
                        exon_t e = sg->node[GTF_acc_id].e;
                        uint32_t acc_id = sg_bin_sch_node(*sr_sg, e, &hit);
                        _insert(acc_id, sr_sg->site[acc_site_id].exon_id, sr_sg->site[acc_site_id].exon_n, sr_sg->site[acc_site_id].exon_m, uint32_t)
                        if (e.end == CHR_MAX_END) _insert(acc_id, sr_sg->v.pre_id, sr_sg->v.pre_n, sr_sg->v.pre_m, uint32_t)
                    }
                }    
                for (j = 0; j < sr_sg->site[don_site_id].exon_n; ++j) {
                    uint32_t don_id = sr_sg->site[don_site_id].exon_id[j];
                    for (k = 0; k < sr_sg->site[acc_site_id].exon_n; ++k) {
                        uint32_t acc_id = sr_sg->site[acc_site_id].exon_id[k];
                        // set next/pre
                        if (sr_sg->node[don_id].next_n >= sr_sg->node[don_id].next_m) _realloc(sr_sg->node[don_id].next_id, sr_sg->node[don_id].next_m, uint32_t)
                            sr_sg->node[don_id].next_id[sr_sg->node[don_id].next_n++] = acc_id;
                        if (sr_sg->node[acc_id].pre_n >= sr_sg->node[acc_id].pre_m) _realloc(sr_sg->node[acc_id].pre_id, sr_sg->node[acc_id].pre_m, uint32_t)
                            sr_sg->node[acc_id].pre_id[sr_sg->node[acc_id].pre_n++] = don_id;

                        sr_sg->path_map[don_id][acc_id] = 1;
                    }
                }
            }
            i++;
        }
    }
    // NEW-END
    // OLD
    /*last_sg_i = 0;
    for (i = 0; i < sj_n; ++i) {
        // for sorted sj_group
        while (last_sg_i < sg_g.SG_n && (sj_group[i].tid > sg_g.SG[last_sg_i]->tid || sj_group[i].don > sg_g.SG[last_sg_i]->end))
            last_sg_i++;
        if (last_sg_i >= sg_g.SG_n) break;

        sg_i = last_sg_i;
        while (sg_i < sg_g.SG_n && sj_group[i].tid == sg_g.SG[sg_i]->tid && sj_group[i].don >= sg_g.SG[sg_i]->start && sj_group[i].acc <= sg_g.SG[sg_i]->end) {
            SG *sg = sg_g.SG[sg_i], *sr_sg = sr_sg_g->SG[sg_i];
            if (sr_sg->node_n == 0) break;
            // 4.0. (don_site, acc_site) => (don_site_id, acc_site_id)
            don_site_id = sg_bin_sch_site(*sr_sg, sj_group[i].don, &hit); if (hit == 0) goto NEXT2;
            acc_site_id = sg_bin_sch_site(*sr_sg, sj_group[i].acc, &hit); if (hit == 0) goto NEXT2;
            // XXX 4.1. update edge(sj_group[i])
            sg_update_edge_pred(sr_sg, sj_group[i], don_site_id, acc_site_id);
            // 4.2. update node()
            uint32_t GTF_don_site_id = sg_bin_sch_site(*sg, sj_group[i].don, &hit);
            for (j = 0; j < sg->site[GTF_don_site_id].exon_n; ++j) {
                uint32_t GTF_don_id = sg->site[GTF_don_site_id].exon_id[j];
                if (node_map[sg_i][GTF_don_id] == 7) {
                    exon_t e = sg->node[GTF_don_id].e;
                    uint32_t don_id = sg_bin_sch_node(*sr_sg, e, &hit);
                    _insert(don_id, sr_sg->site[don_site_id].exon_id, sr_sg->site[don_site_id].exon_n, sr_sg->site[don_site_id].exon_m, uint32_t)
                        if (e.start == 0) _insert(don_id, sr_sg->v.next_id, sr_sg->v.next_n, sr_sg->v.next_m, uint32_t)
                }
            }
            uint32_t GTF_acc_site_id = sg_bin_sch_site(*sg, sj_group[i].acc, &hit);
            for (j = 0; j < sg->site[GTF_acc_site_id].exon_n; ++j) {
                uint32_t GTF_acc_id = sg->site[GTF_acc_site_id].exon_id[j];
                if (node_map[sg_i][GTF_acc_id] == 7) {
                    exon_t e = sg->node[GTF_acc_id].e;
                    uint32_t acc_id = sg_bin_sch_node(*sr_sg, e, &hit);
                    _insert(acc_id, sr_sg->site[acc_site_id].exon_id, sr_sg->site[acc_site_id].exon_n, sr_sg->site[acc_site_id].exon_m, uint32_t)
                        if (e.end == CHR_MAX_END) _insert(acc_id, sr_sg->v.pre_id, sr_sg->v.pre_n, sr_sg->v.pre_m, uint32_t)
                }
            }    
            for (j = 0; j < sr_sg->site[don_site_id].exon_n; ++j) {
                uint32_t don_id = sr_sg->site[don_site_id].exon_id[j];
                for (k = 0; k < sr_sg->site[acc_site_id].exon_n; ++k) {
                    uint32_t acc_id = sr_sg->site[acc_site_id].exon_id[k];
                    // set next/pre
                    if (sr_sg->node[don_id].next_n >= sr_sg->node[don_id].next_m) _realloc(sr_sg->node[don_id].next_id, sr_sg->node[don_id].next_m, uint32_t)
                        sr_sg->node[don_id].next_id[sr_sg->node[don_id].next_n++] = acc_id;
                    if (sr_sg->node[acc_id].pre_n >= sr_sg->node[acc_id].pre_m) _realloc(sr_sg->node[acc_id].pre_id, sr_sg->node[acc_id].pre_m, uint32_t)
                        sr_sg->node[acc_id].pre_id[sr_sg->node[acc_id].pre_n++] = don_id;

                    sr_sg->path_map[don_id][acc_id] = 1;
                }
            }
            break;
NEXT2: sg_i++;
        }
    }*/
    // OLD-END

    for (i = 0; i < sg_g.SG_n; ++i) free(node_map[i]); free(node_map);
    int m;
    for (m = 0; m < sr_sg_g->SG_n; ++m) {
        SG *sr_sg = sr_sg_g->SG[m];
        for (i = 0; i < sr_sg->node_n; ++i) {
            for (j = i+2; j < sr_sg->node_n; ++j) {
                if (sr_sg->path_map[i][j] == 1) continue;
                for (k = i+1; k < j; ++k) {
                    if (sr_sg->path_map[i][k] > 0 && sr_sg->path_map[k][j] > 0)
                        sr_sg->path_map[i][j] = 2;
                }
            }
        }
        cal_pre_domn(sr_sg); cal_post_domn(sr_sg);
    }
    return 0;
}

// sr_sg_g.SG_n == sg_g.SG_n
SG_group *predict_SpliceGraph(SG_group sg_g, FILE *sj_p, chr_name_t *cname, int no_novel_sj)
{
    SG_group *sr_sg_g = sg_init_group(sg_g.SG_n);
    sj_t *sj_group = (sj_t*)_err_malloc(10000 * sizeof(sj_t));
    int sj_n, sj_m = 10000;
    sj_n = read_sj_group(sj_p, cname, &sj_group, sj_m);
    predict_SpliceGraph_core(sg_g, sj_group, sj_n, sr_sg_g, no_novel_sj);
    free(sj_group);
    return sr_sg_g;
}
/****************************************************************/

int pred_sg(int argc, char *argv[])
{
    return 0;
}
