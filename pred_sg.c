#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "gtf.h"
#include "build_sg.h"

extern const char PROG[20];
int pred_sg_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s pred_sg [option] <in.gtf> <in.sj> > out.sg\n\n", PROG);
    err_printf("Options:\n\n");
    err_printf("         -n --novel-sj             allow novel splice-junction in the ASM. [False]\n");
    err_printf("         -f --gtf-format  [STR]    file format of GTF input: \n");
    err_printf("                                     plain GTF file(GTF) or binary splice-graph file(SG). [SG]\n");
    err_printf("         -o --out-sg      [STR]    file name of output splice-graph. [stdout]\n");
    err_printf("\n");
    return 1;
}

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
    uint32_t GTF_don_site_id, GTF_acc_site_id;
    uint8_t **node_map = (uint8_t**)_err_malloc(sg_g.SG_n * sizeof(uint8_t*));
    for (i = 0; i < sg_g.SG_n; ++i) {
        sr_sg_g->SG[i]->tid = sg_g.SG[i]->tid;
        sr_sg_g->SG[i]->is_rev = sg_g.SG[i]->is_rev;
        // update v_start & v_end
        sg_update_node(sr_sg_g->SG[i], (exon_t){sg_g.SG[i]->tid, sg_g.SG[i]->is_rev, 0, 0});
        sg_update_node(sr_sg_g->SG[i], (exon_t){sg_g.SG[i]->tid, sg_g.SG[i]->is_rev, CHR_MAX_END, CHR_MAX_END});

        node_map[i] = (uint8_t*)_err_calloc(sg_g.SG[i]->node_n, sizeof(uint8_t));
        for (j = 0; j < sg_g.SG[i]->node_n; ++j) {
            if (sg_g.SG[i]->node[j].e.start == 0) node_map[i][j] |= 2;
            if (sg_g.SG[i]->node[j].e.end == CHR_MAX_END) node_map[i][j] |= 1;
        }
    }
    int last_sg_i=0, sg_i;
    
    i = 0, last_sg_i = 0;
    while (i < sj_n && last_sg_i < sg_g.SG_n) {
        if (sg_g.SG[last_sg_i]->node_n <= 3 || sg_g.SG[last_sg_i]->start == CHR_MAX_END || sg_g.SG[last_sg_i]->end == 0) { ++last_sg_i; continue; }
        int comp_res = comp_sj_sg(sj_group[i], *(sg_g.SG[last_sg_i]));
        if (comp_res < 0) { i++; continue; }
        else if (comp_res > 0) { last_sg_i++; continue; }
        else { 
            for (sg_i = last_sg_i; sg_i < sg_g.SG_n; ++sg_i) {
                if (comp_sj_sg(sj_group[i], *(sg_g.SG[sg_i])) < 0) break;
                SG *sg = sg_g.SG[sg_i], *sr_sg = sr_sg_g->SG[sg_i];
                // 0. search site/edge: (GTF_don_site_id, GTF_acc_site_id) => GTF_edge_id
                GTF_don_site_id = sg_bin_sch_site(sg->don_site, sg->don_site_n, sj_group[i].don, &hit); if (hit == 0) continue;
                GTF_acc_site_id = sg_bin_sch_site(sg->acc_site, sg->acc_site_n, sj_group[i].acc, &hit); if (hit == 0) continue;
                sg_bin_sch_edge(*sg, GTF_don_site_id, GTF_acc_site_id, &hit); if (hit == 0 && no_novel_sj == 1) continue;
                // 1. update node & site
                //    1.1 map[exon] = 1
                for (j = 0; j < sg->don_site[GTF_don_site_id].exon_n; ++j) {
                    uint32_t GTF_don_id = sg->don_site[GTF_don_site_id].exon_id[j];
                    node_map[sg_i][GTF_don_id] |= 1;
                    if (node_map[sg_i][GTF_don_id] == 3) { // update node
                        node_map[sg_i][GTF_don_id] |= 4;
                        exon_t e = sg->node[GTF_don_id].e;
                        sg_update_node(sr_sg, e);
                        // 1.2 update site(don_site, e.start)
                        sg_update_site(sr_sg, sj_group[i].don, DON_SITE_F);
                        if (e.start != 0) sg_update_site(sr_sg, e.start-1, ACC_SITE_F);
                    }
                }
                for (j = 0; j < sg->acc_site[GTF_acc_site_id].exon_n; ++j) {
                    uint32_t GTF_acc_id = sg->acc_site[GTF_acc_site_id].exon_id[j];
                    node_map[sg_i][GTF_acc_id] |= 2;
                    if (node_map[sg_i][GTF_acc_id] == 3) { // update node
                        exon_t e = sg->node[GTF_acc_id].e;
                        node_map[sg_i][GTF_acc_id] |= 4;
                        sg_update_node(sr_sg, e);
                        // 1.3 update site(acc_site, e.end)
                        sg_update_site(sr_sg, sj_group[i].acc, ACC_SITE_F);
                        if (e.end != CHR_MAX_END) sg_update_site(sr_sg, e.end+1, DON_SITE_F);
                    }
                }
                break;
            }
            i++;
        }
    }

    // 2. alloc mem for node&site
    for (i = 0; i < sr_sg_g->SG_n; ++i) {
        sg_init_node(sr_sg_g->SG[i]); sg_init_site(sr_sg_g->SG[i]);
        if (sr_sg_g->SG[i]->node_n == 0) continue;
        sr_sg_g->SG[i]->path_map = (uint8_t**)_err_malloc(sr_sg_g->SG[i]->node_n * sizeof(uint8_t*));
        for (j = 0; j < sr_sg_g->SG[i]->node_n; ++j) sr_sg_g->SG[i]->path_map[j] = (uint8_t*)_err_calloc(sr_sg_g->SG[i]->node_n, sizeof(uint8_t));
    }
    // 3. updage_edge(edge[edge_id])
    i = 0, last_sg_i = 0;
    while (i < sj_n && last_sg_i < sr_sg_g->SG_n) {
        if (sr_sg_g->SG[last_sg_i]->node_n <= 3 || sr_sg_g->SG[last_sg_i]->start == CHR_MAX_END || sr_sg_g->SG[last_sg_i]->end == 0) { ++last_sg_i; continue; }
        int comp_res = comp_sj_sg(sj_group[i], *(sr_sg_g->SG[last_sg_i]));
        if (comp_res < 0) { ++i; continue; }
        else if (comp_res > 0) { ++last_sg_i; continue; }
        else {
            for (sg_i = last_sg_i; sg_i < sg_g.SG_n; ++sg_i) {
                SG *sg = sg_g.SG[sg_i], *sr_sg = sr_sg_g->SG[sg_i];
                if (sr_sg->node_n <= 3 || (sr_sg->don_site_n+sr_sg->acc_site_n) <= 1 || sr_sg->start == CHR_MAX_END || sr_sg->end == 0) continue;
                if (comp_sj_sg(sj_group[i], *sr_sg) < 0) break;
                // 3.0. (don_site, acc_site) => (don_site_id, acc_site_id)
                don_site_id = sg_bin_sch_site(sr_sg->don_site, sr_sg->don_site_n, sj_group[i].don, &hit); if (hit == 0) continue;
                acc_site_id = sg_bin_sch_site(sr_sg->acc_site, sr_sg->acc_site_n, sj_group[i].acc, &hit); if (hit == 0) continue;
                GTF_don_site_id = sg_bin_sch_site(sg->don_site, sg->don_site_n, sj_group[i].don, &hit);
                GTF_acc_site_id = sg_bin_sch_site(sg->acc_site, sg->acc_site_n, sj_group[i].acc, &hit);
                sg_bin_sch_edge(*sg, GTF_don_site_id, GTF_acc_site_id, &hit); if (hit == 0 && no_novel_sj == 1) continue;
                // XXX 3.1. update edge(sj_group[i])
                sg_update_edge_pred(sr_sg, sj_group[i], don_site_id, acc_site_id);
                // 3.2. update node()
                for (j = 0; j < sg->don_site[GTF_don_site_id].exon_n; ++j) {
                    uint32_t GTF_don_id = sg->don_site[GTF_don_site_id].exon_id[j];
                    if (node_map[sg_i][GTF_don_id] == 7) {
                        exon_t e = sg->node[GTF_don_id].e;
                        uint32_t don_id = sg_bin_sch_node(*sr_sg, e, &hit); if (hit == 0) err_fatal_simple("Can not hit node. (pred_sg 1)\n");
                        _insert(don_id, sr_sg->don_site[don_site_id].exon_id, sr_sg->don_site[don_site_id].exon_n, sr_sg->don_site[don_site_id].exon_m, uint32_t)
                        // update v_start
                        if (e.start == 0) {
                            _insert(don_id, sr_sg->node[0].next_id, sr_sg->node[0].next_n, sr_sg->node[0].next_m, uint32_t)
                            _insert(0, sr_sg->node[don_id].pre_id, sr_sg->node[don_id].pre_n, sr_sg->node[don_id].pre_m, uint32_t)
                        }
                    }
                }
                for (j = 0; j < sg->acc_site[GTF_acc_site_id].exon_n; ++j) {
                    uint32_t GTF_acc_id = sg->acc_site[GTF_acc_site_id].exon_id[j];
                    if (node_map[sg_i][GTF_acc_id] == 7) {
                        exon_t e = sg->node[GTF_acc_id].e;
                        uint32_t acc_id = sg_bin_sch_node(*sr_sg, e, &hit);
                        _insert(acc_id, sr_sg->acc_site[acc_site_id].exon_id, sr_sg->acc_site[acc_site_id].exon_n, sr_sg->acc_site[acc_site_id].exon_m, uint32_t)
                        // update v_end
                        if (e.end == CHR_MAX_END) {
                            _insert(acc_id, sr_sg->node[sr_sg->node_n-1].pre_id, sr_sg->node[sr_sg->node_n-1].pre_n, sr_sg->node[sr_sg->node_n-1].pre_m, uint32_t)
                            _insert((uint32_t)sr_sg->node_n-1, sr_sg->node[acc_id].next_id, sr_sg->node[acc_id].next_n, sr_sg->node[acc_id].next_m, uint32_t)
                        }
                    }
                }    
                for (j = 0; j < sr_sg->don_site[don_site_id].exon_n; ++j) {
                    uint32_t don_id = sr_sg->don_site[don_site_id].exon_id[j];
                    for (k = 0; k < sr_sg->acc_site[acc_site_id].exon_n; ++k) {
                        uint32_t acc_id = sr_sg->acc_site[acc_site_id].exon_id[k];
                        // set next/pre
                        if (sr_sg->node[don_id].next_n == sr_sg->node[don_id].next_m) _realloc(sr_sg->node[don_id].next_id, sr_sg->node[don_id].next_m, uint32_t)
                            sr_sg->node[don_id].next_id[sr_sg->node[don_id].next_n++] = acc_id;
                        if (sr_sg->node[acc_id].pre_n == sr_sg->node[acc_id].pre_m) _realloc(sr_sg->node[acc_id].pre_id, sr_sg->node[acc_id].pre_m, uint32_t)
                            sr_sg->node[acc_id].pre_id[sr_sg->node[acc_id].pre_n++] = don_id;

                        sr_sg->path_map[don_id][acc_id] = 1;
                    }
                }
                break;
            }
            i++;
        }
    }
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

const struct option pred_sg_long_opt [] = {
    { "novel-sj", 0, NULL, 'n' },
    { "gtf-format", 1, NULL, 'f' }, 
    { "out-sg", 1, NULL, 'o' },

    { 0, 0, 0, 0}
};

int pred_sg(int argc, char *argv[])
{
    int c;
    int no_novel_sj = 1, SG_format = 1; char out_file[1024]="-";
    while ((c = getopt_long(argc, argv, "nf:o:", pred_sg_long_opt, NULL)) >= 0) {
        switch (c) {
            case 'n': no_novel_sj = 0; break;
            case 'f': if (strcmp(optarg, "GTF")==0) SG_format = 0; 
                      else if (strcmp(optarg, "SG")==0) SG_format = 1; 
                      else return pred_sg_usage();
                  break;
            case 'o': strcpy(out_file, optarg); break;
            default: err_printf("Error: unknown optin: %s.\n", optarg);
                     return pred_sg_usage();
        }
    }
    if (argc - optind != 2) return pred_sg_usage();
    FILE *sj_fp;

    sj_fp = xopen(argv[optind+1], "r");

    chr_name_t *cname = chr_name_init();
    // build splice-graph with GTF
    SG_group *sg_g;
    if (SG_format) sg_g = sg_restore(argv[optind]);
    else  {
        FILE *gtf_fp = xopen(argv[optind], "r");
        sg_g = construct_SpliceGraph(gtf_fp, cname);
        err_fclose(gtf_fp);
    }
    // predict splice-graph with GTF-based splice-graph and splice-junciton
    SG_group *sr_sg_g = predict_SpliceGraph(*sg_g, sj_fp, cname, no_novel_sj);

    // dump predicted splice-graph to file
    sg_dump(*sg_g, out_file);

    sg_free_group(sg_g); sg_free_group(sr_sg_g);
    err_fclose(sj_fp); chr_name_free(cname);
    return 0;
}
