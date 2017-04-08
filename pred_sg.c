#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "gtf.h"
#include "build_sg.h"
#include "bam2sj.h"
#include "kstring.h"

extern const char PROG[20];
int pred_sg_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s pred_sg [option] <in.gtf> <in.sj>\n\n", PROG);
    err_printf("Options:\n\n");
    err_printf("         -n --novel-sj             allow novel splice-junction in the ASM. [False]\n");
    err_printf("         -N --novel-com            allow novel splice-junction in the ASM. [False]\n");
    err_printf("         -m --use-multi            use both uniq- and multi-mapped reads in the bam input.[False (uniq only)]\n");
    err_printf("         -o --out-prefix  [STR]    prefix of output splice-graph file. [in.sj]\n");
    err_printf("\n");
    return 1;
}


sg_para *sg_init_para(void)
{
    sg_para *sgp = (sg_para*)_err_malloc(sizeof(sg_para));
    sgp->rep_n = NULL; sgp->in_name = NULL;
    sgp->sam_n = 0; sgp->tol_rep_n = 0; sgp->BAM_input=1; 
    sgp->no_novel_sj = 1; sgp->no_novel_com = 1; sgp->only_novel = 0;
    sgp->use_multi = 0;
    sgp->merge_out = 0;
    return sgp;
}

// ':' separates samples, ',' separates replicates
int sg_par_input(sg_para *sgp, char *in)
{
    ks_tokaux_t aux1, aux2; char *p1, *p2;
    kstring_t *s1=(kstring_t*)_err_calloc(1, sizeof(kstring_t)),  *s2=(kstring_t*)_err_calloc(1, sizeof(kstring_t));
    int sam_n = 0, rep_n = 0;
    for (p1 = kstrtok(in, ":", &aux1); p1; p1 = kstrtok(0, 0, &aux1)) {
        if (p1 != NULL) {
            sam_n++;
            kputsn(p1, aux1.p-p1, s1);
            for (p2 = kstrtok(s1->s, ",", &aux2); p2; p2 = kstrtok(0, 0, &aux2))
                rep_n++;
            free(s1->s); ks_release(s1);
        }
    }
    sgp->rep_n = (int*)_err_malloc(sam_n * sizeof(int));
    sgp->in_name = (char**)_err_malloc(rep_n * sizeof(char*));
    sgp->tol_rep_n = rep_n;
    int i = 0;
    for (p1 = kstrtok(in, ":", &aux1); p1; p1 = kstrtok(0, 0, &aux1)) {
        if (p1 != NULL) {
            kputsn(p1, aux1.p-p1, s1);
            for (p2 = kstrtok(s1->s, ",", &aux2); p2; p2 = kstrtok(0, 0, &aux2)) {
                kputsn(p2, aux2.p-p2, s2);
                sgp->in_name[i++] = strdup(s2->s);
                free(s2->s); ks_release(s2);
                sgp->rep_n[sgp->sam_n]++;
            }
            free(s1->s); ks_release(s1);
            sgp->sam_n++;
        }
    }
    free(s1); free(s2);
    return sgp->tol_rep_n;
}

void sg_free_para(sg_para *sgp)
{
    if (sgp->in_name != NULL) {
        int i;
        for (i = 0; i < sgp->tol_rep_n; ++i)
            free(sgp->in_name[i]);
        free(sgp->in_name);
    }
    if (sgp->rep_n != NULL) free(sgp->rep_n);
    free(sgp);
}

/****************************************************************
 * predict splice graph with splice-junctions (short-read data) *
 * based on reference-splice-graph                              *
 ****************************************************************/
int sg_update_edge_pred(SG *sg, sj_t sj, uint32_t don_site_id, uint32_t acc_site_id)
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

int predict_SpliceGraph_core(SG_group sg_g, sj_t *sj_group, int sj_n, SG_group *sr_sg_g, sg_para *sgp)
{
    int no_novel_sj = sgp->no_novel_sj, no_novel_com = sgp->no_novel_com;
    int i, j, k, hit; uint32_t don_site_id, acc_site_id;
    uint32_t GTF_don_site_id, GTF_acc_site_id;
    uint8_t **node_map = (uint8_t**)_err_malloc(sg_g.SG_n * sizeof(uint8_t*));
    for (i = 0; i < sg_g.SG_n; ++i) {
        sr_sg_g->SG[i]->tid = sg_g.SG[i]->tid;
        sr_sg_g->SG[i]->is_rev = sg_g.SG[i]->is_rev;
        // update v_start & v_end
        sg_update_node(sr_sg_g->SG[i], (exon_t){sg_g.SG[i]->tid, sg_g.SG[i]->is_rev, 0, 0}, 0, 0);
        sg_update_node(sr_sg_g->SG[i], (exon_t){sg_g.SG[i]->tid, sg_g.SG[i]->is_rev, MAX_SITE, MAX_SITE}, MAX_SITE, MAX_SITE);

        node_map[i] = (uint8_t*)_err_calloc(sg_g.SG[i]->node_n, sizeof(uint8_t));
        for (j = 0; j < sg_g.SG[i]->node_n; ++j) {
            if (sg_g.SG[i]->node[j].is_init == 1) node_map[i][j] |= 2;
            if (sg_g.SG[i]->node[j].is_termi == 1) node_map[i][j] |= 1;
        }
    }
    int last_sg_i=0, sg_i;
    
    i = 0, last_sg_i = 0;
    while (i < sj_n && last_sg_i < sg_g.SG_n) {
        if (sg_g.SG[last_sg_i]->node_n <= 3 || sg_g.SG[last_sg_i]->start == MAX_SITE || sg_g.SG[last_sg_i]->end == 0) { ++last_sg_i; continue; }
        int comp_res = comp_sj_sg(sj_group[i], *(sg_g.SG[last_sg_i]));
        if (comp_res < 0) { i++; continue; }
        else if (comp_res > 0) { last_sg_i++; continue; }
        else { 
            for (sg_i = last_sg_i; sg_i < sg_g.SG_n; ++sg_i) {
                if (comp_sj_sg(sj_group[i], *(sg_g.SG[sg_i])) < 0) break;
                SG *sg = sg_g.SG[sg_i], *sr_sg = sr_sg_g->SG[sg_i];
                SGsite *don_site = sg->don_site, *acc_site = sg->acc_site; int32_t acc_n = sg->acc_site_n, don_n = sg->don_site_n;
                SGnode *node = sg->node;
                // 0. search site/edge: (GTF_don_site_id, GTF_acc_site_id) => GTF_edge_id
                GTF_don_site_id = sg_bin_sch_site(don_site, don_n, sj_group[i].don, &hit); if (hit == 0) continue;
                GTF_acc_site_id = sg_bin_sch_site(acc_site, acc_n, sj_group[i].acc, &hit); if (hit == 0) continue;
                sg_bin_sch_edge(sg, GTF_don_site_id, GTF_acc_site_id, &hit); if (hit == 0 && no_novel_sj == 1) continue;
                // 1. update node & site
                //    1.1 map[exon] = 1
                for (j = 0; j < don_site[GTF_don_site_id].exon_n; ++j) {
                    uint32_t GTF_don_id = don_site[GTF_don_site_id].exon_id[j];
                    node_map[sg_i][GTF_don_id] |= 1;
                    if (node_map[sg_i][GTF_don_id] == 3) { // update node
                        node_map[sg_i][GTF_don_id] |= 4;
                        exon_t node_e = node[GTF_don_id].node_e;
                        int32_t start = node[GTF_don_id].start, end = node[GTF_don_id].end;
                        sg_update_node(sr_sg, node_e, start, end);
                        // 1.2 update site(don_site, e.start)
                        sg_update_site(sr_sg, sj_group[i].don, DON_SITE_F);
                        if (node_e.start != 0) sg_update_site(sr_sg, node_e.start-1, ACC_SITE_F);
                    }
                }
                for (j = 0; j < acc_site[GTF_acc_site_id].exon_n; ++j) {
                    uint32_t GTF_acc_id = acc_site[GTF_acc_site_id].exon_id[j];
                    node_map[sg_i][GTF_acc_id] |= 2;
                    if (node_map[sg_i][GTF_acc_id] == 3) { // update node
                        exon_t node_e = node[GTF_acc_id].node_e;
                        int32_t start = node[GTF_acc_id].start, end = node[GTF_acc_id].end;
                        node_map[sg_i][GTF_acc_id] |= 4;
                        sg_update_node(sr_sg, node_e, start, end);
                        // 1.3 update site(acc_site, e.end)
                        sg_update_site(sr_sg, sj_group[i].acc, ACC_SITE_F);
                        if (node_e.end != MAX_SITE) sg_update_site(sr_sg, node_e.end+1, DON_SITE_F);
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
    }
    // 3. updage_edge(edge[edge_id])
    i = 0, last_sg_i = 0;
    while (i < sj_n && last_sg_i < sr_sg_g->SG_n) {
        if (sr_sg_g->SG[last_sg_i]->node_n <= 3 || sr_sg_g->SG[last_sg_i]->start == MAX_SITE || sr_sg_g->SG[last_sg_i]->end == 0) { ++last_sg_i; continue; }
        int comp_res = comp_sj_sg(sj_group[i], *(sr_sg_g->SG[last_sg_i]));
        if (comp_res < 0) { ++i; continue; }
        else if (comp_res > 0) { ++last_sg_i; continue; }
        else {
            for (sg_i = last_sg_i; sg_i < sg_g.SG_n; ++sg_i) {
                SG *sg = sg_g.SG[sg_i], *sr_sg = sr_sg_g->SG[sg_i];
                SGsite *acc_site = sr_sg->acc_site, *don_site = sr_sg->don_site; int32_t acc_n = sr_sg->acc_site_n, don_n = sr_sg->don_site_n;
                SGnode *node = sr_sg->node; int32_t node_n = sr_sg->node_n;
                SGnode *sg_node = sg->node; SGsite *sg_acc = sg->acc_site, *sg_don = sg->don_site;
                int32_t sg_don_n = sg->don_site_n, sg_acc_n = sg->acc_site_n;
                if (node_n <= 3 || (don_n+acc_n) <= 1 || sr_sg->start == MAX_SITE || sr_sg->end == 0) continue;
                if (comp_sj_sg(sj_group[i], *sr_sg) < 0) break;
                // 3.0. (don_site, acc_site) => (don_site_id, acc_site_id)
                don_site_id = sg_bin_sch_site(don_site, don_n, sj_group[i].don, &hit); if (hit == 0) continue;
                acc_site_id = sg_bin_sch_site(acc_site, acc_n, sj_group[i].acc, &hit); if (hit == 0) continue;
                GTF_don_site_id = _err_sg_bin_sch_site(sg_don, sg_don_n, sj_group[i].don, &hit);
                GTF_acc_site_id = _err_sg_bin_sch_site(sg_acc, sg_acc_n, sj_group[i].acc, &hit);
                sg_bin_sch_edge(sg, GTF_don_site_id, GTF_acc_site_id, &hit); 
                if (hit == 0 && no_novel_sj == 1) continue;
                else if (hit == 0) sj_group[i].is_anno = 0;

                // 3.1. update edge(sj_group[i])
                sg_update_edge_pred(sr_sg, sj_group[i], don_site_id, acc_site_id);
                // 3.2. update node()
                for (j = 0; j < sg_don[GTF_don_site_id].exon_n; ++j) {
                    uint32_t GTF_don_id = sg_don[GTF_don_site_id].exon_id[j];
                    if (node_map[sg_i][GTF_don_id] == 7) {
                        exon_t node_e = sg_node[GTF_don_id].node_e;
                        uint32_t don_id = _err_sg_bin_sch_node(sr_sg, node_e, &hit);
                        _insert(don_id, don_site[don_site_id].exon_id, don_site[don_site_id].exon_n, don_site[don_site_id].exon_m, uint32_t)
                        // update v_start
                        if (sg_node[GTF_don_id].is_init == 1) {
                            _insert(don_id, node[0].next_id, node[0].next_n, node[0].next_m, uint32_t)
                            _insert(0, node[don_id].pre_id, node[don_id].pre_n, node[don_id].pre_m, uint32_t)
                            node[don_id].is_init = 1;
                        }
                    }
                }
                for (j = 0; j < sg_acc[GTF_acc_site_id].exon_n; ++j) {
                    uint32_t GTF_acc_id = sg_acc[GTF_acc_site_id].exon_id[j];
                    if (node_map[sg_i][GTF_acc_id] == 7) {
                        exon_t node_e = sg_node[GTF_acc_id].node_e;
                        uint32_t acc_id = _err_sg_bin_sch_node(sr_sg, node_e, &hit);
                        _insert(acc_id, acc_site[acc_site_id].exon_id, acc_site[acc_site_id].exon_n, acc_site[acc_site_id].exon_m, uint32_t)
                        // update v_end
                        if (sg_node[GTF_acc_id].is_termi == 1) {
                            _insert(acc_id, node[node_n-1].pre_id, node[node_n-1].pre_n, node[node_n-1].pre_m, uint32_t)
                            _insert((uint32_t)node_n-1, node[acc_id].next_id, node[acc_id].next_n, node[acc_id].next_m, uint32_t)
                            node[acc_id].is_termi = 1;
                        }
                    }
                }    
                for (j = 0; j < don_site[don_site_id].exon_n; ++j) {
                    uint32_t don_id = don_site[don_site_id].exon_id[j];
                    exon_t don_e = node[don_id].node_e;
                    uint32_t GTF_don_id = _err_sg_bin_sch_node(sg, don_e, &hit);
                    for (k = 0; k < acc_site[acc_site_id].exon_n; ++k) {
                        uint32_t acc_id = acc_site[acc_site_id].exon_id[k];
                        exon_t acc_e = node[acc_id].node_e;
                        uint32_t GTF_acc_id = _err_sg_bin_sch_node(sg, acc_e, &hit);
                        // set next/pre
                        if (no_novel_com) {
                            // if (GTF_don_id.next == GTF_acc_id)
                            int m;
                            for (m = 0; m < sg_node[GTF_don_id].next_n; ++m) {
                                if (GTF_acc_id == sg_node[GTF_don_id].next_id[m]) {
                                    if (node[don_id].next_n == node[don_id].next_m) _realloc(node[don_id].next_id, node[don_id].next_m, uint32_t)
                                        node[don_id].next_id[node[don_id].next_n++] = acc_id;
                                    if (node[acc_id].pre_n == node[acc_id].pre_m) _realloc(node[acc_id].pre_id, node[acc_id].pre_m, uint32_t)
                                        node[acc_id].pre_id[node[acc_id].pre_n++] = don_id;
                                    node[don_id].e_site_id = don_site_id;
                                    node[acc_id].s_site_id = acc_site_id;
                                    break;
                                }
                            }
                        } else {
                            if (node[don_id].next_n == node[don_id].next_m) _realloc(node[don_id].next_id, node[don_id].next_m, uint32_t)
                                node[don_id].next_id[node[don_id].next_n++] = acc_id;
                            if (node[acc_id].pre_n == node[acc_id].pre_m) _realloc(node[acc_id].pre_id, node[acc_id].pre_m, uint32_t)
                                node[acc_id].pre_id[node[acc_id].pre_n++] = don_id;
                            node[don_id].e_site_id = don_site_id;
                            node[acc_id].s_site_id = acc_site_id;
                        }
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
        cal_pre_domn(sr_sg); cal_post_domn(sr_sg);
    }
    return 0;
}

// sr_sg_g.SG_n == sg_g.SG_n
SG_group *predict_SpliceGraph(SG_group sg_g, sj_t *sj_group, int sj_n, sg_para *sgp)
{
    print_format_time(stderr); err_printf("[%s] predicting splice-graph with splice-junction and GTF-SG ...\n", __func__);
    SG_group *sr_sg_g = sg_init_group(sg_g.SG_n);

    int i;
    if (sg_g.cname->chr_n > sr_sg_g->cname->chr_m) {
        sr_sg_g->cname->chr_name = (char**)_err_realloc(sr_sg_g->cname->chr_name, sg_g.cname->chr_n * sizeof(char*));
        for (i = sr_sg_g->cname->chr_m; i < sg_g.cname->chr_n; ++i) sr_sg_g->cname->chr_name[i] = (char*)_err_malloc(100 * sizeof(char));
        sr_sg_g->cname->chr_m = sg_g.cname->chr_n;
    }
    sr_sg_g->cname->chr_n = sg_g.cname->chr_n;
    for (i = 0; i < sr_sg_g->cname->chr_n; ++i) strcpy(sr_sg_g->cname->chr_name[i], sg_g.cname->chr_name[i]);

    predict_SpliceGraph_core(sg_g, sj_group, sj_n, sr_sg_g, sgp);
    free(sj_group);
    print_format_time(stderr); err_printf("[%s] predicting splice-graph with splice-junction and GTF-SG done!\n", __func__);
    return sr_sg_g;
}
/****************************************************************/

const struct option pred_sg_long_opt [] = {
    { "novel-sj", 0, NULL, 'n' },
    { "novel-com", 0, NULL, 'N' },
    { "use-multi", 0, NULL, 'm' },
    { "out-prefix", 1, NULL, 'o' },

    { 0, 0, 0, 0}
};

int pred_sg(int argc, char *argv[])
{
    int c;
    sg_para *sgp = sg_init_para(); char out_prefix[1024]="";
    while ((c = getopt_long(argc, argv, "nNmo:", pred_sg_long_opt, NULL)) >= 0) {
        switch (c) {
            case 'n': sgp->no_novel_sj=0, sgp->no_novel_com=0; break;
            case 'N': sgp->no_novel_com = 0; break;
            case 'o': strcpy(out_prefix, optarg); break;
            default: err_printf("Error: unknown optin: %s.\n", optarg);
                     return pred_sg_usage();
        }
    }
    if (argc - optind != 2) return pred_sg_usage();

    gzFile genome_fp = gzopen(argv[optind], "r");
    if (genome_fp == NULL) err_fatal(__func__, "Can not open genome file. %s\n", argv[optind]);

    // build splice-graph with GTF
    SG_group *sg_g;
    FILE *gtf_fp = xopen(argv[optind+1], "r");
    chr_name_t *cname = chr_name_init();
    sg_g = construct_SpliceGraph(gtf_fp, cname);
    err_fclose(gtf_fp);  chr_name_free(cname);

    // get splice-junction
    /* based on .sj file
    FILE *sj_fp = xopen(argv[optind+2], "r");
    sj_t *sj_group = (sj_t*)_err_malloc(10000 * sizeof(sj_t));
    int sj_n, sj_m = 10000;
    sj_n = read_sj_group(sj_fp, sg_g->cname, &sj_group, sj_m);

    SG_group *sr_sg_g = predict_SpliceGraph(*sg_g, sj_group, sj_n, no_novel_sj, no_novel_com);
    err_fclose(sj_fp); sg_free_group(sg_g);
    */
    // based on .bam file
    samFile *in; bam_hdr_t *h; bam1_t *b;
    if ((in = sam_open(argv[optind+2], "rb")) == NULL) err_fatal_core(__func__, "Cannot open \"%s\"\n", argv[optind+2]);
    if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", argv[optind+2]);
    b = bam_init1(); 
    sj_t *sj_group = (sj_t*)_err_malloc(10000 * sizeof(sj_t)); int sj_m = 10000;
    int seq_n = 0, seq_m; kseq_t *seq = kseq_load_genome(genome_fp, &seq_n, &seq_m);
    int sj_n = bam2sj_core(in, h, b, seq, seq_n, &sj_group, sj_m); 
    // predict splice-graph with GTF-based splice-graph and splice-junciton
    SG_group *sr_sg_g = predict_SpliceGraph(*sg_g, sj_group, sj_n, sgp);

    // dump predicted splice-graph to file
    if (strlen(out_prefix) == 0) strcpy(out_prefix, argv[optind+2]);
    int i; for (i = 0; i < seq_n; ++i) {
        free(seq[i].name.s); free(seq[i].seq.s);
    } free(seq);
    sg_free_group(sg_g); sg_free_group(sr_sg_g); gzclose(genome_fp); sg_free_para(sgp);
    return 0;
}
