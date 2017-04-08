#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include "utils.h"
#include "gtf.h"
#include "build_sg.h"
#include "pred_sg.h"
#include "bam2sj.h"
#include "kstring.h"

extern const char PROG[20];
int pred_asm_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s asm [option] <ref.fa> <in.gtf> <in.bam/sj>\n\n", PROG);
    err_printf("Note:    for multi-sample and multi-replicate should be this format: \n");
    err_printf("             \"SAM1-REP1,REP2,REP3;SAM2-REP1,REP2,REP3\"\n");
    err_printf("         use \':\' to separate samples, \',\' to separate replicates.\n\n");
    err_printf("Options:\n\n");
    err_printf("         -n --novel-sj             allow novel splice-junction in the ASM. [False]\n");
    err_printf("         -N --novel-com            allow novel combination of known exons in the ASM. [False]\n");
    err_printf("         -s --sj-file              input with splice-junction file instead of BAM file. [False]\n");
    err_printf("                                   with splice-junction input, the .cnt output will have no count information.\n");
    err_printf("         -m --use-multi            use both uniq- and multi-mapped reads in the bam input.[False (uniq only)]\n");
    err_printf("         -o --output      [STR]    prefix of file name of output ASM & COUNT. [in.bam/sj]\n");
    err_printf("                                   prefix.ASM & prefix.JCNT & prefix.ECNT\n");
	err_printf("\n");
	return 1;
}

/*****************************
 *       generate ASM        *
 *****************************/
SGasm *sg_init_asm(uint32_t sg_id, uint32_t v_start, uint32_t v_end)
{
    SGasm *sg_asm = (SGasm*)_err_malloc(sizeof(SGasm));
    sg_asm->SG_id = sg_id;
    sg_asm->v_start = v_start; sg_asm->v_end = v_end;
    sg_asm->start = MAX_SITE; sg_asm->end = 0;
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

SGasm_group *sg_realloc_asm_group(SGasm_group *asm_g)
{
    asm_g->sg_asm_m <<= 1;
    asm_g->sg_asm = (SGasm**)_err_realloc(asm_g->sg_asm, asm_g->sg_asm_m * sizeof(SGasm*));
    int i;
    for (i = asm_g->sg_asm_m >> 1; i < asm_g->sg_asm_m; ++i) asm_g->sg_asm[i] = sg_init_asm(0, 0, 0);
    return asm_g;
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

void sg_update_asm_edge(SG *sg, SGasm *sg_asm, uint32_t pre_id, uint32_t cur_id)
{
    int hit = 0;
    if (sg->node[pre_id].node_e.end == 0 || sg->node[cur_id].node_e.start == MAX_SITE) return;
    uint32_t pre_site_id = _err_sg_bin_sch_site(sg->don_site, sg->don_site_n, sg->node[pre_id].node_e.end+1, &hit);
    uint32_t cur_site_id = _err_sg_bin_sch_site(sg->acc_site, sg->acc_site_n, sg->node[cur_id].node_e.start-1, &hit);
    uint32_t edge_i = _err_sg_bin_sch_edge(sg, pre_site_id, cur_site_id, &hit); 
    _bin_insert(edge_i, sg_asm->edge_id, sg_asm->edge_n, sg_asm->edge_m, uint32_t)
}

void sg_update_asm(SG *sg, SGasm *sg_asm, uint32_t pre_id, uint32_t cur_id)
{
    _bin_insert(cur_id, sg_asm->node_id, sg_asm->node_n, sg_asm->node_m, uint32_t)
    if (sg->node[cur_id].node_e.start < sg_asm->start) sg_asm->start = sg->node[cur_id].node_e.start;
    if (sg->node[cur_id].node_e.end > sg_asm->end) sg_asm->end = sg->node[cur_id].node_e.end;
    sg->node[cur_id].is_asm = 1;
    sg_update_asm_edge(sg, sg_asm, pre_id, cur_id);
}

void sub_splice_graph(SG *sg, int **node_visit, SGasm *sg_asm, uint32_t cur_id, uint32_t e_id)
{
    if ((*node_visit)[cur_id] == 1) return; else (*node_visit)[cur_id] = 1;

    if (cur_id == e_id) return;
    int i;
    for (i = 0; i < sg->node[cur_id].next_n; ++i) {
        if (sg->node[cur_id].next_id[i] == e_id) {
            sg_update_asm_edge(sg, sg_asm, cur_id, e_id);
            continue;
        } else {
            sg_update_asm(sg, sg_asm, cur_id, sg->node[cur_id].next_id[i]); 
            sub_splice_graph(sg, node_visit, sg_asm, sg->node[cur_id].next_id[i], e_id);
        }
    }
}

int sg_asm_group_add(SGasm_group *asm_g, SGasm *sg_asm)
{
    if (asm_g->sg_asm_n == asm_g->sg_asm_m) sg_realloc_asm_group(asm_g);
    SGasm *a = asm_g->sg_asm[asm_g->sg_asm_n];

    a->SG_id = sg_asm->SG_id;
    a->v_start = sg_asm->v_start; a->v_end = sg_asm->v_end;
    if (sg_asm->node_n > a->node_m) a->node_id = (uint32_t*)_err_realloc(a->node_id, sg_asm->node_n * sizeof(uint32_t));
    if (sg_asm->edge_n > a->edge_m) a->edge_id = (uint32_t*)_err_realloc(a->edge_id, sg_asm->edge_n * sizeof(uint32_t));
    a->node_n = a->node_m = sg_asm->node_n; a->edge_n = a->edge_m = sg_asm->edge_n;
    int i;
    for (i = 0; i < sg_asm->node_n; ++i) a->node_id[i] = sg_asm->node_id[i];
    for (i = 0; i < sg_asm->edge_n; ++i) a->edge_id[i] = sg_asm->edge_id[i];
    a->start = sg_asm->start; a->end = sg_asm->end;

    asm_g->sg_asm_n++;
    return asm_g->sg_asm_n;
}

int check_novel_asm(SG *sg, SGasm *sg_asm)
{
    int i;
    for (i = 0; i < sg_asm->edge_n; ++i) {
        uint32_t eid = sg_asm->edge_id[i];
        if (sg->edge[eid].is_anno == 0) return 1;
    }
    return 0;
}

int check_uniq_asm(SG *sg, SGasm *sg_asm)
{
    int i;
    for (i = 0; i < sg_asm->edge_n; ++i) {
        uint32_t eid = sg_asm->edge_id[i];
        if (sg->edge[eid].uniq_c == 0) return 0;
    }
    return 1;
}

SGasm_group *gen_asm(SG_group *sg_g, sg_para *sgp)
{
    print_format_time(stderr); err_printf("[%s] generating alternative-splice-module with predicted SJ-SG ...\n", __func__);
    int only_novel = sgp->only_novel, use_multi = sgp->use_multi;
    int entry_n, exit_n; uint32_t *entry, *exit;
    int sg_i;
    SGasm_group *asm_g = sg_init_asm_group();
    for (sg_i = 0; sg_i < sg_g->SG_n; ++sg_i) {
        SG *sg = sg_g->SG[sg_i];
        cal_cand_node(*sg, &entry, &exit, &entry_n, &exit_n);
        if (entry_n == 0 || exit_n == 0) goto END;

        int i, j, hit;
        for (i = 0; i < entry_n; ++i) {
            hit = 0;
            for (j = 0; j < exit_n; ++j) {
                int post_domn_n = sg->node[entry[i]].post_domn_n;
                int pre_domn_n = sg->node[exit[j]].pre_domn_n;
                if (post_domn_n > 1 && pre_domn_n > 1 
                        && sg->node[entry[i]].post_domn[1] == exit[j] 
                        && sg->node[exit[j]].pre_domn[1] == entry[i]) {
                    SGasm *sg_asm = sg_init_asm(sg_i, entry[i], exit[j]);
                    int *node_visit = (int*)_err_calloc(sg->node_n, sizeof(int));
                    sub_splice_graph(sg, &node_visit, sg_asm, entry[i], exit[j]);
                    free(node_visit);
                    if ((only_novel == 0 || check_novel_asm(sg, sg_asm) == 1)
                    && (use_multi == 1 || check_uniq_asm(sg, sg_asm) == 1))
                        sg_asm_group_add(asm_g, sg_asm);
                    sg_free_asm(sg_asm);
                    hit = 1; break;
                }
            }
            if (hit) continue;
        }
        END: free(entry); free(exit);
    }
    print_format_time(stderr); err_printf("[%s] generating alternative-splice-module with predicted SJ-SG done!\n", __func__);
    return asm_g;
}

int cal_asm_exon_cnt(SG_group *sg_g, samFile *in, bam_hdr_t *h, bam1_t *b)
{
    print_format_time(stderr); err_printf("[%s] calculating read count for AS exons ...\n", __func__);
    int i, j, last_sg_i = 0;
    while (sam_read1(in, h, b) >= 0) {
        int bam_start = b->core.pos+1;
        int bam_end = b->core.pos+bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
        uint8_t is_uniq = bam_is_uniq_NH(b);

        for (i = last_sg_i; i < sg_g->SG_n; ++i) {
            SG *sg = sg_g->SG[i];
            int32_t tid = sg->tid, start = sg->start, end = sg->end;
            if (start == MAX_SITE || end == 0 || tid < b->core.tid || (tid == b->core.tid && end <= bam_start)) {
                if (i == last_sg_i) last_sg_i++; continue;
            } else if (tid > b->core.tid || (tid == b->core.tid && start >= bam_end)) break;
            else {
                for (j = 0; j < sg->node_n; ++j) {
                    if (sg->node[j].is_asm) {
                        if (sg->node[j].start <= bam_start && sg->node[j].end >= bam_end) {
                            if (is_uniq) sg->node[j].uniq_c++;
                            else sg->node[j].multi_c++;
                        }
                    }
                }
            }
        }
        if (last_sg_i == sg_g->SG_n) break;
    }
    print_format_time(stderr); err_printf("[%s] calculating read count for AS exons done!\n", __func__);
    return 0;
}

int asm_output(char *in_fn, char *prefix, SG_group *sg_g, SGasm_group *asm_g, sg_para *sgp)
{
    int i, j, out_n=3;
    char suf[3][10] = { ".ASM", ".JCNT", ".ECNT" };
    char suff[20] = "";
    if (sgp->use_multi==1) strcat(suff, ".multi");
    if (sgp->no_novel_sj==1) strcat(suff, ".anno");
    if (sgp->only_novel==1) strcat(suff, ".novel");
    char **out_fn = (char**)_err_malloc(sizeof(char*) * out_n);
    if (strlen(prefix) == 0) {
        for (i = 0; i < out_n; ++i) {
            out_fn[i] = (char*)_err_malloc(strlen(in_fn)+10); strcpy(out_fn[i], in_fn); strcat(out_fn[i], suff); strcat(out_fn[i], suf[i]);
        }
    } else {
        for (i = 0; i < out_n; ++i) {
            out_fn[i] = (char*)_err_malloc(strlen(prefix)+10); strcpy(out_fn[i], prefix); strcat(out_fn[i], suff); strcat(out_fn[i], suf[i]);
        }
    }

    FILE **out_fp = (FILE**)_err_malloc(sizeof(FILE*) * out_n);
    for (i = 0; i < out_n; ++i)
        out_fp[i] = xopen(out_fn[i], "w");

    chr_name_t *cname = sg_g->cname;
    fprintf(out_fp[0], "ASM_ID\tSG_ID\tSTRAND\tCHR\tSTART_NODE\tEND_NODE\tTOTAL_NODES_NUM\tUCSC_POS\n");
    fprintf(out_fp[1], "ASM_ID\tSG_ID\tSJ_ID\tSTRAND\tCHR\tINTRON_START\tINTRON_END\tUNIQ_READ_COUNT\tMULTI_READ_COUNT\n");
    fprintf(out_fp[2], "ASM_ID\tSG_ID\tEXON_ID\tSTRAND\tCHR\tEXON_START\tEXON_END\tUNIQ_READ_COUNT\tMULTI_READ_COUNT\n");
    
    for (i = 0; i < asm_g->sg_asm_n; ++i) {
        SGasm *a = asm_g->sg_asm[i];
        int sg_i = a->SG_id; SG *sg = sg_g->SG[sg_i]; 
        SGnode *node = sg->node; SGsite *acc_site = sg->acc_site; SGsite *don_site = sg->don_site; SGedge *edge = sg->edge;
        int start, end; uint32_t v_s = a->v_start, v_e = a->v_end;

        if (node[v_s].node_e.start == 0) start = sg->start-100; else start = node[v_s].node_e.start;
        if (node[v_e].node_e.end == MAX_SITE) end = sg->end+100; else end = node[v_e].node_e.end;

        fprintf(out_fp[0], "%d\t%d\t%c\t%s\t(%d,%d)\t(%d,%d)\t%d\t%s:%d-%d\n", i, sg_i, "+-"[sg->is_rev], cname->chr_name[sg->tid], node[a->v_start].node_e.start, node[a->v_start].node_e.end, node[a->v_end].node_e.start, node[a->v_end].node_e.end, a->node_n, cname->chr_name[sg->tid], start, end);
        for (j = 0; j < a->edge_n; ++j) {
            uint32_t e_id = a->edge_id[j];
            fprintf(out_fp[1], "%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\n", i, sg_i, j, "+-"[sg->is_rev], cname->chr_name[sg->tid], don_site[edge[e_id].don_site_id].site, acc_site[edge[e_id].acc_site_id].site,  edge[e_id].uniq_c, edge[e_id].multi_c);
        }
        for (j = 0; j < a->node_n; ++j) {
            uint32_t asm_n_id = a->node_id[j];
            fprintf(out_fp[2], "%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\n", i, sg_i, j, "+-"[sg->is_rev], cname->chr_name[sg->tid], node[asm_n_id].start, node[asm_n_id].end, node[asm_n_id].uniq_c, node[asm_n_id].multi_c);
        }
    }
    
    // for SE
    /*for (i = 0; i < asm_g->sg_asm_n; ++i) {
        if (asm_g->sg_asm[i]->node_n != 1) continue;
        SGasm *sg_asm = asm_g->sg_asm[i];
        sg_i = sg_asm->SG_id; SG *sg = sg_g->SG[sg_i]; 
        SGnode *node = sg->node; SGsite *acc_site = sg->acc_site; SGsite *don_site = sg->don_site; SGedge *edge = sg->edge;
        int start, end; uint32_t v_s = sg_asm->v_start, v_e = sg_asm->v_end;

        if (node[v_s].node_e.start == 0) start = sg->start-100; else start = node[v_s].node_e.start;
        if (node[v_e].node_e.end == MAX_SITE) end = sg->end+100; else end = node[v_e].node_e.end;

        int ij_cnt, ej_cnt, e_cnt;
        uint32_t i1_id = sg_asm->edge_id[0], i2_id = sg_asm->edge_id[1], ej_id = sg_asm->edge_id[2];
        ij_cnt = edge[i1_id].uniq_c + edge[i2_id].uniq_c;
        ej_cnt = edge[ej_id].uniq_c;
        e_cnt = node[sg_asm->node_id[0]].uniq_c;

        fprintf(jcnt_out, "%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\n", i, sg_i, "+-"[sg->is_rev], cname->chr_name[sg->tid], node[sg_asm->node_id[0]].start, node[sg_asm->node_id[0]].end,  ij_cnt+ej_cnt, e_cnt);
    }*/
    // END of SE
    for (i = 0; i < out_n; ++i) {
        free(out_fn[i]); err_fclose(out_fp[i]);
    }
    free(out_fn); free(out_fp);
    return 0;
}

/*****************************/
const struct option asm_long_opt [] = {
    { "novel-sj", 0, NULL, 'n' },
    { "novel-com", 0, NULL, 'N' },
    { "sj-file", 0, NULL, 's' },
    { "use-multi", 0, NULL, 'm' },
    { "output", 1, NULL, 'o' },

    { 0, 0, 0, 0 }
};

int pred_asm(int argc, char *argv[])
{
    int c; char out_fn[1024]="";
    sg_para *sgp = sg_init_para();
	while ((c = getopt_long(argc, argv, "nNsmo:", asm_long_opt, NULL)) >= 0) {
        switch (c) {
            case 'n': sgp->no_novel_sj=0, sgp->no_novel_com=0; break;
            case 'N': sgp->no_novel_com = 0; break;
            case 's': sgp->BAM_input= 0; break;
            case 'm': sgp->use_multi = 1; break;
            case 'o': strcpy(out_fn, optarg); break;
            default: err_printf("Error: unknown option: %s.\n", optarg);
                     return pred_asm_usage();
        }
    }
    if (argc - optind != 3) return pred_asm_usage();

    gzFile genome_fp = gzopen(argv[optind], "r");
    if (genome_fp == NULL) err_fatal(__func__, "Can not open genome file. %s\n", argv[optind]);
    int seq_n = 0, seq_m; kseq_t *seq = kseq_load_genome(genome_fp, &seq_n, &seq_m);

    // parse input name
    if (sg_par_input(sgp, argv[optind+2]) <= 0) return pred_asm_usage();
    
    chr_name_t *cname = chr_name_init();
    // set cname if input is BAM
    samFile *in; bam_hdr_t *h; bam1_t *b;
    if (sgp->BAM_input) {
        if ((in = sam_open(sgp->in_name[0], "rb")) == NULL) err_fatal_core(__func__, "Cannot open \"%s\"\n", sgp->in_name[0]);
        if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", sgp->in_name[0]);
        bam_set_cname(h, cname);
        bam_hdr_destroy(h); sam_close(in);
    }
    // build splice-graph with GTF
    SG_group *sg_g;
    FILE *gtf_fp = xopen(argv[optind+1], "r");
    sg_g = construct_SpliceGraph(gtf_fp, cname);
    err_fclose(gtf_fp); chr_name_free(cname);

    int i;
    SG_group **sr_sg_g_rep = (SG_group**)_err_malloc(sgp->tol_rep_n * sizeof(SG_group*));
    SGasm_group **asm_g_rep = (SGasm_group**)_err_malloc(sgp->tol_rep_n * sizeof(SGasm_group*));
    for (i = 0; i < sgp->tol_rep_n; ++i) {
        SG_group *sr_sg_g = sr_sg_g_rep[i]; SGasm_group *asm_g = asm_g_rep[i];
        char *in_name = sgp->in_name[i];
        // get splice-junction
        int sj_n, sj_m; sj_t *sj_group;
        if (sgp->BAM_input) { // based on .bam file
            b = bam_init1(); 
            if ((in = sam_open(in_name, "rb")) == NULL) err_fatal_core(__func__, "Cannot open \"%s\"\n", in_name);
            if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", in_name);
            sj_group = (sj_t*)_err_malloc(10000 * sizeof(sj_t)); sj_m = 10000;
            // FIXME bam2itv.tmp
            sj_n = bam2sj_core(in, h, b, seq, seq_n, &sj_group, sj_m);
            bam_destroy1(b); bam_hdr_destroy(h); sam_close(in);
        } else  { // based on .sj file
            FILE *sj_fp = xopen(in_name, "r");
            sj_group = (sj_t*)_err_malloc(10000 * sizeof(sj_t)); sj_m = 10000;
            sj_n = read_sj_group(sj_fp, sg_g->cname, &sj_group, sj_m);
            err_fclose(sj_fp);
        }
        // predict splice-graph with GTF-based splice-graph and splice-junciton
        sr_sg_g = predict_SpliceGraph(*sg_g, sj_group, sj_n, sgp);

        // generate ASM with short-read splice-graph
        asm_g = gen_asm(sr_sg_g, sgp);

        // calculate number of reads falling into exon-body
        if (sgp->BAM_input) {
            samFile *in; bam_hdr_t *h; bam1_t *b;
            in = sam_open(in_name, "rb"); h = sam_hdr_read(in); b = bam_init1(); 
            cal_asm_exon_cnt(sr_sg_g, in, h, b);
            bam_destroy1(b); bam_hdr_destroy(h); sam_close(in);
        }
        asm_output(in_name, out_fn, sr_sg_g, asm_g, sgp);
    }
    sg_free_group(sg_g);
    for (i = 0; i < sgp->tol_rep_n; ++i) {
        sg_free_group(sr_sg_g_rep[i]); sg_free_asm_group(asm_g_rep[i]);
    } free(sr_sg_g_rep); free(asm_g_rep);

    // output to one file
    gzclose(genome_fp); sg_free_para(sgp);
    for (i = 0; i < seq_n; ++i) {
        free(seq[i].name.s); free(seq[i].seq.s);
    } free(seq);
    return 0;
}
