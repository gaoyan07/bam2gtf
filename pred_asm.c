#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include "gtf.h"
#include "build_sg.h"
#include "pred_sg.h"
#include "bam2sj.h"

extern const char PROG[20];
int pred_asm_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s asm [option] <in.gtf> <in.bam/sj>\n\n", PROG);
    err_printf("Options:\n\n");
    err_printf("         -n --novel-sj             allow novel splice-junction in the ASM. [False]\n");
    err_printf("         -s --sj-file              input with splice-junction file instead of BAM file. [false]\n");
    err_printf("                                     with splice-junction input, the .cnt output will have no count information.\n");
    err_printf("         -f --gtf-format  [STR]    file format of GTF input: \n");
    err_printf("                                     plain GTF file(GTF) or binary splice-graph file(SG). [SG]\n");
    err_printf("         -o --output      [STR]    prefix of file name of output ASM & COUNT. [in.bam/sj]\n");
    err_printf("                                     prefix.asm & prefix.cnt\n");
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
    uint32_t pre_site_id = sg_bin_sch_site(sg->don_site, sg->don_site_n, sg->node[pre_id].node_e.end+1, &hit); if (hit==0) err_fatal_core(__func__, "Can not hit site: (%d).(3)\n", sg->node[pre_id].node_e.end+1);
    uint32_t cur_site_id = sg_bin_sch_site(sg->acc_site, sg->acc_site_n, sg->node[cur_id].node_e.start-1, &hit); if (hit==0) err_fatal_simple("Can not hit sitnode_e.(3)\n");
    uint32_t edge_i = sg_bin_sch_edge(sg, pre_site_id, cur_site_id, &hit); 
    if (hit == 0) err_fatal_core(__func__, "Can not hit edgnode_e.(%d,%d) (3)\n", sg->node[pre_id].node_e.end, sg->node[cur_id].node_e.start);
    _insert(edge_i, sg_asm->edge_id, sg_asm->edge_n, sg_asm->edge_m, uint32_t)
}

void sg_update_asm(SG *sg, SGasm *sg_asm, uint32_t pre_id, uint32_t cur_id)
{
    _insert(cur_id, sg_asm->node_id, sg_asm->node_n, sg_asm->node_m, uint32_t)
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

SGasm_group *gen_asm(SG_group *sg_g)
{
    err_printf("[%s] generating alternative-splice-module with predicted SJ-SG ...\n", __func__);
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
                    sg_asm_group_add(asm_g, sg_asm);
                    sg_free_asm(sg_asm);
                    hit = 1; break;
                }
            }
            if (hit) continue;
        }
        END: free(entry); free(exit);
    }
    err_printf("[%s] generating alternative-splice-module with predicted SJ-SG done!\n", __func__);
    return asm_g;
}

int cal_asm_exon_cnt(SG_group *sg_g, samFile *in, bam_hdr_t *h, bam1_t *b)
{
    err_printf("[%s] calculating read count for AS exons ...\n", __func__);
    int i, j, last_sg_i = 0;
    while (sam_read1(in, h, b) >= 0) {
        int bam_start = b->core.pos+1;
        int bam_end = b->core.pos+bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
        /*
        for (i = last_asm_i; i < asm_g->sg_asm_n; ++i) {
            SGasm *a = asm_g->sg_asm[i];
            SG *sg = sg_g->SG[a->SG_id];
            int32_t tid = sg->tid, start = a->start, end = a->end;
            if (tid < b->core.tid || (tid == b->core.tid && end <= bam_start)) {
                if (i == last_asm_i) last_asm_i++;
                continue; 
            } else if (tid > b->core.tid || (tid == b->core.tid && start >= bam_end)) break;
            else {
                for (j = 0; j < a->node_n; ++j) {
                    int node_i = a->node_id[j];
                    if (sg->node[node_i].start <= bam_start && sg->node[node_i].end >= bam_end) 
                        sg->node[node_i].cnt++;
                }
            }
        }
        if (last_asm_i == asm_g->sg_asm_n) break;
        */

        for (i = last_sg_i; i < sg_g->SG_n; ++i) {
            SG *sg = sg_g->SG[i];
            int32_t tid = sg->tid, start = sg->start, end = sg->end;
            if (start == MAX_SITE || end == 0 || tid < b->core.tid || (tid == b->core.tid && end <= bam_start)) {
                if (i == last_sg_i) last_sg_i++; continue;
            } else if (tid > b->core.tid || (tid == b->core.tid && start >= bam_end)) break;
            else {
                for (j = 0; j < sg->node_n; ++j) {
                    if (sg->node[j].is_asm) {
                        if (sg->node[j].start <= bam_start && sg->node[j].end >= bam_end) 
                        sg->node[j].cnt++;
                    }
                }
            }
        }
        if (last_sg_i == sg_g->SG_n) break;
    }
    err_printf("[%s] calculating read count for AS exons done!\n", __func__);
    return 0;
}

/*****************************/
const struct option asm_long_opt [] = {
    { "novel-sj", 0, NULL, 'n' },
    { "gtf-format", 1, NULL, 'f' },
    { "sj-file", 0, NULL, 's' },
    { "output", 1, NULL, 'o' },

    { 0, 0, 0, 0 }
};

int pred_asm(int argc, char *argv[])
{
    int c;
    int no_novel_sj=1, SG_format=1, BAM_format=1; char out_fn[1024]="";
	while ((c = getopt_long(argc, argv, "nsf:o:", asm_long_opt, NULL)) >= 0) {
        switch (c) {
            case 'n': no_novel_sj = 0; break;
            case 'f': if (strcmp(optarg, "GTF")==0) SG_format = 0; 
                      else if (strcmp(optarg, "SG")==0) SG_format = 1; 
                      else return pred_asm_usage();
                  break;
            case 's': BAM_format = 0; break;
            case 'o': strcpy(out_fn, optarg); break;
            default: err_printf("Error: unknown option: %s.\n", optarg);
                     return pred_asm_usage();
        }
    }
    if (argc - optind != 2) return pred_asm_usage();

    // FIXME GTF.tid != BAM.tid
    // build splice-graph with GTF
    SG_group *sg_g;
    if (SG_format) sg_g = sg_restore(argv[optind]);
    else  {
        FILE *gtf_fp = xopen(argv[optind], "r");
        chr_name_t *cname = chr_name_init();
        sg_g = construct_SpliceGraph(gtf_fp, cname);
        err_fclose(gtf_fp); chr_name_free(cname);
    }

    // get splice-junction
    int sj_n, sj_m; sj_t *sj_group;
    if (BAM_format) { // based on .bam file
        samFile *in; bam_hdr_t *h; bam1_t *b;
        if ((in = sam_open(argv[optind+1], "rb")) == NULL) err_fatal_core(__func__, "Cannot open \"%s\"\n", argv[optind+1]);
        if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", argv[optind+1]);
        b = bam_init1(); 
        sj_group = (sj_t*)_err_malloc(10000 * sizeof(sj_t)); sj_m = 10000;
        sj_n = bam2sj_core(in, h, b, &sj_group, sj_m);
        bam_destroy1(b); bam_hdr_destroy(h); sam_close(in);
    } else  { // based on .sj file
        FILE *sj_fp = xopen(argv[optind+1], "r");
        sj_group = (sj_t*)_err_malloc(10000 * sizeof(sj_t)); sj_m = 10000;
        sj_n = read_sj_group(sj_fp, sg_g->cname, &sj_group, sj_m);
        err_fclose(sj_fp);
    }
    // predict splice-graph with GTF-based splice-graph and splice-junciton
    SG_group *sr_sg_g = predict_SpliceGraph(*sg_g, sj_group, sj_n, no_novel_sj);
    sg_free_group(sg_g);

    // generate ASM with short-read splice-graph
    SGasm_group *asm_g = gen_asm(sr_sg_g);

    // calculate number of reads falling into exon-body
    if (BAM_format) {
        samFile *in; bam_hdr_t *h; bam1_t *b;
        in = sam_open(argv[optind+1], "rb"); h = sam_hdr_read(in); b = bam_init1(); 
        cal_asm_exon_cnt(sr_sg_g, in, h, b);
        bam_destroy1(b); bam_hdr_destroy(h); sam_close(in);
    }

    // output
    char *asm_fn, *cnt_fn;
    if (strlen(out_fn) == 0) {
        asm_fn = (char*)_err_malloc(strlen(argv[optind+1])+10);
        cnt_fn = (char*)_err_malloc(strlen(argv[optind+1])+10);
        strcpy(asm_fn, argv[optind+1]); strcat(asm_fn, ".asm");
        strcpy(cnt_fn, argv[optind+1]); strcat(cnt_fn, ".cnt");
    } else {
        asm_fn = (char*)_err_malloc(strlen(out_fn)+10);
        cnt_fn = (char*)_err_malloc(strlen(out_fn)+10);
        strcpy(asm_fn, argv[optind+1]); strcat(asm_fn, ".asm");
        strcpy(cnt_fn, argv[optind+1]); strcat(cnt_fn, ".cnt");
    }
    FILE *asm_out = xopen(asm_fn, "w"); FILE *cnt_out = xopen(cnt_fn, "w");
    chr_name_t *cname = sr_sg_g->cname;
    fprintf(asm_out, "ASM_ID\tSG_ID\tSTRAND\tCHR\tSTART_NODE\tEND_NODE\tTOTAL_NODES_NUM\tUCSC_POS\n");
    fprintf(cnt_out, "ASM_ID\tSG_ID\tEXON_ID\tSTRAND\tCHR\tEXON_START\tEXON_END\tREAD_COUNT\n");
    int i, j, sg_i;
    for (i = 0; i < asm_g->sg_asm_n; ++i) {
        sg_i = asm_g->sg_asm[i]->SG_id;
        int start, end;
        uint32_t v_s = asm_g->sg_asm[i]->v_start, v_e = asm_g->sg_asm[i]->v_end;
        if (sr_sg_g->SG[sg_i]->node[v_s].node_e.start == 0) start = sr_sg_g->SG[sg_i]->start-100; else start = sr_sg_g->SG[sg_i]->node[v_s].node_e.start;
        if (sr_sg_g->SG[sg_i]->node[v_e].node_e.end == MAX_SITE) end = sr_sg_g->SG[sg_i]->end+100; else end = sr_sg_g->SG[sg_i]->node[v_e].node_e.end;

        fprintf(asm_out, "%d\t%d\t%c\t%s\t(%d,%d)\t(%d,%d)\t%d\t%s:%d-%d\n", i+1, sg_i, "+-"[sr_sg_g->SG[sg_i]->is_rev], cname->chr_name[sr_sg_g->SG[sg_i]->tid], sr_sg_g->SG[sg_i]->node[asm_g->sg_asm[i]->v_start].node_e.start, sr_sg_g->SG[sg_i]->node[asm_g->sg_asm[i]->v_start].node_e.end, sr_sg_g->SG[sg_i]->node[asm_g->sg_asm[i]->v_end].node_e.start,sr_sg_g->SG[sg_i]->node[asm_g->sg_asm[i]->v_end].node_e.end, asm_g->sg_asm[i]->node_n, cname->chr_name[sr_sg_g->SG[sg_i]->tid], start, end);
        for (j = 0; j < asm_g->sg_asm[i]->node_n; ++j) {
            fprintf(cnt_out, "%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\n", i+1, sg_i, j, "+-"[sr_sg_g->SG[sg_i]->is_rev], cname->chr_name[sr_sg_g->SG[sg_i]->tid], sr_sg_g->SG[sg_i]->node[asm_g->sg_asm[i]->node_id[j]].start, sr_sg_g->SG[sg_i]->node[asm_g->sg_asm[i]->node_id[j]].end, sr_sg_g->SG[sg_i]->node[asm_g->sg_asm[i]->node_id[j]].cnt);
        }
    }
    sg_free_group(sr_sg_g); sg_free_asm_group(asm_g); err_fclose(asm_out); err_fclose(cnt_out); free(asm_fn); free(cnt_fn);
    return 0;
}
