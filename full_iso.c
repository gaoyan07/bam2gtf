/*
 * full_iso.c
 * use long read data to infer full-length transcript
 */
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <libgen.h>
#include <string.h>
#include <pthread.h>
#include "iso.h"
#include "utils.h"
#include "gtf.h"
#include "splice_graph.h"
#include "update_sg.h"
#include "iso.h"
#include "parse_bam.h"
#include "kstring.h"

const struct option full_iso_long_opt [] = {
    { "thread", 1, NULL, 't' },

    { "novel-sj", 0, NULL, 'n' },
    { "novel-exon", 0, NULL, 'N' },
    { "anchor-len", 1, NULL, 'a' },
    { "intron-len", 1, NULL, 'i' },
    { "genome-file", 1, NULL, 'g' },

    { "only-novel", 0, NULL, 'l' },

    { "exon-thres", 1, NULL, 'T' },
    { "asm-exon-thres", 1, NULL, 'e' },
    { "iso-cnt-thres", 1, NULL, 'C' },

    { "junc-cnt", 1, NULL, 'c' },
    { "novel-jcnt", 1, NULL, 'v' },
    { "edge-wei", 1, NULL, 'w' },

    { "output", 1, NULL, 'o' },

    { 0, 0, 0, 0 }
};

extern const char PROG[20];

uint8_t bit_table16[65536];

extern void gen_bit_table16(void);

int full_iso_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s full-iso [option] <in.gtf> <in.bam> <-o out_dir>\n\n", PROG);

    err_printf("Options:\n\n");

    err_printf("         -t --thread               number of threads to use. [1]\n");
    err_printf("\n");

    err_printf("         -n --novel-sj             allow novel splice-junction to be added in SG. [False]\n");
    err_printf("         -N --novel-exon           allow novel exon/splice-site to be added in SG. [False]\n");
    err_printf("         -a --anchor-len  [INT]    minimum anchor length for junction read. [%d].\n", ANCHOR_MIN_LEN);
    err_printf("         -i --intron-len  [INT]    minimum intron length for junction read. [%d]\n", INTRON_MIN_LEN);
    err_printf("         -g --genome-file [STR]    genome.fa. Use genome sequence to classify intron-motif. \n");
    err_printf("                                   If no genome file is give, intron-motif will be set as 0(non-canonical) [None]\n");
    err_printf("\n");

    err_printf("         -l --only-novel           only output isoform with novel-junctions. [False]\n");
    err_printf("\n");

    err_printf("         -e --iso-exon-thres [INT] maximum number of exons for SG to generate candidate isoforms. [%d]\n", ISO_EXON_MAX); 
    err_printf("         -C --iso-cnt-thres  [INT] maximum number of isoform count to retain SG. [%d]\n", ISO_CNT_MAX); 
    err_printf("\n");

    err_printf("         -c --junc-cnt    [INT]    minimum number of read count for junction. [%d]\n", JUNC_READ_CNT_MIN); 
    err_printf("         -v --novel-jcnt  [INT]    minimum number of read count for novel-junction. [%d]\n", NOVEL_JUNC_READ_CNT_MIN); 
    err_printf("         -w --edge-wei    [DOU]    during weight-based candidate isoform generation, ignore edge in\n");
    err_printf("                                   splicing-graph whose weight is less than specified value. [%.2f]\n", ISO_EDGE_MIN_WEI);
    err_printf("\n");

    err_printf("         -o --output      [STR]    directory output .IsoMatrix & .IsoExon\n");
	err_printf("\n");
	return 1;
}

typedef struct {
    int tid;
    SG_group *sg_g;
    sg_para *sgp;
    FILE **out_fp;
} full_iso_aux_t;

int REP_I;
pthread_rwlock_t RWLOCK;

FILE **full_iso_output(sg_para *sgp, char *out_dir)
{
    char gtf_suf[20] = { ".IsoGTF" };
    char *out_fn = (char*)_err_malloc(strlen(out_dir) + strlen(sgp->in_name[0]) + 30);
    strcpy(out_fn, out_dir); strcat(out_fn, "/"); strcat(out_fn, basename(sgp->in_name[0])); strcat(out_fn, gtf_suf);
    FILE **out_fp = (FILE**)_err_malloc(sizeof(FILE*));
    out_fp[0] = xopen(out_fn, "w");
    free(out_fn);
    return out_fp;
}

void filter_W(double **W, int node_n, int junc_min_cnt) {
    int i, j;
    for (i = 1; i < node_n-2; ++i) {
        for (j = i+1; j < node_n-1; ++j) {
            if (W[i][j] < junc_min_cnt) W[i][j] = 0;
        }
    }
}

void gen_full_iso(SG *sg, char **cname,  double **W, read_exon_map *M, sg_para *sgp, int *iso_i) {
    int i, j;
    uint8_t **con_matrix = (uint8_t**)_err_calloc(sg->node_n, sizeof(uint8_t*)); // [1]
    for (i = 0; i < sg->node_n; ++i) {
        con_matrix[i] = (uint8_t*)_err_calloc(sg->node_n, sizeof(uint8_t));
        for (j = 0; j < sg->node[i].next_n; ++j) {
            con_matrix[i][sg->node[i].next_id[j]] = 3;
        }
    }
    // for internal-terminal exon
    // add_pseu_wei(sg, W, con_matrix);
    int sg_node_n = sg_travl(sg, 0, sg->node_n-1);

    int entry=0, exit=sg->node_n-1;
    int iso_n=0;

    if (sg_node_n <= sgp->asm_exon_max) {
        iso_n = bias_flow_gen_full_iso(sg, cname, W, con_matrix, entry, exit, sgp);
    }
    iso_i += iso_n;
    for (i = 0; i < sg->node_n; ++i) free(con_matrix[i]); free(con_matrix);
}

int full_iso_core(SG_group *sg_g, sg_para *sgp, bam_aux_t *bam_aux)
{
    err_func_format_printf(__func__, "generate candidate isoform and read-isoform compatible matrix for each gene ...\n");
    int sg_i, iso_i; SG *sg; 
    int i, bundle_n;

    iso_i = 0;
    for (sg_i = 0; sg_i < sg_g->SG_n; ++sg_i) {
        sg = sg_g->SG[sg_i];
        // if (sg->node_n - 2 < 3) continue; // ignore simple genes XXX

        // 0. read bam bundle for each gene
        // 1. generate read-exon compatible array
        // 2. update SG with bamBundle
        // OR (optional) only known transcript
        double **W = (double**)_err_calloc(sg->node_n, sizeof(double*));
        // M will be needed in HeaviestHundling, not in enum or flow-decom
        read_exon_map *M = bam_sg_cmptb(bam_aux, W, &bundle_n, sg, sg_g->cname->chr_name, sgp); 

        filter_W(W, sg->node_n, sgp->junc_cnt_min);
        // 3. flow network decomposition OR enumerate all possible isoforms
        if (bundle_n > 0) gen_full_iso(sg, sg_g->cname->chr_name, W, M, sgp, &iso_i);
        
        // free variables
            for (i = 0; i < sg->node_n; ++i) {
                free(W[i]);
            } free(W);
            for (i = 0; i < bundle_n; ++i) {
                free(M[i].map);
            } free(M);
    }
    err_func_format_printf(__func__, "generate candidate isoform and read-isoform compatible matrix for each gene done!\n");
    return 0;
}

/*****************************/
int full_iso(int argc, char *argv[])
{
    int c, i; char ref_fn[1024]="", out_dir[1024]="", *p;
    sg_para *sgp = sg_init_para();
    while ((c = getopt_long(argc, argv, "t:nNa:i:g:le:C:c:v:w:o:", full_iso_long_opt, NULL)) >= 0) {
        switch (c) {
            case 't': sgp->n_threads = atoi(optarg); break;

            case 'n': sgp->no_novel_sj=0; break;
            case 'N': sgp->no_novel_exon=0; break;
            case 'a': sgp->anchor_len[0] = strtol(optarg, &p, 10);
                      if (*p != 0) sgp->anchor_len[1] = strtol(p+1, &p, 10); else return full_iso_usage();
                      if (*p != 0) sgp->anchor_len[2] = strtol(p+1, &p, 10); else return full_iso_usage();
                      if (*p != 0) sgp->anchor_len[3] = strtol(p+1, &p, 10); else return full_iso_usage();
                      if (*p != 0) sgp->anchor_len[4] = strtol(p+1, &p, 10); else return full_iso_usage();
                      break;
            case 'i': sgp->intron_len = atoi(optarg); break;
            case 'g': strcpy(ref_fn, optarg); break;

            case 'l': sgp->only_novel = 1, sgp->no_novel_sj=0; break; // XXX only novel


            case 'e': sgp->asm_exon_max = atoi(optarg); break;
            case 'C': sgp->iso_cnt_max = atoll(optarg); break;

            case 'c': sgp->junc_cnt_min = atoi(optarg); break;
            case 'v': sgp->novel_junc_cnt_min = atoi(optarg); break;
            case 'w': sgp->rm_edge = 1, sgp->edge_wt = atof(optarg); break;

            case 'o': strcpy(out_dir, optarg); break;
            default: err_printf("Error: unknown option: %s.\n", optarg); return full_iso_usage();
        }
    }
    if (argc - optind != 2) return full_iso_usage();
    if (strlen(out_dir) == 0) {
        err_printf("Please specify output directory with \"-o\" option.\n");
        return full_iso_usage();
    }
    sgp->read_type = SING_T;

    int seq_n = 0, seq_m; kseq_t *seq=0;
    if (strlen(ref_fn) != 0) {
        gzFile genome_fp = gzopen(ref_fn, "r");
        if (genome_fp == NULL) { err_fatal(__func__, "Can not open genome file. %s\n", ref_fn); }
        seq = kseq_load_genome(genome_fp, &seq_n, &seq_m);
        err_gzclose(genome_fp); 
    }
    // 0. parse input bam file names
    bam_aux_t **bam_aux = sg_par_input(sgp, argv[optind+1]);
    sgp->fp_n = 1;
    if (sgp->tot_rep_n <= 0) return full_iso_usage();

    // 1. set cname --- 1 thread
    chr_name_t *cname = chr_name_init();
    bam_set_cname(bam_aux[0]->h, cname);
    // 2. build splice-graph --- 1 thread (Optional in future, infer exon-intron boundaries by bam records)
    FILE *gtf_fp = xopen(argv[optind], "r");
    // build from GTF file // XXX OR from bam file
    SG_group *sg_g = construct_SpliceGraph(gtf_fp, argv[optind], cname);
    err_fclose(gtf_fp); chr_name_free(cname);

    // 3. core process
    sgp->out_fp = full_iso_output(sgp, out_dir); gen_bit_table16();
    full_iso_core(sg_g, sgp, bam_aux[0]);

    sg_free_group(sg_g);
    if (seq_n > 0) {
        for (i = 0; i < seq_n; ++i) { free(seq[i].name.s), free(seq[i].seq.s); }
        free(seq);
    }
    for (i = 0; i < sgp->tot_rep_n; ++i) bam_aux_destroy(bam_aux[i]); 
    free(bam_aux); sg_free_para(sgp);
    return 0;
}
