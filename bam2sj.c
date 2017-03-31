#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "htslib/htslib/sam.h"
#include "bam2gtf.h"
#include "utils.h"
#include "gtf.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)
extern const char PROG[20];
const int intron_motif_n = 6;
const char intron_motif[6][10] = {
    "GTAG", "CTAC", 
    "GCAG", "CTGC", 
    "ATAC", "GTAT"
};
const int intron_motif_strand[6] = {
    1, 2, 1, 2, 1, 2
};

int bam2sj_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s bam2sj [option] <genome.fa> <in.bam> > out.sj\n\n", PROG);
    err_printf("Note:    in.bam should be sorted in advance\n\n");
    err_printf("Options:\n\n");
    err_printf("         -g --gtf-anno    [INT]    GTF annotation file. [NULL]\n");
    err_printf("         -m --use-multi            use both uniq- and multi-mapped reads in the bam input.[false (uniq only)]\n");
    //err_printf("         -s --source      [STR]    source field in GTF, program, database or project name. [NONE]\n");
	err_printf("\n");
	return 1;
}

const struct option bam2sj_long_opt [] = {
    { "gtf-anno", 1, NULL, 'g' },
    { "use-multi", 0, NULL, 'm' },

    { 0, 0, 0, 0}
};

int add_sj(sj_t **sj, int *sj_n, int *sj_m, int32_t tid, int32_t don, int32_t acc, uint8_t strand, uint8_t motif_i, uint8_t is_uniq)
{
    if (*sj_n == *sj_m) _realloc(*sj, *sj_m, sj_t)
    (*sj)[*sj_n].tid = tid;
    (*sj)[*sj_n].don = don;
    (*sj)[*sj_n].acc = acc;
    (*sj)[*sj_n].strand = strand;
    (*sj)[*sj_n].motif = motif_i;
    (*sj)[*sj_n].uniq_c = is_uniq; 
    (*sj)[*sj_n].multi_c = 1-is_uniq;
    (*sj_n)++;
    return 0;
}

uint8_t bam_is_uniq_NH(bam1_t *b)
{
    uint8_t *p = bam_aux_get(b, "NH");
    if (p == 0) {
        err_printf("No \"NH\" tag.\n");
        return 0;
    }
    return (bam_aux2i(p) == 1);
}

uint8_t intr_deri_str(kseq_t *seq, int seq_n, int32_t tid, int32_t start, int32_t end, uint8_t *motif_i)
{
    *motif_i = 0;
    if (tid >= seq_n) err_fatal(__func__, "unknown tid: %d\n", tid); 
    char intron[10]="";
    strncpy(intron, seq[tid].seq.s+start-1, 2);
    strncpy(intron+2, seq[tid].seq.s+end-2, 2);
    int i;
    for (i = 0; i < intron_motif_n; ++i) {
        if (strcmp(intron, intron_motif[i]) == 0) {
            *motif_i = i+1;
            return intron_motif_strand[i];
        }
    }
    return 0;
}

int gen_sj(bam1_t *b, kseq_t *seq, int seq_n, sj_t **sj, int *sj_m, int use_multi)
{
    if (bam_unmap(b)) return 0;
    uint32_t n_cigar, *c;
    if ((n_cigar = b->core.n_cigar) < 3) return 0;
    c = bam_get_cigar(b);

    int32_t tid = b->core.tid, end = b->core.pos;/*1-base*/
    uint8_t strand, motif_i, is_uniq; 
    is_uniq = bam_is_uniq_NH(b);
    if (is_uniq == 0 && use_multi == 0) return 0;
    
    uint32_t i;
    int sj_n = 0;

    for (i = 0; i < n_cigar-1; ++i) {
        int l = bam_cigar_oplen(c[i]);
        switch (bam_cigar_op(c[i])) {
            case BAM_CREF_SKIP: // N(0 1)
                if (l >= INTRON_MIN_LEN) {
                    strand = intr_deri_str(seq, seq_n, tid, end+1, end+l, &motif_i);
                    add_sj(sj, &sj_n, sj_m, tid, end+1, end+l, strand, motif_i, is_uniq);
                }
                end += l;
                break;
            case BAM_CDEL : // D(0 1)
                end += l;
                break;
            case BAM_CMATCH: // 1 1
            case BAM_CEQUAL:
            case BAM_CDIFF:
                end += l;
                break;
            case BAM_CINS: // 1 0
            case BAM_CSOFT_CLIP:
            case BAM_CHARD_CLIP:
                break;
            case BAM_CPAD: // 0 0
            case BAM_CBACK:
                break;
            default:
                err_printf("Error: unknown cigar type: %d.\n", bam_cigar_op(c[i]));
                break;
        }
    }

    return sj_n;
}

int sj_sch_group(sj_t *SJ, int SJ_n, sj_t sj, int *hit)
{
    *hit = 0;
    if (SJ_n == 0) return 0;

    int i;
    int32_t tid, don, acc;
    for (i = SJ_n-1; i >= 0; i--) {
        tid = SJ[i].tid, don = SJ[i].don, acc = SJ[i].acc;
        if (tid == sj.tid &&  don == sj.don && acc == sj.acc) { *hit = 1; return i; }
        else if (tid < sj.tid || don < sj.don || (don == sj.don && acc < sj.acc)) // SJ[i] < sj
            return i+1; 
    }
    return 0;
}

int sj_update_group(sj_t **SJ_group, int *SJ_n, int *SJ_m, sj_t *sj, int sj_n)
{
    if (sj_n + *SJ_n > *SJ_m) _realloc(*SJ_group, *SJ_m, sj_t)
    int i, hit=0;
    for (i = 0; i < sj_n; ++i) {
        int sj_i = sj_sch_group(*SJ_group, *SJ_n, sj[i], &hit);
        if (hit == 0) {
            if ((*SJ_n)++ >= *SJ_m) _realloc(*SJ_group, *SJ_m, sj_t)
            // memmove
            if (sj_i <= *SJ_n-2)
                memmove(*SJ_group+sj_i+1, *SJ_group+sj_i, (*SJ_n-sj_i-1) * sizeof(sj_t));
            // set sj
            (*SJ_group)[sj_i].tid = sj[i].tid;
            (*SJ_group)[sj_i].don = sj[i].don;
            (*SJ_group)[sj_i].acc = sj[i].acc;
            (*SJ_group)[sj_i].strand = sj[i].strand;
            (*SJ_group)[sj_i].motif = sj[i].motif;
            (*SJ_group)[sj_i].uniq_c = sj[i].uniq_c; 
            (*SJ_group)[sj_i].multi_c = sj[i].multi_c;
        } else {
            (*SJ_group)[sj_i].uniq_c += sj[i].uniq_c;
            (*SJ_group)[sj_i].multi_c += sj[i].multi_c;
            if ((*SJ_group)[sj_i].strand != sj[i].strand)
                (*SJ_group)[sj_i].strand = 0; // undifined
        }
    }
    return 0;
}

kseq_t *kseq_load_genome(gzFile genome_fp, int *_seq_n, int *_seq_m)
{
    int seq_n = 0, seq_m = 30;
    kseq_t *kseq = kseq_init(genome_fp), *seq = (kseq_t*)_err_malloc(30 * sizeof(kseq_t));

    print_format_time(stderr); err_printf("[%s] loading genome fasta file ...\n", __func__);
    while (kseq_read(kseq) >= 0) {
        kseq_copy(seq+seq_n, *kseq);
        seq_n++;
        if (seq_n == seq_m) {
            seq_m <<= 1;
            seq = (kseq_t*)_err_realloc(seq, seq_m * sizeof(kseq_t));
        }
    }
    print_format_time(stderr); err_printf("[%s] loading genome fasta file done!\n", __func__);
    kseq_destroy(kseq);
    *_seq_n = seq_n; *_seq_m = seq_m;
    return seq;
}

int bam2sj_core(samFile *in, bam_hdr_t *h, bam1_t *b, gzFile genome_fp, sj_t **SJ_group, int SJ_m, int use_multi)
{
    int seq_n = 0, seq_m; kseq_t *seq = kseq_load_genome(genome_fp, &seq_n, &seq_m);
    print_format_time(stderr); err_printf("[%s] generating splice-junction with BAM file ...\n", __func__);
    int SJ_n = 0, sj_m = 1; sj_t *sj = (sj_t*)_err_malloc(sizeof(sj_t));
    while (sam_read1(in, h, b) >= 0) {
        int sj_n = gen_sj(b, seq, seq_n, &sj, &sj_m, use_multi);
        if (sj_n > 0) sj_update_group(SJ_group, &SJ_n, &SJ_m, sj, sj_n);
    }
    free(sj);
    print_format_time(stderr); err_printf("[%s] generating splice-junction with BAM file done!\n", __func__);
    int i;
    for (i = 0; i < seq_n; ++i) {
        free(seq[i].name.s); free(seq[i].seq.s);
    } free(seq);
    return SJ_n;
}

void print_sj(sj_t *sj_group, int sj_n, FILE *out, char **cname)
{
    int i;
    fprintf(out, "###STRAND 0:undefined, 1:+, 2:-\n");
    fprintf(out, "###MOTIF 0:non-canonical, 1:GT/AG, 2:CT/AC, 3:GC/AG, 4:CT/GC, 5:AT/AC, 6:GT/AT\n");
    fprintf(out, "#CHR\tSTART\tEND\tSTRAND\tUNIQ_C\tMULTI_C\tMOTIF\n");
    for (i = 0; i < sj_n; ++i) {
        sj_t sj = sj_group[i];
        fprintf(out, "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", cname[sj.tid], sj.don, sj.acc, sj.strand, sj.uniq_c, sj.multi_c, sj.motif);
    }
}

int bam2sj(int argc, char *argv[])
{
    int c, use_multi = 0;
    FILE *gtf_fp=NULL; // TODO gtf anno

	while ((c = getopt_long(argc, argv, "g:", bam2sj_long_opt, NULL)) >= 0) {
        switch (c) {
            case 'g': gtf_fp = xopen(optarg, "r");
            case 'm': use_multi = 1;
            default: err_printf("Error: unknown option: %s.\n", optarg); 
                     return bam2sj_usage();
        }
    }
    if (argc - optind != 2) return bam2sj_usage();

    gzFile genome_fp = gzopen(argv[optind], "r");
    if (genome_fp == NULL) { err_fatal(__func__, "Can not open genome file. %s\n", argv[optind]); }

    samFile *in; bam_hdr_t *h; bam1_t *b;

    if ((in = sam_open(argv[optind+1], "rb")) == NULL) err_fatal_core(__func__, "Cannot open \"%s\"\n", argv[optind+1]);
    if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", argv[optind+1]);
    b = bam_init1(); 

    sj_t *sj_group = (sj_t*)_err_malloc(10000 * sizeof(sj_t)); int sj_m = 10000;

    int sj_n = bam2sj_core(in, h, b, genome_fp, &sj_group, sj_m, use_multi);

    print_sj(sj_group, sj_n, stdout, h->target_name);

    bam_destroy1(b); sam_close(in); bam_hdr_destroy(h); 
    free(sj_group); gzclose(genome_fp); 
    if (gtf_fp != NULL) err_fclose(gtf_fp);
    return 0;
}
