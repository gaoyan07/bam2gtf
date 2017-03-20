#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "htslib/htslib/sam.h"
#include "bam2gtf.h"
#include "utils.h"
#include "gtf.h"

extern const char PROG[20];
int bam2sj_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s bam2sj [option] <in.bam> > out.sj\n\n", PROG);
    err_printf("Note:    in.bam should be sorted in advance\n\n");
    err_printf("Options:\n\n");
    err_printf("         -g --gtf-anno    [INT]    GTF annotation file. [NULL]\n");
    //err_printf("         -s --source      [STR]    source field in GTF, program, database or project name. [NONE]\n");
	err_printf("\n");
	return 1;
}

const struct option bam2sj_long_opt [] = {
    { "gtf-anno", 1, NULL, 'g' },

    { 0, 0, 0, 0}
};

int add_sj(sj_t **sj, int *sj_n, int *sj_m, int32_t tid, int32_t don, int32_t acc, uint8_t is_rev, uint8_t is_uniq)
{
    if (*sj_n == *sj_m) _realloc(*sj, *sj_m, sj_t)
    (*sj)[*sj_n].tid = tid;
    (*sj)[*sj_n].don = don;
    (*sj)[*sj_n].acc = acc;
    (*sj)[*sj_n].strand = (is_rev==0?1:2);
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

int gen_sj(bam1_t *b, sj_t **sj, int *sj_m)
{
    if (bam_unmap(b)) return 0;
    uint32_t n_cigar, *c;
    if ((n_cigar = b->core.n_cigar) < 3) return 0;
    c = bam_get_cigar(b);

    int32_t tid = b->core.tid, end = b->core.pos;/*1-base*/
    uint8_t *p, is_rev, is_uniq; 
    p = bam_aux_get(b, "XS");
    if (p == 0) is_rev = bam_is_rev(b); else is_rev = ((bam_aux2A(p)=='+')?0 : 1);
    is_uniq = bam_is_uniq_NH(b);
    
    uint32_t i;
    int sj_n = 0;

    for (i = 0; i < n_cigar-1; ++i) {
        int l = bam_cigar_oplen(c[i]);
        switch (bam_cigar_op(c[i])) {
            case BAM_CREF_SKIP: // N(0 1)
                if (l >= INTRON_MIN_LEN) add_sj(sj, &sj_n, sj_m, tid, end+1, end+l, is_rev, is_uniq);
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
            (*SJ_group)[sj_i].uniq_c = sj[i].uniq_c; 
            (*SJ_group)[sj_i].multi_c = sj[i].multi_c;
        } else {
            (*SJ_group)[sj_i].uniq_c += sj[i].uniq_c;
            (*SJ_group)[sj_i].multi_c += sj[i].multi_c;
            if ((*SJ_group)[sj_i].strand != sj[i].strand)
                (*SJ_group)[sj_i].strand = 0;
        }
    }
    return 0;
}

int bam2sj_core(samFile *in, bam_hdr_t *h, bam1_t *b, sj_t **SJ_group, int SJ_m)
{
    err_printf("[%s] generating splice-junction with BAM file ...\n", __func__);
    int SJ_n = 0, sj_m = 1; sj_t *sj = (sj_t*)_err_malloc(sizeof(sj_t));
    while (sam_read1(in, h, b) >= 0) {
        int sj_n = gen_sj(b, &sj, &sj_m);
        if (sj_n > 0) sj_update_group(SJ_group, &SJ_n, &SJ_m, sj, sj_n);
    }
    free(sj);
    err_printf("[%s] generating splice-junction with BAM file done!\n", __func__);
    return SJ_n;
}

void print_sj(sj_t *sj_group, int sj_n, FILE *out)
{
    int i;
    for (i = 0; i < sj_n; ++i) {
        sj_t sj = sj_group[i];
        fprintf(out, "%d\t%d\t%d\t%d\t%d\t%d\n", sj.tid, sj.don, sj.acc, sj.strand, sj.uniq_c, sj.multi_c);
    }
}

int bam2sj(int argc, char *argv[])
{
    int c;
    FILE *gtf_fp=NULL;

	while ((c = getopt_long(argc, argv, "g:", bam2sj_long_opt, NULL)) >= 0) {
        switch (c) {
            case 'g': gtf_fp = xopen(optarg, "r");
            default: err_printf("Error: unknown option: %s.\n", optarg); 
                     return bam2sj_usage();
        }
    }
    if (argc - optind != 1) return bam2sj_usage();

    samFile *in; bam_hdr_t *h; bam1_t *b;

    if ((in = sam_open(argv[optind], "rb")) == NULL) err_fatal_core(__func__, "Cannot open \"%s\"\n", argv[optind]);
    if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", argv[optind]);
    b = bam_init1(); 
    sj_t *sj_group = (sj_t*)_err_malloc(10000 * sizeof(sj_t)); int sj_m = 10000;

    int sj_n = bam2sj_core(in, h, b, &sj_group, sj_m);
    bam_destroy1(b); bam_hdr_destroy(h); sam_close(in);

    print_sj(sj_group, sj_n, stdout);

    free(sj_group);
    if (gtf_fp != NULL) err_fclose(gtf_fp);
    return 0;
}
