#ifndef _BAM_SJ_H
#define _BAM_SJ_H
#include <stdlib.h>
#include "htslib/htslib/sam.h"
#include "gtf.h"
#include "kseq.h"
#include "build_sg.h"

KSEQ_INIT(gzFile, gzread)

static inline uint8_t bam_is_uniq_NH(bam1_t *b)
{
    uint8_t *p = bam_aux_get(b, "NH");
    if (p == 0) {
        err_printf("No \"NH\" tag.\n");
        return 0;
    }
    return (bam_aux2i(p) == 1);
}

#define bam_is_prop(b) (((b)->core.flag&BAM_FPROPER_PAIR) != 0)

kseq_t *kseq_load_genome(gzFile genome_fp, int *_seq_n, int *_seq_m);
int bam2cnt_core(samFile *in, bam_hdr_t *h, bam1_t *b, kseq_t *seq, int seq_n, sj_t **SJ_group, int SJ_m, SG_group *sg_g, sg_para *sgp);
int bam2sj_core(samFile *in, bam_hdr_t *h, bam1_t *b, kseq_t *seq, int seq_n, sj_t **SJ_group, int SJ_m, sg_para *sgp);
int bam2sj(int argc, char *argv[]);

#endif
