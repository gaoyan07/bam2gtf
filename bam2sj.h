#ifndef _BAM_SJ_H
#define _BAM_SJ_H
#include <stdlib.h>
#include "htslib/htslib/sam.h"
#include "gtf.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

uint8_t bam_is_uniq_NH(bam1_t *b);
kseq_t *kseq_load_genome(gzFile genome_fp, int *_seq_n, int *_seq_m);
int bam2sj_core(samFile *in, bam_hdr_t *h, bam1_t *b, kseq_t *seq, int seq_n, sj_t **SJ_group, int SJ_m);
int bam2sj(int argc, char *argv[]);

#endif
