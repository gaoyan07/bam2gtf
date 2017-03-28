#ifndef _BAM_SJ_H
#define _BAM_SJ_H
#include <stdlib.h>
#include "htslib/htslib/sam.h"
#include "gtf.h"


uint8_t bam_is_uniq_NH(bam1_t *b);
int bam2sj_core(samFile *in, bam_hdr_t *h, bam1_t *b, gzFile genome_fp, sj_t **SJ_group, int SJ_m);
int bam2sj(int argc, char *argv[]);

#endif
