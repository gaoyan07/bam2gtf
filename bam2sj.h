#ifndef _BAM_SJ_H
#define _BAM_SJ_H
#include "htslib/htslib/sam.h"
#include "gtf.h"

int bam2sj_core(samFile *in, bam_hdr_t *h, bam1_t *b, sj_t **SJ_group, int SJ_m);
int bam2sj(int argc, char *argv[]);

#endif
