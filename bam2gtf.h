#ifndef _BAM2GTF_H
#define _BAM2GTF_H
#include "htslib/htslib/sam.h"
#include "gtf.h"

#define bam_unmap(b) ((b)->core.flag & BAM_FUNMAP)

int read_bam_trans(samFile *in, bam_hdr_t *h, bam1_t *b, int exon_min, int intron_len, read_trans_t *T);
int read_intron_group(intron_group_t *I, FILE *fp);
int read_anno_trans1(read_trans_t *T, FILE *fp);

int bam2gtf(int argc, char *argv[]);

#endif
