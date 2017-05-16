#ifndef _BAM_SJ_H
#define _BAM_SJ_H
#include <stdlib.h>
#include "gtf.h"
#include "kseq.h"
#include "build_sg.h"

KSEQ_INIT(gzFile, gzread)


#define bam_is_prop(b) (((b)->core.flag&BAM_FPROPER_PAIR) != 0)

kseq_t *kseq_load_genome(gzFile genome_fp, int *_seq_n, int *_seq_m);
int parse_bam_record(samFile *in, bam_hdr_t *h, bam1_t *b, kseq_t *seq, int seq_n, SG_group *sg_g, int *sg_ad_idx, ad_t **AD_group, int *AD_n, int AD_m, sj_t **SJ_group, int *SJ_n, int SJ_m, sg_para *sgp);
int bam2cnt_core(samFile *in, bam_hdr_t *h, bam1_t *b, kseq_t *seq, int seq_n, sj_t **SJ_group, int SJ_m, sg_para *sgp);
int bam2sj_core(samFile *in, bam_hdr_t *h, bam1_t *b, kseq_t *seq, int seq_n, sj_t **SJ_group, int SJ_m, sg_para *sgp);
int bam2sj(int argc, char *argv[]);
void free_sj_group(sj_t *sj_g, int sj_n);
void free_ad_group(ad_t *ad_g, int ad_n);
uint8_t bam_is_uniq_NH(bam1_t *b);

sj_t *generate_SpliceJunction(sg_para* sgp, kseq_t *seq, int seq_n, int *sj_group_n);

#endif
