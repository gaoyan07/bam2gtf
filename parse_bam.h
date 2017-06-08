#ifndef _BAM_SJ_H
#define _BAM_SJ_H
#include <stdlib.h>
#include "htslib/htslib/sam.h"
#include "gtf.h"
#include "kseq.h"
#include "splice_graph.h"

KSEQ_INIT(gzFile, gzread)


#define bam_is_prop(b) (((b)->core.flag&BAM_FPROPER_PAIR) != 0)


typedef struct {
    char fn[1024];
    hts_idx_t *idx;
    samFile *in;
    bam_hdr_t *h;

    bam1_t *b;
    hts_itr_t *itr;
} bam_aux_t;

int ad_sim_comp(ad_t *ad1, ad_t *ad2);
int ad_comp(ad_t *ad1, ad_t *ad2);
ad_t *ad_init(int n);
void ad_copy(ad_t *dest, ad_t *src);
int parse_bam_record1(bam1_t *b, ad_t *ad, sg_para *sgp);

kseq_t *kseq_load_genome(gzFile genome_fp, int *_seq_n, int *_seq_m);
int parse_bam_record(samFile *in, bam_hdr_t *h, bam1_t *b, kseq_t *seq, int seq_n, SG_group *sg_g, int *sg_ad_idx, ad_t **AD_group, int *AD_n, int AD_m, sj_t **SJ_group, int *SJ_n, int SJ_m, sg_para *sgp);
int bam2cnt_core(samFile *in, bam_hdr_t *h, bam1_t *b, kseq_t *seq, int seq_n, sj_t **SJ_group, int SJ_m, sg_para *sgp);
int bam2sj_core(samFile *in, bam_hdr_t *h, bam1_t *b, kseq_t *seq, int seq_n, sj_t **SJ_group, int SJ_m, sg_para *sgp);
int bam2sj(int argc, char *argv[]);
void free_sj_group(sj_t *sj_g, int sj_n);
void free_ad_group(ad_t *ad_g, int ad_n);
uint8_t bam_is_uniq_NH(bam1_t *b);

bam_aux_t *bam_aux_init();
void bam_aux_destroy(bam_aux_t *aux);

bam_aux_t **sg_par_input(sg_para *sgp, char *in);
bam_aux_t **sg_par_input_list(sg_para *sgp, const char *list);


sj_t *generate_SpliceJunction(sg_para* sgp, kseq_t *seq, int seq_n, int *sj_group_n);

#define err_sam_open(in, fn) { if ((in = sam_open(fn, "rb")) == NULL) err_fatal(__func__, "fail to open \"%s\"\n", fn); }
#define err_sam_hdr_read(h, in, fn) { if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "fail to read header for \"%s\"\n", fn); }
#define err_sam_idx_load(idx, in, fn) { if ((idx = sam_index_load(in, fn)) == NULL) err_fatal(__func__, "fail to load the BAM index for \"%s\"\n", fn); }

#endif
