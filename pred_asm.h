#ifndef _ASM_H
#define _ASM_H
#include "build_sg.h"

#define ISO_EXON_MAX 50
#define ISO_CNT_MAX 10
#define ISO_READ_CNT_MIN 1

void sg_free_asm_group(SGasm_group *asm_g);
SGasm_group *gen_asm(SG_group *sg_g, sg_para *sgp);
int cal_asm_exon_cnt(SG_group *sg_g, samFile *in, bam_hdr_t *h, bam1_t *b);
int asm_output(char *fn, char *prefix, SG_group *sg_g, SGasm_group *asm_g, sg_para *sgp);
int pred_asm(int argc, char *argv[]);

#endif
