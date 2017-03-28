#ifndef _ASM_H
#define _ASM_H
#include "build_sg.h"

void sg_free_asm_group(SGasm_group *asm_g);
SGasm_group *gen_asm(SG_group *sg_g);
int cal_asm_exon_cnt(SG_group *sg_g, samFile *in, bam_hdr_t *h, bam1_t *b);
int asm_output(char *fn, char *prefix, SG_group *sg_g, SGasm_group *asm_g);
int pred_asm(int argc, char *argv[]);

#endif
