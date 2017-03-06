#ifndef _PRED_SG_H
#define _PRED_SG_H

#include "build_sg.h"

SG_group *predict_SpliceGraph(SG_group sg_g, FILE *sj_p, chr_name_t *cname, int no_novel_sj);

#endif
