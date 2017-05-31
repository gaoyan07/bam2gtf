#ifndef _PRED_SG_H
#define _PRED_SG_H

#include "splice_graph.h"

int pred_sg(int argc, char *argv[]);
SG_group *predict_SpliceGraph(SG_group sg_g, sj_t *sj_group, int sj_n, sg_para *sgp);

#endif
