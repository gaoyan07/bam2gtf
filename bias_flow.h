#ifndef _BIAS_FLOW_H_
#define _BIAS_FLOW_H_
#include "cand_iso.h"
#include "splice_graph.h"

cmptb_map_t **bias_flow_gen_cand_iso(SG *sg, double ***weight, int src, int sink, int rep_n, int *iso_n, sg_para *sgp);

#endif
