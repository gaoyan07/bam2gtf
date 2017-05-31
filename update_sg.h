#ifndef _UPDATE_SG_H
#define _UPDATE_SG_H

int update_SpliceGraph(SG_group *sg_g, sj_t *sj_group, int sj_n, sg_para *sgp);
int add_novel_sg_edge(SG *sg, gec_t *exon_id, gec_t exon_n, sg_para *sgp);

#endif
