#include <stdio.h>
#include <stdlib.h>
#include "gtf.h"
#include "splice_graph.h"
#include "utils.h"

int add_novel_sg_edge(SG *sg, double **wei_matrix, gec_t *exon_id, gec_t exon_n, sg_para *sgp) {
    if (sgp->no_novel_sj) return 0;
    int i, hit; gec_t don_id, acc_id;
    SGnode *node = sg->node; int edge_id;
    for (i = 0; i < exon_n-1; ++i) {
        don_id = exon_id[i], acc_id = exon_id[i+1];
        edge_id =sg_bin_sch_edge(sg, don_id, acc_id, &hit);
        if (wei_matrix[don_id][acc_id] >= sgp->junc_cnt_min && hit == 0) {
            uint8_t is_anno = 0, is_rev = sg->is_rev;

            _bin_insert(acc_id, node[don_id].next_id, node[don_id].next_n, node[don_id].next_m, gec_t)
            _bin_insert(don_id, node[acc_id].pre_id, node[acc_id].pre_n, node[acc_id].pre_m, gec_t)
            sg_add_edge(sg->edge, edge_id, (sg->edge_n), (sg->edge_m), don_id, acc_id, is_rev, is_anno)
        }
    }
    return 0;
}
