#include <stdio.h>
#include <stdlib.h>
#include "gtf.h"
#include "cand_iso.h"
#include "splice_graph.h"
#include "utils.h"
#include "kdq.h"
#include "kstring.h"

void cal_flow_bias_recur(double *bias, uint8_t *node_visit, SG *sg, double **W, int cur_id, int src, int sink) {
    if (node_visit[cur_id-src] == 1) return; else node_visit[cur_id-src] = 1;
    if (cur_id == sink) return;

    int i; double wei_in=0, wei_out=0;
    for (i = 0; i < sg->node[cur_id].next_n; ++i) {
        cal_flow_bias_recur(bias, node_visit, sg, W, sg->node[cur_id].next_id[i], src, sink);
        wei_out += W[cur_id][sg->node[cur_id].next_id[i]];
    }
    if (cur_id == src) return;
    // cal bias
    for (i = 0; i < sg->node[cur_id].pre_n; ++i) {
        wei_in += W[sg->node[cur_id].pre_id[i]][cur_id];
    }
    if (wei_in == 0) bias[cur_id-src] = 0;
    else bias[cur_id-src] = wei_out/wei_in;
}

// cal bias factor for each node
double *cal_flow_bias(SG *sg, double **W, int src, int sink) {
    uint8_t *node_visit = (uint8_t*)_err_calloc(sink-src+1, sizeof(uint8_t));
    double *bias = (double*)_err_calloc(sink-src+1, sizeof(double));
    cal_flow_bias_recur(bias, node_visit, sg, W, src, src, sink);
    free(node_visit);
    return bias;
}

int heaviest_in_edge(SG *sg, double **W, int cur_id, double min_w) {
    int i, m = -1;
    for (i = 0; i < sg->node[cur_id].pre_n; ++i) {
        int don_id = sg->node[cur_id].pre_id[i];
        if (W[don_id][cur_id] > min_w) {
            min_w = W[don_id][cur_id];
            m = don_id;
        }
    }
    return m;
}

int heaviest_out_edge(SG *sg, double **W, int cur_id, double min_w) {
    int i, m = -1;
    for (i = 0; i < sg->node[cur_id].next_n; ++i) {
        int acc_id = sg->node[cur_id].next_id[i];
        if (W[cur_id][acc_id] > min_w) {
            min_w = W[cur_id][acc_id];
            m = acc_id;
        }
    }
    return m;
}

gec_t heaviest_path(SG *sg, double **W, double *bias, int src, int sink, int edge_s_id, int edge_e_id, gec_t *node_id, double *cap, double *bv, double min_w) {
    int max_i = src, max_j = sink, i, j; double w = min_w;
    SGedge *edge = sg->edge; int don_id, acc_id;
    for (i = edge_s_id; i <= edge_e_id; ++i) {
        don_id = edge[i].don_id, acc_id = edge[i].acc_id;
        if (W[don_id][acc_id] > w) {
            w = W[don_id][acc_id];
            max_i = don_id, max_j = acc_id;
        }
    }
    if (w == min_w) return 0;

    gec_t l = 0;
    i = max_i;
    while (1) {
        i = heaviest_in_edge(sg, W, i, min_w);
        if (i < 0) break;
        node_id[l++] = i;
    }
    //if (node_id[l-1] != src) err_fatal_simple("0 edge in heaviest path.(1)\n");
    // reverse
    for (i = 0; i < l/2; ++i) node_id[i] = node_id[l-1-i];
    node_id[l++] = max_i, node_id[l++] = max_j;
    j = max_j;
    while (1) {
        j = heaviest_out_edge(sg, W, j, min_w);
        if (j < 0) break;
        node_id[l++] = j;
    }
    //if (node_id[l-1] != sink) err_fatal_simple("0 edge in heaviest path.(2)\n");
    for (i = 0; i < l-1; ++i) {
        cap[i] = W[node_id[i]][node_id[i+1]];
        bv[i] = bias[node_id[i]-src];
    }

    return l;
}

void normalize_bias_factor(double *bias, int src, int sink) {
    int i; double b=1.0;
    bias[src-src] = 1.0;
    for (i = src+1; i < sink; ++i) {
        bias[i-src] = b * bias[i-src];
        b = bias[i-src];
    }
    //for (i = src+1; i < sink; ++i) printf("%f\n", bias[i-src]);
}

double bias_flow(double *cap, double *bias, int src, int sink) {
    int i, min_i; double min_f;
    min_i = src; min_f = cap[src-src];
    for (i = src+1; i < sink; ++i) {
        if (cap[i-src] / bias[i-src] < min_f) {
            min_f = cap[i-src] / bias[i-src];
            min_i = i;
        }
    }
    for (i = src; i < sink; ++i)
        cap[i-src] -= (min_f * bias[i-src]);
    return min_f;
}

void recal_flow_cap(double **W, double *cap, gec_t *node_id, int l) {
    int i;
    for (i = 0; i < l-1; ++i) {
        W[node_id[i]][node_id[i+1]] = cap[i];
    }
}

void bias_flow_iso_core(SG *sg, double **W, int src, int sink, int map_n, cmptb_map_t **iso_map, int *iso_i, int iso_max, sg_para *sgp) {
    double *bias = cal_flow_bias(sg, W, src, sink);
    gec_t *node_id = (gec_t*)_err_malloc((sink-src+1) * sizeof(gec_t));
    double *capacities = (double*)_err_malloc((sink-src+1) * sizeof(double));
    double *bv =  (double*)_err_malloc((sink-src+1) * sizeof(double));

    // cal edge s_id/e_id
    int edge_s_id = _err_sg_bin_sch_edge(sg, src, sg->node[src].next_id[0]);
    int edge_e_id = _err_sg_bin_sch_edge(sg, sg->node[sink].pre_id[sg->node[sink].pre_n-1], sink);

    gec_t l = heaviest_path(sg, W, bias, src, sink, edge_s_id, edge_e_id, node_id, capacities, bv, sgp->edge_wt);
    while (l > 0) {
        int s = 0, t = l-1;
        normalize_bias_factor(bv, s, t);
        bias_flow(capacities, bv, s, t);
        recal_flow_cap(W, capacities, node_id, l);
        // node_id, l => iso_exon_map
        cmptb_map_t *iso_m = gen_iso_exon_map(node_id, l, map_n);
        insert_iso_exon_map(iso_map, iso_i, map_n, iso_m);
        if (*iso_i == iso_max) break;
        free(iso_m);

        // next path
        l = heaviest_path(sg, W, bias, src, sink, edge_s_id, edge_e_id, node_id, capacities, bv, sgp->edge_wt);
    }
    free(bias); free(node_id); free(capacities); free(bv);
}

int bias_flow_gen_cand_iso(SG *sg, double ***W, int src, int sink, int rep_n, cmptb_map_t **iso_map, int map_n, sg_para *sgp) {
    int i, iso_max = sgp->iso_cnt_max, iso_i=0;
    for (i = 0; i < iso_max; ++i) iso_map[i] = (cmptb_map_t*)_err_calloc(map_n, sizeof(cmptb_map_t));

    for (i = 0; i < rep_n; ++i) {
        bias_flow_iso_core(sg, W[i], src, sink, map_n, iso_map, &iso_i, iso_max, sgp);
    }

    return iso_i;
}
