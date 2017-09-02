#include <stdio.h>
#include <stdlib.h>
#include "gtf.h"
#include "iso.h"
#include "splice_graph.h"
#include "utils.h"
#include "kdq.h"
#include "kstring.h"

KDQ_INIT(int)
#define kdq_gec_t kdq_t(int)

void cal_flow_bias_recur(double *bias, uint8_t *node_visit, SG *sg, double **W, uint8_t **con_matrix, int cur_id, int src, int sink) {
    if (node_visit[cur_id-src] == 1) return; else node_visit[cur_id-src] = 1;
    if (cur_id == sink) return;

    int i; double wei_in=0, wei_out=0;
    for (i = 0; i < sg->node[cur_id].next_n; ++i) {
        if (is_con_matrix(con_matrix, cur_id, sg->node[cur_id].next_id[i])) {
            cal_flow_bias_recur(bias, node_visit, sg, W, con_matrix, sg->node[cur_id].next_id[i], src, sink);
            wei_out += W[cur_id][sg->node[cur_id].next_id[i]];
        }
    }
    if (cur_id == src) return;
    // cal bias
    for (i = 0; i < sg->node[cur_id].pre_n; ++i) {
        if (is_con_matrix(con_matrix, sg->node[cur_id].pre_id[i], cur_id)) {
            wei_in += W[sg->node[cur_id].pre_id[i]][cur_id];
        }
    }
    if (wei_in == 0) bias[cur_id-src] = 0;
    else bias[cur_id-src] = wei_out/wei_in;
}

// cal bias factor for each node
double *cal_flow_bias(SG *sg, double **W, uint8_t **con_matrix, int src, int sink) {
    uint8_t *node_visit = (uint8_t*)_err_calloc(sink-src+1, sizeof(uint8_t));
    double *bias = (double*)_err_calloc(sink-src+1, sizeof(double));
    cal_flow_bias_recur(bias, node_visit, sg, W, con_matrix, src, src, sink);
    free(node_visit);
    return bias;
}

int heaviest_in_edge(SG *sg, double **W, uint8_t **con_matrix, int cur_id, double min_w, int src) {
    if (cur_id == src) return -1;
    int i, m = -1;
    for (i = 0; i < sg->node[cur_id].pre_n; ++i) {
        int don_id = sg->node[cur_id].pre_id[i];
        if (is_con_matrix(con_matrix, don_id, cur_id) && W[don_id][cur_id] > min_w) {
            min_w = W[don_id][cur_id];
            m = don_id;
        }
    }
    return m;
}

int heaviest_out_edge(SG *sg, double **W, uint8_t **con_matrix, int cur_id, double min_w, int sink) {
    if (cur_id == sink) return -1;
    int i, m = -1;
    for (i = 0; i < sg->node[cur_id].next_n; ++i) {
        int acc_id = sg->node[cur_id].next_id[i];
        if (is_con_matrix(con_matrix, cur_id, acc_id) && W[cur_id][acc_id] > min_w) {
            min_w = W[cur_id][acc_id];
            m = acc_id;
        }
    }
    return m;
}

void heaviest_edge_recur(SG *sg, double **W, uint8_t **con_matrix, uint8_t *node_visit, int cur_id, int src, int sink, double *w, int *max_i, int *max_j) {
    if (node_visit[cur_id-src] == 1) return; else node_visit[cur_id-src]  = 1;
    if (cur_id == sink) return;

    int i, acc_id;
    for (i = 0; i < sg->node[cur_id].next_n; ++i) {
        acc_id = sg->node[cur_id].next_id[i];
        if (is_con_matrix(con_matrix, cur_id, acc_id)) {
            heaviest_edge_recur(sg, W, con_matrix, node_visit, acc_id, src, sink, w, max_i, max_j);
            if (W[cur_id][acc_id] > *w) {
                *w = W[cur_id][acc_id];
                *max_i = cur_id, *max_j = acc_id;
            }
        }
    }
}

double heaviest_edge(SG *sg, double **W, uint8_t **con_matrix, double min_w, int src, int sink, int *max_i, int *max_j) {
    double w = min_w;
    uint8_t *node_visit = (uint8_t*)_err_calloc(sink-src+1, sizeof(uint8_t));
    heaviest_edge_recur(sg, W, con_matrix, node_visit, src, src, sink, &w, max_i, max_j);
    free(node_visit);
    return w;
}

gec_t heaviest_path(SG *sg, double **W, uint8_t **con_matrix, double *bias, int src, int sink, gec_t *node_id, double *cap, double *bv, double min_w) {
    int max_i = src, max_j = sink, i, j; double w = min_w;
    // SG traversal
    w = heaviest_edge(sg, W, con_matrix, w, src, sink, &max_i, &max_j);

    if (w == min_w) return 0;

    gec_t l = 0;
    i = max_i;
    while (1) {
        i = heaviest_in_edge(sg, W, con_matrix, i, min_w, src);
        if (i < 0) break;
        node_id[l++] = i;
    }
    // reverse
    int tmp;
    for (i = 0; i < l/2; ++i) {
        tmp = node_id[i];
        node_id[i] = node_id[l-1-i];
        node_id[l-1-i] = tmp;
    }
    node_id[l++] = max_i;
    // check isoform integrity
    int id = node_id[0], hit = 0;
    if (id == src) {
        hit = 1;
    } else {
        //if (src != 0) 
        //    err_printf("Non-0 src\n"); // because remaining edge weight < min_w
        for (i = 0; i < sg->node[id].pre_n; ++i) {
            if (sg->node[id].pre_id[i] == src) {
                hit = 1; break;
            }
        }
    }
    if (hit == 0) return 0;

    node_id[l++] = max_j;
    j = max_j;
    while (1) {
        j = heaviest_out_edge(sg, W, con_matrix, j, min_w, sink);
        if (j < 0) break;
        node_id[l++] = j;
    }

    // check isoform integrity
    id = node_id[l-1]; hit = 0;
    if (id == sink) {
        hit = 1;
    } else {
        //if (sink != sg->node_n-1) 
        //    err_printf("Non-(n-1) sink\n"); // because remaining edge weight < min_w
        for (i = 0; i < sg->node[id].next_n; ++i) {
            if (sg->node[id].next_id[i] == sink) {
                hit = 1; break;
            }
        }
    }
    if (hit == 0) return 0;

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
    //for (i = src+1; i < sink; ++i) err_printf("%f\t", bias[i-src]);
    //err_printf("\n");
}

double bias_flow(double *cap, double *bias, int src, int sink) {
    int i; double min_f;
    min_f = cap[src-src];
    for (i = src+1; i < sink; ++i) {
        if (cap[i-src] / bias[i-src] < min_f) {
            min_f = cap[i-src] / bias[i-src];
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

int path_filter(gec_t *id, gec_t l, gec_t sg_node_n) {
    int n = 0, i;
    for (i = 0; i < l; ++i) {
        if (id[i] != 0 && id[i] != sg_node_n-1) n++;
    }
    return (n>=2 ? 1 : 0);
}

void bias_flow_iso_core(SG *sg, double **W, uint8_t **con_matrix, int src, int sink, int map_n, cmptb_map_t **iso_map, int **iso_se, int *iso_n, int iso_max, sg_para *sgp) {
    double *bias = cal_flow_bias(sg, W, con_matrix, src, sink);
    gec_t *node_id = (gec_t*)_err_malloc((sink-src+1) * sizeof(gec_t));
    double *capacities = (double*)_err_malloc((sink-src+1) * sizeof(double));
    double *bv =  (double*)_err_malloc((sink-src+1) * sizeof(double));

    gec_t l;
    //err_printf("%d %d\t%d %d\n", src, sink, sg->node[src].start, sg->node[sink].end);

    while (1) {
        if ((l = heaviest_path(sg, W, con_matrix, bias, src, sink, node_id, capacities, bv, sgp->edge_wt)) <= 0) break;

        int s = 0, t = l-1;
        //err_printf("%d\n", l);
        normalize_bias_factor(bv, s, t);
        bias_flow(capacities, bv, s, t);
        recal_flow_cap(W, capacities, node_id, l);
        // node_id, l => iso_exon_map
        if (path_filter(node_id, l, sg->node_n) == 0) continue;
        int *se = (int*)_err_malloc(2 * sizeof(int));
        cmptb_map_t *iso_m = gen_iso_exon_map(node_id, l, map_n, sg->node_n, se);
        insert_iso_exon_map(iso_map, iso_se, iso_n, map_n, iso_m, se);
        free(iso_m); free(se);
        if (*iso_n == iso_max) break; // XXX top iso_max, weight based
    }
    free(bias); free(node_id); free(capacities); free(bv);
}

int bias_flow_gen_cand_asm(SG *sg, double **rep_W, uint8_t **con_matrix, int src, int sink, cmptb_map_t **iso_map, int **iso_se, int map_n, sg_para *sgp) {
    int i, iso_max = sgp->iso_cnt_max, iso_n=0;
    // XXX add annotation iso to iso_map firstly
    for (i = 0; i < iso_max; ++i) {
        iso_map[i] = (cmptb_map_t*)_err_calloc(map_n, sizeof(cmptb_map_t));
        iso_se[i] = (int*)_err_calloc(2, sizeof(int));
    } 

    bias_flow_iso_core(sg, rep_W, con_matrix, src, sink, map_n, iso_map, iso_se, &iso_n, iso_max, sgp);

    return iso_n;
}

void enum_gen_cand_asm_core(cmptb_map_t **iso_map, int **iso_se, int *iso_n, int iso_max, int map_n, SG *sg, uint8_t **con_matrix, int cur_id, int src, int sink, uint8_t *node_visit, gec_t *path, gec_t *path_idx) {
    int next_id, i;
    SGnode *node = sg->node;

    node_visit[cur_id-src] = 1;
    path[*path_idx] = cur_id;
    (*path_idx)++;

    if (cur_id == sink) { // src-sink path
        if (path_filter(path, *path_idx, sg->node_n)) {
            if (*iso_n < iso_max) {
                int *se = (int*)_err_malloc(2 * sizeof(int));
                cmptb_map_t *iso_m = gen_iso_exon_map(path, *path_idx, map_n, sg->node_n, se);
                insert_iso_exon_map(iso_map, iso_se, iso_n, map_n, iso_m, se);
                free(iso_m); free(se);
            } else {
                (*iso_n)++;
            }
        }
    } else { // recursively all next-node
        for (i = 0; i < node[cur_id].next_n; ++i) {
            next_id = node[cur_id].next_id[i];
            if (!node_visit[next_id-src] && is_con_matrix(con_matrix, cur_id, next_id)) {
                enum_gen_cand_asm_core(iso_map, iso_se, iso_n, iso_max, map_n, sg, con_matrix, next_id, src, sink, node_visit, path, path_idx);
            }
        }
    }
    // remove current node and mark it as unvisited
    (*path_idx)--;
    node_visit[cur_id-src] = 0;
}

int enum_gen_cand_asm(SG *sg, uint8_t **con_matrix, int src, int sink, cmptb_map_t **iso_map, int **iso_se, int map_n, sg_para *sgp) {
    int i, iso_n = 0, iso_max = sgp->iso_cnt_max;
    gec_t *path = (gec_t*)_err_calloc(sink-src+1, sizeof(gec_t)); gec_t path_idx = 0;
    uint8_t *node_visit = (uint8_t*)_err_calloc(sink-src+1, sizeof(uint8_t));

    for (i = 0; i < iso_max; ++i) {
        iso_map[i] = (cmptb_map_t*)_err_calloc(map_n, sizeof(cmptb_map_t));
        iso_se[i] = (int*)_err_calloc(2, sizeof(int));
    }

    enum_gen_cand_asm_core(iso_map, iso_se, &iso_n, iso_max, map_n,  sg, con_matrix, src, src, sink, node_visit, path, &path_idx);
    if (iso_n > iso_max) iso_n = 0;
    free(node_visit); free(path);
    return iso_n;
}


int bias_flow_full_iso_core(SG *sg, char **cname, double **W, uint8_t **con_matrix, int src, int sink, int iso_max, sg_para *sgp) {
    char source[20] = "gtools";
    double *bias = cal_flow_bias(sg, W, con_matrix, src, sink);
    gec_t *node_id = (gec_t*)_err_malloc((sink-src+1) * sizeof(gec_t));
    double *capacities = (double*)_err_malloc((sink-src+1) * sizeof(double));
    double *bv =  (double*)_err_malloc((sink-src+1) * sizeof(double));

    int iso_n = 0;
    gec_t l;

    while (1) {
        if ((l = heaviest_path(sg, W, con_matrix, bias, src, sink, node_id, capacities, bv, sgp->edge_wt)) <= 0) break;
        gtf_print_trans(sgp->out_fp[0], source, sg->gene_name, sg->gene_id, cname[sg->tid], "+-"[sg->is_rev], sg, node_id, l, iso_n);

        int s = 0, t = l-1;
        normalize_bias_factor(bv, s, t);
        bias_flow(capacities, bv, s, t);
        recal_flow_cap(W, capacities, node_id, l);
        if (path_filter(node_id, l, sg->node_n) == 0) continue;
        if (++iso_n > iso_max) {
            iso_n = 0;
            break; // XXX top iso_max, weight based
        }
    }
    free(bias); free(node_id); free(capacities); free(bv);
    return iso_n;
}

int bias_flow_gen_full_iso(SG *sg, char **cname, double **W, uint8_t **con_matrix, int src, int sink, sg_para *sgp) {
    // XXX add annotation iso to iso_map firstly
    int iso_n = bias_flow_full_iso_core(sg, cname, W, con_matrix, src, sink, sgp->iso_cnt_max, sgp);

    return iso_n;
}

int heaviest_bundling_gen_asm(SG *sg, uint8_t **con_matrix, int src, int sink, cmptb_map_t **iso_map, int **iso_se, int map_n, sg_para *sgp) {
    return 0;
}


int heaviest_bundling_gen_full_iso(SG *sg, char **cname, double **W, uint8_t **con_matrix, int src, int sink, sg_para *sgp) {
    
    return 0;
}
