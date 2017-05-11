#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include "pred_asm.h"
#include "utils.h"
#include "gtf.h"
#include "build_sg.h"
#include "update_sg.h"
#include "parse_bam.h"
#include "kstring.h"
#include "kdq.h"

KDQ_INIT(int)
#define kdq_gec_t kdq_t(int)

const struct option asm_long_opt [] = {
    { "novel-sj", 0, NULL, 'n' },
    { "novel-com", 0, NULL, 'N' },
    { "proper-pair", 1, NULL, 'p' },
    { "anchor-len", 1, NULL, 'a' },
    { "intron-len", 1, NULL, 'i' },
    { "genome-file", 1, NULL, 'g' },

    { "edge-wei", 1, NULL, 'w' },
    { "only-novel", 0, NULL, 'l' },
    { "use-multi", 0, NULL, 'm' },

    { "asm-exon", 1, NULL, 'e' },
    { "iso-cnt", 1, NULL, 'C' },
    { "read-cnt", 1, NULL, 'c' },

    { "output", 1, NULL, 'o' },

    { 0, 0, 0, 0 }
};

extern const char PROG[20];
int pred_asm_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s asm [option] <in.gtf> <in.bam/sj>\n\n", PROG);
    err_printf("Note:    for multi-sample and multi-replicate should be this format: \n");
    err_printf("             \"SAM1-REP1,REP2,REP3;SAM2-REP1,REP2,REP3\"\n");
    err_printf("         use \':\' to separate samples, \',\' to separate replicates.\n\n");
    err_printf("Options:\n\n");
    err_printf("         -n --novel-sj             allow novel splice-junction in the ASM. [False]\n");
    err_printf("         -N --novel-com            allow novel combination of known exons in the ASM. [False]\n");
    err_printf("         -p --prop-pair            set -p to force to filter out reads mapped in improper pair. [False]\n");
    err_printf("         -a --anchor-len  [INT]    minimum anchor length for junction read. [%d].\n", ANCHOR_MIN_LEN);
    err_printf("         -i --intron-len  [INT]    minimum intron length for junction read. [%d]\n", INTRON_MIN_LEN);
    err_printf("         -g --genome-file [STR]    genome.fa. Use genome sequence to classify intron-motif. \n");
    err_printf("                                   If no genome file is give, intron-motif will be set as 0(non-canonical) [None]\n");
    err_printf("\n");
    err_printf("         -w --edge-wei    [INT]    remove edge in splice-graph whose weight is less than specified value. [%d]\n", 0);
    err_printf("         -l --only-novel           only output ASM/ASE with novel-junctions. [False]\n");

    err_printf("         -m --use-multi            use both uniq- and multi-mapped reads in the bam input.[False (uniq only)]\n");
    err_printf("\n");
    err_printf("         -e --iso-exon    [INT]    maximum number of exons for ASM to generate candidate isoforms. [%d]\n", ISO_EXON_MAX); 
    err_printf("         -C --iso-cnt     [INT]    maximum number of isoform count to keep ASM to candidate isoforms. [%d]\n", ISO_CNT_MAX); 
    err_printf("         -c --read-cnt    [INT]    minimum number of read count for candidate isoforms. [%d]\n", ISO_READ_CNT_MIN); 
    err_printf("         -o --output      [STR]    prefix of file name of output ASM & COUNT. [in.bam/sj]\n");
    err_printf("                                   prefix.ASM & prefix.JCNT & prefix.ECNT\n");
	err_printf("\n");
	return 1;
}

// init and free
SGasm *sg_init_asm(int sg_id, int v_start, int v_end)
{
    SGasm *sg_asm = (SGasm*)_err_malloc(sizeof(SGasm));
    sg_asm->SG_id = sg_id;
    sg_asm->v_start = v_start; sg_asm->v_end = v_end;
    sg_asm->start = MAX_SITE; sg_asm->end = 0;
    sg_asm->node_n = 0; sg_asm->node_m = 1;
    sg_asm->node_id = (gec_t*)_err_malloc(sizeof(gec_t));
    sg_asm->edge_n = 0; sg_asm->edge_m = 1;
    sg_asm->edge_id = (int*)_err_malloc(sizeof(int));
    return sg_asm;
}

SGasm_group *sg_init_asm_group(void)
{
    SGasm_group *asm_g = (SGasm_group*)_err_malloc(sizeof(SGasm_group));
    asm_g->sg_asm_n = 0; asm_g->sg_asm_m = 1;
    asm_g->sg_asm = (SGasm**)_err_malloc(sizeof(SGasm*));
    asm_g->sg_asm[0] = sg_init_asm(0, 0, 0);
    return asm_g;
}

SGasm_group *sg_realloc_asm_group(SGasm_group *asm_g)
{
    asm_g->sg_asm_m <<= 1;
    asm_g->sg_asm = (SGasm**)_err_realloc(asm_g->sg_asm, asm_g->sg_asm_m * sizeof(SGasm*));
    int i;
    for (i = asm_g->sg_asm_m >> 1; i < asm_g->sg_asm_m; ++i) asm_g->sg_asm[i] = sg_init_asm(0, 0, 0);
    return asm_g;
}

void sg_free_asm(SGasm *sg_asm)
{
    free(sg_asm->node_id); free(sg_asm->edge_id);
    free(sg_asm);
}

void sg_free_asm_group(SGasm_group *asm_g)
{
    int i; for (i = 0; i < asm_g->sg_asm_m; ++i) sg_free_asm(asm_g->sg_asm[i]);
    free(asm_g->sg_asm); free(asm_g);
}

SGiso *sg_init_iso(int asm_id, int sg_id, int v_start, int v_end)
{
    SGiso *sg_iso = (SGiso*)_err_malloc(sizeof(SGiso));
    sg_iso->ASM_id = asm_id; sg_iso->SG_id = sg_id;
    sg_iso->v_start = v_start; sg_iso->v_end = v_end;
    return sg_iso;
}

void sg_free_iso(SGiso *sg_iso)
{
    free(sg_iso->node_id); 
    free(sg_iso->uniq_sj_c); free(sg_iso->uniq_tot_c); 
#ifndef _RMATS_
    free(sg_iso->multi_sj_c); free(sg_iso->multi_tot_c);
#endif
    free(sg_iso);
}

/*****************************
 *       generate ASM        *
 *****************************/
void cal_cand_node(SG sg, int **entry, int **exit, int *entry_n, int *exit_n)
{
    int i, n1=0, n2=0;
    for (i = 0; i < sg.node_n; ++i) {
        if (sg.node[i].next_n > 1) n1++;
        if (sg.node[i].pre_n > 1) n2++;
    }
    *entry_n = n1, *exit_n = n2;
    *entry = (int*)_err_malloc(n1 * sizeof(int));
    *exit = (int*)_err_malloc(n2 * sizeof(int));

    n1 = 0, n2 = 0;
    for (i = 0; i < sg.node_n; ++i) {
        if (sg.node[i].next_n > 1) (*entry)[n1++] = sg.node[i].node_id;
        if (sg.node[i].pre_n > 1) (*exit)[n2++] = sg.node[i].node_id;;
    }
}

void sg_update_asm_edge(SG *sg, SGasm *sg_asm, int pre_id, int cur_id)
{
    if (sg->node[pre_id].node_e.end == 0 || sg->node[cur_id].node_e.start == MAX_SITE) return;
    int pre_site_id = _err_sg_bin_sch_site(sg->don_site, sg->don_site_n, sg->node[pre_id].node_e.end+1);
    int cur_site_id = _err_sg_bin_sch_site(sg->acc_site, sg->acc_site_n, sg->node[cur_id].node_e.start-1);
    int edge_i = _err_sg_bin_sch_edge(sg, pre_site_id, cur_site_id); 
    _bin_insert(edge_i, sg_asm->edge_id, sg_asm->edge_n, sg_asm->edge_m, int)
}

void sg_update_asm(SG *sg, SGasm *sg_asm, int pre_id, int cur_id)
{
    _bin_insert(cur_id, sg_asm->node_id, sg_asm->node_n, sg_asm->node_m, gec_t)
    if (sg->node[cur_id].node_e.start < sg_asm->start) sg_asm->start = sg->node[cur_id].node_e.start;
    if (sg->node[cur_id].node_e.end > sg_asm->end) sg_asm->end = sg->node[cur_id].node_e.end;
    sg->node[cur_id].is_asm = 1;
    sg_update_asm_edge(sg, sg_asm, pre_id, cur_id);
}

void sub_splice_graph(SG *sg, uint64_t **node_visit, SGasm *sg_asm, int cur_id, int e_id)
{
    if ((*node_visit)[cur_id] == 1) return; else (*node_visit)[cur_id] = 1;

    if (cur_id == e_id) return;
    int i;
    for (i = 0; i < sg->node[cur_id].next_n; ++i) {
        if (sg->node[cur_id].next_id[i] == e_id) {
            sg_update_asm_edge(sg, sg_asm, cur_id, e_id);
        } else {
            sg_update_asm(sg, sg_asm, cur_id, sg->node[cur_id].next_id[i]); 
            sub_splice_graph(sg, node_visit, sg_asm, sg->node[cur_id].next_id[i], e_id);
        }
    }
}

void sg_asm_cal_node_n(SG *sg, uint64_t *node_visit, uint64_t *iso_n, gec_t **node_n, gec_t **next_id_idx, int cur_id, int s_id, int e_id)
{
    if (node_visit[cur_id] > 0) return;
    else node_visit[cur_id] = 1;
    if (cur_id == e_id) {
        node_n[cur_id] = (gec_t*)_err_malloc(sizeof(gec_t));
        node_n[cur_id][0] = 1; 
        return;
    }

    int i, next_id;
    // calculate iso_n
    for (i = 0; i < sg->node[cur_id].next_n; ++i) {
        next_id = sg->node[cur_id].next_id[i];
        sg_asm_cal_node_n(sg, node_visit, iso_n, node_n, next_id_idx, next_id, s_id, e_id);
    }
    // calculate node_n for each iso
    // calculate next_id index for each node
    node_n[cur_id] = (gec_t*)_err_malloc(iso_n[cur_id] * sizeof(gec_t));
    next_id_idx[cur_id] = (gec_t*)_err_malloc(iso_n[cur_id] * sizeof(gec_t));
    int start_i = 0; uint64_t j, k, m=0;
    for (i = 0; i < sg->node[cur_id].next_n; ++i) {
        next_id = sg->node[cur_id].next_id[i];
        for (j = start_i, k = 0; k < iso_n[next_id]; ++j, ++k) {
            node_n[cur_id][j] = node_n[next_id][k] + 1;
            next_id_idx[cur_id][m++] = next_id;
        }
        start_i += iso_n[next_id];
    }
}

void sg_asm_cal_iso_n(SG *sg, gec_t *node_visit, gec_t *tot_exon_n, uint64_t *iso_n, gec_t cur_id, gec_t s_id, gec_t e_id)
{
    if (node_visit[cur_id] > 0) return;
    else {
        node_visit[cur_id] = 1;
        (*tot_exon_n)++;
    }
    if (cur_id == e_id) {
        iso_n[cur_id] = 1;
        return; 
    }
    iso_n[cur_id] = 0;
    int i, next_id;
    // calculate iso_n
    for (i = 0; i < sg->node[cur_id].next_n; ++i) {
        next_id = sg->node[cur_id].next_id[i];
        sg_asm_cal_iso_n(sg, node_visit, tot_exon_n, iso_n, next_id, s_id, e_id);
        iso_n[cur_id] += iso_n[next_id];
    }
}

// non-recursive
/*void sg_asm_cal_iso_n(SG *sg, int *iso_n, gec_t **id_idx, int s_id, int e_id)
{
    kdq_gec_t *node_q = kdq_init_int();
    int cur_id, pre_id, i, *v; int cur_n; 
    SGnode *node = sg->node;
    // cal iso_n for each node
    kdq_push_int(node_q, e_id);
    iso_n[e_id] = 1;
    while ((v = kdq_shift_int(node_q)) != 0) {
        cur_id = *v;
        cur_n = iso_n[cur_id];
        for (i = 0; i < node[cur_id].pre_n; ++i) {
            pre_id = node[cur_id].pre_id[i];
            if (iso_n[pre_id] == 0) {
                kdq_push_int(node_q, pre_id);
                iso_n[pre_id] = cur_n;
            } else {
                iso_n[pre_id] += cur_n;
            }
        }
    }
    kdq_destroy_int(node_q);
}*/

void sg_iso_group_init_malloc(SGiso *iso, uint64_t iso_n, gec_t *node_n)
{
    uint64_t i, tot=0;
    for (i = 0; i < iso_n; ++i) tot += node_n[i];
    err_printf("iso: %lld\ttotal exon: %lld\n", (long long)iso_n, (long long)tot);
    iso->node_id = (gec_t*)_err_malloc(tot * sizeof(gec_t));
    iso->uniq_sj_c = (int*)_err_calloc(iso_n, sizeof(int));
    iso->uniq_tot_c = (int*)_err_calloc(iso_n, sizeof(int));
#ifndef _RMATS_
    iso->multi_sj_c = (int*)_err_calloc(iso_n, sizeof(int));
    iso->multi_tot_c = (int*)_err_calloc(iso_n, sizeof(int));
#endif
}

void sg_iso_fill_col(SGiso *iso, uint64_t *count, uint64_t *iso_n, gec_t *node_n, gec_t **next_id_idx, gec_t s_id)
{
    gec_t last_id, exon_i; 
    uint64_t iso_i, id_i, last_i=0;
    for (iso_i = 0; iso_i < iso_n[s_id]; ++iso_i) {
        iso->node_id[last_i] = s_id;
        for (exon_i = 1; exon_i < node_n[iso_i]; ++exon_i) {
            id_i = exon_i + last_i;
            last_id = iso->node_id[id_i-1];
            iso->node_id[id_i] = next_id_idx[last_id][count[last_id] % iso_n[last_id]];
            ++count[last_id];
        }
        last_i += node_n[iso_i];
    }
}

int sg_asm_filt(SG *sg, gec_t *node_visit, uint64_t *iso_n, gec_t s_id, gec_t e_id, sg_para *sgp)
{
    gec_t tot_exon_n = 0;
    sg_asm_cal_iso_n(sg, node_visit, &tot_exon_n, iso_n, s_id, s_id, e_id);
    if (tot_exon_n > sgp->asm_exon_max || iso_n[s_id] > sgp->iso_cnt_max) return -1;
    else return tot_exon_n;
}

int sg_asm_gen_iso(SG *sg, SGiso *iso, uint64_t *iso_n, gec_t **node_n, gec_t s_id, gec_t e_id)
{
    // 1st round: calculate iso_n and node_n
    //            build counter index
    uint64_t *node_visit = (uint64_t*)_err_calloc(sg->node_n, sizeof(uint64_t));
    gec_t **next_id_idx = (gec_t**)_err_malloc(sg->node_n * sizeof(gec_t*));
    sg_asm_cal_node_n(sg, node_visit, iso_n, node_n, next_id_idx, s_id, s_id, e_id);
    // 1st round: calcuate node_n for each node
    // 2nd round: malloc iso[iso_n * node_n]
    sg_iso_group_init_malloc(iso, iso_n[s_id], node_n[s_id]); // 23G (229,173,153 iso)
    // 3rd round: fill node_id by column
    memset(node_visit, 0, sg->node_n * sizeof(uint64_t));
    sg_iso_fill_col(iso, node_visit, iso_n, node_n[s_id], next_id_idx, s_id);
    int i;
    for (i = 0; i < sg->node_n; ++i)
        if (node_visit[i] > 0) free(next_id_idx[i]);
    free(node_visit); free(next_id_idx);
    return 0;
}

int sg_asm_group_add(SGasm_group *asm_g, SGasm *sg_asm)
{
    if (asm_g->sg_asm_n == asm_g->sg_asm_m) sg_realloc_asm_group(asm_g);
    SGasm *a = asm_g->sg_asm[asm_g->sg_asm_n];

    a->SG_id = sg_asm->SG_id;
    a->v_start = sg_asm->v_start; a->v_end = sg_asm->v_end;
    if (sg_asm->node_n > a->node_m) a->node_id = (gec_t*)_err_realloc(a->node_id, sg_asm->node_n * sizeof(gec_t));
    if (sg_asm->edge_n > a->edge_m) a->edge_id = (int*)_err_realloc(a->edge_id, sg_asm->edge_n * sizeof(int));
    a->node_n = a->node_m = sg_asm->node_n; a->edge_n = a->edge_m = sg_asm->edge_n;
    int i;
    for (i = 0; i < sg_asm->node_n; ++i) a->node_id[i] = sg_asm->node_id[i];
    for (i = 0; i < sg_asm->edge_n; ++i) a->edge_id[i] = sg_asm->edge_id[i];
    a->start = sg_asm->start; a->end = sg_asm->end;

    asm_g->sg_asm_n++;
    return asm_g->sg_asm_n;
}

int check_novel_asm(SG *sg, SGasm *sg_asm)
{
    int i;
    for (i = 0; i < sg_asm->edge_n; ++i) {
        int eid = sg_asm->edge_id[i];
        if (sg->edge[eid].is_anno == 0) return 1;
    }
    return 0;
}

int check_uniq_asm(SG *sg, SGasm *sg_asm)
{
    int i;
    for (i = 0; i < sg_asm->edge_n; ++i) {
        int eid = sg_asm->edge_id[i];
        if (sg->edge[eid].uniq_c == 0) return 0;
    }
    return 1;
}

SGasm_group *gen_asm(SG_group *sg_g, sg_para *sgp)
{
    err_func_format_printf(__func__, "generating alternative-splice-module with predicted SJ-SG ...\n");
    int only_novel = sgp->only_novel, use_multi = sgp->use_multi;
    int entry_n, exit_n; int *entry, *exit;
    int sg_i;
    SGasm_group *asm_g = sg_init_asm_group();
    for (sg_i = 0; sg_i < sg_g->SG_n; ++sg_i) {
        SG *sg = sg_g->SG[sg_i];
        cal_cand_node(*sg, &entry, &exit, &entry_n, &exit_n);
        if (entry_n == 0 || exit_n == 0) goto END;

        int i, j, hit;
        for (i = 0; i < entry_n; ++i) {
            hit = 0;
            for (j = 0; j < exit_n; ++j) {
                int post_domn_n = sg->node[entry[i]].post_domn_n;
                int pre_domn_n = sg->node[exit[j]].pre_domn_n;
                if (post_domn_n > 1 && pre_domn_n > 1 
                        && sg->node[entry[i]].post_domn[1] == exit[j] 
                        && sg->node[exit[j]].pre_domn[1] == entry[i]) {
                    SGasm *sg_asm = sg_init_asm(sg_i, entry[i], exit[j]);
                    uint64_t *node_visit = (uint64_t*)_err_calloc(sg->node_n, sizeof(uint64_t));
                    sub_splice_graph(sg, &node_visit, sg_asm, entry[i], exit[j]);
                    free(node_visit);
                    if ((only_novel == 0 || check_novel_asm(sg, sg_asm) == 1)
                    && (use_multi == 1 || check_uniq_asm(sg, sg_asm) == 1))
                        sg_asm_group_add(asm_g, sg_asm);
                    sg_free_asm(sg_asm);
                    hit = 1; break;
                }
            }
            if (hit) continue;
        }
        END: free(entry); free(exit);
    }
    err_func_format_printf(__func__, "generating alternative-splice-module with predicted SJ-SG done!\n");
    return asm_g;
}

int comp_ad_iso(ad_t *ad, SGiso *iso, SG *sg)
{
    if (ad->tid < sg->tid) return -1;
    else if (ad->tid > sg->tid) return 1;
    else {
        // for sj-ad
        // if (ad->end <= sg->node[iso->v_start].end) return -1;
        // else if (ad->start >= sg->node[iso->v_end].start) return 1;
        // else return 0; // fully fall in OR share at one junction
        // for all-ad
        if (ad->end < sg->node[iso->v_start].start) return -1;
        else if (ad->start > sg->node[iso->v_end].end) return 1;
        else return 0; // fully fall in OR overlap 
    }
}

int comp_ad_sg(ad_t *ad, SG *sg)
{
    if (ad->tid < sg->tid) return -1;
    else if (ad->tid > sg->tid) return 1;
    else {
        if (ad->end < sg->start) return -1;
        else if (ad->start > sg->end) return 1;
        else return 0; // fully fall in OR share at one junction
    }
}

// return value:
// 0: not consistent
//  fully fall in iso: overlap with intronic region
//  NOT fully fall in iso: no shared junction
int ad_bin_sch_node(ad_t *ad, int *ad_i, SGnode *node, gec_t *node_id, int node_n, int *node_i)
{
    if (ad->start >= node[node_id[0]].start) {
        *ad_i = 0;
        //int end = ad->exon_end[0];
        int pos = ad->start;
        int left = 0, right = node_n-1, mid;
        int mid_s, mid_e;
        if (right == -1) return 0;
        while (left <= right) {
            mid = (left + right) >> 1;
            mid_s = node[node_id[mid]].start;
            mid_e = node[node_id[mid]].end;
            if (pos >= mid_s && pos <= mid_e) {
                *node_i = mid;
                return 1;
            } else if (pos < mid_s) 
                right = mid-1;
            else // pos > mid_e
                left = mid+1;
        }
    } else {
        *node_i = 0;
        // int end = node[node_id[0]].end;
        int pos = node[node_id[0]].start;
        int left = 0, right = ad->intv_n-1, mid;
        int mid_s, mid_e;
        if (right == -1) return 0;
        while (left <= right) {
            mid = (left + right) >> 1;
            mid_s = (mid == 0 ? ad->start : ad->intr_end[mid-1]+1);
            mid_e = ad->exon_end[mid];
            if (pos >= mid_s && pos <= mid_e) {
                *ad_i = mid;
                return 1;
            } else if (pos < mid_s) 
                right = mid-1;
            else // end > mid_e
                left = mid+1;
        } 
    }
    return 0;
}

int check_consis_ad_core(ad_t *ad, int _ad_i, SGnode *node, gec_t *node_id, int node_n, int node_i)
{
    int ad_i, n_i;
    int ad_pos; SGnode n;
    for (ad_i = _ad_i, n_i = node_i; ad_i < ad->intv_n && n_i < node_n; ++ad_i, ++n_i) {
        n = node[node_id[n_i]];
//#ifdef _RMATS_
//        if (n_i != 0) {
//#endif
        if (ad_i != 0) { // check exon start
            ad_pos = ad->intr_end[ad_i-1] + 1;
            if (ad_pos != n.start) return 0;
        } else {
            ad_pos = ad->start;
            if (ad_pos < n.start) return 0; 
        }
//#ifdef _RMATS_
//        }
//        if (n_i != node_n-1) {
//#endif
        ad_pos = ad->exon_end[ad_i];
        if (ad_i != ad->intv_n-1) { // check exon end
            if (ad_pos != n.end) return 0;
        } else {
            if (ad_pos > n.end) return 0;
        }
//#ifdef _RMATS_
//        }
//#endif
    }
    return 1;
}

int check_consis_ad(ad_t *ad, SGnode *node, gec_t *node_id, gec_t node_n)
{
    int ad_i, node_i;
    if (ad_bin_sch_node(ad, &ad_i, node, node_id, node_n, &node_i) == 0) return 0;
    return check_consis_ad_core(ad, ad_i, node, node_id, node_n, node_i);
}

int sg_iso_ad_cnt(SG *sg, ad_t *ad_g, int ad_n, int ad_i, SGiso *iso, uint64_t iso_n, gec_t *node_n, uint8_t **iso_ad_map, uint64_t iso_ad_m, sg_para *sgp, uint64_t *iso_filt_n)
{
    *iso_filt_n = 0;
    if (ad_i == 0) return 0;
    else ad_i--;
    uint64_t iso_ad_n = 0, last_iso_ad_map_i = 0, iso_i;
    SGnode *node = sg->node;
    do {
        ad_t *ad = ad_g+ad_i;
        if (ad->tid > sg->tid || (ad->tid == sg->tid && ad->start >= node[iso->v_end].end)) break;
        else if (ad->tid < sg->tid) { ad_i++; continue; }

        int comp_res = comp_ad_iso(ad, iso, sg);
        if (comp_res > 0) break;
        else if (comp_res == 0) { //
            uint64_t last_i=0, hit=0;
            for (iso_i = 0; iso_i < iso_n; ++iso_i) {
                // check if ad is consistent with iso->node
                if (check_consis_ad(ad, node, iso->node_id+last_i, node_n[iso_i])) {
                    hit = 1;
                    if (iso_ad_n == iso_ad_m) {
                        iso_ad_m <<= 1;
                        *iso_ad_map = (uint8_t*)_err_realloc(*iso_ad_map, iso_ad_m * iso_n * sizeof(uint8_t));
                        memset((*iso_ad_map)+iso_ad_n*iso_n, 0, iso_ad_n * iso_n * sizeof(uint8_t));
                    }
                    (*iso_ad_map)[last_iso_ad_map_i + iso_i] = 1;
#ifndef _RMATS_
                    if (ad->is_uniq) {
#endif
                        if (ad->is_splice) iso->uniq_sj_c[iso_i]++;
                        iso->uniq_tot_c[iso_i]++;
#ifndef _RMATS_
                    } else {
                        if (ad->is_splice) iso->multi_sj_c[iso_i]++;
                        iso->multi_tot_c[iso_i]++;
                    }
#endif
                }
                last_i += node_n[iso_i];
            }
            if (hit) {
                iso_ad_n++;
                last_iso_ad_map_i += iso_n;
            }
        }
        ad_i++;
    } while (ad_i < ad_n);
    for (iso_i = 0; iso_i < iso_n; ++iso_i) if (iso->uniq_tot_c[iso_i] >= sgp->iso_read_cnt_min) ++(*iso_filt_n);
    return iso_ad_n;
}

void sg_iso_eb_cnt(SG *sg, SGiso *iso, uint64_t iso_n, gec_t *node_n)
{
    SGnode *node = sg->node;
    gec_t n_i; uint64_t iso_i, last_i = 0;
    for (iso_i = 0; iso_i < iso_n; ++iso_i) {
        for (n_i = 0; n_i < node_n[iso_i]; ++n_i) {
            iso->uniq_tot_c[iso_i] += node[iso->node_id[last_i+n_i]].uniq_c;
#ifndef _RMATS_
            iso->multi_tot_c[iso_i] += node[iso->node_id[last_i+n_i]].multi_c;
#endif
        }
        last_i += node_n[iso_i];
    }
}

void sg_per_iso_output(FILE **out_fp, SG *sg, chr_name_t *cname, SGiso *iso, uint64_t iso_n, uint64_t iso_filt_n, gec_t *node_n, uint8_t *iso_ad_map, uint64_t iso_ad_n, gec_t *node_visit, int tot_exon_n, sg_para *sgp)
{
    int asm_i = iso->ASM_id; uint64_t iso_i, ISO_i, last_i = 0;
    SGnode *node = sg->node;
    /*
    int vs = iso->v_start, ve = iso->v_end;
    int sg_i = sg->SG_id;


    for (iso_i = 0; iso_i < iso_n; ++iso_i) {
        if (iso->uniq_tot_c[iso_i] < sgp->iso_read_cnt_min) continue;
        fprintf(out_fp[0], "%lld\t%d\t%d\t%c\t%s\t%d\t%d,%d\t%d,%d\t%d\t%d\n", (long long)iso_i, asm_i, sg_i, "+-"[sg->is_rev], cname->chr_name[sg->tid], node_n[iso_i], node[vs].start, node[vs].end, node[ve].start, node[ve].end, iso->uniq_sj_c[iso_i], iso->uniq_tot_c[iso_i]);
        fprintf(out_fp[1], "%lld\t%d\t%d\t%c\t%s\t%d", (long long)iso_i, asm_i, sg_i, "+-"[sg->is_rev], cname->chr_name[sg->tid], node_n[iso_i]);
        int k;
        for (k = 0; k < node_n[iso_i]; ++k) {
            fprintf(out_fp[1], "\t%d,%d", node[iso->node_id[last_i+k]].start,node[iso->node_id[last_i+k]].end);
        }
        fprintf(out_fp[1], "\n");
        last_i += node_n[iso_i];
    }*/
    // iso-read map
    // ASM/ISO header
    fprintf(out_fp[0], "ASM#%d\t%lld\t%lld\n", asm_i, (long long)iso_ad_n, (long long)iso_filt_n);
    //fprintf(out_fp[0], "ISO-READ-MAP");
    //for (iso_i=0, ISO_i=0; iso_i < iso_n; ++iso_i, ++ISO_i) {
    //    if (iso->uniq_tot_c[iso_i] < sgp->iso_read_cnt_min) continue;
    //    fprintf(out_fp[0], "\tISO#%lld", (long long)ISO_i);
    //}
    //fprintf(out_fp[0], "\n");
    uint64_t iso_ad_i, last_iso_ad_map_i = 0;
    for (iso_ad_i = 0; iso_ad_i < iso_ad_n; ++iso_ad_i) {
        //fprintf(out_fp[0], "READ#%lld\t", (long long)iso_ad_i);
        for (iso_i=0, ISO_i=0; iso_i < iso_n; ++iso_i) {
            if (iso->uniq_tot_c[iso_i] < sgp->iso_read_cnt_min) continue;
            if (ISO_i == 0) fprintf(out_fp[0], "%d", iso_ad_map[last_iso_ad_map_i+iso_i]);
            else  fprintf(out_fp[0], "\t%d", iso_ad_map[last_iso_ad_map_i+iso_i]);
            ISO_i++;
        }
        fprintf(out_fp[0], "\n");
        last_iso_ad_map_i += iso_n;
    }
    
    // remove virtual start/end nodes
    int v_n=0;
    int tmp1 = node_visit[0], tmp2 = node_visit[sg->node_n-1];
    if (node_visit[0] == 1) {
        node_visit[0] = 0;
        v_n++;
    }
    if (node_visit[sg->node_n-1] == 1) {
        node_visit[sg->node_n-1] = 0;
        v_n++;
    }

    // iso exon coordinates
    // ASM exon heder
    fprintf(out_fp[1], "ASM#%d\t%d\t%lld\t%c\t%s\n", asm_i, tot_exon_n-v_n, (long long)ISO_i, "+-"[sg->is_rev], cname->chr_name[sg->tid]);
    int i; gec_t exon_i = 0;
    //fprintf(out_fp[1], "EXON\t");
    for (i = 0; i < sg->node_n; ++i) {
        if (node_visit[i]) { 
            if (exon_i == 0) fprintf(out_fp[1], "%d,%d", node[i].start, node[i].end);
            else fprintf(out_fp[1], "\t%d,%d", node[i].start, node[i].end);
            node_visit[i] = ++exon_i;
        }
    } fprintf(out_fp[1], "\n");
    for (iso_i = 0; iso_i < iso_n; ++iso_i) {
        //fprintf(out_fp[1], "ISO#%lld\t%d", (long long)iso_i, node_n[iso_i]-v_n);
        if (iso->uniq_tot_c[iso_i] >= sgp->iso_read_cnt_min) {
            fprintf(out_fp[1], "%d", node_n[iso_i]-v_n);
            for (i = 0; i < node_n[iso_i]; ++i) {
                if (node_visit[iso->node_id[last_i+i]] != 0)
                    fprintf(out_fp[1], "\t%d", node_visit[iso->node_id[last_i+i]]-1);
            }
            fprintf(out_fp[1], "\n");
        }
        last_i += node_n[iso_i];
    }
    node_visit[0] = tmp1, node_visit[sg->node_n-1] = tmp2;
}

int gen_asm_iso(SG_group *sg_g, int *sg_ad_idx, ad_t *ad_group, int ad_n, sg_para *sgp, FILE **out_fp, int f_n)
{
    err_func_format_printf(__func__, "generating candidate isoforms of alternative-splice-module ...\n");
    int entry_n, exit_n; int *entry, *exit;
    int sg_i, asm_i=0;

    for (sg_i = 0; sg_i < sg_g->SG_n; ++sg_i) {
        SG *sg = sg_g->SG[sg_i];
        gec_t *node_visit = (gec_t*)_err_malloc(sg->node_n * sizeof(gec_t));
        uint64_t *iso_n = (uint64_t*)_err_malloc(sg->node_n * sizeof(uint64_t));
        gec_t **node_n = (gec_t**)_err_malloc(sg->node_n * sizeof(gec_t*));
        cal_cand_node(*sg, &entry, &exit, &entry_n, &exit_n);
        if (entry_n == 0 || exit_n == 0) goto END;

        int i, j, hit;
        for (i = 0; i < entry_n; ++i) {
            hit = 0;
            for (j = 0; j < exit_n; ++j) {
                int post_domn_n = sg->node[entry[i]].post_domn_n;
                int pre_domn_n = sg->node[exit[j]].pre_domn_n;
                if (post_domn_n > 1 && pre_domn_n > 1 
                        && sg->node[entry[i]].post_domn[1] == exit[j] 
                        && sg->node[exit[j]].pre_domn[1] == entry[i]) {
                    memset(node_visit, 0, sg->node_n * sizeof(gec_t));
                    int tot_exon_n;
                    if ((tot_exon_n = sg_asm_filt(sg, node_visit, iso_n, entry[i], exit[j], sgp)) < 0) continue;

                    err_func_format_printf(__func__, "SG: %d\tASM: %d\tVS: %d\tVE: %d\t", sg_i, asm_i, entry[i], exit[j]);
                    // generate isoforms
                    SGiso *iso = sg_init_iso(asm_i, sg_i, entry[i], exit[j]);
                    sg_asm_gen_iso(sg, iso, iso_n, node_n, entry[i], exit[j]);
                    // calculate iso cnt with sg_ad_idx
                    uint64_t iso_ad_m = 100000; uint8_t *iso_ad_map = (uint8_t*)_err_calloc(iso_ad_m * iso_n[entry[i]], sizeof(uint8_t));
                    uint64_t iso_filt_n, iso_ad_n;
                    iso_ad_n = sg_iso_ad_cnt(sg, ad_group, ad_n, sg_ad_idx[sg_i], iso, iso_n[entry[i]], node_n[entry[i]], &iso_ad_map, iso_ad_m, sgp, &iso_filt_n);
                    //sg_iso_eb_cnt(sg, iso, iso_n[entry[i]], node_n[entry[i]]); // exon-body-cnt
                    // output read count and isoform exons
                    // XXX check novel, check uniq
                    if (iso_ad_n > 0 && iso_filt_n > 0) {
                        sg_per_iso_output(out_fp, sg, sg_g->cname, iso, iso_n[entry[i]], iso_filt_n, node_n[entry[i]], iso_ad_map, iso_ad_n, node_visit, tot_exon_n, sgp);
                        asm_i++;
                    }
                    free(iso_ad_map);
                    sg_free_iso(iso);
                    int n; for (n = 0; n < sg->node_n; ++n) {
                        if (node_visit[n]) free(node_n[n]);
                    }

                    hit = 1; break;
                }
            }
            if (hit) continue;
        }
END: free(entry); free(exit); free(node_visit); free(iso_n); free(node_n);
    }
    int i;
    for (i = 0; i < f_n; ++i) err_fclose(out_fp[i]);
    free(out_fp);
    err_func_format_printf(__func__, "generating candidate isoforms of alternative-splice-module done!\n");
    return asm_i;
}

int cal_asm_exon_cnt(SG_group *sg_g, samFile *in, bam_hdr_t *h, bam1_t *b)
{
    err_func_format_printf(__func__, "calculating read count for AS exons ...\n");
    int i, j, last_sg_i = 0;
    while (sam_read1(in, h, b) >= 0) {
        int bam_start = b->core.pos+1;
        int bam_end = b->core.pos+bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
        uint8_t is_uniq = bam_is_uniq_NH(b);

        for (i = last_sg_i; i < sg_g->SG_n; ++i) {
            SG *sg = sg_g->SG[i];
            int tid = sg->tid, start = sg->start, end = sg->end;
            if (start == MAX_SITE || end == 0 || tid < b->core.tid || (tid == b->core.tid && end <= bam_start)) {
                if (i == last_sg_i) last_sg_i++; continue;
            } else if (tid > b->core.tid || (tid == b->core.tid && start >= bam_end)) break;
            else {
                for (j = 0; j < sg->node_n; ++j) {
                    if (sg->node[j].is_asm) {
                        if (sg->node[j].start <= bam_start && sg->node[j].end >= bam_end) {
                            if (is_uniq) sg->node[j].uniq_c++;
                            else sg->node[j].multi_c++;
                        }
                    }
                }
            }
        }
        if (last_sg_i == sg_g->SG_n) break;
    }
    err_func_format_printf(__func__, "calculating read count for AS exons done!\n");
    return 0;
}

int asm_output(char *in_fn, char *prefix, SG_group *sg_g, SGasm_group *asm_g, sg_para *sgp)
{
    int i, j, out_n=3;
    char suf[3][10] = { ".ASM", ".JCNT", ".ECNT" };
    char suff[20] = "";
    if (sgp->use_multi==1) strcat(suff, ".multi");
    if (sgp->no_novel_sj==1) strcat(suff, ".anno");
    if (sgp->only_novel==1) strcat(suff, ".novel");
    char **out_fn = (char**)_err_malloc(sizeof(char*) * out_n);
    if (strlen(prefix) == 0) {
        for (i = 0; i < out_n; ++i) {
            out_fn[i] = (char*)_err_malloc(strlen(in_fn)+30); strcpy(out_fn[i], in_fn); strcat(out_fn[i], suff); strcat(out_fn[i], suf[i]);
        }
    } else {
        for (i = 0; i < out_n; ++i) {
            out_fn[i] = (char*)_err_malloc(strlen(prefix)+30); strcpy(out_fn[i], prefix); strcat(out_fn[i], suff); strcat(out_fn[i], suf[i]);
        }
    }

    FILE **out_fp = (FILE**)_err_malloc(sizeof(FILE*) * out_n);
    for (i = 0; i < out_n; ++i)
        out_fp[i] = xopen(out_fn[i], "w");

    chr_name_t *cname = sg_g->cname;
    fprintf(out_fp[0], "ASM_ID\tSG_ID\tSTRAND\tCHR\tSTART_NODE\tEND_NODE\tTOTAL_NODES_NUM\tUCSC_POS\n");
    fprintf(out_fp[1], "ASM_ID\tSG_ID\tSJ_ID\tSTRAND\tCHR\tINTRON_START\tINTRON_END\tUNIQ_READ_COUNT\n");
    fprintf(out_fp[2], "ASM_ID\tSG_ID\tEXON_ID\tSTRAND\tCHR\tEXON_START\tEXON_END\tUNIQ_READ_COUNT\n");
    
    for (i = 0; i < asm_g->sg_asm_n; ++i) {
        SGasm *a = asm_g->sg_asm[i];
        int sg_i = a->SG_id; SG *sg = sg_g->SG[sg_i]; 
        SGnode *node = sg->node; SGsite *acc_site = sg->acc_site; SGsite *don_site = sg->don_site; SGedge *edge = sg->edge;
        int start, end; int v_s = a->v_start, v_e = a->v_end;

        if (node[v_s].node_e.start == 0) start = sg->start-100; else start = node[v_s].node_e.start;
        if (node[v_e].node_e.end == MAX_SITE) end = sg->end+100; else end = node[v_e].node_e.end;

        fprintf(out_fp[0], "%d\t%d\t%c\t%s\t(%d,%d)\t(%d,%d)\t%d\t%s:%d-%d\n", i, sg_i, "+-"[sg->is_rev], cname->chr_name[sg->tid], node[a->v_start].node_e.start, node[a->v_start].node_e.end, node[a->v_end].node_e.start, node[a->v_end].node_e.end, a->node_n, cname->chr_name[sg->tid], start, end);
        for (j = 0; j < a->edge_n; ++j) {
            int e_id = a->edge_id[j];
            fprintf(out_fp[1], "%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\n", i, sg_i, j, "+-"[sg->is_rev], cname->chr_name[sg->tid], don_site[edge[e_id].don_site_id].site, acc_site[edge[e_id].acc_site_id].site,  edge[e_id].uniq_c);
        }
        for (j = 0; j < a->node_n; ++j) {
            int asm_n_id = a->node_id[j];
            fprintf(out_fp[2], "%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\n", i, sg_i, j, "+-"[sg->is_rev], cname->chr_name[sg->tid], node[asm_n_id].start, node[asm_n_id].end, node[asm_n_id].uniq_c);
        }
    }
    
    for (i = 0; i < out_n; ++i) {
        free(out_fn[i]); err_fclose(out_fp[i]);
    }
    free(out_fn); free(out_fp);
    return 0;
}

FILE **iso_header_output(char *in_fn, char *prefix, int *fn)
{
    int i, out_n=2;
    char suf[2][20] = { ".IsoMatrix" , ".IsoExon" };
    char suff[20] = "";
    //if (sgp->use_multi==1) strcat(suff, ".multi");
    //if (sgp->no_novel_sj==1) strcat(suff, ".anno");
    //if (sgp->only_novel==1) strcat(suff, ".novel");
    char **out_fn = (char**)_err_malloc(sizeof(char*) * out_n);
    if (strlen(prefix) == 0) {
        for (i = 0; i < out_n; ++i) {
            out_fn[i] = (char*)_err_malloc(strlen(in_fn)+30); strcpy(out_fn[i], in_fn); strcat(out_fn[i], suff); strcat(out_fn[i], suf[i]);
        }
    } else {
        for (i = 0; i < out_n; ++i) {
            out_fn[i] = (char*)_err_malloc(strlen(prefix)+30); strcpy(out_fn[i], prefix); strcat(out_fn[i], suff); strcat(out_fn[i], suf[i]);
        }
    }

    FILE **out_fp = (FILE**)_err_malloc(sizeof(FILE*) * out_n);
    for (i = 0; i < out_n; ++i) out_fp[i] = xopen(out_fn[i], "w");
    //fprintf(out_fp[0], "ISO_ID\tASM_ID\tSG_ID\tSTRAND\tCHR\tEXON_CNT\tSTART_EXON\tEND_EXON\tUNIQ_SJ_CNT\tUNIQ_TOT_COUNT\n");
    //fprintf(out_fp[1], "ISO_ID\tASM_ID\tSG_ID\tSTRAND\tCHR\tEXON_CNT\tISO_EXON\n");
    *fn = out_n;
    for (i = 0; i < out_n; ++i) free(out_fn[i]); free(out_fn);
    return out_fp;
}

/*****************************/
int pred_asm(int argc, char *argv[])
{
    int c, i; char out_fn[1024]="", ref_fn[1024]="", *p;
    sg_para *sgp = sg_init_para();
	while ((c = getopt_long(argc, argv, "nNpa:i:g:w:lme:c:C:o:M:U:A:", asm_long_opt, NULL)) >= 0) {
        switch (c) {
            case 'n': sgp->no_novel_sj=0, sgp->no_novel_com=0; break;
            case 'N': sgp->no_novel_com = 0; break;
            case 'l': sgp->only_novel = 1, sgp->no_novel_sj=0, sgp->no_novel_com=0; break;
            case 'm': sgp->use_multi = 1; break;
            case 'M': sgp->merge_out = 1; break;
            case 'e': sgp->asm_exon_max = atoi(optarg); break;
            case 'C': sgp->iso_cnt_max = atoll(optarg); break;
            case 'c': sgp->iso_read_cnt_min = atoi(optarg); break;
            case 'p': sgp->read_type = PAIR_T; break;
            case 'a': sgp->anchor_len[0] = strtol(optarg, &p, 10);
                      if (*p != 0) sgp->anchor_len[1] = strtol(p+1, &p, 10); else return pred_asm_usage();
                      if (*p != 0) sgp->anchor_len[2] = strtol(p+1, &p, 10); else return pred_asm_usage();
                      if (*p != 0) sgp->anchor_len[3] = strtol(p+1, &p, 10); else return pred_asm_usage();
                      if (*p != 0) sgp->anchor_len[4] = strtol(p+1, &p, 10); else return pred_asm_usage();
                      break;
            case 'U': sgp->uniq_min[0] = strtol(optarg, &p, 10);
                      if (*p != 0) sgp->uniq_min[1] = strtol(p+1, &p, 10); else return pred_asm_usage();
                      if (*p != 0) sgp->uniq_min[2] = strtol(p+1, &p, 10); else return pred_asm_usage();
                      if (*p != 0) sgp->uniq_min[3] = strtol(p+1, &p, 10); else return pred_asm_usage();
                      if (*p != 0) sgp->uniq_min[4] = strtol(p+1, &p, 10); else return pred_asm_usage();
                      break; 
            case 'A': sgp->all_min[0] = strtol(optarg, &p, 10);
                      if (*p != 0) sgp->all_min[1] = strtol(p+1, &p, 10); else return pred_asm_usage();
                      if (*p != 0) sgp->all_min[2] = strtol(p+1, &p, 10); else return pred_asm_usage();
                      if (*p != 0) sgp->all_min[3] = strtol(p+1, &p, 10); else return pred_asm_usage();
                      if (*p != 0) sgp->all_min[4] = strtol(p+1, &p, 10); else return pred_asm_usage();
                      break;
            case 'i': sgp->intron_len = atoi(optarg); break;
            case 'w': sgp->rm_edge = 1, sgp->edge_wt = atoi(optarg); break;
            case 'g': strcpy(ref_fn, optarg); break;
            case 'o': strcpy(out_fn, optarg); break;
            default: err_printf("Error: unknown option: %s.\n", optarg); return pred_asm_usage();
        }
    }
    if (argc - optind != 2) return pred_asm_usage();

    int seq_n = 0, seq_m; kseq_t *seq;
    if (strlen(ref_fn) != 0) {
        gzFile genome_fp = gzopen(ref_fn, "r");
        if (genome_fp == NULL) { err_fatal(__func__, "Can not open genome file. %s\n", ref_fn); }
        seq = kseq_load_genome(genome_fp, &seq_n, &seq_m);
        err_gzclose(genome_fp); 
    }
    // parse input name
    if (sg_par_input(sgp, argv[optind+1]) <= 0) return pred_asm_usage();
    
    // set cname
    chr_name_t *cname = chr_name_init();
    samFile *in; bam_hdr_t *h; bam1_t *b;
    if ((in = sam_open(sgp->in_name[0], "rb")) == NULL) err_fatal_core(__func__, "Cannot open \"%s\"\n", sgp->in_name[0]);
    if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", sgp->in_name[0]);
    bam_set_cname(h, cname);
    bam_hdr_destroy(h); sam_close(in);
    // build splice-graph
    FILE *gtf_fp = xopen(argv[optind], "r");
    SG_group *sg_g = construct_SpliceGraph(gtf_fp, cname);
    err_fclose(gtf_fp); chr_name_free(cname);

    SGasm_group **asm_g_rep = (SGasm_group**)_err_malloc(sgp->tol_rep_n * sizeof(SGasm_group*));
    for (i = 0; i < sgp->tol_rep_n; ++i) {
        char *in_name = sgp->in_name[i];
        // get splice-junction
        int sj_n, sj_m; sj_t *sj_group; int ad_n, ad_m; ad_t *ad_group;
        b = bam_init1(); 
        if ((in = sam_open(in_name, "rb")) == NULL) err_fatal_core(__func__, "Cannot open \"%s\"\n", in_name);
        if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", in_name);
        sj_m = 10000; sj_group = (sj_t*)_err_malloc(sj_m * sizeof(sj_t));
        ad_m = 10000; ad_group = (ad_t*)_err_malloc(ad_m * sizeof(ad_t));
        int *sg_ad_idx = (int*)_err_calloc(sg_g->SG_n, sizeof(int));
        // calculate number of junction-reads and exon-body reads
        //sj_n = bam2cnt_core(in, h, b, seq, seq_n, &sj_group, sj_m, sg_g, sgp);
        sj_n = parse_bam_record(in, h, b, seq, seq_n, sg_g, sg_ad_idx, &ad_group, &ad_n, ad_m, &sj_group, &sj_n, sj_m, sgp);
        bam_destroy1(b); bam_hdr_destroy(h); sam_close(in);

        // update edge weight and add novel edge for GTF-SG
        update_SpliceGraph(sg_g, sj_group, sj_n, sgp);
        free(sj_group); // free sj

        int fn; 
        FILE **out_fp = iso_header_output(in_name, out_fn, &fn);
        gen_asm_iso(sg_g, sg_ad_idx, ad_group, ad_n, sgp, out_fp, fn);

        free_ad_group(ad_group, ad_n); // free ad
        free(sg_ad_idx);
        // output asm
        //asm_output(in_name, out_fn, sg_g, asm_g, sgp);
        //asm_g_rep[i] = asm_g; 
    }
    for (i = 0; i < sgp->tol_rep_n; ++i) {
        //sg_free_asm_group(asm_g_rep[i]); 
    }
    free(asm_g_rep); sg_free_group(sg_g);

    // output to one file
    sg_free_para(sgp);
    if (seq_n > 0) {
        for (i = 0; i < seq_n; ++i)  
            free(seq[i].name.s), free(seq[i].seq.s);
        free(seq);
    }
    return 0;
}
