#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "utils.h"
#include "gtf.h"
#include "parse_bam.h"
#include "build_sg.h"
#include "update_sg.h"
#include "pred_asm.h"

extern const char PROG[20];
int asm2ase_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s ase [option] <in.gtf> <in.bam>\n\n", PROG);
    err_printf("Note:    for multi-sample and multi-replicate, input bam should be in this format: \n");
    err_printf("             \"SAM1-REP1,REP2,REP3:SAM2-REP1,REP2,REP3\"\n");
    err_printf("         use \':\' to separate samples, \',\' to separate replicates.\n\n");
    err_printf("Options:\n\n");
    err_printf("         -n --novel-sj             allow novel splice-junction in the ASM. [False]\n");
    err_printf("         -N --novel-com            allow novel combination of known sites with known junctions in the ASM. [False]\n");
    err_printf("         -p --prop-pair            set -p to force to filter out reads mapped in improper pair. [False]\n");
    err_printf("         -a --anchor-len  [INT]    minimum anchor length for junction read. [%d].\n", ANCHOR_MIN_LEN);
    err_printf("         -i --intron-len  [INT]    minimum intron length for junction read. [%d]\n", INTRON_MIN_LEN);
    err_printf("         -g --genome-file [STR]    genome.fa. Use genome sequence to classify intron-motif. \n");
    err_printf("                                   If no genome file is give, intron-motif will be set as 0(non-canonical) [None]\n");
    err_printf("         -l --only-novel           only output ASM/ASE with novel-junctions. [False]\n");
    err_printf("         -m --use-multi            use both uniq- and multi-mapped reads in the bam input.[False (uniq only)]\n");
    err_printf("         -o --output      [STR]    prefix of file name of output ASM & COUNT & ASE. [in.bam/sj]\n");
    err_printf("                                   prefix.ASM & prefix.J/ECNT & prefix.SE/A5SS/A3SS/MXE/RI\n");
    err_printf("         -M --merge                merge multi-sample/replicate output\n");
    err_printf("\n");
    return 1;
}

#define add_asm_se(ase, up_e, se_e, down_e, up_c, down_c, both_c, skip_c, body_c, asm_i, sg_i) { \
    int _i, _flag=0;                         \
    for (_i = ase->se_n-1; _i >= 0; --_i) { \
        if (ase->se[_i].sg_i != sg_i) break;  \
        if (ase->se[_i].up == up_e && ase->se[_i].se == se_e && ase->se[_i].down == down_e) {  \
            _flag = 1; break;    \
        }   \
    }   \
    if (_flag == 0) { \
        if (ase->se_n == ase->se_m) _realloc(ase->se, ase->se_m, SE_t)  \
        ase->se[ase->se_n].up = up_e;       \
        ase->se[ase->se_n].se = se_e;       \
        ase->se[ase->se_n].down = down_e;   \
        ase->se[ase->se_n].up_c = up_c;   \
        ase->se[ase->se_n].down_c = down_c;   \
        ase->se[ase->se_n].ud_both_c = both_c;   \
        ase->se[ase->se_n].skip_c = skip_c;   \
        ase->se[ase->se_n].body_c = body_c;   \
        ase->se[ase->se_n].asm_i = asm_i;   \
        ase->se[ase->se_n].sg_i = sg_i;     \
        ase->se_n++;    \
    }   \
}

// XXX optimization
void asm2se(SG *sg, SGasm *a, ASE_t *ase, int asm_i, int sg_i, int use_multi, int only_novel)
{
    int i, j, k, hit; int pre_don_site_i, next_acc_site_i;
    SGnode *n = sg->node, cur, pre, next; SGedge *ed = sg->edge;
    for (i = 0; i < a->node_n; ++i) {
        cur = n[a->node_id[i]];
        for (j = 0; j < cur.pre_n; ++j) {
            if (cur.pre_id[j] == 0) continue;
            pre = n[cur.pre_id[j]];
            pre_don_site_i = pre.e_site_id;
            for (k = 0; k < cur.next_n; ++k) {
                if (cur.next_id[k] == (int)sg->node_n-1) continue;
                next = n[cur.next_id[k]];
                next_acc_site_i = next.s_site_id;
                int ej_id = sg_bin_sch_edge(sg, pre_don_site_i, next_acc_site_i, &hit);

                if (hit == 1) {
                    int ij1_id = _err_sg_bin_sch_edge(sg, pre.e_site_id, cur.s_site_id);
                    int ij2_id = _err_sg_bin_sch_edge(sg, cur.e_site_id, next.s_site_id);
                    if ((use_multi == 1 || (sg->edge[ej_id].uniq_c > 0 && sg->edge[ij1_id].uniq_c > 0 && sg->edge[ij2_id].uniq_c > 0)) 
                    && (only_novel == 0 || ed[ej_id].is_anno == 0 || ed[ij1_id].is_anno == 0 || ed[ij2_id].is_anno == 0)) {
                        int up_c=0, down_c=0, both_c=0, skip_c=0, body_c=0;
                        // XXX count read cnt
                        // transform se to 2 isoform
                        // invoke: iso_cal_cnt(ad)

                        int up = pre.node_id, se = cur.node_id, down = next.node_id;
                        add_asm_se(ase, up, se, down, up_c, down_c, both_c, skip_c, body_c, asm_i, sg_i)
                    }
                }
            }
        }
    }
}

#define add_asm_a5ss(ase, shor_e, long_e, down_e, shor_c, lon_c, pj_c, pl_both_c, body_c, asm_i, sg_i) {  \
    int _i, _flag=0;                           \
    for (_i = ase->a5ss_n-1; _i >= 0; --_i) { \
        if (ase->a5ss[_i].sg_i != sg_i) break;  \
        if (ase->a5ss[_i].shor == shor_e && ase->a5ss[_i].lon == long_e && ase->a5ss[_i].down == down_e) {  \
            _flag = 1; break;    \
        }   \
    }   \
    if (_flag == 0) { \
        if (ase->a5ss_n == ase->a5ss_m) _realloc(ase->a5ss, ase->a5ss_m, A5SS_t)  \
        ase->a5ss[ase->a5ss_n].shor = shor_e;  \
        ase->a5ss[ase->a5ss_n].lon = long_e;    \
        ase->a5ss[ase->a5ss_n].down = down_e;   \
        ase->a5ss[ase->a5ss_n].lon_c = lon_c;    \
        ase->a5ss[ase->a5ss_n].shor_c = shor_c;   \
        ase->a5ss[ase->a5ss_n].pj_c = pj_c;   \
        ase->a5ss[ase->a5ss_n].pl_both_c = pl_both_c;   \
        ase->a5ss[ase->a5ss_n].body_c = body_c;   \
        ase->a5ss[ase->a5ss_n].asm_i = asm_i;   \
        ase->a5ss[ase->a5ss_n].sg_i = sg_i;     \
        ase->a5ss_n++;  \
    }   \
}

void asm2a5ss(SG *sg, SGasm *a, ASE_t *ase, int asm_i, int sg_i, int use_multi, int only_novel)
{
    int i, j, k; uint8_t is_rev = sg->is_rev; 
    SGnode *n = sg->node, cur; SGedge *ed = sg->edge;
    if (is_rev) {
        SGnode next1, next2;
        int next1_s, next1_e, next2_s, next2_e;
        for (i = 0; i < a->node_n; ++i) {
            cur = n[a->node_id[i]];
            for (j = 0; j < cur.next_n-1; ++j) {
                if (cur.next_id[j] == (int)sg->node_n-1) continue;
                next1 = n[cur.next_id[j]];
                next1_s = next1.start; next1_e = next1.end;
                for (k = j+1; k < cur.next_n; ++k) {
                    if (cur.next_id[k] == (int)sg->node_n-1) continue;
                    next2 = n[cur.next_id[k]];
                    next2_s = next2.start; next2_e = next2.end;
                    if (next1_e == next2_e && next1_s != next2_s) {
                        int sj1_id = _err_sg_bin_sch_edge(sg, cur.e_site_id, next1.s_site_id);
                        int sj2_id = _err_sg_bin_sch_edge(sg, cur.e_site_id, next2.s_site_id);
                        if ((use_multi == 1 || (ed[sj1_id].uniq_c > 0 && ed[sj2_id].uniq_c > 0)) 
                        && (only_novel == 0 || ed[sj1_id].is_anno == 0 || ed[sj2_id].is_anno == 0)) {
                                int lon_c=0, shor_c=0, pj_c=0, pl_both_c=0, body_c=0;
                                // XXX count read cnt

                                int up = cur.node_id, lon = next1.node_id, shor = next2.node_id;
                                add_asm_a5ss(ase, shor, lon, up, shor_c, lon_c, pj_c, pl_both_c, body_c, asm_i, sg_i)
                        }
                    }
                }
            }
        }
    } else {
        SGnode pre1, pre2;
        int pre1_s, pre1_e, pre2_s, pre2_e;
        for (i = 0; i < a->node_n; ++i) {
            cur = n[a->node_id[i]];
            for (j = 0; j < cur.pre_n-1; ++j) {
                if (cur.pre_id[j] == 0) continue;
                pre1 = n[cur.pre_id[j]];
                pre1_s = pre1.start; pre1_e = pre1.end;
                for (k = j+1; k < cur.pre_n; ++k) {
                    if (cur.pre_id[k] == 0) continue;
                    pre2 = n[cur.pre_id[k]];
                    pre2_s = pre2.start; pre2_e = pre2.end;
                    if (pre1_s == pre2_s && pre1_e != pre2_e) {
                        int sj1_id = _err_sg_bin_sch_edge(sg, pre1.e_site_id, cur.s_site_id);
                        int sj2_id = _err_sg_bin_sch_edge(sg, pre2.e_site_id, cur.s_site_id);
                        if ((use_multi == 1 || (ed[sj1_id].uniq_c > 0 && ed[sj2_id].uniq_c > 0))
                        && (only_novel == 0 || ed[sj1_id].is_anno == 0 || ed[sj2_id].is_anno == 0)) {
                            int lon_c=0,  shor_c=0, pj_c=0, pl_both_c=0, body_c=0;
                            // XXX count read cnt

                            int shor = pre1.node_id, lon = pre2.node_id, down = cur.node_id;
                            add_asm_a5ss(ase, shor, lon, down, shor_c, lon_c, pj_c, pl_both_c, body_c, asm_i, sg_i)
                        }
                    }
                }
            }
        }
    }
}

#define add_asm_a3ss(ase, up_e, long_e, shor_e, lon_c, shor_c, pj_c, lp_both_c, body_c, asm_i, sg_i) {  \
    int _i, _flag=0;                           \
    for (_i = ase->a3ss_n-1; _i >= 0; --_i) { \
        if (ase->a3ss[_i].sg_i != sg_i) break;  \
        if (ase->a3ss[_i].shor == shor_e && ase->a3ss[_i].lon == long_e && ase->a3ss[_i].up == up_e) {   \
            _flag = 1; break;    \
        }   \
    }   \
    if (_flag == 0) { \
        if (ase->a3ss_n == ase->a3ss_m) _realloc(ase->a3ss, ase->a3ss_m, A3SS_t)  \
        ase->a3ss[ase->a3ss_n].up = up_e;        \
        ase->a3ss[ase->a3ss_n].lon = long_e;    \
        ase->a3ss[ase->a3ss_n].shor = shor_e;  \
        ase->a3ss[ase->a3ss_n].lon_c = lon_c;    \
        ase->a3ss[ase->a3ss_n].shor_c = shor_c;  \
        ase->a3ss[ase->a3ss_n].pj_c = pj_c;    \
        ase->a3ss[ase->a3ss_n].lp_both_c = lp_both_c;  \
        ase->a3ss[ase->a3ss_n].body_c = body_c;  \
        ase->a3ss[ase->a3ss_n].asm_i = asm_i;     \
        ase->a3ss[ase->a3ss_n].sg_i = sg_i;     \
        ase->a3ss_n++;  \
    }   \
}

void asm2a3ss(SG *sg, SGasm *a, ASE_t *ase, int asm_i, int sg_i, int use_multi, int only_novel)
{
    int i, j, k; uint8_t is_rev = sg->is_rev;
    SGnode *n = sg->node, cur; SGedge *ed = sg->edge;
    if (is_rev) {
        SGnode pre1, pre2;
        int pre1_s, pre1_e, pre2_s, pre2_e;
        for (i = 0; i < a->node_n; ++i) {
            cur = n[a->node_id[i]];
            for (j = 0; j < cur.pre_n-1; ++j) {
                if (cur.pre_id[j] == 0) continue;
                pre1 = n[cur.pre_id[j]];
                pre1_s = pre1.start; pre1_e = pre1.end;
                for (k = j+1; k < cur.pre_n; ++k) {
                    if (cur.pre_id[k] == 0) continue;
                    pre2 = n[cur.pre_id[k]];
                    pre2_s = pre2.start; pre2_e = pre2.end;
                    if (pre1_s == pre2_s && pre1_e != pre2_e) {
                        int sj1_id = _err_sg_bin_sch_edge(sg, pre1.e_site_id, cur.s_site_id);
                        int sj2_id = _err_sg_bin_sch_edge(sg, pre2.e_site_id, cur.s_site_id);
                        if ((use_multi == 1 || (ed[sj1_id].uniq_c > 0 && ed[sj2_id].uniq_c > 0))
                        && (only_novel == 0 || ed[sj1_id].is_anno == 0 || ed[sj2_id].is_anno == 0)) {
                            int lon_c=0, shor_c=0, pj_c=0, lp_both_c=0, body_c=0;
                            // XXX count read cnt

                            int shor = pre1.node_id, lon = pre2.node_id, down = cur.node_id;
                            add_asm_a3ss(ase, down, lon, shor, lon_c, shor_c, pj_c, lp_both_c, body_c, asm_i, sg_i)
                        }
                    }
                }
            }
        }
    } else {
        SGnode next1, next2;
        int next1_s, next1_e, next2_s, next2_e;
        for (i = 0; i < a->node_n; ++i) {
            cur = n[a->node_id[i]];
            for (j = 0; j < cur.next_n-1; ++j) {
                if (cur.next_id[j] == (int)sg->node_n-1) continue;
                next1 = n[cur.next_id[j]];
                next1_s = next1.start; next1_e = next1.end;
                for (k = j+1; k < cur.next_n; ++k) {
                    if (cur.next_id[k] == (int)sg->node_n-1) continue;
                    next2 = n[cur.next_id[k]];
                    next2_s = next2.start; next2_e = next2.end;
                    if (next1_e == next2_e && next1_s != next2_s) {
                        int sj1_id = _err_sg_bin_sch_edge(sg, cur.e_site_id, next1.s_site_id);
                        int sj2_id = _err_sg_bin_sch_edge(sg, cur.e_site_id, next2.s_site_id);
                        if ((use_multi == 1 || (ed[sj1_id].uniq_c > 0 && ed[sj2_id].uniq_c > 0))
                        && (only_novel == 0 || ed[sj1_id].is_anno == 0 || ed[sj2_id].is_anno == 0)) {
                            int lon_c=0, shor_c=0, pj_c=0, lp_both_c=0, body_c=0;
                            // XXX count read cnt

                            int up = cur.node_id, lon = next1.node_id, shor = next2.node_id;
                            add_asm_a3ss(ase, up, lon, shor, lon_c, shor_c, pj_c, lp_both_c, body_c, asm_i, sg_i)
                        }
                    }
                }
            }
        }
    }
}

#define add_asm_mxe(ase, up_e, fir_e, sec_e, down_e, fir_up_c, fir_down_c, fir_both_c, sec_up_c, sec_down_c, sec_both_c, fir_body_c, sec_body_c, asm_i, sg_i) {  \
    int _i, _flag=0;                           \
    for (_i = ase->mxe_n-1; _i >= 0; --_i) { \
        if (ase->mxe[_i].sg_i != sg_i) break;  \
        if (ase->mxe[_i].up == up_e && ase->mxe[_i].fir == fir_e && ase->mxe[_i].sec == sec_e && ase->mxe[_i].down == down_e) {   \
            _flag = 1; break;    \
        }   \
    }   \
    if (_flag == 0) { \
        if (ase->mxe_n == ase->mxe_m) _realloc(ase->mxe, ase->mxe_m, MXE_t)  \
        ase->mxe[ase->mxe_n].up = up_e;     \
        ase->mxe[ase->mxe_n].fir = fir_e;   \
        ase->mxe[ase->mxe_n].sec = sec_e;   \
        ase->mxe[ase->mxe_n].down = down_e; \
        ase->mxe[ase->mxe_n].fir_up_c = fir_up_c;     \
        ase->mxe[ase->mxe_n].fir_down_c = fir_down_c;     \
        ase->mxe[ase->mxe_n].fir_both_c = fir_both_c;     \
        ase->mxe[ase->mxe_n].sec_up_c = sec_up_c;     \
        ase->mxe[ase->mxe_n].sec_down_c = sec_down_c;     \
        ase->mxe[ase->mxe_n].sec_both_c = sec_both_c;     \
        ase->mxe[ase->mxe_n].fir_body_c = fir_body_c;     \
        ase->mxe[ase->mxe_n].sec_body_c = sec_body_c;     \
        ase->mxe[ase->mxe_n].asm_i = asm_i; \
        ase->mxe[ase->mxe_n].sg_i = sg_i; \
        ase->mxe_n++;   \
    }   \
}

void asm2mxe(SG *sg, SGasm *a, ASE_t *ase, int asm_i, int sg_i, int use_multi, int only_novel)
{
    int i, j, k, m, n, l;
    SGnode *node = sg->node, mx1, mx2, pre, next; SGedge *ed = sg->edge;
    for (i = 0; i < a->node_n-1; ++i) {
        mx1 = node[a->node_id[i]];
        for (j = i+1; j < a->node_n; ++j) {
            mx2 = node[a->node_id[j]];
            if (mx2.start <= mx1.end) continue;
            for (k = 0; k < mx1.pre_n; ++k) {
                if (mx1.pre_id[k] == 0) continue;
                for (m = 0; m < mx2.pre_n; ++m) {
                    if (mx2.pre_id[m] == 0) continue;
                    if (mx1.pre_id[k] == mx2.pre_id[m]) {
                        pre = node[mx1.pre_id[k]];
                        for (n = 0; n < mx1.next_n; ++n) {
                            if (mx1.next_id[n] == (int)sg->node_n-1) continue;
                            for (l = 0; l < mx2.next_n; ++l) {
                                if (mx2.next_id[l] == (int)sg->node_n-1) continue;
                                if (mx1.next_id[n] == mx2.next_id[l]) {
                                    next = node[mx1.next_id[n]];
                                    int sj1_id = _err_sg_bin_sch_edge(sg, pre.e_site_id, mx1.s_site_id);
                                    int sj2_id = _err_sg_bin_sch_edge(sg, pre.e_site_id, mx2.s_site_id);
                                    int sj3_id = _err_sg_bin_sch_edge(sg, mx1.e_site_id, next.s_site_id);
                                    int sj4_id = _err_sg_bin_sch_edge(sg, mx2.e_site_id, next.s_site_id);
                                    if ((use_multi == 1 || (ed[sj1_id].uniq_c > 0 && ed[sj2_id].uniq_c > 0 && ed[sj3_id].uniq_c > 0 && ed[sj4_id].uniq_c > 0))
                                    && (only_novel == 0 ||  ed[sj1_id].is_anno == 0 || ed[sj2_id].is_anno == 0 || ed[sj3_id].is_anno == 0 || ed[sj4_id].is_anno == 0)){
                                        int fir_up_c=0,fir_down_c=0,fir_both_c=0,sec_up_c=0,sec_down_c=0,sec_both_c=0, fir_body_c=0, sec_body_c=0;
                                        // XXX count read cnt

                                        int up = pre.node_id, fir = mx1.node_id, sec = mx2.node_id, down = next.node_id;
                                        add_asm_mxe(ase, up, fir, sec, down, fir_up_c, fir_down_c, fir_both_c, sec_up_c, sec_down_c, sec_both_c, fir_body_c, sec_body_c, asm_i, sg_i)
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

#define add_asm_ri(ase, up_e, down_e, in_e, ej_c, pj1_c, pj2_c, pj_both_c, body_c, asm_i, sg_i) {  \
    int _i, _flag=0;                           \
    for (_i = ase->ri_n-1; _i >= 0; --_i) { \
        if (ase->ri[_i].sg_i != sg_i) break;  \
        if (ase->ri[_i].down == down_e && ase->ri[_i].up == up_e && ase->ri[_i].in == in_e) {   \
            _flag = 1; break;    \
        }   \
    }   \
    if (_flag == 0) { \
        if (ase->ri_n == ase->ri_m) _realloc(ase->ri, ase->ri_m, RI_t)  \
        ase->ri[ase->ri_n].up = up_e;     \
        ase->ri[ase->ri_n].down = down_e; \
        ase->ri[ase->ri_n].in = in_e; \
        ase->ri[ase->ri_n].ej_c = ej_c; \
        ase->ri[ase->ri_n].pj1_c = pj1_c; \
        ase->ri[ase->ri_n].pj2_c = pj2_c; \
        ase->ri[ase->ri_n].pj_both_c = pj_both_c; \
        ase->ri[ase->ri_n].body_c = body_c; \
        ase->ri[ase->ri_n].asm_i = asm_i; \
        ase->ri[ase->ri_n].sg_i = sg_i; \
        ase->ri_n++;    \
    }   \
}

void asm2ri(SG *sg, SGasm *a, ASE_t *ase, int asm_i, int sg_i, int use_multi, int only_novel)
{
    int i, j, k;
    SGnode *n = sg->node, pre, ri, next; SGedge *ed = sg->edge;
    for (i = 0; i < a->node_n-2; ++i) {
        pre = n[a->node_id[i]];
        for (j = i + 1; j < a->node_n-1; ++j) {
            ri = n[a->node_id[j]];
            if (pre.start == ri.start) {
                for (k = 0; k < pre.next_n; ++k) {
                    next = n[pre.next_id[k]];
                    if (ri.end == next.end) {
                        int sj_id = _err_sg_bin_sch_edge(sg, pre.e_site_id, next.s_site_id);
                        if ((use_multi == 1 || (ed[sj_id].uniq_c > 0))
                        && (only_novel == 0 || ed[sj_id].is_anno == 0)) {
                            int ej_c=0, pj1_c=0, pj2_c=0, pj_both_c=0, body_c=0;
                            // XXX count read cnt

                            int up = pre.node_id, down = next.node_id, in = ri.node_id;
                            add_asm_ri(ase, up, down, in, ej_c, pj1_c, pj2_c, pj_both_c, body_c, asm_i, sg_i)
                        }
                    }
                }
            }
        }
    }
}

ASE_t *ase_init(void)
{
    ASE_t *ase = (ASE_t*)_err_malloc(sizeof(ASE_t));
    ase->se_m = 1000; ase->se_n = 0; ase->se = (SE_t*)_err_malloc(sizeof(SE_t) * 1000);
    ase->a5ss_m = 1000; ase->a5ss_n = 0; ase->a5ss = (A5SS_t*)_err_malloc(sizeof(A5SS_t) * 1000);
    ase->a3ss_m = 1000; ase->a3ss_n = 0; ase->a3ss = (A3SS_t*)_err_malloc(sizeof(A3SS_t) * 1000);
    ase->mxe_m = 1000; ase->mxe_n = 0; ase->mxe = (MXE_t*)_err_malloc(sizeof(MXE_t) * 1000);
    ase->ri_m = 1000; ase->ri_n = 0; ase->ri = (RI_t*)_err_malloc(sizeof(RI_t) * 1000);
    return ase;
}

void ase_free(ASE_t *ase)
{
    free(ase->se); free(ase->a5ss); free(ase->a3ss); free(ase->mxe); free(ase->ri);
    free(ase);
}

void ase_output(char *in_fn, char *prefix, SG_group *sg_g,  ASE_t *ase, sg_para *sgp)
{
    int i, out_n=10;
    char suf[10][10] = { ".SE", ".A5SS", ".A3SS", ".MXE", ".RI",
        ".SE_CNT", ".A5SS_CNT", ".A3SS_CNT", ".MXE_CNT", ".RI_CNT"};
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
    for (i = 0; i < out_n; ++i) out_fp[i] = xopen(out_fn[i], "w");

    chr_name_t *cname = sg_g->cname;
    // head-line
    {
        // SE/A5SS/A3SS/MXE/RI
        fprintf(out_fp[0], "ID\tASM_ID\tSG_ID\tSTRAND\tCHR\tSE_ES\tSE_EE\tUP_ES\tUP_EE\tDOWN_ES\tDOWN_EE\n");
        fprintf(out_fp[1], "ID\tASM_ID\tSG_ID\tSTRAND\tCHR\tSHORT_ES\tSHORT_EE\tLONG_ES\tLONG_EE\tDOWN_ES\tDOWN_EE\n");
        fprintf(out_fp[2], "ID\tASM_ID\tSG_ID\tSTRAND\tCHR\tUP_ES\tUP_EE\tLOND_ES\tLONG_EE\tSHORT_ES\tSHORT_EE\n");
        fprintf(out_fp[3], "ID\tASM_ID\tSG_ID\tSTRAND\tCHR\tUP_ES\tUP_EE\tFIRST_ES\tFIRST_EE\tSECOND_ES\tSECOND_EE\tDOWN_ES\tDOWN_EE\n");
        fprintf(out_fp[4], "ID\tASM_ID\tSG_ID\tSTRAND\tCHR\tUP_ES\tUP_EE\tDOWN_ES\tDOWN_EE\n");
        // SE/A5SS/A3SS/MXE/RI_CNT
        fprintf(out_fp[5], "ID\tSTRAND\tCHR\tIJ1_UCNT\tIJ2_UCNT\tBOTH_UCNT\tEJ_UCNT\n");
        fprintf(out_fp[6], "ID\tSTRAND\tCHR\tLO_UCNT\tSH_UCNT\n");
        fprintf(out_fp[7], "ID\tSTRAND\tCHR\tLO_UCNT\tSH_UCNT\n");
        fprintf(out_fp[8], "ID\tSTRAND\tCHR\tFIR1_UCNT\tFIR2_UCNT\tFIR_BOTH_UCNT\tSEC1_UCNT\tSEC2_UCNT\tSEC_BOTH_UCNT\n");
        fprintf(out_fp[9], "ID\tSTRAND\tCHR\tRE_UCNT\tSJ_UCNT\n");
    }

    { // SE
        for (i = 0; i < ase->se_n; ++i) {
            int sg_i = ase->se[i].sg_i, asm_i = ase->se[i].asm_i;
            SG *sg = sg_g->SG[sg_i]; SE_t e=ase->se[i]; SGnode *n = sg->node;
            fprintf(out_fp[0], "%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", i, asm_i, sg_i, "+-"[sg->is_rev], cname->chr_name[sg->tid], n[e.se].start, n[e.se].end, n[e.up].start, n[e.up].end, n[e.down].start, n[e.down].end);
            fprintf(out_fp[5], "%d\t%c\t%s\t%d\t%d\t%d\t%d\n", i, "+-"[sg->is_rev], cname->chr_name[sg->tid], e.up_c, e.down_c, e.ud_both_c, e.skip_c);
        }
    }
    { // A5SS
        for (i = 0; i < ase->a5ss_n; ++i) {
            int sg_i = ase->a5ss[i].sg_i, asm_i = ase->a5ss[i].asm_i;
            SG *sg = sg_g->SG[sg_i]; A5SS_t e = ase->a5ss[i]; SGnode *n = sg->node;
            fprintf(out_fp[1], "%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", i, asm_i, sg_i, "+-"[sg->is_rev], cname->chr_name[sg->tid], n[e.shor].start, n[e.shor].end, n[e.lon].start, n[e.lon].end, n[e.down].start, n[e.down].end);
            fprintf(out_fp[6], "%d\t%c\t%s\t%d\t%d\n", i, "+-"[sg->is_rev], cname->chr_name[sg->tid], e.lon_c, e.shor_c);
        }
    }
    { // A3SS
        for (i = 0; i < ase->a3ss_n; ++i) {
            int sg_i = ase->a3ss[i].sg_i, asm_i = ase->a3ss[i].asm_i;
            SG *sg = sg_g->SG[sg_i]; A3SS_t e = ase->a3ss[i]; SGnode *n = sg->node;
            fprintf(out_fp[2], "%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", i, asm_i, sg_i, "+-"[sg->is_rev], cname->chr_name[sg->tid], n[e.up].start, n[e.up].end, n[e.lon].start, n[e.lon].end, n[e.shor].start, n[e.shor].end);
            fprintf(out_fp[7], "%d\t%c\t%s\t%d\t%d\n", i, "+-"[sg->is_rev], cname->chr_name[sg->tid], e.lon_c, e.shor_c);
        }
    }
    { // MXE
        for (i = 0; i < ase->mxe_n; ++i) {
            int sg_i = ase->mxe[i].sg_i, asm_i = ase->mxe[i].asm_i;
            SG *sg = sg_g->SG[sg_i]; MXE_t e = ase->mxe[i]; SGnode *n = sg->node;
            fprintf(out_fp[3], "%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", i, asm_i, sg_i, "+-"[sg->is_rev], cname->chr_name[sg->tid], n[e.up].start, n[e.up].end, n[e.fir].start, n[e.fir].end, n[e.sec].start, n[e.sec].end, n[e.down].start, n[e.down].end);
            fprintf(out_fp[8], "%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", i, "+-"[sg->is_rev], cname->chr_name[sg->tid], e.fir_up_c, e.fir_down_c, e.fir_both_c, e.sec_up_c, e.sec_down_c, e.sec_both_c);
        }
    }
    { // RI
        for (i = 0; i < ase->ri_n; ++i) {
            int sg_i = ase->ri[i].sg_i, asm_i = ase->ri[i].asm_i;
            SG *sg = sg_g->SG[sg_i]; RI_t e = ase->ri[i]; SGnode *n = sg->node;
            fprintf(out_fp[4], "%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\n", i, asm_i, sg_i, "+-"[sg->is_rev], cname->chr_name[sg->tid], n[e.up].start, n[e.up].end, n[e.down].start, n[e.down].end);
            fprintf(out_fp[9], "%d\t%c\t%s\t%d\t%d\n", i, "+-"[sg->is_rev], cname->chr_name[sg->tid], n[e.in].uniq_c-n[e.up].uniq_c-n[e.down].uniq_c, e.ej_c);
        }
    }

    for (i = 0; i < out_n; ++i) {
        free(out_fn[i]); err_fclose(out_fp[i]);
    } free(out_fn); free(out_fp);
}

/*void ase_merge_output(char *prefix, SG_group *sg_g_rep, ASE_t *ase_rep, sg_para *sgp)
{
    int i;
    int out_n=10;
    char suf[10][10] = { ".SE", ".A5SS", ".A3SS", ".MXE", ".RI",
        ".SE_CNT", ".A5SS_CNT", ".A3SS_CNT", ".MXE_CNT", ".RI_CNT"};
    char suff[20] = "";
    if (sgp->use_multi==1) strcat(suff, ".multi");
    if (sgp->no_novel_sj==1) strcat(suff, ".anno");
    if (sgp->only_novel==1) strcat(suff, ".novel");
    char **out_fn = (char**)_err_malloc(sizeof(char*) * out_n);

    for (i = 0; i < out_n; ++i) {
        out_fn[i] = (char*)_err_malloc(strlen(prefix)+10); strcpy(out_fn[i], prefix); strcat(out_fn[i], suff); strcat(out_fn[i], suf[i]);
    }

    FILE **out_fp = (FILE**)_err_malloc(sizeof(FILE*) * out_n);
    for (i = 0; i < out_n; ++i) out_fp[i] = xopen(out_fn[i], "w");

    chr_name_t *cname = sg_g_rep[0].cname;
    // head-line
    // SE/A5SS/A3SS/MXE/RI
    fprintf(out_fp[0], "ID\tASM_ID\tSG_ID\tSTRAND\tCHR\tSE_ES\tSE_EE\tUP_ES\tUP_EE\tDOWN_ES\tDOWN_EE\n");
    fprintf(out_fp[1], "ID\tASM_ID\tSG_ID\tSTRAND\tCHR\tSHORT_ES\tSHORT_EE\tLONG_ES\tLONG_EE\tDOWN_ES\tDOWN_EE\n");
    fprintf(out_fp[2], "ID\tASM_ID\tSG_ID\tSTRAND\tCHR\tUP_ES\tUP_EE\tLOND_ES\tLONG_EE\tSHORT_ES\tSHORT_EE\n");
    fprintf(out_fp[3], "ID\tASM_ID\tSG_ID\tSTRAND\tCHR\tUP_ES\tUP_EE\tFIRST_ES\tFIRST_EE\tSECOND_ES\tSECOND_EE\tDOWN_ES\tDOWN_EE\n");
    fprintf(out_fp[4], "ID\tASM_ID\tSG_ID\tSTRAND\tCHR\tUP_ES\tUP_EE\tDOWN_ES\tDOWN_EE\n");
    // SE/A5SS/A3SS/MXE/RI_CNT
    fprintf(out_fp[5], "ID\tSTRAND\tCHR\tIJ1_UCNT\tIJ2_UCNT\tEJ_UCNT\tIJ1_MCNT\tIJ2_MCNT\tEJ_MCNT\n");
    fprintf(out_fp[6], "ID\tSTRAND\tCHR\tLO_UCNT\tSH_UCNT\tLO_MCNT\tSH_MCNT\n");
    fprintf(out_fp[7], "ID\tSTRAND\tCHR\tLO_UCNT\tSH_UCNT\tLO_MCNT\tSH_MCNT\n");
    fprintf(out_fp[8], "ID\tSTRAND\tCHR\tFIR1_UCNT\tFIR2_UCNT\tSEC1_UCNT\tSEC2_UCNT\tFIR1_MCNT\tFIR2_MCNT\tSEC1_MCNT\tSEC2_MCNT\n");
    fprintf(out_fp[9], "ID\tSTRAND\tCHR\tRE_UCNT\tSJ_UCNT\tRE_MCNT\tSJ_MCNT\n");

    // SE
    int *cur_i = (int*)_err_calloc(sgp->tol_rep_n, sizeof(int));
    int min_rep;
    while (1) {
        // get min-ase
        if (min_rep = (get_min_se(ase_rep, sg_g_rep, sgp->tol_rep_n, cur_i)) >= 0) {
            SE_t e = ase_rep[min_rep].se[cur_i[min_rep]];
            int sg_i = e.sg_i, asm_i = e.asm_i;
            SG *sg = sg_g_rep[min_rep].SG[sg_i]; SGnode *n = sg->node; SGedge *ed = sg->edge;
            // output SE
            fprintf(out_fp[0], "%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", i, asm_i, sg_i, "+-"[sg->is_rev], cname->chr_name[sg->tid], n[e.se].start, n[e.se].end, n[e.up].start, n[e.up].end, n[e.down].start, n[e.down].end);
            // output SE_CNT
        }
    }

    for (i = 0; i < ase->se_n; ++i) {
        int sg_i = ase->se[i].sg_i, asm_i = ase->se[i].asm_i;
        SG *sg = sg_g->SG[sg_i]; SE_t e=ase->se[i]; SGnode *n = sg->node; SGedge *ed = sg->edge;
        fprintf(out_fp[0], "%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", i, asm_i, sg_i, "+-"[sg->is_rev], cname->chr_name[sg->tid], n[e.se].start, n[e.se].end, n[e.up].start, n[e.up].end, n[e.down].start, n[e.down].end);
        // up->se, se->down, up->down
        int ij1_id = _err_sg_bin_sch_edge(sg, n[e.up].e_site_id, n[e.se].s_site_id, &hit);
        int ij2_id = _err_sg_bin_sch_edge(sg, n[e.se].e_site_id, n[e.down].s_site_id, &hit);
        int ej_id = _err_sg_bin_sch_edge(sg, n[e.up].e_site_id, n[e.down].s_site_id, &hit);
        fprintf(out_fp[5], "%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", i, "+-"[sg->is_rev], cname->chr_name[sg->tid], ed[ij1_id].uniq_c, ed[ij2_id].uniq_c, ed[ej_id].uniq_c, ed[ij1_id].multi_c, ed[ij2_id].multi_c, ed[ej_id].multi_c);
    }

    free(cur_i);
    for (i = 0; i < out_n; ++i) {
        free(out_fn[i]); err_fclose(out_fp[i]);
    } free(out_fn); free(out_fp);
}*/

int asm2ase_core(SG_group *sg_g, SGasm_group *asm_g, ASE_t *ase, sg_para *sgp)
{
    err_func_format_printf(__func__, "generating alternative splice-events from alternative splice-modules ...\n");
    int i, use_multi = sgp->use_multi, only_novel = sgp->only_novel;
    for (i = 0; i < asm_g->sg_asm_n; ++i) {
        SGasm *sg_asm = asm_g->sg_asm[i];
        int sg_i = sg_asm->SG_id; SG *sg = sg_g->SG[sg_i]; 
        asm2se(sg, sg_asm, ase, i, sg_i, use_multi, only_novel);       // SE
        asm2a5ss(sg, sg_asm, ase, i, sg_i, use_multi, only_novel);     // A5SS 
        asm2a3ss(sg, sg_asm, ase, i, sg_i, use_multi, only_novel);     // A3SS
        asm2mxe(sg, sg_asm, ase, i, sg_i, use_multi, only_novel);      // MXE
        asm2ri(sg, sg_asm, ase, i, sg_i, use_multi, only_novel);     // RI
    }
    err_func_format_printf(__func__, "generating alternative splice-events from alternative splice-modules done!\n");
    return ase->se_n+ase->a5ss_n+ase->a3ss_n+ase->mxe_n+ase->ri_n;
}

const struct option se_long_opt [] = {
    { "novel-sj", 0, NULL, 'n' },
    { "novel-com", 0, NULL, 'N' },
    { "proper-pair", 1, NULL, 'p' },
    { "anchor-len", 1, NULL, 'a' },
    { "intron-len", 1, NULL, 'i' },
    { "genome-file", 1, NULL, 'g' },
    { "only-novel", 0, NULL, 'l' },
    { "use-multi", 0, NULL, 'm' },
    { "output", 1, NULL, 'o' },
    { "merge", 0, NULL, 'M' },

    { 0, 0, 0, 0 }
};

int asm2ase(int argc, char *argv[])
{
    // same to pred_asm START
    int c, i; char out_pre[1024]="", ref_fn[1024]="", *p;
    sg_para *sgp = sg_init_para();
	while ((c = getopt_long(argc, argv, "nNlmMp:g:a:o:", se_long_opt, NULL)) >= 0) {
        switch (c) {
            case 'n': sgp->no_novel_sj=0, sgp->no_novel_com=0; break;
            case 'N': sgp->no_novel_com = 0; break;
            case 'l': sgp->only_novel = 1, sgp->no_novel_sj=0, sgp->no_novel_com=0; break;
            case 'm': sgp->use_multi = 1; break;
            case 'M': sgp->merge_out = 1; break;
            case 'p': sgp->read_type = PAIR_T; break;
            case 'a': sgp->anchor_len[0] = strtol(optarg, &p, 10);
                      if (*p != 0) sgp->anchor_len[1] = strtol(p+1, &p, 10); else return asm2ase_usage();
                      if (*p != 0) sgp->anchor_len[2] = strtol(p+1, &p, 10); else return asm2ase_usage();
                      if (*p != 0) sgp->anchor_len[3] = strtol(p+1, &p, 10); else return asm2ase_usage();
                      if (*p != 0) sgp->anchor_len[4] = strtol(p+1, &p, 10); else return asm2ase_usage();
                      break;
            case 'U': sgp->uniq_min[0] = strtol(optarg, &p, 10);
                      if (*p != 0) sgp->uniq_min[1] = strtol(p+1, &p, 10); else return asm2ase_usage();
                      if (*p != 0) sgp->uniq_min[2] = strtol(p+1, &p, 10); else return asm2ase_usage();
                      if (*p != 0) sgp->uniq_min[3] = strtol(p+1, &p, 10); else return asm2ase_usage();
                      if (*p != 0) sgp->uniq_min[4] = strtol(p+1, &p, 10); else return asm2ase_usage();
                      break; 
            case 'A': sgp->all_min[0] = strtol(optarg, &p, 10);
                      if (*p != 0) sgp->all_min[1] = strtol(p+1, &p, 10); else return asm2ase_usage();
                      if (*p != 0) sgp->all_min[2] = strtol(p+1, &p, 10); else return asm2ase_usage();
                      if (*p != 0) sgp->all_min[3] = strtol(p+1, &p, 10); else return asm2ase_usage();
                      if (*p != 0) sgp->all_min[4] = strtol(p+1, &p, 10); else return asm2ase_usage();
                      break;
            case 'i': sgp->intron_len = atoi(optarg); break;
            case 'g': strcpy(ref_fn, optarg); break;
            case 'o': strcpy(out_pre, optarg); break;
            default: err_printf("Error: unknown option: %s.\n", optarg); return asm2ase_usage();
        }
    }
    if (argc - optind != 2) return asm2ase_usage();

    int seq_n = 0, seq_m; kseq_t *seq;
    if (strlen(ref_fn) != 0) {
        gzFile genome_fp = gzopen(ref_fn, "r");
        if (genome_fp == NULL) { err_fatal(__func__, "Can not open genome file. %s\n", ref_fn); }
        seq = kseq_load_genome(genome_fp, &seq_n, &seq_m);
        err_gzclose(genome_fp); 
    }

    // parse input name
    if (sg_par_input(sgp, argv[optind+1]) <= 0) return asm2ase_usage();
    
    // set cname
    chr_name_t *cname = chr_name_init();
    samFile *in; bam_hdr_t *h; bam1_t *b;
    if ((in = sam_open(sgp->in_name[0], "rb")) == NULL) err_fatal_core(__func__, "Cannot open \"%s\"\n", sgp->in_name[0]);
    if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", sgp->in_name[0]);
    bam_set_cname(h, cname);
    bam_hdr_destroy(h); sam_close(in);
    // build splice-graph with GTF
    FILE *gtf_fp = xopen(argv[optind], "r");
    SG_group *sg_g = construct_SpliceGraph(gtf_fp, cname);
    err_fclose(gtf_fp); chr_name_free(cname);

    SGasm_group **asm_g_rep = (SGasm_group**)_err_malloc(sgp->tol_rep_n * sizeof(SGasm_group*));
    ASE_t **ase_rep = (ASE_t**)_err_malloc(sgp->tol_rep_n * sizeof(ASE_t*));
    for (i = 0; i < sgp->tol_rep_n; ++i) {
        char *in_name = sgp->in_name[i];
        // get splice-junction
        int sj_n, sj_m; sj_t *sj_group; int ad_n, ad_m; ad_t *ad_group;
        // based on .bam file
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

        // generate ASM with short-read splice-graph
        SGasm_group *asm_g = gen_asm(sg_g, sgp);

        // same to pred_asm END
        ASE_t *ase = ase_init();
        asm2ase_core(sg_g, asm_g, ase, sgp);

        free_ad_group(ad_group, ad_n); // free ad
        free(sg_ad_idx);

        // output asm
        asm_output(in_name, out_pre, sg_g, asm_g, sgp);
        // output ase
        if (sgp->merge_out == 0) ase_output(in_name, out_pre, sg_g, ase, sgp);
        asm_g_rep[i] = asm_g; ase_rep[i] = ase;
    }
    //if (sgp->merge_out == 1) ase_merge_output(out_pre, sr_sg_g_rep, ase_rep, sgp);

    for (i = 0; i < sgp->tol_rep_n; ++i) {
        sg_free_asm_group(asm_g_rep[i]); ase_free(ase_rep[i]);
    } free(asm_g_rep); free(ase_rep);
    sg_free_group(sg_g);

    // output to one file
    sg_free_para(sgp);
    if (seq_n > 0) {
        for (i = 0; i < seq_n; ++i) 
            free(seq[i].name.s), free(seq[i].seq.s);
        free(seq);
    }
    return 0;
}
