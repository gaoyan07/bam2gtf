#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "utils.h"
#include "gtf.h"
#include "build_sg.h"
#include "pred_asm.h"
#include "pred_sg.h"
#include "bam2sj.h"

extern const char PROG[20];
int pred_se_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s ase [option] <ref.fa> <in.gtf> <in.bam/sj>\n\n", PROG);
    err_printf("Options:\n\n");
    err_printf("         -n --novel-sj             allow novel splice-junction in the ASM. [False]\n");
    err_printf("         -N --novel-com            allow novel combination of known exons in the ASM. [False]\n");
    err_printf("         -s --sj-file              input with splice-junction file instead of BAM file. [false]\n");
    err_printf("                                   with splice-junction input, the .cnt output will have no count information.\n");
    err_printf("         -m --use-multi            use both uniq- and multi-mapped reads in the bam input.[false (uniq only)]\n");
    err_printf("         -o --output      [STR]    prefix of file name of output ASM & COUNT & ASE. [in.bam/sj]\n");
    err_printf("                                   prefix.ASM & prefix.J/ECNT & prefix.SE/A5SS/A3SS/MXE/RI\n");
    err_printf("\n");
    return 1;
}

#define add_asm_se(ase, up_e, se_e, down_e, asm_i, sg_i) { \
    int _i, flag=0;                         \
    for (_i = ase->se_n-1; _i >= 0; --_i) { \
        if (ase->se[_i].sg_i != sg_i) break;  \
        if (ase->se[_i].up == up_e && ase->se[_i].se == se_e && ase->se[_i].down == down_e) {  \
            flag = 1; break;    \
        }   \
    }   \
    if (flag == 0) { \
        if (ase->se_n == ase->se_m) _realloc(ase->se, ase->se_m, SE_t)  \
        ase->se[ase->se_n].up = up_e;       \
        ase->se[ase->se_n].se = se_e;       \
        ase->se[ase->se_n].down = down_e;   \
        ase->se[ase->se_n].asm_i = asm_i;   \
        ase->se[ase->se_n].sg_i = sg_i;     \
        ase->se_n++;    \
    }   \
}
// XXX optimization
void asm2se(SG *sg, SGasm *a, ASE_t *ase, int asm_i, int sg_i)
{
    int i, j, k, hit; uint32_t pre_site_i, next_site_i;
    SGnode *n = sg->node, cur_n, pre, next;
    for (i = 0; i < a->node_n; ++i) {
        cur_n = n[a->node_id[i]];
        for (j = 0; j < cur_n.pre_n; ++j) {
            if (cur_n.pre_id[j] == 0) continue;
            pre = n[cur_n.pre_id[j]];
            pre_site_i = _err_sg_bin_sch_site(sg->don_site, sg->don_site_n, pre.end+1, &hit); 
            for (k = 0; k < cur_n.next_n; ++k) {
                if (cur_n.next_id[k] == (uint32_t)sg->node_n-1) continue;
                next = n[cur_n.next_id[k]];
                next_site_i = _err_sg_bin_sch_site(sg->acc_site, sg->acc_site_n, next.start-1, &hit);
                sg_bin_sch_edge(sg, pre_site_i, next_site_i, &hit);

                if (hit == 1) {
                    uint32_t up, se, down;
                    up = pre.node_id;
                    se = cur_n.node_id;
                    down = next.node_id;
                    add_asm_se(ase, up, se, down, asm_i, sg_i)
                }
            }
        }
    }
}

#define add_asm_a5ss(ase, short_e, long_e, down_e, asm_i, sg_i) {  \
    int _i, flag=0;                           \
    for (_i = ase->a5ss_n-1; _i >= 0; --_i) { \
        if (ase->a5ss[_i].sg_i != sg_i) break;  \
        if (ase->a5ss[_i].shor == short_e && ase->a5ss[_i].lon == long_e && ase->a5ss[_i].down == down_e) {  \
            flag = 1; break;    \
        }   \
    }   \
    if (flag == 0) { \
        if (ase->a5ss_n == ase->a5ss_m) _realloc(ase->a5ss, ase->a5ss_m, A5SS_t)  \
        ase->a5ss[ase->a5ss_n].shor = short_e;  \
        ase->a5ss[ase->a5ss_n].lon = long_e;    \
        ase->a5ss[ase->a5ss_n].down = down_e;   \
        ase->a5ss[ase->a5ss_n].asm_i = asm_i;   \
        ase->a5ss[ase->a5ss_n].sg_i = sg_i;     \
        ase->a5ss_n++;  \
    }   \
}

void asm2a5ss(SG *sg, SGasm *a, ASE_t *ase, int asm_i, int sg_i)
{
    int i, j, k; uint8_t is_rev = sg->is_rev; 
    SGnode *n = sg->node, cur_n;
    if (is_rev) {
        SGnode next1, next2;
        int32_t next1_s, next1_e, next2_s, next2_e;
        for (i = 0; i < a->node_n; ++i) {
            cur_n = n[a->node_id[i]];
            for (j = 0; j < cur_n.next_n-1; ++j) {
                if (cur_n.next_id[j] == (uint32_t)sg->node_n-1) continue;
                next1 = n[cur_n.next_id[j]];
                next1_s = next1.start; next1_e = next1.end;
                for (k = j+1; k < cur_n.next_n; ++k) {
                    if (cur_n.next_id[k] == (uint32_t)sg->node_n-1) continue;
                    next2 = n[cur_n.next_id[k]];
                    next2_s = next2.start; next2_e = next2.end;

                    if (next1_e == next2_e && next1_s != next2_s) {
                        uint32_t up, lon, shor;
                        up = cur_n.node_id;
                        lon = next1.node_id;
                        shor = next2.node_id;
                        add_asm_a5ss(ase, shor, lon, up, asm_i, sg_i)
                    }
                }
            }
        }
    } else {
        SGnode pre1, pre2;
        int32_t pre1_s, pre1_e, pre2_s, pre2_e;
        for (i = 0; i < a->node_n; ++i) {
            cur_n = n[a->node_id[i]];
            for (j = 0; j < cur_n.pre_n-1; ++j) {
                if (cur_n.pre_id[j] == 0) continue;
                pre1 = n[cur_n.pre_id[j]];
                pre1_s = pre1.start; pre1_e = pre1.end;
                for (k = j+1; k < cur_n.pre_n; ++k) {
                    if (cur_n.pre_id[k] == 0) continue;
                    pre2 = n[cur_n.pre_id[k]];
                    pre2_s = pre2.start; pre2_e = pre2.end;

                    if (pre1_s == pre2_s && pre1_e != pre2_e) {
                        uint32_t shor, lon, down;
                        shor = pre1.node_id;
                        lon = pre2.node_id;
                        down = cur_n.node_id;
                        add_asm_a5ss(ase, shor, lon, down, asm_i, sg_i)
                    }
                }
            }
        }
    }
}

#define add_asm_a3ss(ase, up_e, long_e, short_e, asm_i, sg_i) {  \
    int _i, flag=0;                           \
    for (_i = ase->a3ss_n-1; _i >= 0; --_i) { \
        if (ase->a3ss[_i].sg_i != sg_i) break;  \
        if (ase->a3ss[_i].shor == short_e && ase->a3ss[_i].lon == long_e && ase->a3ss[_i].up == up_e) {   \
            flag = 1; break;    \
        }   \
    }   \
    if (flag == 0) { \
        if (ase->a3ss_n == ase->a3ss_m) _realloc(ase->a3ss, ase->a3ss_m, A3SS_t)  \
        ase->a3ss[ase->a3ss_n].up = up_e;        \
        ase->a3ss[ase->a3ss_n].lon = long_e;    \
        ase->a3ss[ase->a3ss_n].shor = short_e;  \
        ase->a3ss[ase->a3ss_n].asm_i = asm_i;     \
        ase->a3ss[ase->a3ss_n].sg_i = sg_i;     \
        ase->a3ss_n++;  \
    }   \
}

void asm2a3ss(SG *sg, SGasm *a, ASE_t *ase, int asm_i, int sg_i)
{
    int i, j, k; uint8_t is_rev = sg->is_rev;
    SGnode *n = sg->node, cur_n;
    if (is_rev) {
        SGnode pre1, pre2;
        int32_t pre1_s, pre1_e, pre2_s, pre2_e;
        for (i = 0; i < a->node_n; ++i) {
            cur_n = n[a->node_id[i]];
            for (j = 0; j < cur_n.pre_n-1; ++j) {
                if (cur_n.pre_id[j] == 0) continue;
                pre1 = n[cur_n.pre_id[j]];
                pre1_s = pre1.start; pre1_e = pre1.end;
                for (k = j+1; k < cur_n.pre_n; ++k) {
                    if (cur_n.pre_id[k] == 0) continue;
                    pre2 = n[cur_n.pre_id[k]];
                    pre2_s = pre2.start; pre2_e = pre2.end;

                    if (pre1_s == pre2_s && pre1_e != pre2_e) {
                        uint32_t shor, lon, down;
                        shor = pre1.node_id;
                        lon = pre2.node_id;
                        down = cur_n.node_id;
                        add_asm_a3ss(ase, down, lon, shor, asm_i, sg_i)
                    }
                }
            }
        }
    } else {
        SGnode next1, next2;
        int32_t next1_s, next1_e, next2_s, next2_e;
        for (i = 0; i < a->node_n; ++i) {
            cur_n = n[a->node_id[i]];
            for (j = 0; j < cur_n.next_n-1; ++j) {
                if (cur_n.next_id[j] == (uint32_t)sg->node_n-1) continue;
                next1 = n[cur_n.next_id[j]];
                next1_s = next1.start; next1_e = next1.end;
                for (k = j+1; k < cur_n.next_n; ++k) {
                    if (cur_n.next_id[k] == (uint32_t)sg->node_n-1) continue;
                    next2 = n[cur_n.next_id[k]];
                    next2_s = next2.start; next2_e = next2.end;

                    if (next1_e == next2_e && next1_s != next2_s) {
                        uint32_t up, lon, shor;
                        up = cur_n.node_id;
                        lon = next1.node_id;
                        shor = next2.node_id;
                        add_asm_a3ss(ase, up, lon, shor, asm_i, sg_i)
                    }
                }
            }
        }
    }
}

#define add_asm_mxe(ase, up_e, fir_e, sec_e, down_e, asm_i, sg_i) {  \
    int _i, flag=0;                           \
    for (_i = ase->mxe_n-1; _i >= 0; --_i) { \
        if (ase->mxe[_i].sg_i != sg_i) break;  \
        if (ase->mxe[_i].up == up_e && ase->mxe[_i].fir == fir_e && ase->mxe[_i].sec == sec_e && ase->mxe[_i].down == down_e) {   \
            flag = 1; break;    \
        }   \
    }   \
    if (flag == 0) { \
        if (ase->mxe_n == ase->mxe_m) _realloc(ase->mxe, ase->mxe_m, MXE_t)  \
        ase->mxe[ase->mxe_n].up = up_e;     \
        ase->mxe[ase->mxe_n].fir = fir_e;   \
        ase->mxe[ase->mxe_n].sec = sec_e;   \
        ase->mxe[ase->mxe_n].down = down_e; \
        ase->mxe[ase->mxe_n].asm_i = asm_i; \
        ase->mxe[ase->mxe_n].sg_i = sg_i; \
        ase->mxe_n++;   \
    }   \
}

void asm2mxe(SG *sg, SGasm *a, ASE_t *ase, int asm_i, int sg_i)
{
    int i, j, k, m, n, l;
    SGnode *node = sg->node, mx1, mx2, pre, next;
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
                            if (mx1.next_id[n] == (uint32_t)sg->node_n-1) continue;
                            for (l = 0; l < mx2.next_n; ++l) {
                                if (mx2.next_id[l] == (uint32_t)sg->node_n-1) continue;
                                if (mx1.next_id[n] == mx2.next_id[l]) {
                                    next = node[mx1.next_id[n]];
                                    uint32_t up, fir, sec, down;
                                    up = pre.node_id;
                                    fir = mx1.node_id;
                                    sec = mx2.node_id;
                                    down = next.node_id;
                                    add_asm_mxe(ase, up, fir, sec, down, asm_i, sg_i)
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

#define add_asm_ri(ase, up_e, down_e, in_e, asm_i, sg_i) {  \
    int _i, flag=0;                           \
    for (_i = ase->ri_n-1; _i >= 0; --_i) { \
        if (ase->ri[_i].sg_i != sg_i) break;  \
        if (ase->ri[_i].down == down_e && ase->ri[_i].up == up_e && ase->ri[_i].in == in_e) {   \
            flag = 1; break;    \
        }   \
    }   \
    if (flag == 0) { \
        if (ase->ri_n == ase->ri_m) _realloc(ase->ri, ase->ri_m, RI_t)  \
        ase->ri[ase->ri_n].up = up_e;     \
        ase->ri[ase->ri_n].down = down_e; \
        ase->ri[ase->ri_n].in = in_e; \
        ase->ri[ase->ri_n].asm_i = asm_i; \
        ase->ri[ase->ri_n].sg_i = sg_i; \
        ase->ri_n++;    \
    }   \
}

void asm2ri(SGnode *n, SGasm *a, ASE_t *ase, int asm_i, int sg_i)
{
    int i, j, k;
    SGnode pre, ri, next;
    for (i = 0; i < a->node_n-2; ++i) {
        pre = n[a->node_id[i]];
        for (j = i + 1; j < a->node_n-1; ++j) {
            ri = n[a->node_id[j]];
            if (pre.start == ri.start) {
                for (k = 0; k < pre.next_n; ++k) {
                    next = n[pre.next_id[k]];
                    if (ri.end == next.end) {
                        uint32_t up, down, in;
                        up = pre.node_id;
                        down = next.node_id;
                        in = ri.node_id;
                        add_asm_ri(ase, up, down, in, asm_i, sg_i)
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

void ase_output(char *in_fn, char *prefix, SG_group *sg_g,  ASE_t *ase)
{
    int i, out_n=10;
    char suf[10][10] = { ".SE", ".A5SS", ".A3SS", ".MXE", ".RI",
                         ".SE_CNT", ".A5SS_CNT", ".A3SS_CNT", ".MXE_CNT", ".RI_CNT"};
    char **out_fn = (char**)_err_malloc(sizeof(char*) * out_n);
    if (strlen(prefix) == 0) {
        for (i = 0; i < out_n; ++i) {
            out_fn[i] = (char*)_err_malloc(strlen(in_fn)+10); strcpy(out_fn[i], in_fn); strcat(out_fn[i], suf[i]);
        }
    } else {
        for (i = 0; i < out_n; ++i) {
            out_fn[i] = (char*)_err_malloc(strlen(prefix)+10); strcpy(out_fn[i], prefix); strcat(out_fn[i], suf[i]);
        }
    }
    FILE **out_fp = (FILE**)_err_malloc(sizeof(FILE*) * out_n);
    for (i = 0; i < out_n; ++i)
        out_fp[i] = xopen(out_fn[i], "w");

    chr_name_t *cname = sg_g->cname;
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


    int hit;
    for (i = 0; i < ase->se_n; ++i) {
        int sg_i = ase->se[i].sg_i, asm_i = ase->se[i].asm_i;
        SG *sg = sg_g->SG[sg_i]; SE_t e=ase->se[i]; SGnode *n = sg->node; SGedge *ed = sg->edge;
        fprintf(out_fp[0], "%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", i, asm_i, sg_i, "+-"[sg->is_rev], cname->chr_name[sg->tid], n[e.se].start, n[e.se].end, n[e.up].start, n[e.up].end, n[e.down].start, n[e.down].end);
        // up->se, se->down, up->down
        uint32_t ij1_id = _err_sg_bin_sch_edge(sg, n[e.up].e_site_id, n[e.se].s_site_id, &hit);
        uint32_t ij2_id = _err_sg_bin_sch_edge(sg, n[e.se].e_site_id, n[e.down].s_site_id, &hit);
        uint32_t ej_id = _err_sg_bin_sch_edge(sg, n[e.up].e_site_id, n[e.down].s_site_id, &hit);
        fprintf(out_fp[5], "%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", i, "+-"[sg->is_rev], cname->chr_name[sg->tid], ed[ij1_id].uniq_c, ed[ij2_id].uniq_c, ed[ej_id].uniq_c, ed[ij1_id].multi_c, ed[ij2_id].multi_c, ed[ej_id].multi_c);
    }
    for (i = 0; i < ase->a5ss_n; ++i) {
        int sg_i = ase->a5ss[i].sg_i, asm_i = ase->a5ss[i].asm_i;
        SG *sg = sg_g->SG[sg_i]; A5SS_t e = ase->a5ss[i]; SGnode *n = sg->node; SGedge *ed = sg->edge;
        fprintf(out_fp[1], "%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", i, asm_i, sg_i, "+-"[sg->is_rev], cname->chr_name[sg->tid], n[e.shor].start, n[e.shor].end, n[e.lon].start, n[e.lon].end, n[e.down].start, n[e.down].end);
        // lon -> down, shor -> down
        uint32_t lon_id, shor_id;
        if (sg->is_rev == 0) { // forward
            lon_id = _err_sg_bin_sch_edge(sg, n[e.lon].e_site_id, n[e.down].s_site_id, &hit); 
            shor_id = _err_sg_bin_sch_edge(sg, n[e.shor].e_site_id, n[e.down].s_site_id, &hit); 
        } else { // reverse
            lon_id = _err_sg_bin_sch_edge(sg, n[e.down].e_site_id, n[e.lon].s_site_id, &hit); 
            shor_id = _err_sg_bin_sch_edge(sg, n[e.down].e_site_id, n[e.shor].s_site_id, &hit); 
        }
        fprintf(out_fp[6], "%d\t%c\t%s\t%d\t%d\t%d\t%d\n", i, "+-"[sg->is_rev], cname->chr_name[sg->tid], ed[lon_id].uniq_c, ed[shor_id].uniq_c, ed[lon_id].multi_c, ed[shor_id].multi_c);
    }
    for (i = 0; i < ase->a3ss_n; ++i) {
        int sg_i = ase->a3ss[i].sg_i, asm_i = ase->a3ss[i].asm_i;
        SG *sg = sg_g->SG[sg_i]; A3SS_t e = ase->a3ss[i]; SGnode *n = sg->node; SGedge *ed = sg->edge;
        fprintf(out_fp[2], "%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", i, asm_i, sg_i, "+-"[sg->is_rev], cname->chr_name[sg->tid], n[e.up].start, n[e.up].end, n[e.lon].start, n[e.lon].end, n[e.shor].start, n[e.shor].end);
        // up -> lon, up -> shor
        uint32_t lon_id, shor_id;
        if (sg->is_rev == 0) { // forward
            lon_id = _err_sg_bin_sch_edge(sg, n[e.up].e_site_id, n[e.lon].s_site_id, &hit); 
            shor_id = _err_sg_bin_sch_edge(sg, n[e.up].e_site_id, n[e.shor].s_site_id, &hit); 
        } else { // reverse
            lon_id = _err_sg_bin_sch_edge(sg, n[e.lon].e_site_id, n[e.up].s_site_id, &hit); 
            shor_id = _err_sg_bin_sch_edge(sg, n[e.shor].e_site_id, n[e.up].s_site_id, &hit); 
        }
        fprintf(out_fp[7], "%d\t%c\t%s\t%d\t%d\t%d\t%d\n", i, "+-"[sg->is_rev], cname->chr_name[sg->tid], ed[lon_id].uniq_c, ed[shor_id].uniq_c, ed[lon_id].multi_c, ed[shor_id].multi_c);
    }
    for (i = 0; i < ase->mxe_n; ++i) {
        int sg_i = ase->mxe[i].sg_i, asm_i = ase->mxe[i].asm_i;
        SG *sg = sg_g->SG[sg_i]; MXE_t e = ase->mxe[i]; SGnode *n = sg->node; SGedge *ed = sg->edge;
        fprintf(out_fp[3], "%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", i, asm_i, sg_i, "+-"[sg->is_rev], cname->chr_name[sg->tid], n[e.up].start, n[e.up].end, n[e.fir].start, n[e.fir].end, n[e.sec].start, n[e.sec].end, n[e.down].start, n[e.down].end);
        // up -> fir, fir -> down
        // up -> sec, sec -> down
        uint32_t fir1_id = _err_sg_bin_sch_edge(sg, n[e.up].e_site_id, n[e.fir].s_site_id, &hit);
        uint32_t fir2_id = _err_sg_bin_sch_edge(sg, n[e.fir].e_site_id, n[e.down].s_site_id, &hit);
        uint32_t sec1_id = _err_sg_bin_sch_edge(sg, n[e.up].e_site_id, n[e.sec].s_site_id, &hit);
        uint32_t sec2_id = _err_sg_bin_sch_edge(sg, n[e.sec].e_site_id, n[e.down].s_site_id, &hit);
        fprintf(out_fp[8], "%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", i, "+-"[sg->is_rev], cname->chr_name[sg->tid], ed[fir1_id].uniq_c, ed[fir2_id].uniq_c, ed[sec1_id].uniq_c, ed[sec2_id].uniq_c, ed[fir1_id].multi_c, ed[fir2_id].multi_c, ed[sec1_id].multi_c, ed[sec2_id].multi_c);
    }
    for (i = 0; i < ase->ri_n; ++i) {
        int sg_i = ase->ri[i].sg_i, asm_i = ase->ri[i].asm_i;
        SG *sg = sg_g->SG[sg_i]; RI_t e = ase->ri[i]; SGnode *n = sg->node; SGedge *ed = sg->edge;
        fprintf(out_fp[4], "%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\n", i, asm_i, sg_i, "+-"[sg->is_rev], cname->chr_name[sg->tid], n[e.up].start, n[e.up].end, n[e.down].start, n[e.down].end);
        // exon_cnt, up -> down
        uint32_t sj_id = _err_sg_bin_sch_edge(sg, n[e.up].e_site_id, n[e.down].s_site_id, &hit);
        fprintf(out_fp[9], "%d\t%c\t%s\t%d\t%d\t%d\t%d\n", i, "+-"[sg->is_rev], cname->chr_name[sg->tid], n[e.in].uniq_c-n[e.up].uniq_c-n[e.down].uniq_c, ed[sj_id].uniq_c, n[e.in].multi_c-n[e.up].multi_c-n[e.down].multi_c, ed[sj_id].multi_c);
    }

    for (i = 0; i < out_n; ++i) {
        free(out_fn[i]); err_fclose(out_fp[i]);
    } free(out_fn); free(out_fp);
}

int asm2ase(SG_group *sg_g, SGasm_group *asm_g, ASE_t *ase)
{
    int i;
    print_format_time(stderr); err_printf("[%s] generating alternative splice-events from alternative splice-modules ...\n", __func__);
    for (i = 0; i < asm_g->sg_asm_n; ++i) {
        SGasm *sg_asm = asm_g->sg_asm[i];
        int sg_i = sg_asm->SG_id; SG *sg = sg_g->SG[sg_i]; 
        SGnode *node = sg->node;
        asm2se(sg, sg_asm, ase, i, sg_i);       // SE
        asm2a5ss(sg, sg_asm, ase, i, sg_i);     // A5SS 
        asm2a3ss(sg, sg_asm, ase, i, sg_i);     // A3SS
        asm2mxe(sg, sg_asm, ase, i, sg_i);      // MXE
        asm2ri(node, sg_asm, ase, i, sg_i);     // RI
    }
    print_format_time(stderr); err_printf("[%s] generating alternative splice-events from alternative splice-modules done\n", __func__);
    return ase->se_n+ase->a5ss_n+ase->a3ss_n+ase->mxe_n+ase->ri_n;
}

const struct option se_long_opt [] = {
    { "novel-sj", 0, NULL, 'n' },
    { "novel-com", 0, NULL, 'N' },
    { "sj-file", 0, NULL, 's' },
    { "use-multi", 0, NULL, 'm' },
    { "output", 1, NULL, 'o' },

    { 0, 0, 0, 0 }
};

int pred_ase(int argc, char *argv[])
{
    // same to pred_asm START
    int c, no_novel_sj=1, no_novel_com=1, BAM_format=1, use_multi=0; char out_fn[1024]="";
    while ((c = getopt_long(argc, argv, "nNmso:", se_long_opt, NULL)) >= 0) {
        switch (c) {
            case 'n': no_novel_sj=0, no_novel_com=0; break;
            case 'N': no_novel_com = 0; break;
            case 's': BAM_format = 0; break;
            case 'm': use_multi = 1; break;
            case 'o': strcpy(out_fn, optarg); break;
            default: err_printf("Error: unknown option: %s.\n", optarg);
                     return pred_se_usage();
        }
    }
    if (argc - optind != 3) return pred_se_usage();

    gzFile genome_fp = gzopen(argv[optind], "r");
    if (genome_fp == NULL) { err_fatal(__func__, "Can not open genome file. %s\n", argv[optind]); }
    chr_name_t *cname = chr_name_init();
    // set cname if input is BAM
    samFile *in; bam_hdr_t *h; bam1_t *b;
    if (BAM_format) {
        if ((in = sam_open(argv[optind+2], "rb")) == NULL) err_fatal_core(__func__, "Cannot open \"%s\"\n", argv[optind+2]);
        if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", argv[optind+2]);
        bam_set_cname(h, cname);
    }
    // build splice-graph with GTF
    SG_group *sg_g;
    FILE *gtf_fp = xopen(argv[optind+1], "r");
    sg_g = construct_SpliceGraph(gtf_fp, cname);
    err_fclose(gtf_fp); chr_name_free(cname);

    // get splice-junction
    int sj_n, sj_m; sj_t *sj_group;
    if (BAM_format) { // based on .bam file
        b = bam_init1(); 
        sj_group = (sj_t*)_err_malloc(10000 * sizeof(sj_t)); sj_m = 10000;
        sj_n = bam2sj_core(in, h, b, genome_fp, &sj_group, sj_m, use_multi);
        bam_destroy1(b); bam_hdr_destroy(h); sam_close(in);
    } else  { // based on .sj file
        FILE *sj_fp = xopen(argv[optind+2], "r");
        sj_group = (sj_t*)_err_malloc(10000 * sizeof(sj_t)); sj_m = 10000;
        sj_n = read_sj_group(sj_fp, sg_g->cname, &sj_group, sj_m);
        err_fclose(sj_fp);
    }
    // predict splice-graph with GTF-based splice-graph and splice-junciton
    SG_group *sr_sg_g = predict_SpliceGraph(*sg_g, sj_group, sj_n, no_novel_sj, no_novel_com);
    sg_free_group(sg_g);

    // generate ASM with short-read splice-graph
    SGasm_group *asm_g = gen_asm(sr_sg_g);

    // calculate number of reads falling into exon-body
    if (BAM_format) {
        samFile *in; bam_hdr_t *h; bam1_t *b;
        in = sam_open(argv[optind+2], "rb"); h = sam_hdr_read(in); b = bam_init1(); 
        cal_asm_exon_cnt(sr_sg_g, in, h, b);
        bam_destroy1(b); bam_hdr_destroy(h); sam_close(in);
    }
    // same to pred_asm END

    ASE_t *ase = ase_init();
    asm2ase(sr_sg_g, asm_g, ase);

    // output
    asm_output(argv[optind+2], out_fn, sr_sg_g, asm_g);
    ase_output(argv[optind+2], out_fn, sr_sg_g, ase);

    sg_free_group(sr_sg_g); sg_free_asm_group(asm_g);
    gzclose(genome_fp); ase_free(ase);
    return 0;
}
