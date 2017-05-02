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
    err_printf("         -l --only-novel           only output ASM/ASE with novel-junctions. [False]\n");

    err_printf("         -m --use-multi            use both uniq- and multi-mapped reads in the bam input.[False (uniq only)]\n");
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
    sg_asm->node_id = (int*)_err_malloc(sizeof(int));
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
    sg_iso->start = MAX_SITE; sg_iso->end = 0;
    sg_iso->iso_n = 0;
    return sg_iso;
}

SGiso *sg_init_tmp_iso(int asm_id, int sg_id, int v_start, int v_end)
{
    SGiso *sg_iso = (SGiso*)_err_malloc(sizeof(SGiso));
    sg_iso->ASM_id = asm_id; sg_iso->SG_id = sg_id;
    sg_iso->v_start = v_start; sg_iso->v_end = v_end;
    sg_iso->start = MAX_SITE; sg_iso->end = 0;
    sg_iso->iso_n = 0; sg_iso->iso_m = 10;
    sg_iso->node_n = (int*)_err_calloc(10, sizeof(int));
    sg_iso->node_id = (int**)_err_malloc(10 * sizeof(int));
    sg_iso->edge_n = (int*)_err_calloc(10, sizeof(int));
    sg_iso->edge_id = (int**)_err_malloc(10 * sizeof(int));
    return sg_iso;
}

void sg_free_iso(SGiso *sg_iso)
{
    int i;
    for (i = 0; i < sg_iso->iso_n; ++i) {
        free(sg_iso->node_id[i]);
        free(sg_iso->edge_id[i]);
    }
    free(sg_iso->node_id); free(sg_iso->edge_id);
    free(sg_iso->node_n); free(sg_iso->edge_n);
    free(sg_iso);
}

SGiso_group *sg_init_tmp_iso_group(int asm_id, int sg_id, int v_start, int v_end)
{
    SGiso_group *iso_g = (SGiso_group*)_err_malloc(sizeof(SGiso_group));
    iso_g->sg_iso_n = v_end-v_start+1; iso_g->sg_iso_m = v_end-v_start+1;
    iso_g->sg_iso = (SGiso**)_err_malloc(iso_g->sg_iso_m * sizeof(SGiso*));
    int i; for (i = 0; i < iso_g->sg_iso_m; ++i) iso_g->sg_iso[i] = sg_init_tmp_iso(asm_id, sg_id, v_start, v_end);
    int end = v_end-v_start;
    SGiso *s_iso = iso_g->sg_iso[0], *e_iso = iso_g->sg_iso[end];
    s_iso->iso_n = 1;
    s_iso->node_n[0] = 1; iso_g->sg_iso[0]->edge_n[0] = 0;
    s_iso->node_id[0] = (int*)_err_malloc(sizeof(int));
    s_iso->node_id[0][0] = v_start;
    e_iso->iso_n = 1;
    e_iso->node_n[0] = 1; iso_g->sg_iso[end]->edge_n[0] = 0;
    e_iso->node_id[0] = (int*)_err_malloc(sizeof(int));
    e_iso->node_id[0][0] = v_end;
    return iso_g;
}

SGiso_group *sg_init_iso_group(void)
{
    SGiso_group *iso_g = (SGiso_group*)_err_malloc(sizeof(SGiso_group));
    iso_g->sg_iso_n = 0; iso_g->sg_iso_m = 1;
    iso_g->sg_iso = (SGiso**)_err_malloc(sizeof(SGiso*));
    iso_g->sg_iso[0] = sg_init_iso(0, 0, 0, 0);
    return iso_g;
}

SGiso_group *sg_realloc_iso_group(SGiso_group *iso_g)
{
    iso_g->sg_iso_m <<= 1;
    iso_g->sg_iso = (SGiso**)_err_realloc(iso_g->sg_iso, iso_g->sg_iso_m * sizeof(SGiso*));
    int i;
    for (i = iso_g->sg_iso_m >> 1; i < iso_g->sg_iso_m; ++i) iso_g->sg_iso[i] = sg_init_iso(0, 0, 0, 0);
    return iso_g;
}

void sg_free_iso_group(SGiso_group *iso_g)
{
    int i; for (i = 0; i < iso_g->sg_iso_m; ++i) sg_free_iso(iso_g->sg_iso[i]);
    free(iso_g->sg_iso); free(iso_g);
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
    _bin_insert(cur_id, sg_asm->node_id, sg_asm->node_n, sg_asm->node_m, int)
    if (sg->node[cur_id].node_e.start < sg_asm->start) sg_asm->start = sg->node[cur_id].node_e.start;
    if (sg->node[cur_id].node_e.end > sg_asm->end) sg_asm->end = sg->node[cur_id].node_e.end;
    sg->node[cur_id].is_asm = 1;
    sg_update_asm_edge(sg, sg_asm, pre_id, cur_id);
}

void sub_splice_graph(SG *sg, int **node_visit, SGasm *sg_asm, int cur_id, int e_id)
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

void sg_update_iso_edge(SG *sg, SGiso *iso, SGiso *next_iso, int cur_id, int next_id)
{
    int i, j;
    if (sg->node[cur_id].node_e.end == 0 || sg->node[next_id].node_e.start == MAX_SITE) {
        for (i = 0; i < next_iso->iso_n; ++i) {
            // add edge_id
            iso->edge_n[iso->iso_n+i] = next_iso->edge_n[i];
            iso->edge_id[iso->iso_n+i] = (int*)_err_malloc((next_iso->edge_n[i]) * sizeof(int));
            for (j = 0; j < next_iso->edge_n[i]; ++j)
                iso->edge_id[iso->iso_n+i][j] = next_iso->edge_id[i][j];
        }
    } else {
        int cur_site_id = _err_sg_bin_sch_site(sg->don_site, sg->don_site_n, sg->node[cur_id].node_e.end+1);
        int next_site_id = _err_sg_bin_sch_site(sg->acc_site, sg->acc_site_n, sg->node[next_id].node_e.start-1);
        int edge_id = _err_sg_bin_sch_edge(sg, cur_site_id, next_site_id); 

        for (i = 0; i < next_iso->iso_n; ++i) {
            // add edge_id
            iso->edge_n[iso->iso_n+i] = next_iso->edge_n[i]+1;
            iso->edge_id[iso->iso_n+i] = (int*)_err_malloc((next_iso->edge_n[i]+1) * sizeof(int));
            iso->edge_id[iso->iso_n+i][0] = edge_id;
            for (j = 0; j < next_iso->edge_n[i]; ++j)
                iso->edge_id[iso->iso_n+i][j+1] = next_iso->edge_id[i][j];
        }
    }

}

void sg_update_iso(SG *sg, SGiso *iso, SGiso *next_iso, int cur_id, int next_id)
{
    if (iso->iso_n + next_iso->iso_n > iso->iso_m) {
        iso->iso_m = iso->iso_n + next_iso->iso_n;
        iso->node_n = (int*)_err_realloc(iso->node_n, iso->iso_m * sizeof(int));
        iso->node_id = (int**)_err_realloc(iso->node_id, iso->iso_m * sizeof(int*));
        iso->edge_n = (int*)_err_realloc(iso->edge_n, iso->iso_m * sizeof(int));
        iso->edge_id = (int**)_err_realloc(iso->edge_id, iso->iso_m * sizeof(int*));
    }
    // add node
    int i, j;
    for (i = 0; i < next_iso->iso_n; ++i) {
        // add cur_id
        iso->node_n[iso->iso_n+i] = next_iso->node_n[i]+1;
        iso->node_id[iso->iso_n+i] = (int*)_err_malloc((next_iso->node_n[i]+1) * sizeof(int));
        iso->node_id[iso->iso_n+i][0] = cur_id;
        for (j = 0; j < next_iso->node_n[i]; ++j)
            iso->node_id[iso->iso_n+i][j+1] = next_iso->node_id[i][j];
    }
    // add edge
    sg_update_iso_edge(sg, iso, next_iso, cur_id, next_id);
    
    iso->iso_n += next_iso->iso_n;
}

int sg_asm_gen_iso(SG *sg, int **node_visit, SGiso_group *iso_g, int cur_id, int v_id, int e_id)
{
    SGiso *iso = iso_g->sg_iso[cur_id-v_id];
    if ((*node_visit)[cur_id] == 1) {
        return iso->iso_n;
    } else (*node_visit)[cur_id] = 1;
    if (cur_id == e_id) return 1;

    iso->iso_n = 0;
    int iso_n=0, i;
    for (i = 0; i < sg->node[cur_id].next_n; ++i) {
        int next_id = sg->node[cur_id].next_id[i];
        iso_n = sg_asm_gen_iso(sg, node_visit, iso_g, next_id, v_id, e_id);
        SGiso *next_iso = iso_g->sg_iso[next_id-v_id];

        sg_update_iso(sg, iso, next_iso, cur_id, next_id);
    }

    return iso_n;
}

int sg_asm_group_add(SGasm_group *asm_g, SGasm *sg_asm)
{
    if (asm_g->sg_asm_n == asm_g->sg_asm_m) sg_realloc_asm_group(asm_g);
    SGasm *a = asm_g->sg_asm[asm_g->sg_asm_n];

    a->SG_id = sg_asm->SG_id;
    a->v_start = sg_asm->v_start; a->v_end = sg_asm->v_end;
    if (sg_asm->node_n > a->node_m) a->node_id = (int*)_err_realloc(a->node_id, sg_asm->node_n * sizeof(int));
    if (sg_asm->edge_n > a->edge_m) a->edge_id = (int*)_err_realloc(a->edge_id, sg_asm->edge_n * sizeof(int));
    a->node_n = a->node_m = sg_asm->node_n; a->edge_n = a->edge_m = sg_asm->edge_n;
    int i;
    for (i = 0; i < sg_asm->node_n; ++i) a->node_id[i] = sg_asm->node_id[i];
    for (i = 0; i < sg_asm->edge_n; ++i) a->edge_id[i] = sg_asm->edge_id[i];
    a->start = sg_asm->start; a->end = sg_asm->end;

    asm_g->sg_asm_n++;
    return asm_g->sg_asm_n;
}

int sg_iso_group_add(SGiso_group *iso_g, SGiso *sg_iso)
{
    if (iso_g->sg_iso_n == iso_g->sg_iso_m) sg_realloc_iso_group(iso_g);
    SGiso *a = iso_g->sg_iso[iso_g->sg_iso_n];
    int i, j;

    a->ASM_id = sg_iso->ASM_id; a->SG_id = sg_iso->SG_id;
    a->v_start = sg_iso->v_start; a->v_end = sg_iso->v_end;
    a->start = sg_iso->start; a->end = sg_iso->end;
    a->iso_n = sg_iso->iso_n;
    for (i = 0; i < a->iso_n; ++i) {
        a->node_n = (int*)_err_malloc(a->iso_n * sizeof(int));
        a->node_id = (int**)_err_malloc(a->iso_n * sizeof(int*));
        a->node_n[i] = sg_iso->node_n[i];
        for (j = 0; j < a->node_n[i]; ++j) {
            a->node_id[i] = (int*)_err_malloc(a->node_n[i] * sizeof(int));
            a->node_id[i][j] = sg_iso->node_id[i][j];
        }
        a->edge_n = (int*)_err_malloc(a->iso_n * sizeof(int));
        a->edge_id = (int**)_err_malloc(a->iso_n * sizeof(int*));
        a->edge_n[i] = sg_iso->edge_n[i];
        for (j = 0; j < a->edge_n[i]; ++j) {
            a->edge_id[i] = (int*)_err_malloc(a->edge_n[i] * sizeof(int));
            a->edge_id[i][j] = sg_iso->edge_id[i][j];
        }
    }

    iso_g->sg_iso_n++;
    return iso_g->sg_iso_n;
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

int check_novel_iso(SG *sg, SGiso *sg_iso)
{
    int i, j; 
    for (i = 0; i < sg_iso->iso_n; ++i) {
        for (j = 0; j < sg_iso->edge_n[i]; ++j) {
            int eid = sg_iso->edge_id[i][j];
            if (sg->edge[eid].is_anno == 0) return 1;
        }
    }
    return 0;
}

int check_uniq_iso(SG *sg, SGiso *sg_iso)
{
    int i, j;
    for (i = 0; i < sg_iso->iso_n; ++i) {
        for (j = 0; j < sg_iso->node_n[i]; ++j) {
            int eid = sg_iso->edge_id[i][j];
            if (sg->edge[eid].uniq_c == 0) return 0;
        }
    }
    return 1;
}

SGasm_group *gen_asm(SG_group *sg_g, sg_para *sgp)
{
    print_format_time(stderr); err_printf("[%s] generating alternative-splice-module with predicted SJ-SG ...\n", __func__);
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
                    int *node_visit = (int*)_err_calloc(sg->node_n, sizeof(int));
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
    print_format_time(stderr); err_printf("[%s] generating alternative-splice-module with predicted SJ-SG done!\n", __func__);
    return asm_g;
}

SGiso_group *gen_asm_iso(SG_group *sg_g, sg_para *sgp)
{
    print_format_time(stderr); err_printf("[%s] generating candidate isoforms of alternative-splice-module ...\n", __func__);
    int only_novel = sgp->only_novel, use_multi = sgp->use_multi;
    int entry_n, exit_n; int *entry, *exit;
    int sg_i, asm_i=0;
    SGiso_group *iso_g = sg_init_iso_group();

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
                    int *node_visit = (int*)_err_calloc(sg->node_n, sizeof(int));
                    SGiso_group *iso_tmp_g = sg_init_tmp_iso_group(asm_i, sg_i, entry[i], exit[j]);
                    //sub_splice_graph(sg, &node_visit, sg_asm, entry[i], exit[j]);
                    sg_asm_gen_iso(sg, &node_visit, iso_tmp_g, entry[i], entry[i], exit[j]);
                    free(node_visit);
                    SGiso *iso = iso_tmp_g->sg_iso[0];
                    if ((only_novel == 0 || check_novel_iso(sg, iso) == 1)
                    && (use_multi == 1 || check_uniq_iso(sg, iso) == 1)) {
                        sg_iso_group_add(iso_g, iso);
                        asm_i++;
                    }
                    sg_free_iso_group(iso_tmp_g);
                    hit = 1; break;
                }
            }
            if (hit) continue;
        }
        END: free(entry); free(exit);
    }
    print_format_time(stderr); err_printf("[%s] generating candidate isoforms of alternative-splice-module done!\n", __func__);
    return iso_g;
}

int cal_asm_exon_cnt(SG_group *sg_g, samFile *in, bam_hdr_t *h, bam1_t *b)
{
    print_format_time(stderr); err_printf("[%s] calculating read count for AS exons ...\n", __func__);
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
    print_format_time(stderr); err_printf("[%s] calculating read count for AS exons done!\n", __func__);
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

int iso_output(char *in_fn, char *prefix, SG_group *sg_g, SGiso_group *iso_g, sg_para *sgp)
{
    int i, j, out_n=1;
    char suf[3][10] = { ".ISO" };
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
    fprintf(out_fp[0], "ISO_ID\tASM_ID\tSG_ID\tSTRAND\tCHR\tSTART_NODE\tEND_NODE\tTOTAL_NODES_NUM\tUNIQ_READ_COUNT\tMULTI_READ_COUNT\tUCSC_POS\n");
    
    int iso_i = 0;
    for (i = 0; i < iso_g->sg_iso_n; ++i) {
        SGiso *a = iso_g->sg_iso[i];
        int asm_i = a->ASM_id, sg_i = a->SG_id; SG *sg = sg_g->SG[sg_i]; 
        SGnode *node = sg->node;
        int start, end; int v_s = a->v_start, v_e = a->v_end;

        if (node[v_s].node_e.start == 0) start = sg->start-100; else start = node[v_s].node_e.start;
        if (node[v_e].node_e.end == MAX_SITE) end = sg->end+100; else end = node[v_e].node_e.end;
        for (j = 0; j < a->iso_n; ++j) {
            fprintf(out_fp[0], "%d\t%d\t%d\t%c\t%s\t(%d,%d)\t(%d,%d)\t%d\t%d\t%d\t%s:%d-%d\n", iso_i, asm_i, sg_i, "+-"[sg->is_rev], cname->chr_name[sg->tid], node[a->v_start].node_e.start, node[a->v_start].node_e.end, node[a->v_end].node_e.start, node[a->v_end].node_e.end, a->node_n[j], 0, 0, cname->chr_name[sg->tid], start, end);
            iso_i++;
        }
    }
    
    for (i = 0; i < out_n; ++i) {
        free(out_fn[i]); err_fclose(out_fp[i]);
    }
    free(out_fn); free(out_fp);
    return 0;
}

/*****************************/
const struct option asm_long_opt [] = {
    { "novel-sj", 0, NULL, 'n' },
    { "novel-com", 0, NULL, 'N' },
    { "proper-pair", 1, NULL, 'p' },
    { "anchor-len", 1, NULL, 'a' },
    { "intron-len", 1, NULL, 'i' },
    { "genome-file", 1, NULL, 'g' },
    { "use-multi", 0, NULL, 'm' },
    { "output", 1, NULL, 'o' },

    { 0, 0, 0, 0 }
};

int pred_asm(int argc, char *argv[])
{
    int c, i; char out_fn[1024]="", ref_fn[1024]="", *p;
    sg_para *sgp = sg_init_para();
	while ((c = getopt_long(argc, argv, "nNlmMpa:U:A:i:g:G:o:", asm_long_opt, NULL)) >= 0) {
        switch (c) {
            case 'n': sgp->no_novel_sj=0, sgp->no_novel_com=0; break;
            case 'N': sgp->no_novel_com = 0; break;
            case 'l': sgp->only_novel = 1, sgp->no_novel_sj=0, sgp->no_novel_com=0; break;
            case 'm': sgp->use_multi = 1; break;
            case 'M': sgp->merge_out = 1; break;
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
        gzclose(genome_fp); 
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
    SGiso_group **iso_g_rep = (SGiso_group**)_err_malloc(sgp->tol_rep_n * sizeof(SGiso_group*));
    for (i = 0; i < sgp->tol_rep_n; ++i) {
        char *in_name = sgp->in_name[i];
        // get splice-junction
        int sj_n, sj_m; sj_t *sj_group; int ad_n, ad_m; ad_t *ad_group;
        b = bam_init1(); 
        if ((in = sam_open(in_name, "rb")) == NULL) err_fatal_core(__func__, "Cannot open \"%s\"\n", in_name);
        if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", in_name);
        sj_m = 10000; sj_group = (sj_t*)_err_malloc(sj_m * sizeof(sj_t));
        ad_m = 10000; ad_group = (ad_t*)_err_malloc(ad_m * sizeof(ad_t));
        // calculate number of junction-reads and exon-body reads
        //sj_n = bam2cnt_core(in, h, b, seq, seq_n, &sj_group, sj_m, sg_g, sgp);
        sj_n = parse_bam_record(in, h, b, seq, seq_n, &ad_group, &ad_n, ad_m, &sj_group, &sj_n, sj_m, sg_g, sgp);
        bam_destroy1(b); bam_hdr_destroy(h); sam_close(in);

        // update edge weight and add novel edge for GTF-SG
        update_SpliceGraph(sg_g, sj_group, sj_n, sgp);

        free_sj_group(sj_group, sj_n); // free sj

        // generate ASM with short-read splice-graph
        //SGasm_group *asm_g = gen_asm(sg_g, sgp);
        SGiso_group *iso_g = gen_asm_iso(sg_g, sgp);

        free_ad_group(ad_group, ad_n); // free ad
        // output asm
        //asm_output(in_name, out_fn, sg_g, asm_g, sgp);
        iso_output(in_name, out_fn, sg_g, iso_g, sgp);
        //asm_g_rep[i] = asm_g; 
        iso_g_rep[i] = iso_g;
    }
    for (i = 0; i < sgp->tol_rep_n; ++i) {
        //sg_free_asm_group(asm_g_rep[i]); 
        sg_free_iso_group(iso_g_rep[i]); 
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
