#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include "gtf.h"
#include "build_sg.h"
#include "pred_sg.h"

extern const char PROG[20];
int pred_asm_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s asm [option] <in.gtf> <in.sj> > out.asm\n\n", PROG);
    err_printf("Options:\n\n");
    err_printf("         -n --novel-sj             allow novel splice-junction in the ASM. [False]\n");
    err_printf("         -f --gtf-format  [STR]    file format of GTF input: plain GTF file(GTF) or binary splice-graph file(SG). [SG]\n");
	err_printf("\n");
	return 1;
}

/*****************************
 *       generate ASM        *
 *****************************/
SGasm *sg_init_asm(uint32_t sg_id, uint32_t v_start, uint32_t v_end)
{
    SGasm *sg_asm = (SGasm*)_err_malloc(sizeof(SGasm));
    sg_asm->SG_id = sg_id;
    sg_asm->v_start = v_start; sg_asm->v_end = v_end;
    sg_asm->node_n = 0; sg_asm->node_m = 1;
    sg_asm->node_id = (uint32_t*)_err_malloc(sizeof(uint32_t));
    sg_asm->edge_n = 0; sg_asm->edge_m = 1;
    sg_asm->edge_id = (uint32_t*)_err_malloc(sizeof(uint32_t));
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
    free(asm_g->sg_asm);
    free(asm_g);
}

// XXX sort edge by interval
/*int sg_max_edge(SG *sg)
{
    int i, j;
    for (i = 1; i < sg->edge_n; ++i) {
        for (j = 0; j < i; ++j) {
            if ([j] > [i])
                [i] is NOT max
        }
    }
    return 0;
}*/

void cal_cand_node(SG sg, uint32_t **entry, uint32_t **exit, int *entry_n, int *exit_n)
{
    int i, n1=0, n2=0;
    for (i = 0; i < sg.node_n; ++i) {
        if (sg.node[i].next_n > 1) n1++;
        if (sg.node[i].pre_n > 1) n2++;
    }
    *entry_n = n1, *exit_n = n2;
    *entry = (uint32_t*)_err_malloc(n1 * sizeof(uint32_t));
    *exit = (uint32_t*)_err_malloc(n2 * sizeof(uint32_t));

    n1 = 0, n2 = 0;
    for (i = 0; i < sg.node_n; ++i) {
        if (sg.node[i].next_n > 1) (*entry)[n1++] = sg.node[i].node_id;
        if (sg.node[i].pre_n > 1) (*exit)[n2++] = sg.node[i].node_id;;
    }
}

void sg_update_asm_edge(SG sg, SGasm *sg_asm, uint32_t pre_id, uint32_t cur_id)
{
    int hit = 0;
    if (sg.node[pre_id].e.end == 0 || sg.node[cur_id].e.start == CHR_MAX_END) return;
    uint32_t pre_site_id = sg_bin_sch_site(sg.don_site, sg.don_site_n, sg.node[pre_id].e.end+1, &hit); if (hit==0) err_fatal_core(__func__, "Can not hit site: (%d).(3)\n", sg.node[pre_id].e.end+1);
    uint32_t cur_site_id = sg_bin_sch_site(sg.acc_site, sg.acc_site_n, sg.node[cur_id].e.start-1, &hit); if (hit==0) err_fatal_simple("Can not hit site.(3)\n");
    uint32_t edge_i = sg_bin_sch_edge(sg, pre_site_id, cur_site_id, &hit); 
    if (hit == 0) err_fatal_core(__func__, "Can not hit edge.(%d,%d) (3)\n", sg.node[pre_id].e.end, sg.node[cur_id].e.start);
    _insert(edge_i, sg_asm->edge_id, sg_asm->edge_n, sg_asm->edge_m, uint32_t)
}

void sg_update_asm(SG sg, SGasm *sg_asm, uint32_t pre_id, uint32_t cur_id)
{
    _insert(cur_id, sg_asm->node_id, sg_asm->node_n, sg_asm->node_m, uint32_t)
    sg_update_asm_edge(sg, sg_asm, pre_id, cur_id);
}


void sub_splice_graph(SG sg, SGasm *sg_asm, uint32_t cur_id, uint32_t e_id)
{
    //if (sg.node[cur_id].next_n == 0) err_fatal_simple("No next node!\n");

    if (cur_id == e_id) return;
    int i;
    for (i = 0; i < sg.node[cur_id].next_n; ++i) {
        if (sg.node[cur_id].next_id[i] == e_id) {
            sg_update_asm_edge(sg, sg_asm, cur_id, e_id);
            continue;
        } else {
            sg_update_asm(sg, sg_asm, cur_id, sg.node[cur_id].next_id[i]); 
            sub_splice_graph(sg, sg_asm, sg.node[cur_id].next_id[i], e_id);
        }
    }
}

int sg_asm_group_add(SGasm_group *asm_g, SGasm *sg_asm)
{
    if (asm_g->sg_asm_n == asm_g->sg_asm_m) sg_realloc_asm_group(asm_g);
    SGasm *a = asm_g->sg_asm[asm_g->sg_asm_n];

    a->SG_id = sg_asm->SG_id;
    a->v_start = sg_asm->v_start; a->v_end = sg_asm->v_end;
    if (sg_asm->node_n > a->node_m) a->node_id = (uint32_t*)_err_realloc(a->node_id, sg_asm->node_n * sizeof(uint32_t));
    if (sg_asm->edge_n > a->edge_m) a->edge_id = (uint32_t*)_err_realloc(a->edge_id, sg_asm->edge_n * sizeof(uint32_t));
    a->node_n = a->node_m = sg_asm->node_n; a->edge_n = a->edge_m = sg_asm->edge_n;
    int i;
    for (i = 0; i < sg_asm->node_n; ++i) a->node_id[i] = sg_asm->node_id[i];
    for (i = 0; i < sg_asm->edge_n; ++i) a->edge_id[i] = sg_asm->edge_id[i];

    asm_g->sg_asm_n++;
    return asm_g->sg_asm_n;
}

SGasm_group *gen_ASM(SG_group sg_g)
{
    int entry_n, exit_n; uint32_t *entry, *exit;
    int sg_i;
    SGasm_group *asm_g = sg_init_asm_group();
    for (sg_i = 0; sg_i < sg_g.SG_n; ++sg_i) {
        SG sg = *(sg_g.SG[sg_i]);
        cal_cand_node(sg, &entry, &exit, &entry_n, &exit_n);
        if (entry_n == 0 || exit_n == 0) goto END;

        int i, j, hit;
        for (i = 0; i < entry_n; ++i) {
            hit = 0;
            for (j = 0; j < exit_n; ++j) {
                int post_domn_n = sg.node[entry[i]].post_domn_n;
                int pre_domn_n = sg.node[exit[j]].pre_domn_n;
                if (post_domn_n > 1 && pre_domn_n > 1 
                        && sg.node[entry[i]].post_domn[1] == exit[j] 
                        && sg.node[exit[j]].pre_domn[1] == entry[i]) {
                    SGasm *sg_asm = sg_init_asm(sg_i, entry[i], exit[j]);
                    sub_splice_graph(sg, sg_asm, entry[i], exit[j]);
                    sg_asm_group_add(asm_g, sg_asm);
                    //asm_sg = rm_max_edge(asm_sg);
                    //get_ASM(asm_sg, asm_g);
                    sg_free_asm(sg_asm);
                    hit = 1; break;
                }
            }
            if (hit) continue;
        }
END: free(entry); free(exit);
    }
    return asm_g;
}

/*****************************/
const struct option asm_long_opt [] = {
    { "novel-sj", 0, NULL, 'n' },
    { "gtf-format", 1, NULL, 'f' },

    { 0, 0, 0, 0 }
};

int pred_asm(int argc, char *argv[])
{
    int c;
    int no_novel_sj=1, SG_format=1;
	while ((c = getopt_long(argc, argv, "nf:", asm_long_opt, NULL)) >= 0) {
        switch (c) {
            case 'n': no_novel_sj = 0; break;
            case 'f': if (strcmp(optarg, "GTF")==0) SG_format = 0; 
                      else if (strcmp(optarg, "SG")==0) SG_format = 1; 
                      else return pred_asm_usage();
                  break;
            default: err_printf("Error: unknown option: %s.\n", optarg);
                     return pred_asm_usage();
        }
    }
    
    if (argc - optind != 2) return pred_asm_usage();

    FILE *sj_fp;

    sj_fp = xopen(argv[optind+1], "r");

    chr_name_t *cname = chr_name_init();
    // build splice-graph with GTF
    SG_group *sg_g;
    if (SG_format) sg_g = sg_restore(argv[optind]);
    else  {
        FILE *gtf_fp = xopen(argv[optind], "r");
        sg_g = construct_SpliceGraph(gtf_fp, cname);
        err_fclose(gtf_fp);
    }
    // predict splice-graph with GTF-based splice-graph and splice-junciton
    SG_group *sr_sg_g = predict_SpliceGraph(*sg_g, sj_fp, cname, no_novel_sj);

    // predict ASM with short-read splice-graph
    SGasm_group *asm_g = gen_ASM(*sr_sg_g);
    int i;
    printf("%d\n", asm_g->sg_asm_n);
    for (i = 0; i < asm_g->sg_asm_n; ++i) {
        int sg_i = asm_g->sg_asm[i]->SG_id;
        int start, end;
        uint32_t v_s = asm_g->sg_asm[i]->v_start, v_e = asm_g->sg_asm[i]->v_end;
        if (sr_sg_g->SG[sg_i]->node[v_s].e.start == 0) start = sr_sg_g->SG[sg_i]->start-100; else start = sr_sg_g->SG[sg_i]->node[v_s].e.start;
        if (sr_sg_g->SG[sg_i]->node[v_e].e.end == CHR_MAX_END) end = sr_sg_g->SG[sg_i]->end+100; else end = sr_sg_g->SG[sg_i]->node[v_e].e.end;

        printf("%d(%d):\t%c%s: (%d,%d)-(%d,%d) %s:%d-%d\n", i+1, sg_i, "+-"[sr_sg_g->SG[sg_i]->is_rev], cname->chr_name[sr_sg_g->SG[sg_i]->tid], sr_sg_g->SG[sg_i]->node[asm_g->sg_asm[i]->v_start].e.start, sr_sg_g->SG[sg_i]->node[asm_g->sg_asm[i]->v_start].e.end, sr_sg_g->SG[sg_i]->node[asm_g->sg_asm[i]->v_end].e.start,sr_sg_g->SG[sg_i]->node[asm_g->sg_asm[i]->v_end].e.end, cname->chr_name[sr_sg_g->SG[sg_i]->tid], start, end);
    }

    sg_free_group(sg_g); sg_free_group(sr_sg_g); sg_free_asm_group(asm_g);
    err_fclose(sj_fp); chr_name_free(cname);
    return 0;
}
