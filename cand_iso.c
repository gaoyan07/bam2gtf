#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <pthread.h>
#include "cand_iso.h"
#include "utils.h"
#include "gtf.h"
#include "splice_graph.h"
#include "update_sg.h"
#include "bias_flow.h"
#include "parse_bam.h"
#include "kstring.h"
#include "kdq.h"

KDQ_INIT(int)
#define kdq_gec_t kdq_t(int)

const struct option asm_long_opt [] = {
    { "thread", 1, NULL, 't' },

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
int cand_iso_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s iso [option] <in.gtf> <in.bam/sj>\n\n", PROG);
    err_printf("Note:    for multi-sample and multi-replicate should be this format: \n");
    err_printf("             \"SAM1-REP1,REP2,REP3;SAM2-REP1,REP2,REP3\"\n");
    err_printf("         use \':\' to separate samples, \',\' to separate replicates.\n\n");
    err_printf("Options:\n\n");
    err_printf("         -t --thread               number of threads to use. [1]\n");
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

/*****************************
 *       generate ASM        *
 *****************************/
void cal_cand_node(SG *sg, int **entry, int **exit, int *entry_n, int *exit_n)
{
    int i, n1=0, n2=0;
    for (i = 0; i < sg->node_n; ++i) {
        if (sg->node[i].next_n > 1) n1++;
        if (sg->node[i].pre_n > 1) n2++;
    }
    *entry_n = n1, *exit_n = n2;
    *entry = (int*)_err_malloc(n1 * sizeof(int));
    *exit = (int*)_err_malloc(n2 * sizeof(int));

    n1 = 0, n2 = 0;
    for (i = 0; i < sg->node_n; ++i) {
        if (sg->node[i].next_n > 1) (*entry)[n1++] = sg->node[i].node_id;
        if (sg->node[i].pre_n > 1) (*exit)[n2++] = sg->node[i].node_id;;
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

FILE **iso_output(sg_para *sgp, char *prefix)
{
    // .IsoMatrix && .IsoExon
    int i, out_n = sgp->tot_rep_n+1;
    char mat_suf[20] = { ".IsoMatrix" };
    char exon_suf[20] = { ".IsoExon" };
    //char suff[20] = "";
    //if (sgp->use_multi==1) strcat(suff, ".multi");
    //if (sgp->no_novel_sj==1) strcat(suff, ".anno");
    //if (sgp->only_novel==1) strcat(suff, ".novel");
    char **out_fn = (char**)_err_malloc(sizeof(char*) * out_n);
    if (strlen(prefix) == 0) {
        for (i = 0; i < out_n-1; ++i) {
            out_fn[i] = (char*)_err_malloc(strlen(sgp->in_name[i])+30); strcpy(out_fn[i], sgp->in_name[i]); strcat(out_fn[i], mat_suf);
        }
        out_fn[i] = (char*)_err_malloc(strlen(sgp->in_name[i])+30); strcpy(out_fn[i], sgp->in_name[i]); strcat(out_fn[i], exon_suf);
    } else {
        for (i = 0; i < out_n-1; ++i) {
            out_fn[i] = (char*)_err_malloc(strlen(prefix)+30); strcpy(out_fn[i], prefix); strcat(out_fn[i], mat_suf);
        }
        out_fn[i] = (char*)_err_malloc(strlen(prefix)+30); strcpy(out_fn[i], prefix); strcat(out_fn[i], exon_suf);
    }

    FILE **out_fp = (FILE**)_err_malloc(sizeof(FILE*) * out_n);
    for (i = 0; i < out_n; ++i) out_fp[i] = xopen(out_fn[i], "w");
    for (i = 0; i < out_n; ++i) free(out_fn[i]); free(out_fn);

    return out_fp;
}

typedef struct {
    int tid;
    SG_group *sg_g;
    sg_para *sgp;
    FILE **out_fp;
} asm_iso_aux_t;

int REP_I;
pthread_rwlock_t RWLOCK;


int is_cmptb_read_iso(read_exon_map *read_map, cmptb_map_t *iso_exon_map) {
    int start = read_map->map_s, si = read_map->map_si, end = read_map->map_e, ei = read_map->map_ei;
    int i;
    // set 0 for start
    if (((iso_exon_map[start] << si) >> si) != read_map->map[0]) return 0;
    for (i = start+1; i < end; ++i) {
        // compare
        if (iso_exon_map[i] != read_map->map[i-start]) return 0;
    }
    // set 0 for end
    if (((iso_exon_map[end] >> ei) << ei) != read_map->map[end-start]) return 0;
    return 1;
}

// read: exon#0, #2, #4
// map:  10101000
int gen_read_exon_map(read_exon_map *map, gec_t *exon_id, gec_t exon_n) {
    if (exon_n < 1) err_fatal_core(__func__, "Error: exon_n = %d\n", exon_n);
    int i;
    map->map_s = exon_id[0] >> MAP_STEP_N;        // map_s/e: index in iso_exon_map
    map->map_si = (~exon_id[0]) & MAP_STEP_M; 
    map->map_si = exon_id[0] & MAP_STEP_M;

    map->map_e = exon_id[exon_n-1] >> MAP_STEP_N;
    map->map_ei = (~exon_id[exon_n-1]) & MAP_STEP_M;

    map->map = (cmptb_map_t*)_err_calloc((map->map_e-map->map_s+1), sizeof(cmptb_map_t));

    map->map[0] |= (0x1ULL << map->map_si);
    map->map[map->map_e-map->map_s] |= (0x1ULL << map->map_ei);

    for (i = 1; i < exon_n-1; ++i) {
        map->map[(exon_id[i] >> MAP_STEP_N) - map->map_s] |= (0x1ULL << (~exon_id[i] & MAP_STEP_M));
    }
    return 1;
}

cmptb_map_t *gen_iso_exon_map(gec_t *exon_id, gec_t exon_n, int map_n) {
    if (exon_n < 1) err_fatal_core(__func__, "Error: exon_n = %d\n", exon_n);
    cmptb_map_t *map = (cmptb_map_t*)_err_calloc(map_n, sizeof(cmptb_map_t));
    int i;
    for (i = 0; i < exon_n; ++i)
        map[exon_id[i] >> MAP_STEP_N] |= (0x1ULL << (~exon_id[i] & MAP_STEP_M));
    return map;
}

void insert_iso_exon_map(cmptb_map_t **iso_map, int *iso_i, int map_n, cmptb_map_t *iso_m) {
    int i, j, mis=0;
    for (i = 0; i < *iso_i; ++i) {
        mis = 0;
        for (j = 0; j < map_n; ++j) {
            if (iso_map[i][j] != iso_m[j]) {
                mis = 1; break;
            }
        }
        if (mis == 0) {
            return;
        }
    }                           
    for (i = 0; i < map_n; ++i) {
        iso_map[(*iso_i)++][i] = iso_m[i];          
    }
}
/*int cand_iso_core(SG_group *sg_g, sg_para *sgp, char *out_fn)
{
    FILE **out_fp = iso_output(sgp, out_fn);
    err_func_format_printf(__func__, "generating candidate isoforms of alternative-splice-module ...\n");
    int i, t = sgp->n_threads;
    if (t < 1) t = 1;
    asm_iso_aux_t *aux = (asm_iso_aux_t*)_err_calloc(t, sizeof(asm_iso_aux_t));
    for (i = 0; i < t; ++i) {
        aux[i].tid = i;
        aux[i].sg_g = sg_g;
        aux[i].sgp = sgp;
        aux[i].out_fp = out_fp;
    }

    REP_I = 0;
    pthread_t *tid = (pthread_t*)_err_calloc(t, sizeof(pthread_t));
    pthread_attr_t attr; pthread_attr_init(&attr); pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for (i = 0; i < t; ++i) {
        pthread_create(&tid[i], &attr, thread_generate_ASMISO, aux+i);
    }
    for (i = 0; i < t; ++i) pthread_join(tid[i], 0);

    free(aux); free(tid);
    err_func_format_printf(__func__, "generating candidate isoforms of alternative-splice-module done!\n");
    return 0;
}*/
void read_exon_map_free(read_exon_map *m) {
    free(m->map); free(m);
}

void read_exon_map_copy(read_exon_map *m, read_exon_map *m1) {
    m->weight = m1->weight;
    m->map_s = m1->map_s; m->map_si = m1->map_si;
    m->map_e = m1->map_e; m->map_ei = m1->map_ei;
    int i;
    m->map = (cmptb_map_t*)_err_malloc((m->map_e-m->map_s+1) * sizeof(cmptb_map_t));
    for (i = m->map_s; i <= m->map_e; ++i) m->map[i-m->map_s] = m1->map[i-m->map_s];
}

void exon_id_copy(gec_t **dest, gec_t *dn, gec_t *dm, gec_t *src, gec_t sn) {
    if (*dm < sn) {
        *dest = (gec_t*)_err_realloc(*dest, sn * sizeof(gec_t));
        *dm = sn;
    }
    *dn = sn;
    int i;
    for (i = 0; i < sn; ++i)
        (*dest)[i] = src[i];
}

// sort by 1st exon, 2nd exon, ...
int read_exon_map_comp(read_exon_map *m1, read_exon_map *m2) {
    if (m1->map_s != m2->map_s) return (int)(m1->map_s-m2->map_e);
    int i, n = MIN_OF_TWO(m1->map_e, m2->map_e) - m1->map_s;

    for (i = 0; i <= n; ++i) {
        if (m1->map[i] > m2->map[i]) return 1;
        else if (m1->map[i] < m2->map[i]) return -1;
    }
    return (int) (m1->map_e - m2->map_e);
}

void read_exon_map_insert(read_exon_map **m, int *bundle_n, int *bundle_m, read_exon_map *m1) {
    int i, r;
    for (i = *bundle_n-1; i >= 0; --i) {
        r = read_exon_map_comp((*m)+i, m1);
        if (r == 0) {
            (*m)[i].weight++; return;
        } else if (r < 0) break;
    }
    if (*bundle_n == *bundle_m) {
        (*bundle_m) <<= 1;
        (*m) = (read_exon_map*)_err_realloc(m, *bundle_m * sizeof(read_exon_map));
    }
    // insert m1 to (*m)[i+1]
    if (i+1 <= *bundle_n-1) memmove((*m) + i+1 + 1, (*m) + i+1, (*bundle_n - i-1) * sizeof(read_exon_map));
    read_exon_map_copy((*m)+i+1, m1);
    (*bundle_n)++;
}

inline void update_weight_matrix(double **wei_matrix, gec_t *exon_id, gec_t exon_n) {
    int i;
    for (i = 0; i < exon_n - 1; ++i) ++wei_matrix[exon_id[i]][exon_id[i+1]];
}

read_exon_map *bam_sgnode_cmptb(ad_t *ad, SG *sg, gec_t **exon_id, gec_t *exon_n, gec_t *exon_m, uint8_t *cmptb) {
    gec_t node_id; *exon_n = 0; 
    read_exon_map *m = NULL; *cmptb = 0;

    int ad_i, hit, s_site_i, e_site_i, s_node_i, e_node_i;
    int start, end;
    SGsite *dsite = sg->don_site, *asite = sg->acc_site; int dn = sg->don_site_n, an = sg->acc_site_n;
    SGnode *node = sg->node;

    for (ad_i = 0; ad_i < ad->intv_n; ++ad_i) {
        // cal s_n_i
        if (ad_i != 0) {
            start = ad->intr_end[ad_i-1];
            s_site_i = sg_bin_sch_site(asite, an, start, &hit);
            if (hit == 0) goto map_END;
            s_node_i = asite[s_site_i].exon_id[0];
        } else {
            start = ad->start;
            s_site_i = sg_bin_sch_site(asite, an, start, &hit);
            if (hit == 0) {
                if (s_site_i == 0) goto map_END;
                else s_site_i--;
            }
            s_node_i = asite[s_site_i].exon_id[0];
        }
        // cal e_n_i
        if (ad_i != ad->intv_n-1) {
            end = ad->exon_end[ad_i];
            e_site_i = sg_bin_sch_site(dsite, dn, end, &hit);
            if (hit == 0) goto map_END;
            e_node_i = dsite[e_site_i].exon_id[0];
        } else {
            end = ad->end;
            e_site_i = sg_bin_sch_site(dsite, dn, end, &hit);
            if (hit == 0) {
                if (e_site_i == dn) goto map_END;
            }
            e_node_i = dsite[e_site_i].exon_id[0];
        }
        // push exon_id
        if (node[s_node_i].start > start || node[e_node_i].end < end) goto map_END;
        for (node_id = s_node_i; node_id < e_node_i; ++node_id) {
            if (node[node_id].end + 1 != node[node_id+1].start) goto map_END;
            else _sim_insert(node_id, (*exon_id), (*exon_n), (*exon_m), gec_t)
        }
        _sim_insert(e_node_i, (*exon_id), (*exon_n), (*exon_m), gec_t)
    }

    m = (read_exon_map*)_err_malloc(sizeof(read_exon_map)); 
    *cmptb = 1;
    gen_read_exon_map(m, (*exon_id), (*exon_n));
map_END:
    return m;
}

read_exon_map *bam_sg_cmptb(bam_aux_t *bam_aux, double **wei_matrix, int *b_n, SG *sg, sg_para *sgp) {
    samFile *in = bam_aux->in; hts_itr_t *itr = bam_aux->itr; bam1_t *b = bam_aux->b;
    ad_t *ad, *last_ad; gec_t *exon_id, *last_exon_id, exon_n, exon_m, last_exon_n, last_exon_m; uint8_t cmptb, last_cmptb;
    int i, r, bundle_n = 0, bundle_m = 10000;//, N = sg->node_n;
    read_exon_map *m = (read_exon_map*)_err_malloc(bundle_m * sizeof(read_exon_map));
    wei_matrix = (double**)_err_malloc(sg->node_n * sizeof(double*));
    for (i = 0; i < sg->node_n; ++i) wei_matrix[i] = (double*)_err_calloc(sg->node_n, sizeof(double));

    ad = (ad_t*)_err_malloc(sizeof(ad_t)); 
    last_ad = (ad_t*)_err_malloc(sizeof(ad_t));last_ad->start = 0; 
    exon_id = (gec_t*)_err_malloc(4 * sizeof(gec_t)); exon_m = 4, exon_n = 0;
    last_exon_id = (gec_t*)_err_malloc(4 * sizeof(gec_t)); last_exon_m = 4, last_exon_n = 0;
    last_cmptb = 0;
    while ((r = sam_itr_next(in, itr, b)) >= 0) {
        // parse bam record
        if (parse_bam_record1(b, ad, sgp) <= 0) continue;
        if (ad_sim_comp(ad, last_ad) == 0) {
            if (last_cmptb) {
                m[bundle_n-1].weight++;
                // update edge weight and novel edge
                update_weight_matrix(wei_matrix, last_exon_id, last_exon_n);
            }
            continue;
        }
        read_exon_map *m1 = bam_sgnode_cmptb(ad, sg, &exon_id, &exon_n, &exon_m, &cmptb);
        // insert new bundle
        if (cmptb) {
            read_exon_map_insert(&m, &bundle_n, &bundle_m, m1);
            // update edge weight and novel edge
            update_weight_matrix(wei_matrix, exon_id, exon_n);
            add_novel_sg_edge(sg, exon_id, exon_n, sgp);
        }
        read_exon_map_free(m1);
        last_cmptb = cmptb; ad_copy(last_ad, ad); exon_id_copy(&last_exon_id, &last_exon_n, &last_exon_m, exon_id, exon_n);
    }
    if (r < -1) err_func_format_printf("BAM file error. \"%s\"", bam_aux->fn);

    *b_n = bundle_n;
    free_ad_group(ad, 1); free_ad_group(last_ad, 1);
    free(exon_id), free(last_exon_id);
    return m;
}

void read_iso_cmptb(read_exon_map **read_map, cmptb_map_t **iso_map, int rep_n, char *out_fn) {
}

void gen_cand_iso(SG *sg, read_exon_map **M, double ***W,  int rep_n, sg_para *sgp, char *out_fn) {
    int entry_n, exit_n; int *entry, *exit;

    cal_pre_domn(sg), cal_post_domn(sg);
    cal_cand_node(sg, &entry, &exit, &entry_n, &exit_n);
    if (entry_n == 0 || exit_n == 0) goto END;

    int i, j;
    for (i = 0; i < entry_n; ++i) {
        for (j = 0; j < exit_n; ++j) {
            int post_domn_n = sg->node[entry[i]].post_domn_n;
            int pre_domn_n = sg->node[exit[j]].pre_domn_n;
            if (post_domn_n > 1 && pre_domn_n > 1 
                    && sg->node[entry[i]].post_domn[1] == exit[j] 
                    && sg->node[exit[j]].pre_domn[1] == entry[i]) {

                int iso_n;
                cmptb_map_t **iso_map = bias_flow_gen_cand_iso(sg, W, entry[i], exit[j], rep_n, &iso_n, sgp);
                // cal read-iso cmptb matrix
                read_iso_cmptb(M, iso_map, rep_n, out_fn);

                break;
            }
        }
    }
END: free(entry);
}

int cand_iso_core(SG_group *sg_g, sg_para *sgp, bam_aux_t **bam_aux, char *out_fn)
{
    int sg_i, rep_i; SG *sg; 
    char reg[1024]; 
    read_exon_map **M = (read_exon_map**)_err_malloc(sgp->tot_rep_n * sizeof(read_exon_map*));
    double ***W = (double***)_err_malloc(sgp->tot_rep_n * sizeof(double**));
    int *bundle_n = (int*)_err_malloc(sgp->tot_rep_n * sizeof(int));

    for (sg_i = 0; sg_i < sg_g->SG_n; ++sg_i) {
        sg = sg_g->SG[sg_i];

        // 0. read bam bundle for each gene
        for (rep_i = 0; rep_i < sgp->tot_rep_n; ++rep_i) {
            sprintf(reg, "%s:%d-%d", sg_g->cname->chr_name[sg->tid], sg->start, sg->end);
            bam_aux[rep_i]->itr = sam_itr_querys(bam_aux[rep_i]->idx, bam_aux[rep_i]->h, reg);
            // 1. generate read-exon compatible array
            // 2. update SG with bamBundle
            // OR (optional) only known transcript
            M[rep_i] = bam_sg_cmptb(bam_aux[rep_i], W[rep_i], bundle_n+rep_i, sg, sgp);
        }
        // 3. flow network decomposition
        gen_cand_iso(sg, M, W, sgp->tot_rep_n, sgp, out_fn);
    }
    return 0;
}

/*****************************/
int cand_iso(int argc, char *argv[])
{
    int c, i; char out_fn[1024]="", ref_fn[1024]="", *p;
    sg_para *sgp = sg_init_para();
	while ((c = getopt_long(argc, argv, "t:nNpa:i:g:w:lme:c:C:o:M:U:A:", asm_long_opt, NULL)) >= 0) {
        switch (c) {
            case 't': sgp->n_threads = atoi(optarg); break;
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
                      if (*p != 0) sgp->anchor_len[1] = strtol(p+1, &p, 10); else return cand_iso_usage();
                      if (*p != 0) sgp->anchor_len[2] = strtol(p+1, &p, 10); else return cand_iso_usage();
                      if (*p != 0) sgp->anchor_len[3] = strtol(p+1, &p, 10); else return cand_iso_usage();
                      if (*p != 0) sgp->anchor_len[4] = strtol(p+1, &p, 10); else return cand_iso_usage();
                      break;
            case 'U': sgp->uniq_min[0] = strtol(optarg, &p, 10);
                      if (*p != 0) sgp->uniq_min[1] = strtol(p+1, &p, 10); else return cand_iso_usage();
                      if (*p != 0) sgp->uniq_min[2] = strtol(p+1, &p, 10); else return cand_iso_usage();
                      if (*p != 0) sgp->uniq_min[3] = strtol(p+1, &p, 10); else return cand_iso_usage();
                      if (*p != 0) sgp->uniq_min[4] = strtol(p+1, &p, 10); else return cand_iso_usage();
                      break; 
            case 'A': sgp->all_min[0] = strtol(optarg, &p, 10);
                      if (*p != 0) sgp->all_min[1] = strtol(p+1, &p, 10); else return cand_iso_usage();
                      if (*p != 0) sgp->all_min[2] = strtol(p+1, &p, 10); else return cand_iso_usage();
                      if (*p != 0) sgp->all_min[3] = strtol(p+1, &p, 10); else return cand_iso_usage();
                      if (*p != 0) sgp->all_min[4] = strtol(p+1, &p, 10); else return cand_iso_usage();
                      break;
            case 'i': sgp->intron_len = atoi(optarg); break;
            case 'w': sgp->rm_edge = 1, sgp->edge_wt = atoi(optarg); break;
            case 'g': strcpy(ref_fn, optarg); break;
            case 'o': strcpy(out_fn, optarg); break;
            default: err_printf("Error: unknown option: %s.\n", optarg); return cand_iso_usage();
        }
    }
    if (argc - optind != 2) return cand_iso_usage();

    int seq_n = 0, seq_m; kseq_t *seq;
    if (strlen(ref_fn) != 0) {
        gzFile genome_fp = gzopen(ref_fn, "r");
        if (genome_fp == NULL) { err_fatal(__func__, "Can not open genome file. %s\n", ref_fn); }
        seq = kseq_load_genome(genome_fp, &seq_n, &seq_m);
        err_gzclose(genome_fp); 
    }
    // 0. parse input bam file names
    bam_aux_t **bam_aux = sg_par_input(sgp, argv[optind+1]);
    if (sgp->tot_rep_n <= 0) return cand_iso_usage();
    
    // 1. set cname --- 1 thread
    chr_name_t *cname = chr_name_init();
    bam_set_cname(bam_aux[0]->h, cname);
    // 2. build splice-graph --- 1 thread (Optional in future, infer exon-intron boundaries by bam records)
    // XXX read gene and build SG one-by-one
    FILE *gtf_fp = xopen(argv[optind], "r");
    // build from GTF file
    // XXX OR from bam file
    SG_group *sg_g = construct_SpliceGraph(gtf_fp, cname);
    err_fclose(gtf_fp); chr_name_free(cname);

    // 3. core process
    cand_iso_core(sg_g, sgp, bam_aux, out_fn);

    sg_free_group(sg_g); sg_free_para(sgp);
    if (seq_n > 0) {
        for (i = 0; i < seq_n; ++i) { free(seq[i].name.s), free(seq[i].seq.s); }
        free(seq);
    }
    for (i = 0; i < sgp->tot_rep_n; ++i) bam_aux_destroy(bam_aux[i]); free(bam_aux);
    return 0;
}
