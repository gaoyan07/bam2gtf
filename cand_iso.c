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

    { "list-input", 0, NULL, 'L' },
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

uint8_t bit_table16[65536];

void gen_bit_table16(void) {
    bit_table16[0] = 0; int i;
    for (i = 1; i != 65536; ++i) {
        bit_table16[i] = (i & 1) + bit_table16[i / 2];
    }
}

int cand_iso_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s iso [option] <in.gtf> <in.bam/sj>\n\n", PROG);
    err_printf("Note:    for multi-sample and multi-replicate should be this format: \n");
    err_printf("             \"SAM1-REP1,REP2,REP3;SAM2-REP1,REP2,REP3\"\n");
    err_printf("         use \':\' to separate samples, \',\' to separate replicates.\n\n");
    err_printf("Options:\n\n");
    err_printf("         -t --thread               number of threads to use. [1]\n");
    err_printf("         -L --list-input           use bam list file as input. [False]\n");
    err_printf("                                   refer to \"example.list\"\n");
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
    err_printf("                                   prefix.IsoMatrix & prefix.IsoExon\n");
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

FILE **iso_output(sg_para *sgp, char *prefix)
{
    // .IsoMatrix && .IsoExon
    int i, out_n = sgp->tot_rep_n+1;
    char mat_suf[20] = { ".IsoMatrix" };
    char exon_suf[20] = { ".IsoExon" };
    char **out_fn = (char**)_err_malloc(sizeof(char*) * out_n);

    if (strlen(prefix) != 0) strcat(prefix, ".");
    for (i = 0; i < out_n-1; ++i) {
        out_fn[i] = (char*)_err_malloc(strlen(prefix)+strlen(sgp->in_name[i])+30); strcpy(out_fn[i], prefix), strcat(out_fn[i], sgp->in_name[i]); strcat(out_fn[i], mat_suf);
    }
    if (strlen(prefix) == 0) strcpy(prefix, sgp->in_name[0]);
    out_fn[i] = (char*)_err_malloc(strlen(prefix)+30); strcpy(out_fn[i], prefix); strcat(out_fn[i], exon_suf);
    /* else {
        for (i = 0; i < out_n-1; ++i) {
            out_fn[i] = (char*)_err_malloc(strlen(prefix)+30); strcpy(out_fn[i], prefix); strcat(out_fn[i], mat_suf);
        }
        out_fn[i] = (char*)_err_malloc(strlen(prefix)+30); strcpy(out_fn[i], prefix); strcat(out_fn[i], exon_suf);
    }*/

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
    int i; cmptb_map_t tmp_s, tmp_e, r=1;
    // set 0 for start/end
    tmp_s = iso_exon_map[start];
    tmp_e = iso_exon_map[end];
    iso_exon_map[start] = (iso_exon_map[start] << si) >> si;
    iso_exon_map[end] = (iso_exon_map[end] >> ei) << ei;
    for (i = start; i <= end; ++i) {
        // compare
        if (iso_exon_map[i] != read_map->map[i-start]) {
            r = 0;
            break;
        }
    }
    iso_exon_map[start] = tmp_s;
    iso_exon_map[end] = tmp_e;
    return r;
}

int is_cmptb_exon_iso(int exon_i, cmptb_map_t *iso_exon_map) {
    return (iso_exon_map[exon_i >> MAP_STEP_N] >> (~exon_i & MAP_STEP_M)) & 1;
}

// read: exon#0, #2, #4
// map:  10101000
int gen_read_exon_map(read_exon_map *map, gec_t *exon_id, gec_t exon_n) {
    if (exon_n < 1) err_fatal_core(__func__, "Error: exon_n = %d\n", exon_n);
    int i;
    map->map_s = exon_id[0] >> MAP_STEP_N;        // map_s/e: index in iso_exon_map
    map->map_si = (~exon_id[0]) & MAP_STEP_M; 

    map->map_e = exon_id[exon_n-1] >> MAP_STEP_N;
    map->map_ei = (~exon_id[exon_n-1]) & MAP_STEP_M;

    map->map = (cmptb_map_t*)_err_calloc((map->map_e-map->map_s+1), sizeof(cmptb_map_t));

    map->map[0] |= (0x1ULL << map->map_si);
    map->map[map->map_e-map->map_s] |= (0x1ULL << map->map_ei);

    map->map_si = exon_id[0] & MAP_STEP_M;
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

void read_exon_map_free(read_exon_map *m) {
    free(m->map); free(m);
}

void read_exon_map_copy(read_exon_map *m, read_exon_map *m1) {
    m->weight = m1->weight; m->rlen = m1->rlen;
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
    for (i = 0; i < sn; ++i) (*dest)[i] = src[i];
}

// sort by 1st exon, 2nd exon, ...
int read_exon_map_comp(read_exon_map *m1, read_exon_map *m2) {
    if (m1->rlen != m2->rlen) return (m1->rlen - m2->rlen); // XXX len
    if (m1->map_s != m2->map_s) return (int)(m2->map_s-m1->map_e);
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
        } else if (r > 0) break;
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

void update_weight_matrix(double **wei_matrix, gec_t *exon_id, gec_t exon_n) {
    int i;
    for (i = 0; i < exon_n - 1; ++i) ++wei_matrix[exon_id[i]][exon_id[i+1]];
}

read_exon_map *bam_sgnode_cmptb(ad_t *ad, SG *sg, gec_t **exon_id, gec_t *exon_n, gec_t *exon_m, uint8_t *cmptb) {
    gec_t node_id; *exon_n = 0; 
    read_exon_map *read_map = NULL; *cmptb = 0;

    int ad_i, hit, s_site_i, e_site_i, s_node_i, e_node_i;
    int start, end;
    SGsite *dsite = sg->don_site, *asite = sg->acc_site; int dn = sg->don_site_n, an = sg->acc_site_n;
    SGnode *node = sg->node;

    for (ad_i = 0; ad_i < ad->intv_n; ++ad_i) {
        // cal s_n_i
        if (ad_i != 0) {
            start = ad->intr_end[ad_i-1]+1;
            s_site_i = sg_bin_sch_site(asite, an, start-1, &hit);
            if (hit == 0) goto map_END;
            s_node_i = asite[s_site_i].exon_id[0];
        } else {
            start = ad->start;
            s_site_i = sg_bin_sch_site(asite, an, start-1, &hit);
            if (hit == 0) {
                if (s_site_i == 0) goto map_END;
                else s_site_i--;
            }
            s_node_i = asite[s_site_i].exon_id[0];
        }
        // cal e_n_i
        if (ad_i != ad->intv_n-1) {
            end = ad->exon_end[ad_i];
            e_site_i = sg_bin_sch_site(dsite, dn, end+1, &hit);
            if (hit == 0) goto map_END;
            e_node_i = dsite[e_site_i].exon_id[0];
        } else {
            end = ad->end;
            e_site_i = sg_bin_sch_site(dsite, dn, end+1, &hit);
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

    read_map = (read_exon_map*)_err_malloc(sizeof(read_exon_map)); 
    gen_read_exon_map(read_map, (*exon_id), (*exon_n));
    read_map->weight = 1; read_map->rlen = ad->rlen;
    *cmptb = 1;
map_END:
    return read_map;
}

/*read_exon_map *bam_sg_cmptb(bam_aux_t *bam_aux, double **wei_matrix, int *b_n, SG *sg, char *reg, sg_para *sgp) {
    if (sg->start == 700237)
        printf("debug");
    char fn[1024] = "./data/test2.sort.bam";
    samFile *in;  hts_idx_t *idx; bam_hdr_t *h; bam1_t *b;
    err_sam_open(in, fn); b = bam_init1(); err_sam_hdr_read(h, in, fn); err_sam_idx_load(idx, in, fn);
    hts_itr_t *itr = sam_itr_querys(idx, h, reg);

    ad_t *ad = ad_init(1), *last_ad = ad_init(1); gec_t *exon_id, *last_exon_id, exon_n, exon_m, last_exon_n, last_exon_m; uint8_t cmptb, last_cmptb;
    int i, r, bundle_n = 0, bundle_m = 10000;//, N = sg->node_n;
    read_exon_map *read_map = (read_exon_map*)_err_malloc(bundle_m * sizeof(read_exon_map));
    for (i = 0; i < sg->node_n; ++i) wei_matrix[i] = (double*)_err_calloc(sg->node_n, sizeof(double));

    exon_id = (gec_t*)_err_malloc(4 * sizeof(gec_t)); exon_m = 4, exon_n = 0;
    last_exon_id = (gec_t*)_err_malloc(4 * sizeof(gec_t)); last_exon_m = 4, last_exon_n = 0;
    last_cmptb = 0;
    while ((r = sam_itr_next(in, itr, b)) >= 0) {
        // parse bam record
        if (parse_bam_record1(b, ad, sgp) <= 0) continue;
        if (ad_sim_comp(ad, last_ad) == 0) {
            if (last_cmptb) {
                read_map[bundle_n-1].weight++;
                // update edge weight and novel edge
                update_weight_matrix(wei_matrix, last_exon_id, last_exon_n);
            }
            continue;
        }
        read_exon_map *m1 = bam_sgnode_cmptb(ad, sg, &exon_id, &exon_n, &exon_m, &cmptb);
        // insert new bundle
        if (cmptb) {
            read_exon_map_insert(&read_map, &bundle_n, &bundle_m, m1);
            // update edge weight and novel edge
            update_weight_matrix(wei_matrix, exon_id, exon_n);
            add_novel_sg_edge(sg, exon_id, exon_n, sgp);
            read_exon_map_free(m1);
        }
        last_cmptb = cmptb; exon_id_copy(&last_exon_id, &last_exon_n, &last_exon_m, exon_id, exon_n);
        ad_copy(last_ad, ad); 
    }
    if (r < -1) err_func_format_printf("BAM file error. \"%s\"", bam_aux->fn);

    *b_n = bundle_n; 
    bam_destroy1(b); hts_itr_destroy(itr); hts_idx_destroy(idx); bam_hdr_destroy(h); sam_close(in);
    free_ad_group(ad, 1); free_ad_group(last_ad, 1);
    free(exon_id), free(last_exon_id);
    return read_map;
}*/

read_exon_map *bam_sg_cmptb(bam_aux_t *bam_aux, double **wei_matrix, int *b_n, SG *sg, char *reg, sg_para *sgp) {
    char *fn = bam_aux->fn;
    samFile *in;  hts_idx_t *idx; bam_hdr_t *h; bam1_t *b;
    err_sam_open(in, fn); b = bam_init1(); err_sam_hdr_read(h, in, fn); err_sam_idx_load(idx, in, fn);
    hts_itr_t *itr = sam_itr_querys(idx, h, reg);

    ad_t *ad = ad_init(1), *last_ad = ad_init(1); gec_t *exon_id, *last_exon_id, exon_n, exon_m, last_exon_n, last_exon_m; uint8_t cmptb, last_cmptb;
    int i, r, bundle_n = 0, bundle_m = 10000;//, N = sg->node_n;
    read_exon_map *read_map = (read_exon_map*)_err_malloc(bundle_m * sizeof(read_exon_map));
    for (i = 0; i < sg->node_n; ++i) wei_matrix[i] = (double*)_err_calloc(sg->node_n, sizeof(double));

    exon_id = (gec_t*)_err_malloc(4 * sizeof(gec_t)); exon_m = 4, exon_n = 0;
    last_exon_id = (gec_t*)_err_malloc(4 * sizeof(gec_t)); last_exon_m = 4, last_exon_n = 0;
    last_cmptb = 0;
    while ((r = sam_read1(in, h, b)) >= 0)  {
        if (parse_bam_record1(b, ad, sgp) <= 0) continue;
        if (ad->end < sg->start) continue;
        else if (ad->start > sg->end) break;

        if (ad_sim_comp(ad, last_ad) == 0) {
            if (last_cmptb) {
                read_map[bundle_n-1].weight++;
                // update edge weight and novel edge
                update_weight_matrix(wei_matrix, last_exon_id, last_exon_n);
            }
            continue;
        }
        read_exon_map *m1 = bam_sgnode_cmptb(ad, sg, &exon_id, &exon_n, &exon_m, &cmptb);
        // insert new bundle
        if (cmptb) {
            read_exon_map_insert(&read_map, &bundle_n, &bundle_m, m1);
            // update edge weight and novel edge
            update_weight_matrix(wei_matrix, exon_id, exon_n);
            add_novel_sg_edge(sg, exon_id, exon_n, sgp);
            read_exon_map_free(m1);
        }
        last_cmptb = cmptb; exon_id_copy(&last_exon_id, &last_exon_n, &last_exon_m, exon_id, exon_n);
        ad_copy(last_ad, ad); 
    }
    if (r < -1) err_func_format_printf("BAM file error. \"%s\"", bam_aux->fn);

    *b_n = bundle_n; 
    bam_destroy1(b); hts_itr_destroy(itr); hts_idx_destroy(idx); bam_hdr_destroy(h); sam_close(in);
    free_ad_group(ad, 1); free_ad_group(last_ad, 1);
    free(exon_id), free(last_exon_id);
    return read_map;
}

cmptb_map_t *cmptb_map_union(cmptb_map_t **iso_map, int iso_n, int map_n) {
    cmptb_map_t *um = (cmptb_map_t*)_err_calloc(map_n, sizeof(cmptb_map_t));
    int i, j;
    for (i = 0; i < iso_n; ++i) {
        for (j = 0; j < map_n; ++j) {
            um[j] |= iso_map[i][j];
        }
    }
    return um;
}

int cmptb_map_exon_cnt(cmptb_map_t *union_map, int map_n) {
    int n = 0, i;
    for (i = 0; i < map_n; ++i) {
        cmptb_map_t m = union_map[i];
        n += bit_table16[m & 0xffff] +
             bit_table16[(m >> 16) & 0xffff] + 
             bit_table16[(m >> 32) & 0xffff] + 
             bit_table16[m >> 48];
    }
    return n;
}

void read_iso_cmptb(SG *sg, char **cname, read_exon_map **read_map, int rep_n, int *bundle_n, cmptb_map_t **iso_map, int iso_n, sg_para *sgp, int *asm_i, int src, int sink) {
    int rep_i, b_i, iso_i; read_exon_map *rm; cmptb_map_t *im;
    // .IsoMatrix
    for (rep_i = 0; rep_i < rep_n; ++rep_i) {
        int bn = bundle_n[rep_i];
        // matrix header: ASM_ID, ReadBundle_NUM, Isoform_NUM
        fprintf(sgp->out_fp[rep_i], "ASM#%d\t%d\t%d\n", *asm_i, bn, iso_n);
        for (b_i = 0; b_i < bn; ++b_i) {
            rm = (read_map[rep_i])+b_i;
            // read length
            fprintf(sgp->out_fp[rep_i], "%d\t%d", read_map[rep_i][b_i].rlen, read_map[rep_i][b_i].weight);
            // read-iso compatible matrix
            for (iso_i = 0; iso_i < iso_n; ++iso_i) {
                im = iso_map[iso_i];
                fprintf(sgp->out_fp[rep_i], "\t%d", is_cmptb_read_iso(rm, im));
            } fprintf(sgp->out_fp[rep_i], "\n");
        }
    }
    // .IsoExon
    int exon_i, exon_n, v_n = 0, map_n = 1 + ((sg->node_n-1) >> MAP_STEP_N);
    SGnode *node = sg->node;
    cmptb_map_t *union_map = cmptb_map_union(iso_map, iso_n, map_n);
    exon_n = cmptb_map_exon_cnt(union_map, map_n);

    if (src == 0) { src++; v_n++; }
    if (sink == sg->node_n-1) { sink--; v_n++; }
    // exon header: ASM_ID, Exon_NUM, Isoform_NUM, strand, CHR_Name
    fprintf(sgp->out_fp[rep_n], "ASM#%d\t%d\t%d\t%c\t%s\n", *asm_i, exon_n-v_n, iso_n, "+-"[sg->is_rev], cname[sg->tid]);
    // exon coordinate pair: start, end
    int *exon_index = (int*)_err_calloc(sink-src+1, sizeof(int)); int _exon_i = 0;
    for (exon_i = src; exon_i <= sink; ++exon_i) {
        if (is_cmptb_exon_iso(exon_i, union_map)) {
            fprintf(sgp->out_fp[rep_n], "%d,%d\t", node[exon_i].start, node[exon_i].end);
            exon_index[exon_i-src] = _exon_i++;
        }
    } fprintf(sgp->out_fp[rep_n], "\n");
    for (iso_i = 0; iso_i < iso_n; ++iso_i) {
        im = iso_map[iso_i];
        exon_n = cmptb_map_exon_cnt(im, map_n);
        fprintf(sgp->out_fp[rep_n], "%d", exon_n-v_n);
        for (exon_i = src; exon_i <= sink; ++exon_i) {
            if (is_cmptb_exon_iso(exon_i, im)) {
                fprintf(sgp->out_fp[rep_n], "\t%d", exon_index[exon_i-src]);
            }
        } fprintf(sgp->out_fp[rep_n], "\n");
    }

    free(union_map); free(exon_index);
    (*asm_i)++;
}

void gen_cand_iso(SG *sg, char **cname, read_exon_map **M, double ***W,  int rep_n, int *bundle_n, sg_para *sgp, int *asm_i) {
    int entry_n, exit_n; int *entry, *exit;

    cal_pre_domn(sg), cal_post_domn(sg);
    cal_cand_node(sg, &entry, &exit, &entry_n, &exit_n);
    if (entry_n == 0 || exit_n == 0) {
        free(entry); free(exit); return;
    }

    cmptb_map_t **iso_map = (cmptb_map_t**)_err_malloc(sgp->iso_cnt_max * sizeof(cmptb_map_t*));
    int i, j;
    for (i = 0; i < entry_n; ++i) {
        for (j = 0; j < exit_n; ++j) {
            int post_domn_n = sg->node[entry[i]].post_domn_n;
            int pre_domn_n = sg->node[exit[j]].pre_domn_n;
            if (post_domn_n > 1 && pre_domn_n > 1 
                    && sg->node[entry[i]].post_domn[1] == exit[j] 
                    && sg->node[exit[j]].pre_domn[1] == entry[i]) {

                int iso_n, k, map_n = ((sg->node_n-1) >> MAP_STEP_N) + 1;

                iso_n = bias_flow_gen_cand_iso(sg, W, entry[i], exit[j], rep_n, iso_map, map_n, sgp);
                // cal read-iso cmptb matrix
                read_iso_cmptb(sg, cname, M, rep_n, bundle_n, iso_map, iso_n, sgp, asm_i, entry[i], exit[j]);
                for (k = 0; k < sgp->iso_cnt_max; ++k) free(iso_map[k]);

                break;
            }
        }
    }
    free(iso_map); free(entry); free(exit);
}

int cand_iso_core(SG_group *sg_g, sg_para *sgp, bam_aux_t **bam_aux)
{
    int sg_i, asm_i, rep_i; SG *sg; 
    read_exon_map **M = (read_exon_map**)_err_malloc(sgp->tot_rep_n * sizeof(read_exon_map*));
    double ***W = (double***)_err_malloc(sgp->tot_rep_n * sizeof(double**));
    int *bundle_n = (int*)_err_malloc(sgp->tot_rep_n * sizeof(int));

    asm_i = 0;
    for (sg_i = 0; sg_i < sg_g->SG_n; ++sg_i) {
        sg = sg_g->SG[sg_i];
        if (sg->node_n - 2 < 3) continue;

        // 0. read bam bundle for each gene
        for (rep_i = 0; rep_i < sgp->tot_rep_n; ++rep_i) {
            // 1. generate read-exon compatible array
            // 2. update SG with bamBundle
            // OR (optional) only known transcript
            char reg[1024]; 
            sprintf(reg, "%s:%d-%d", sg_g->cname->chr_name[sg->tid], sg->start, sg->end);
            W[rep_i] = (double**)_err_calloc(sg->node_n, sizeof(double));
            M[rep_i] = bam_sg_cmptb(bam_aux[rep_i], W[rep_i], bundle_n+rep_i, sg, reg, sgp);
        }
        // 3. flow network decomposition
        gen_cand_iso(sg, sg_g->cname->chr_name, M, W, sgp->tot_rep_n, bundle_n, sgp, &asm_i);
        for (rep_i = 0; rep_i < sgp->tot_rep_n; ++rep_i) {
            int i;
            for (i = 0; i < sg->node_n; ++i) free(W[rep_i][i]); free(W[rep_i]); 
            for (i = 0; i < bundle_n[rep_i]; ++i) {
                free(M[rep_i][i].map);
            } free(M[rep_i]);
        }
    }
    free(bundle_n); free(W); free(M);
    return 0;
}

/*****************************/
int cand_iso(int argc, char *argv[])
{
    int c, i; char ref_fn[1024]="", out_pre[1024]="", *p;
    sg_para *sgp = sg_init_para();
    while ((c = getopt_long(argc, argv, "t:LnNpa:i:g:w:lme:c:C:o:M:U:A:", asm_long_opt, NULL)) >= 0) {
        switch (c) {
            case 't': sgp->n_threads = atoi(optarg); break;
            case 'L': sgp->in_list = 1; break;
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
            case 'o': strcpy(out_pre, optarg); break;
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
    bam_aux_t **bam_aux = (sgp->in_list ? sg_par_input_list(sgp, argv[optind+1]) : sg_par_input(sgp, argv[optind+1]));
    if (sgp->tot_rep_n <= 0) return cand_iso_usage();

    // 1. set cname --- 1 thread
    chr_name_t *cname = chr_name_init();
    bam_set_cname(bam_aux[0]->h, cname);
    // 2. build splice-graph --- 1 thread (Optional in future, infer exon-intron boundaries by bam records)
    FILE *gtf_fp = xopen(argv[optind], "r");
    // build from GTF file // XXX OR from bam file
    SG_group *sg_g = construct_SpliceGraph(gtf_fp, cname);
    err_fclose(gtf_fp); chr_name_free(cname);

    // 3. core process
    sgp->out_fp = iso_output(sgp, out_pre); gen_bit_table16();
    cand_iso_core(sg_g, sgp, bam_aux);

    sg_free_group(sg_g);
    if (seq_n > 0) {
        for (i = 0; i < seq_n; ++i) { free(seq[i].name.s), free(seq[i].seq.s); }
        free(seq);
    }
    for (i = 0; i < sgp->tot_rep_n; ++i) bam_aux_destroy(bam_aux[i]); free(bam_aux);
    sg_free_para(sgp);
    return 0;
}
