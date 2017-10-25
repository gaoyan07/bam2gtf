#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <libgen.h>
#include <string.h>
#include <pthread.h>
#include "iso.h"
#include "utils.h"
#include "gtf.h"
#include "splice_graph.h"
#include "update_sg.h"
#include "iso.h"
#include "parse_bam.h"
#include "kstring.h"

const struct option asm_long_opt [] = {
    { "thread", 1, NULL, 't' },

    { "list-input", 0, NULL, 'L' },
    { "novel-sj", 0, NULL, 'n' },
    { "novel-exon", 0, NULL, 'N' },
    { "un-pair", 0, NULL, 'u' },
    { "anchor-len", 1, NULL, 'a' },
    { "intron-len", 1, NULL, 'i' },
    { "genome-file", 1, NULL, 'g' },

    { "edge-wei", 1, NULL, 'w' },
    { "gtf-event", 0, NULL, 'f' },
    { "only-junc", 0, NULL, 'j' },
    { "only-novel", 0, NULL, 'l' },
    { "use-multi", 0, NULL, 'm' },

    { "fully", 0, NULL, 'F' },
    { "exon-thres", 1, NULL, 'T' },
    { "asm-exon-thres", 1, NULL, 'e' },
    { "iso-cnt-thres", 1, NULL, 'C' },

    { "junc-cnt", 1, NULL, 'c' },
    { "novel-jcnt", 1, NULL, 'v' },

    // output specific module
    { "recursive", 0, NULL, 'r' },
    { "module-type", 1, NULL, 'd' },
    { "exon-num", 1, NULL, 'E' },

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

int cand_asm_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s asm [option] <in.gtf> <in.bam/list> -o out_dir\n\n", PROG);
    err_printf("Note:    for multi-sample and multi-replicate, input should be list file(with -L) or in this format: \n");
    err_printf("             \"SAM1-REP1,REP2,REP3;SAM2-REP1,REP2,REP3\"\n");
    err_printf("         use \':\' to separate samples, \',\' to separate replicates.\n\n");

    err_printf("Options:\n\n");

    err_printf("         -t --thread               number of threads to use. [1]\n");
    err_printf("         -L --list-input           use bam list file as input. Refer to \"example.list\". [False]\n");
    err_printf("\n");

    err_printf("         -f --gtf-event            only output gtf annotation event. [False]\n");
    err_printf("         -n --novel-sj             allow novel splice-junction in the ASM.\n");
    err_printf("                                   Novel splice-junction here means novel combination of known splice-sites. [False]\n");
    err_printf("         -N --novel-exon           allow novel exon in the ASM.\n");
    err_printf("                                   Novel exon here means novel splice sites or novel intron-retation. [False]\n");
    err_printf("         -u --un-pair              set -u to use both proper and unproper paired mapped read. [False] (proper paired only)\n");
    err_printf("         -m --use-multi            use both uniq- and multi-mapped reads in the bam file. [False] (uniq-mapped only)\n");
    err_printf("         -a --anchor-len  [INT]    minimum anchor length for junction read. [%d].\n", ANCHOR_MIN_LEN);
    err_printf("         -i --intron-len  [INT]    minimum intron length for junction read. [%d]\n", INTRON_MIN_LEN);
    err_printf("         -g --genome-file [STR]    genome.fa. Use genome sequence to classify intron-motif. \n");
    err_printf("                                   If no genome file is give, intron-motif will be set as 0(non-canonical) [None]\n");
    err_printf("\n");

    err_printf("         -j --only-junc            only count junction-read count. [False]\n");
    err_printf("         -l --only-novel           only output ASM/ASE with novel-junctions. [False]\n");
    err_printf("\n");

    err_printf("         -F --fully                count reads fully inside the module. [False]\n");
    err_printf("         -T --exon-thres  [INT]    maximum number of exon for ASM to enumerate all possible isoforms.\n");
    err_printf("                                   For larger ASM, heuristic method will be applied. [%d]\n", ISO_EXON_ALL_CNT);
    err_printf("         -e --iso-exon-thres [INT] maximum number of exons for ASM to generate candidate isoforms. [%d]\n", ISO_EXON_MAX); 
    err_printf("         -C --iso-cnt-thres  [INT] maximum number of isoform count to retain ASM. [%d]\n", ISO_CNT_MAX); 
    err_printf("\n");

    err_printf("         -c --junc-cnt    [INT]    minimum number of read count for junction. [%d]\n", JUNC_READ_CNT_MIN); 
    err_printf("         -v --novel-jcnt  [INT]    minimum number of read count for novel-junction. [%d]\n", NOVEL_JUNC_READ_CNT_MIN); 
    err_printf("         -w --edge-wei    [DOU]    during weight-based candidate isoform generation, ignore edge in\n");
    err_printf("                                   splicing-graph whose weight is less than specified value. [%.2f]\n", ISO_EDGE_MIN_WEI);
    err_printf("\n");

    err_printf("         -d --module-type [INT]    only output specific type of alternative splice module, skipped exon,\n");
    err_printf("                                   alternative initial/terminal exon, 4 exon module ... Refer to \"module_type.list\". [0]\n");
    err_printf("         -E --exon-num    [INT]    only output alternative splice module with specific exon number. 0 for all. [0]\n");
    err_printf("         -r --recursive            recursively output alternative splice module which are fully embedded\n");
    err_printf("                                   in a larger module. [False]\n");
    err_printf("\n");

    err_printf("         -o --output      [STR]    directory output .IsoMatrix & .IsoExon\n");
    err_printf("         -G --gtf-out     [STR]    file name of formatted GTF\n");
	err_printf("\n");
	return 1;
}

/*****************************
 *       generate ASM        *
 *****************************/
int none_zore_next(SGnode *node, uint8_t **con_matrix) {
    int i, n=0;
    for (i = 0; i < node->next_n; ++i) {
        if (is_con_matrix(con_matrix, node->node_id, node->next_id[i])) n++;
    }
    return n;
}

int none_zore_pre(SGnode *node, uint8_t **con_matrix) {
    int i, n=0;
    for (i = 0; i < node->pre_n; ++i) {
        if (is_con_matrix(con_matrix, node->pre_id[i], node->node_id)) n++;
    }
    return n;
}

void cal_cand_node(SG *sg, int **entry, int **exit, int *entry_n, int *exit_n, uint8_t **con_matrix)
{
    int i, n1=0, n2=0;
    for (i = 0; i < sg->node_n; ++i) {
        if (none_zore_next(sg->node+i, con_matrix) > 1) n1++;
        if (none_zore_pre(sg->node+i, con_matrix) > 1) n2++;
    }
    *entry_n = n1, *exit_n = n2;
    *entry = (int*)_err_malloc(n1 * sizeof(int));
    *exit = (int*)_err_malloc(n2 * sizeof(int));

    n1 = 0, n2 = 0;
    for (i = 0; i < sg->node_n; ++i) {
        if (none_zore_next(sg->node+i, con_matrix) > 1) (*entry)[n1++] = sg->node[i].node_id;
        if (none_zore_pre(sg->node+i, con_matrix) > 1) (*exit)[n2++] = sg->node[i].node_id;;
    }
}

FILE **iso_output(sg_para *sgp, char *out_dir)
{
    // .IsoMatrix && .IsoExon
    int i, out_n = sgp->tot_rep_n+1;
    char mat_suf[20] = { ".IsoMatrix" };
    char exon_suf[20] = { ".IsoExon" };
    char **out_fn = (char**)_err_malloc(sizeof(char*) * out_n);

    for (i = 0; i < out_n-1; ++i) {
         out_fn[i] = (char*)_err_malloc(strlen(out_dir)+strlen(sgp->in_name[i])+30); 
         strcpy(out_fn[i], out_dir); strcat(out_fn[i], "/"); strcat(out_fn[i], basename(sgp->in_name[i])); strcat(out_fn[i], mat_suf);
    }
    out_fn[i] = (char*)_err_malloc(strlen(out_dir)+strlen(sgp->in_name[0])+30); 
    strcpy(out_fn[i], out_dir); strcat(out_fn[i], "/"); strcat(out_fn[i], basename(sgp->in_name[0])); strcat(out_fn[i], exon_suf);
    FILE **out_fp = (FILE**)_err_malloc(sizeof(FILE*) * out_n);
    for (i = 0; i < out_n; ++i) out_fp[i] = xopen(out_fn[i], "w");
    for (i = 0; i < out_n; ++i) free(out_fn[i]); 
    free(out_fn);

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

read_iso_map *read_iso_map_init(int ri_map_m, int map_n) {
    int r_i;
    read_iso_map *ri_map = (read_iso_map*)_err_calloc(ri_map_m, sizeof(read_iso_map));
    for (r_i = 0; r_i < ri_map_m; r_i++) {
        ri_map[r_i].rlen = 0, ri_map[r_i].weight = 0;
        ri_map[r_i].map = (cmptb_map_t*)_err_calloc(map_n, sizeof(cmptb_map_t));
        ri_map[r_i].map_n = map_n;
    }
    return ri_map;
}

int is_cmptb_read_iso(read_exon_map *read_map, cmptb_map_t *iso_exon_map, int *iso_se) {
    int start = read_map->map_s, read_si = read_map->map_si, end = read_map->map_e, read_ei = read_map->map_ei;
    int i, iso_si=iso_se[0], iso_ei=iso_se[1]; cmptb_map_t iso_tmp_s, iso_tmp_e, read_tmp_s, read_tmp_e;
    // set 0 for start/end
    iso_tmp_s = iso_exon_map[start]; iso_tmp_e = iso_exon_map[end];
    read_tmp_s = read_map->map[0]; read_tmp_e = read_map->map[end-start];

    iso_exon_map[start] = (iso_exon_map[start] << read_si) >> read_si;
    iso_exon_map[end] = (iso_exon_map[end] >> read_ei) << read_ei;
    read_map->map[0] = (read_map->map[0] << iso_si) >> iso_si;
    read_map->map[end-start] = (read_map->map[end-start] >> iso_ei) << iso_ei;

    int r = 0;
    for (i = start; i <= end; ++i) {
        if (iso_exon_map[i] == read_map->map[i-start] && iso_exon_map[i] != 0) {
            r = 1;
        } else if (iso_exon_map[i] != read_map->map[i-start]) {
            r = 0;
            break;
        }
    }
    iso_exon_map[start] = iso_tmp_s; iso_exon_map[end] = iso_tmp_e;
    read_map->map[0] =read_tmp_s; read_map->map[end-start] = read_tmp_e;
    return r;
}

int is_fully_cmptb_read_iso(read_exon_map *read_map, cmptb_map_t *iso_exon_map) { //, int *iso_se) { XXX
    int start = read_map->map_s, read_si = read_map->map_si, end = read_map->map_e, read_ei = read_map->map_ei;
    int i; cmptb_map_t iso_tmp_s, iso_tmp_e;
    // set 0 for start/end
    iso_tmp_s = iso_exon_map[start]; iso_tmp_e = iso_exon_map[end];

    iso_exon_map[start] = (iso_exon_map[start] << read_si) >> read_si;
    iso_exon_map[end] = (iso_exon_map[end] >> read_ei) << read_ei;

    int r = 0;
    for (i = start; i <= end; ++i) {
        // compare
        if (iso_exon_map[i] == read_map->map[i-start] && iso_exon_map[i] != 0) {
            r = 1;
        } else if (iso_exon_map[i] != read_map->map[i-start]) {
            r = 0;
            break;
        }
    }
    iso_exon_map[start] = iso_tmp_s; iso_exon_map[end] = iso_tmp_e;
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

cmptb_map_t *gen_iso_exon_map(gec_t *exon_id, gec_t exon_n, int map_n, int sg_node_n, int *iso_se) {
    if (exon_n < 1) err_fatal_core(__func__, "Error: exon_n = %d\n", exon_n);
    cmptb_map_t *map = (cmptb_map_t*)_err_calloc(map_n, sizeof(cmptb_map_t));
    int i;
    for (i = 0; i < exon_n; ++i) {
        if (exon_id[i] != 0 && exon_id[i] != sg_node_n-1) {
            map[exon_id[i] >> MAP_STEP_N] |= (0x1ULL << (~exon_id[i] & MAP_STEP_M));
        }
    }
    iso_se[0] = (exon_id[0] == 0 ? exon_id[1] : exon_id[0]);
    iso_se[1] = MAP_STEP_M - ((exon_id[exon_n-1] == sg_node_n-1 ? exon_id[exon_n-2] : exon_id[exon_n-1]) & MAP_STEP_M);
    return map;
}

// 1. from iso-nor_exon map iso-split_exon map
// 2. from read-split_exon map to read-iso map
read_iso_map *gen_read_iso_map(read_exon_map *split_read_map, cmptb_map_t **split_iso_map, int **split_iso_se, int iso_n, int map_n, int *zero, uint8_t fully) {
    int iso_i; cmptb_map_t *split_im; int *split_se;
    *zero = 0;
    read_iso_map *rim = (read_iso_map*)_err_malloc(sizeof(read_iso_map));
    rim->map = (cmptb_map_t*)_err_calloc(map_n, sizeof(cmptb_map_t));
    rim->weight = split_read_map->weight; rim->rlen = split_read_map->rlen;
    rim->map_n = map_n;
    for (iso_i = 0; iso_i < iso_n; ++iso_i) {
        split_im = split_iso_map[iso_i]; split_se = split_iso_se[iso_i];
        if (fully) {
            if (is_fully_cmptb_read_iso(split_read_map, split_im)) { //, se)) { XXX
                rim->map[iso_i >> MAP_STEP_N] |= (0x1ULL << (~iso_i & MAP_STEP_M));
                *zero = 1;
            }
        } else {
            if (is_cmptb_read_iso(split_read_map, split_im, split_se)) {
                rim->map[iso_i >> MAP_STEP_N] |= (0x1ULL << (~iso_i & MAP_STEP_M));
                *zero = 1;
            }
        }
    }
    return rim;
}

void insert_iso_exon_map(cmptb_map_t **iso_map, int **iso_se, int *iso_i, int map_n, cmptb_map_t *iso_m, int *se) {
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
        iso_map[*iso_i][i] = iso_m[i];          
    } iso_se[*iso_i][0] = se[0]; iso_se[*iso_i][1] = se[1];
    (*iso_i)++;
}

int comp_read_iso_map(read_iso_map *m1, read_iso_map *m2) {
    int i;
    for (i = 0; i < m1->map_n; ++i) {
        if (m1->map[i] > m2->map[i]) return 1;
        else if (m1->map[i] < m2->map[i]) return -1;
    }
    return (m1->rlen - m2->rlen);
}

int map_bin_sch_read_iso(read_iso_map *ri_map, int ri_n, read_iso_map *ri, int *hit) {
    *hit = 0;
    int left = 0, right = ri_n-1, mid;
    if (right == -1) return 0;
    while (left <= right) {
        mid = (left + right) >> 1;
        int comp_res = comp_read_iso_map(ri_map+mid, ri);
        if (comp_res == 0) {
            *hit = 1; return mid;
        } else if (comp_res > 0) {
            if (mid != 0) {
                if (comp_read_iso_map(ri_map+mid-1, ri) < 0) {
                    return mid;
                } else right = mid-1;
            } else {
                return mid;
            }
        } else left = mid+1;
    }
    return ri_n;
}

void read_iso_map_move(read_iso_map *dest, read_iso_map *src, int l) {
    int i, j;
    for (i = l-1; i >= 0; --i) {
        dest[i].rlen = src[i].rlen; src[i].rlen = 0;
        dest[i].weight = src[i].weight; src[i].weight = 0;
        dest[i].map_n = src[i].map_n; src[i].map_n = 0;
        for (j = 0; j < dest[i].map_n; ++j) {
            dest[i].map[j] = src[i].map[j];
        }
    }
}

void read_iso_map_copy(read_iso_map *dest, read_iso_map *src) {
    dest->weight = src->weight;
    dest->rlen = src->rlen;
    dest->map_n = src->map_n;
    int i; 
    for (i = 0; i < dest->map_n; ++i) {
        dest->map[i] = src->map[i];
    }
}

void read_iso_map_free(read_iso_map *m) {
    free(m->map); free(m);
}

void insert_read_iso_map(read_iso_map **ri_map, int *ri_n, int *ri_m, read_iso_map *ri) {
    int i, hit; int map_i = map_bin_sch_read_iso(*ri_map, *ri_n, ri, &hit);
    if (hit == 0) {
        if (*ri_n == *ri_m) {
            _realloc((*ri_map), *ri_m, read_iso_map)
            for (i = *ri_n; i < *ri_m; ++i) {
                (*ri_map)[i].map = (cmptb_map_t*)_err_calloc(ri->map_n, sizeof(cmptb_map_t));
            }
        }
        if (map_i <= *ri_n-1)
            read_iso_map_move((*ri_map)+map_i+1, (*ri_map)+map_i, *ri_n-map_i);
        read_iso_map_copy((*ri_map)+map_i, ri);
        (*ri_n)++;
    } else {
        (*ri_map)[map_i].weight += ri->weight;
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
    if (m1->rlen != m2->rlen) return 1; //(m1->rlen - m2->rlen);
    if (m1->map_s != m2->map_s) return (int)(m1->map_s-m2->map_s);
    if (m1->map_si != m2->map_si) return (int)(m1->map_si-m2->map_si);

    int i, n = MIN_OF_TWO(m1->map_e, m2->map_e) - m1->map_s;

    for (i = 0; i <= n; ++i) {
        if (m1->map[i] != m2->map[i]) return 1;
    }
    return (int) abs(m1->map_e - m2->map_e);
}

void read_exon_map_insert(read_exon_map **m, int *bundle_n, int *last_bi, int *bundle_m, read_exon_map *m1) {
    int i, r;
    for (i = *bundle_n-1; i >= 0; --i) {
        r = read_exon_map_comp((*m)+i, m1);
        if (r == 0) {
            *last_bi = i;
            (*m)[i].weight++; return;
        } else if (r < 0) break;
    }
    if (*bundle_n == *bundle_m) {
        (*bundle_m) <<= 1;
        (*m) = (read_exon_map*)_err_realloc(*m, *bundle_m * sizeof(read_exon_map));
    }
    i = *bundle_n;
    // insert m1 to (*m)[i+1]
    // if (i+1 <= *bundle_n-1) memmove((*m) + i+1 + 1, (*m) + i+1, (*bundle_n - i-1) * sizeof(read_exon_map));
    read_exon_map_copy((*m)+i, m1);
    *last_bi = i;
    (*bundle_n)++;
}

void update_weight_matrix(double **wei_matrix, gec_t *exon_id, gec_t exon_n) {
    int i;
    for (i = 0; i < exon_n - 1; ++i) ++wei_matrix[exon_id[i]][exon_id[i+1]];
}

// sg: split_sg
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

// XXX here nor_wei_matrix is NOT correct, cause it can NOT handle ambiguous junction
void convert_split_nor_wei(double **split_wei_matrix, double ** nor_wei_matrix, SG *split_sg, SG *nor_sg) {
    int i, j, m, n;
    for (i = 1; i < split_sg->node_n-2; ++i) {
        for (j = i + 1; j < split_sg->node_n-1; ++j) {
            if (split_wei_matrix[i][j] > 0) {
                // find corresponding (m,n) for nor_sg
                int don_site = split_sg->node[i].end, acc_site = split_sg->node[j].start;
                for (m = 1; m < nor_sg->node_n; ++m) {
                    if (nor_sg->node[m].end == don_site) {
                        for (n = 0; n < nor_sg->node[m].next_n; ++n) {
                            gec_t next_id = nor_sg->node[m].next_id[n];
                            if (nor_sg->node[next_id].start == acc_site) {
                                nor_wei_matrix[m][next_id] = split_wei_matrix[i][j];
                            }
                        }
                    }
                    if (nor_sg->node[m].start >= don_site) break;
                }
            }
        }
    }
}

// @return: read-split_exon's compatible map
// @outpu: wei_matrix: nor_nodes' weight matrix
read_exon_map *bam_sg_cmptb(bam_aux_t *bam_aux, hts_itr_t *itr, double **nor_wei_matrix, int *b_n, SG *split_sg, SG *nor_sg, sg_para *sgp) {
    samFile *in = bam_aux->in;  bam1_t *b = bam_aux->b;

    ad_t *ad = ad_init(1), *last_ad = ad_init(1); gec_t *exon_id, *last_exon_id, exon_n, exon_m, last_exon_n, last_exon_m; uint8_t cmptb, last_cmptb;
    int i, r, bundle_n = 0, last_bi=-1, bundle_m = 10000;//, N = sg->node_n;

    read_exon_map *split_read_map = (read_exon_map*)_err_malloc(bundle_m * sizeof(read_exon_map));
    double **split_wei_matrix = (double**)_err_calloc(split_sg->node_n, sizeof(double*));
    for (i = 0; i < split_sg->node_n; ++i) split_wei_matrix[i] = (double*)_err_calloc(split_sg->node_n, sizeof(double));
    for (i = 0; i < nor_sg->node_n; ++i) nor_wei_matrix[i] = (double*)_err_calloc(nor_sg->node_n, sizeof(double));

    exon_id = (gec_t*)_err_malloc(4 * sizeof(gec_t)); exon_m = 4, exon_n = 0;
    last_exon_id = (gec_t*)_err_malloc(4 * sizeof(gec_t)); last_exon_m = 4, last_exon_n = 0;
    last_cmptb = 0;
    while ((r = sam_itr_next(in, itr, b)) >= 0)  {
        if (parse_bam_record1(b, ad, sgp) <= 0) continue;
        if (ad->end < split_sg->start) continue;
        if (sgp->only_junc && ad->intv_n <= 1) continue;
        else if (ad->start > split_sg->end) break;

        if (ad_sim_comp(ad, last_ad) == 0) {
            if (last_cmptb) {
                split_read_map[last_bi].weight++;
                // update edge weight and novel edge
                update_weight_matrix(split_wei_matrix, last_exon_id, last_exon_n);
            }
            continue;
        }
        // m1: read-split_exon map
        read_exon_map *m1 = bam_sgnode_cmptb(ad, split_sg, &exon_id, &exon_n, &exon_m, &cmptb);
        // insert new bundle
        if (cmptb) {
            read_exon_map_insert(&split_read_map, &bundle_n, &last_bi, &bundle_m, m1);
            // update edge weight and novel edge
            update_weight_matrix(split_wei_matrix, exon_id, exon_n);
            read_exon_map_free(m1);
        }
        last_cmptb = cmptb; exon_id_copy(&last_exon_id, &last_exon_n, &last_exon_m, exon_id, exon_n);
        ad_copy(last_ad, ad); 
    }
    if (r < -1) err_func_format_printf("BAM file error. \"%s\"", bam_aux->fn);

    // convert split_node weight to nor_node weight
    convert_split_nor_wei(split_wei_matrix, nor_wei_matrix, split_sg, nor_sg);
    if (sgp->no_novel_sj == 0) {
        add_novel_sg_edge(split_sg, split_wei_matrix, sgp);
        add_novel_sg_edge(nor_sg, nor_wei_matrix, sgp);
    }

    *b_n = bundle_n; 
    free_ad_group(ad, 1); free_ad_group(last_ad, 1);
    free(exon_id), free(last_exon_id);
    for (i = 0; i < split_sg->node_n; ++i) free(split_wei_matrix[i]); free(split_wei_matrix);

    return split_read_map;
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

int check_module_type(SG *sg, cmptb_map_t *union_map, cmptb_map_t **iso_map, int src, int sink, int node_n, int iso_n, int type, int *whole_exon_id, int **exon_id, int *iso_exon_n, int *exon_index) {
    int i, j, exon_i, iso_i;
    SGnode *node = sg->node; cmptb_map_t *im;

    int _exon_i = 0;
    for (exon_i = src; exon_i <= sink; ++exon_i) {
        if (is_cmptb_exon_iso(exon_i, union_map)) {
            exon_index[exon_i-src] = _exon_i;
            whole_exon_id[_exon_i++] = exon_i;
        }
    }
    for (iso_i = 0; iso_i < iso_n; ++iso_i) {
        im = iso_map[iso_i];
        _exon_i = 0;
        for (exon_i = src; exon_i <= sink; ++exon_i) {
            if (is_cmptb_exon_iso(exon_i, im)) {
                exon_id[iso_i][_exon_i++] = exon_i;
            }
        }
        iso_exon_n[iso_i] = _exon_i;
    }

    if (type == ALL_TYPE) return 1; // all
    if (type == SE_TYPE) {
        if (iso_n != 2) return 0;
        // inclusion isoform
        if (iso_exon_n[0] != node_n) return 0;
        // first junction
        if (node[exon_id[0][0]].end+1 == node[exon_id[0][1]].start) return 0;
        // no middle junction
        for (i = 1; i < node_n-2; ++i) {
            if (node[exon_id[0][i]].end+1 != node[exon_id[0][i+1]].start) return 0;
        }
        // second junction
        if (node[exon_id[0][node_n-1]].start-1 == node[exon_id[0][node_n-2]].end) return 0;
        // exclusion isoform
        if (iso_exon_n[1] != 2 || exon_id[1][0] != src || exon_id[1][1] != sink) return 0;
    } else if ((type == A5SS_TYPE && sg->is_rev == 0) || (type == A3SS_TYPE && sg->is_rev == 1)) { // TODO
        if (node_n != 4 || iso_n != 2 || iso_exon_n[1] != 3 || iso_exon_n[2] != 3) return 0;
        // same first node
        for (i = 0; i < iso_n; ++i) {
            if (exon_id[i][0] != src) return 0;
        }
        // same last node
        for (i = 0; i < iso_n; ++i) {
            if (exon_id[i][iso_exon_n[i]-1] != sink) return 0;
        } 
        if (node[exon_id[0][1]].start != node[exon_id[1][1]].start) return 0;
    } else if ((type == A3SS_TYPE && sg->is_rev == 0) || (type == A5SS_TYPE && sg->is_rev == 1)) { // TODO
        // if (node_n != 3 || iso_n != 2) return 0;
        // same first node
        for (i = 0; i < iso_n; ++i) {
            if (exon_id[i][0] != src) return 0;
        }
        // same last node
        for (i = 0; i < iso_n; ++i) {
            if (exon_id[i][iso_exon_n[i]-1] != sink) return 0;
        } 
        if (node[exon_id[0][1]].end != node[exon_id[1][1]].end) return 0;
    } else if (type == MXE_TYPE) {
        if (iso_n != 2 || node_n != 4) return 0; // TODO node_n == 4
        if (exon_id[0][0] != src || exon_id[1][0] != src) return 0;
        if (iso_exon_n[0] != 3 || iso_exon_n[1] != 3) return 0;
        if (exon_id[0][2] != sink || exon_id[1][2] != sink) return 0;
        // 4 junction
        for (i = 0; i < node_n-1; ++i) {
            if (node[whole_exon_id[i]].end+1 == node[whole_exon_id[i+1]].start) return 0;
        }        
    } else if (type == RI_TYPE) {
        if (node_n != 5 || iso_n != 2 || iso_exon_n[0] != 4 || iso_exon_n[1] != 3) return 0;
        // same first node
        for (i = 0; i < iso_n; ++i) {
            if (exon_id[i][0] != src) return 0;
        }
        // same last node
        for (i = 0; i < iso_n; ++i) {
            if (exon_id[i][iso_exon_n[i]-1] != sink) return 0;
        } 
        // consecutive first and second ... exon
        if (node[exon_id[0][1]].start != node[exon_id[1][1]].start) return 0;
        if (node[exon_id[0][2]].end != node[exon_id[1][1]].end) return 0;
    } else if (type == _2SE_TYPE) {
        if (iso_n != 2) return 0;
        // inclusion isoform
        if (iso_exon_n[0] != node_n) return 0;
        // first junction
        if (node[exon_id[0][0]].end+1 == node[exon_id[0][1]].start) return 0;
        // 1 middle junction
        int hit=0;
        for (i = 1; i < node_n-2; ++i) {
            if (node[exon_id[0][i]].end+1 != node[exon_id[0][i+1]].start) {
                hit++;
            }
        }
        if (hit != 1) return 0;
        // second junction
        if (node[exon_id[0][node_n-1]].start-1 == node[exon_id[0][node_n-2]].end) return 0;
        // exclusion isoform
        if (iso_exon_n[1] != 2 || exon_id[1][0] != src || exon_id[1][1] != sink) return 0;
    } else if ((type == AIE_TYPE && sg->is_rev == 0) || (type == ATE_TYPE && sg->is_rev == 1)) {
        // different first node
        for (i = 0; i < iso_n-1; ++i) {
            for (j = i+1; j < iso_n; ++j) {
                if (exon_id[i][0] == exon_id[j][0]) return 0;
                if (node[exon_id[i][1]].start != node[exon_id[j][1]].start) return 0;
            }
        }
        // same last node
        for (i = 0; i < iso_n-1; ++i) {
            for (j = i+1; j < iso_n; ++j) {
                if (exon_id[i][iso_exon_n[i]-1] != exon_id[j][iso_exon_n[j]-1]) return 0;
            }
        }
    } else if ((type == ATE_TYPE && sg->is_rev == 0) || (type == AIE_TYPE && sg->is_rev == 1)) {
        // same first node
        for (i = 0; i < iso_n-1; ++i) {
            for (j = i+1; j < iso_n; ++j) {
                if (exon_id[i][0] != exon_id[j][0]) return 0;
            }
        }
        // different last node
        int j;
        for (i = 0; i < iso_n-1; ++i) {
            for (j = i+1; j < iso_n; ++j) {
                if (exon_id[i][iso_exon_n[i]-1] == exon_id[j][iso_exon_n[j]-1]) return 0;
                if (node[exon_id[i][iso_exon_n[i]-2]].end != node[exon_id[j][iso_exon_n[j]-2]].end) return 0;
            }
        }
    }

    return 1;
}

void convert_nor_split_iso_path(gec_t *nor_iso_path, gec_t nor_path_idx, gec_t *split_iso_path, gec_t *split_path_idx, SG *nor_sg, SG *split_sg) {
    int i, j;
    *split_path_idx = 0;
    for (i = 0; i < nor_path_idx; ++i) {
        if (i == 0 && nor_iso_path[i] == 0) {
            split_iso_path[(*split_path_idx)++] = 0;
        } else if (i == nor_path_idx-1 && nor_iso_path[i] == nor_sg->node_n-1) {
            split_iso_path[(*split_path_idx)++] = split_sg->node_n-1;
        } else { // normal exon
            gec_t id = nor_iso_path[i];
            int start = nor_sg->node[id].start, end = nor_sg->node[id].end;
            for (j = 1; j < split_sg->node_n-1; ++j) {
                if (split_sg->node[j].start >= start && split_sg->node[j].end <= end) {
                    split_iso_path[(*split_path_idx)++] = j;
                } else if (split_sg->node[j].end > end) {
                    break;
                }
            }
        }
    }
}

void gen_iso_map(SG *nor_sg, SG *split_sg, gec_t **nor_iso_path, gec_t *nor_path_idx, int iso_n, cmptb_map_t **nor_iso_map, int **nor_iso_se, cmptb_map_t **split_iso_map, int **split_iso_se) {
    int nor_map_n = ((nor_sg->node_n-1) >> MAP_STEP_N) + 1, split_map_n = ((split_sg->node_n-1) >> MAP_STEP_N) + 1;
    int i, iso_i = 0;

    gec_t *split_iso_path = (gec_t*)_err_malloc(split_sg->node_n * sizeof(gec_t)); gec_t split_path_idx;
    for (i = 0; i < iso_n; ++i) {
        nor_iso_map[i] = (cmptb_map_t*)_err_calloc(nor_map_n, sizeof(cmptb_map_t));
        split_iso_map[i] = (cmptb_map_t*)_err_calloc(split_map_n, sizeof(cmptb_map_t));
        nor_iso_se[i] = (int*)_err_calloc(2, sizeof(int));
        split_iso_se[i] = (int*)_err_calloc(2, sizeof(int));

        // for nor_iso_map
        // convert nor_iso_path to nor_iso_map
        int *se = (int*)_err_malloc(2 * sizeof(int));
        cmptb_map_t *iso_m = gen_iso_exon_map(nor_iso_path[i], nor_path_idx[i], nor_map_n, nor_sg->node_n, se);
        insert_iso_exon_map(nor_iso_map, nor_iso_se, &iso_i, nor_map_n, iso_m, se); iso_i--;
        free(iso_m);
        // for split_iso_map
        // first.  convert nor_iso_path to split_iso_path
        // second. convert split_iso_path to split_iso_map
        convert_nor_split_iso_path(nor_iso_path[i], nor_path_idx[i], split_iso_path, &split_path_idx, nor_sg, split_sg);
        iso_m = gen_iso_exon_map(split_iso_path, split_path_idx, split_map_n, split_sg->node_n, se);
        insert_iso_exon_map(split_iso_map, split_iso_se, &iso_i, split_map_n, iso_m, se);
        free(iso_m); free(se);
    }
    free(split_iso_path);
}

// generate read-iso compatible matrix
void read_iso_cmptb(SG *nor_sg, SG *split_sg, uint8_t **nor_con_matrix, char **cname, read_exon_map **split_read_map, int rep_n, int *bundle_n, gec_t **nor_iso_path, gec_t *nor_path_idx,/*cmptb_map_t **nor_iso_map, int **nor_iso_se,*/ int iso_n, sg_para *sgp, int *asm_i, int src, int sink) {
    int **nor_iso_se = (int**)_err_malloc(iso_n * sizeof(int*));
    int **split_iso_se = (int**)_err_malloc(iso_n * sizeof(int*));
    cmptb_map_t **nor_iso_map = (cmptb_map_t**)_err_malloc(iso_n * sizeof(cmptb_map_t*));
    cmptb_map_t **split_iso_map = (cmptb_map_t**)_err_malloc(iso_n * sizeof(cmptb_map_t*));
    gen_iso_map(nor_sg, split_sg, nor_iso_path, nor_path_idx, iso_n, nor_iso_map, nor_iso_se, split_iso_map, split_iso_se);

    int rep_i, b_i, r_i, iso_i; read_exon_map *split_rm;
    // .IsoExon
    int i, exon_i, nor_exon_n, map_n = 1 + ((nor_sg->node_n-1) >> MAP_STEP_N);
    SGnode *nor_node = nor_sg->node;
    cmptb_map_t *nor_union_map = cmptb_map_union(nor_iso_map, iso_n, map_n);
    nor_exon_n = cmptb_map_exon_cnt(nor_union_map, map_n);
    int *whole_exon_id = (int*)_err_calloc(nor_exon_n, sizeof(int));
    int **nor_exon_id = (int**)_err_calloc(iso_n, sizeof(int*));
    for (i = 0; i < iso_n; ++i) nor_exon_id[i] = (int*)_err_calloc(nor_exon_n, sizeof(int));
    int *iso_exon_n = (int*)_err_calloc(iso_n, sizeof(int));
    int *exon_index = (int*)_err_calloc(sink-src+1, sizeof(int));

    if (check_module_type(nor_sg, nor_union_map, nor_iso_map, src, sink, nor_exon_n, iso_n, sgp->module_type, whole_exon_id, nor_exon_id, iso_exon_n, exon_index) == 0) goto ASMEnd;

    // exon header: ASM_ID, Exon_NUM, Isoform_NUM, strand, CHR_Name
    fprintf(sgp->out_fp[rep_n], "ASM#%d\t%d\t%d\t%c\t%s\t%s\t%s\n", *asm_i, nor_exon_n, iso_n, "+-"[nor_sg->is_rev], cname[nor_sg->tid], nor_sg->gene_name, nor_sg->gene_id);
    // exon coordinate pair: start, end
    for (exon_i = src; exon_i <= sink; ++exon_i) {
        if (is_cmptb_exon_iso(exon_i, nor_union_map))
            fprintf(sgp->out_fp[rep_n], "%d,%d\t", nor_node[exon_i].start, nor_node[exon_i].end);
    } fprintf(sgp->out_fp[rep_n], "\n");
    for (iso_i = 0; iso_i < iso_n; ++iso_i) {
        fprintf(sgp->out_fp[rep_n], "%d", iso_exon_n[iso_i]);
        for (exon_i = 0; exon_i < iso_exon_n[iso_i]; ++exon_i) {
            fprintf(sgp->out_fp[rep_n], "\t%d", exon_index[nor_exon_id[iso_i][exon_i]-src]);
        } fprintf(sgp->out_fp[rep_n], "\n");
    }

    // .IsoMatrix
    // calculate up/down pseudo length
    int up_pseu_len = 0, down_pseu_len = 0;
    // longest path from src to 0
    int src_len = (src == 0) ? 0 : (nor_sg->node[src].end - nor_sg->node[src].start + 1);
    up_pseu_len = dag_longest_path_from_src(nor_sg, nor_con_matrix, 0, src) - src_len;
    // longest path track from sink to sg->node_n-1
    int sink_len = (sink == nor_sg->node_n-1) ? 0 : (nor_sg->node[sink].end - nor_sg->node[sink].start + 1);
    down_pseu_len = dag_longest_path_from_sink(nor_sg, nor_con_matrix, sink, nor_sg->node_n-1) - sink_len;

    // for read-iso compatible matrix, use ONLY split_read_map and split_iso_map
    for (rep_i = 0; rep_i < rep_n; ++rep_i) {
        int bn = bundle_n[rep_i];
        int ri_map_n = 0, ri_map_m = 1000, map_n = 1 + ((iso_n-1) >> MAP_STEP_N); // map_n: number of cmptb_map_t per read-iso matrix
        read_iso_map *split_ri_map = read_iso_map_init(ri_map_m, map_n);
        for (b_i = 0; b_i < bn; ++b_i) {
            split_rm = (split_read_map[rep_i])+b_i;
            // gen read-iso map
            int zero;
            // from read-split_exon map to read-nor_exon map, then to read-nor_iso map
            read_iso_map *split_rim = gen_read_iso_map(split_rm, split_iso_map, split_iso_se, iso_n, map_n, &zero, sgp->fully);
            if (zero != 0) {
                // insert read-iso map, update weight
                insert_read_iso_map(&split_ri_map, &ri_map_n, &ri_map_m, split_rim);
            }
            read_iso_map_free(split_rim);
        }
        // matrix header: ASM_ID, ReadBundle_NUM, Isoform_NUM
        fprintf(sgp->out_fp[rep_i], "ASM#%d\t%d\t%d\t%d\t%d\n", *asm_i, ri_map_n, iso_n, up_pseu_len, down_pseu_len);
        // output read-iso compatible matrix
        for (r_i = 0; r_i < ri_map_n; ++r_i) {
            // read length and count
            fprintf(sgp->out_fp[rep_i], "%d\t%d", split_ri_map[r_i].rlen, split_ri_map[r_i].weight);
            // read-iso compatible matrix
            for (iso_i = 0; iso_i < iso_n; ++iso_i) {
                fprintf(sgp->out_fp[rep_i], "\t%ld", 0x1 & (split_ri_map[r_i].map[iso_i >> MAP_STEP_N] >> (~iso_i & MAP_STEP_M)));
            } fprintf(sgp->out_fp[rep_i], "\n");
        }
        for (r_i = 0; r_i < ri_map_m; ++r_i) 
            free(split_ri_map[r_i].map);
        free(split_ri_map);
    }

    (*asm_i)++;

ASMEnd:
    free(nor_union_map); free(exon_index);
    free(iso_exon_n); 
    for (i = 0; i < iso_n; ++i) { 
        free(nor_exon_id[i]); free(nor_iso_map[i]); free(split_iso_map[i]);
        free(nor_iso_se[i]); free(split_iso_se[i]);
    }
    free(nor_exon_id); free(whole_exon_id);
    free(nor_iso_map); free(split_iso_map);
    free(nor_iso_se); free(split_iso_se);
}

void add_pseu_wei(SG *sg, double **W, uint8_t **con_matrix) {
    int i, j, next_id, pre_id; double w;
    for (i = 0; i < sg->node[0].next_n; ++i) {
        next_id = sg->node[0].next_id[i];
        //W[0][next_id] = 10; // XXX pseudo weigh
        w = 0;
        if (is_con_matrix(con_matrix, 0, next_id) && W[0][next_id] == 0) {
            for (j = 0; j < sg->node[next_id].pre_n; ++j) {
                pre_id = sg->node[next_id].pre_id[j];
                if (pre_id != 0 && is_con_matrix(con_matrix, pre_id, next_id)) {
                    w = W[pre_id][next_id];
                    break;
                }

            }
            W[0][next_id] = w * 0.8; // XXX pseudo weigh
        }
    }
    for (i = 0; i < sg->node[sg->node_n-1].pre_n; ++i) {
        pre_id = sg->node[sg->node_n-1].pre_id[i];
        //rep_W[pre_id][sg->node_n-1] = 10; // XXX pseudo weigh
        w = 0;
        if (is_con_matrix(con_matrix, pre_id, sg->node_n-1) && W[pre_id][sg->node_n-1] == 0) {
            for (j = 0; j < sg->node[pre_id].next_n; ++j) {
                next_id = sg->node[pre_id].next_id[j];
                if (next_id != sg->node_n-1 && is_con_matrix(con_matrix, pre_id, next_id)) {
                    w = W[pre_id][next_id];
                }
            }
            W[pre_id][sg->node_n-1] = w * 0.8;
        }
    }
}

void gen_cand_asm(SG *nor_sg, SG *split_sg, gene_t *gene, char **cname, read_exon_map **split_M, double **nor_rep_W, uint8_t **nor_con_matrix, int rep_n, int *bundle_n, sg_para *sgp, int *asm_i) {
    int entry_n, exit_n; int *entry, *exit;

    cal_pre_domn(nor_sg, nor_rep_W, nor_con_matrix, sgp->junc_cnt_min), cal_post_domn(nor_sg, nor_rep_W, nor_con_matrix, sgp->junc_cnt_min);
    cal_cand_node(nor_sg, &entry, &exit, &entry_n, &exit_n, nor_con_matrix);
    if (entry_n == 0 || exit_n == 0) {
        free(entry); free(exit); return;
    }

    // for internal-terminal exon
    add_pseu_wei(nor_sg, nor_rep_W, nor_con_matrix);

    gec_t **nor_iso_path = (gec_t**)_err_malloc(sgp->iso_cnt_max * sizeof(gec_t*));
    gec_t *nor_path_idx = (gec_t*)_err_malloc(sgp->iso_cnt_max * sizeof(gec_t));
    // int **iso_se = (int**)_err_malloc(sgp->iso_cnt_max * sizeof(cmptb_map_t*));
    int last_exit = -1;
    int i, j;
    for (i = 0; i < entry_n; ++i) {
        if (entry[i] < last_exit) continue;
        for (j = 0; j < exit_n; ++j) {
            if (exit[j] < last_exit) continue;
            int post_domn_n = nor_sg->node[entry[i]].post_domn_n;
            int pre_domn_n = nor_sg->node[exit[j]].pre_domn_n;
            if (post_domn_n > 1 && pre_domn_n > 1 
                    && nor_sg->node[entry[i]].post_domn[1] == exit[j] 
                    && nor_sg->node[exit[j]].pre_domn[1] == entry[i]) {

                int iso_n, k; //map_n = ((nor_sg->node_n-1) >> MAP_STEP_N) + 1;

                int asm_node_n = sg_travl_n(nor_sg, entry[i], exit[j], nor_con_matrix);

                if (sgp->exon_num != 0 && sgp->exon_num != asm_node_n) continue;

                iso_n = 0;
                if (sgp->only_gtf) {
                    iso_n = anno_gen_cand_asm(nor_sg, gene, nor_con_matrix, entry[i], exit[j], nor_iso_path, nor_path_idx, sgp);
                } else if (asm_node_n <= sgp->asm_exon_max) {
                    if (asm_node_n <= sgp->exon_thres) { 
                        iso_n = enum_gen_cand_asm(nor_sg, nor_con_matrix, entry[i], exit[j], nor_iso_path, nor_path_idx, sgp);
                    } else {
                        // XXX only when use split_rep_W, bias_flow makes sense
                        // use sum(W) to generate isoform
                        iso_n = bias_flow_gen_cand_asm(nor_sg, nor_rep_W, nor_con_matrix, entry[i], exit[j], nor_iso_path, nor_path_idx, sgp);
                    }
                } else break; 
                // cal read-iso cmptb matrix
                if (iso_n > 1) {
                    read_iso_cmptb(nor_sg, split_sg, nor_con_matrix, cname, split_M, rep_n, bundle_n, nor_iso_path, nor_path_idx, iso_n, sgp, asm_i, entry[i], exit[j]);
                    if (sgp->recur == 0) last_exit = exit[j];
                }
                for (k = 0; k < sgp->iso_cnt_max; ++k) {
                    free(nor_iso_path[k]);
                }
                break;
            }
        }
    }
    free(nor_iso_path); free(nor_path_idx); free(entry); free(exit);
}

void update_rep_W(double **rep_W, double ***W, int rep_n, SG *sg, int junc_cnt_min) {
    int i, j, rep_i;
    for (i = 0; i < sg->node_n-1; ++i) {
        for (j = i+1; j < sg->node_n; ++j) {
            if (sg->node[i].end + 1 != sg->node[j].start) {
                int flag = 0;
                for (rep_i = 0; rep_i < rep_n; ++rep_i) {
                    if (W[rep_i][i][j] >= junc_cnt_min) {
                        flag = 1;
                        break;
                    }
                }
                if (flag == 0) continue;
            }
            for (rep_i = 0; rep_i < rep_n; ++rep_i) {
                rep_W[i][j] += W[rep_i][i][j];
            }
        }
    }
}

int bam_infer_exon(bam_aux_t *bam_aux, hts_itr_t *itr, exon_t **e, int *e_n, int *e_m, int **don, int *don_n, int *don_m, sg_para *sgp) {
    samFile *in = bam_aux->in;  bam1_t *b = bam_aux->b; int r; ad_t *ad = ad_init(1), *last_ad = ad_init(1);
    while ((r = sam_itr_next(in, itr, b)) >= 0)  {
        if (parse_bam_record1(b, ad, sgp) <= 0 || ad_sim_comp(ad, last_ad) == 0) continue;
        // 1. ad => update exonic coordinates
        *e_n = push_exon_coor(e, e_n, e_m, ad);
        // 2. ad => update splice-junction
        push_sj(don, don_n, don_m, ad);

        ad_copy(last_ad, ad); 
    }
    if (r < -1) err_func_format_printf("BAM file error. \"%s\"", bam_aux->fn);

    free_ad_group(ad, 1); free_ad_group(last_ad, 1);
    return *e_n;
}

int cand_asm_core(gene_group_t *gg, int g_n, sg_para *sgp, bam_aux_t **bam_aux)
{
    err_func_format_printf(__func__, "generate candidate isoform and read-isoform compatible matrix for each gene ...\n");
    int g_i, asm_i, rep_i;
    // nor_M: read - split_exon compatible map
    read_exon_map **split_M = (read_exon_map**)_err_malloc(sgp->tot_rep_n * sizeof(read_exon_map*));
    // nor_W: (nor_node1, nor_node2) -> weight, rep_nor_W: summation of nor_W
    double ***nor_W = (double***)_err_malloc(sgp->tot_rep_n * sizeof(double**));
    int *bundle_n = (int*)_err_malloc(sgp->tot_rep_n * sizeof(int));
    int i, hit = 0; char **cname = bam_aux[0]->h->target_name;
    FILE *gtf_fp = sgp->gtf_fp;

    asm_i = 0;
    for (g_i = 0; g_i < g_n; ++g_i) {
        gene_t *gene = gg->g+g_i;
        exon_t *infer_e = NULL; int infer_e_n = 0;

        char reg[1024]; sprintf(reg, "%s:%d-%d", cname[gene->tid], gene->start, gene->end);
        if (sgp->no_novel_exon == 0) {
            exon_t *bam_e = (exon_t*)_err_malloc(sizeof(exon_t)); int bam_e_n = 0, bam_e_m=1;
            int *don=(int*)_err_malloc(sizeof(int)), don_n=0, don_m=1;
            for (rep_i = 0; rep_i < sgp->tot_rep_n; ++rep_i) {
                hts_itr_t *itr = sam_itr_querys(bam_aux[rep_i]->idx, bam_aux[rep_i]->h, reg);
                bam_infer_exon(bam_aux[rep_i], itr, &bam_e, &bam_e_n, &bam_e_m, &don, &don_n, &don_m, sgp);
                hts_itr_destroy(itr); 
            }
            // 3. use splice-junction to split exons
            infer_e = infer_exon_coor(&infer_e_n, bam_e, bam_e_n, don, don_n);
            free(bam_e); free(don);
        }

        // XXX DO NOT modify *gene
        SG *nor_sg = build_SpliceGraph_novel_exon_core(gtf_fp, gene, cname, infer_e, infer_e_n, 0);
        gene_t *split_gene = copy_gene(gene);
        SG *split_sg = build_SpliceGraph_novel_exon_core(gtf_fp, split_gene, cname, infer_e, infer_e_n, 1);  free(infer_e); gene_free(split_gene);

        // rep_W: summation of nor_W
        double **nor_rep_W = (double**)_err_calloc(nor_sg->node_n, sizeof(double*));
        // nor_con_matirx: connectability of nor_nodes
        uint8_t **nor_con_matrix = (uint8_t**)_err_calloc(nor_sg->node_n, sizeof(uint8_t*)); // connect matrix
        for (i = 0; i < nor_sg->node_n; ++i) {
            nor_rep_W[i] = (double*)_err_calloc(nor_sg->node_n, sizeof(double));
            nor_con_matrix[i] = (uint8_t*)_err_calloc(nor_sg->node_n, sizeof(uint8_t));
        }
        hit = 0;
        for (rep_i = 0; rep_i < sgp->tot_rep_n; ++rep_i) {
            // 1. generate read-exon compatible array
            // 2. update SG with bamBundle
            // OR (optional) only known transcript
            nor_W[rep_i] = (double**)_err_calloc(nor_sg->node_n, sizeof(double*));
            hts_itr_t *itr = sam_itr_querys(bam_aux[rep_i]->idx, bam_aux[rep_i]->h, reg);
            split_M[rep_i] = bam_sg_cmptb(bam_aux[rep_i], itr, nor_W[rep_i], bundle_n+rep_i, split_sg, nor_sg, sgp);

            hts_itr_destroy(itr); 
            hit += bundle_n[rep_i];
        }
        update_rep_W(nor_rep_W, nor_W, sgp->tot_rep_n, nor_sg, sgp->junc_cnt_min); // use sum(W) of total reps
        // 3. generate asm
        gen_cand_asm(nor_sg, split_sg, gene, cname, split_M, nor_rep_W, nor_con_matrix, sgp->tot_rep_n, bundle_n, sgp, &asm_i);
        
        // free variables
        for (rep_i = 0; rep_i < sgp->tot_rep_n; ++rep_i) {
            for (i = 0; i < nor_sg->node_n; ++i) {
                free(nor_W[rep_i][i]);
            }
            free(nor_W[rep_i]);
            for (i = 0; i < bundle_n[rep_i]; ++i) {
                free(split_M[rep_i][i].map);
            } free(split_M[rep_i]);
        }
        for (i = 0; i < nor_sg->node_n; ++i) {
            free(nor_rep_W[i]); 
            free(nor_con_matrix[i]);
        } 
        free(nor_rep_W); free(nor_con_matrix);

        sg_free(nor_sg); sg_free(split_sg);
    }

    free(bundle_n); free(nor_W); free(split_M);
    err_func_format_printf(__func__, "generate candidate isoform and read-isoform compatible matrix for each gene done!\n");

    return 0;
}

/*****************************/
int cand_asm(int argc, char *argv[])
{
    int c, i; char ref_fn[1024]="", out_dir[1024]="", out_gtf[1024]="", *p;
    sg_para *sgp = sg_init_para();
    while ((c = getopt_long(argc, argv, "t:LnNuma:i:g:w:fjlFT:e:C:c:v:d:ro:G:", asm_long_opt, NULL)) >= 0) {
        switch (c) {
            case 't': sgp->n_threads = atoi(optarg); break;
            case 'L': sgp->in_list = 1; break;

            case 'f': sgp->only_gtf = 1; break;
            case 'n': sgp->no_novel_sj=0; break;
            case 'N': sgp->no_novel_exon=0; break;
            case 'u': sgp->read_type = SING_T; break;
            case 'm': sgp->use_multi = 1; break;
            case 'a': sgp->anchor_len[0] = strtol(optarg, &p, 10);
                      if (*p != 0) sgp->anchor_len[1] = strtol(p+1, &p, 10); else return cand_asm_usage();
                      if (*p != 0) sgp->anchor_len[2] = strtol(p+1, &p, 10); else return cand_asm_usage();
                      if (*p != 0) sgp->anchor_len[3] = strtol(p+1, &p, 10); else return cand_asm_usage();
                      if (*p != 0) sgp->anchor_len[4] = strtol(p+1, &p, 10); else return cand_asm_usage();
                      break;
            case 'i': sgp->intron_len = atoi(optarg); break;
            case 'g': strcpy(ref_fn, optarg); break;

            case 'w': sgp->rm_edge = 1, sgp->edge_wt = atof(optarg); break;

            case 'j': sgp->only_junc = 1; break;
            case 'l': sgp->only_novel = 1, sgp->no_novel_sj=0; break; // XXX only novel


            case 'F': sgp->fully = 1; break;
            case 'T': sgp->exon_thres = atoi(optarg); break;
            case 'e': sgp->asm_exon_max = atoi(optarg); break;
            case 'C': sgp->iso_cnt_max = atoll(optarg); break;
            case 'c': sgp->junc_cnt_min = atoi(optarg); break;
            case 'v': sgp->novel_junc_cnt_min = atoi(optarg); break;

            case 'd': sgp->module_type = atoi(optarg); break;
            case 'E': sgp->exon_num = atoi(optarg); break;
            case 'r': sgp->recur = 1; break;

            case 'o': strcpy(out_dir, optarg); break;
            case 'G': strcpy(out_gtf, optarg); sgp->gtf_fp = xopen(out_gtf, "w"); break;
            default: err_printf("Error: unknown option: %s.\n", optarg); return cand_asm_usage();
        }
    }
    if (argc - optind != 2) return cand_asm_usage();
    if (strlen(out_dir) == 0) {
        err_printf("Please specify output directory with \"-o\" option.\n");
        return cand_asm_usage();
    }

    int seq_n = 0, seq_m; kseq_t *seq=0;
    if (strlen(ref_fn) != 0) {
        gzFile genome_fp = gzopen(ref_fn, "r");
        if (genome_fp == NULL) { err_fatal(__func__, "Can not open genome file. %s\n", ref_fn); }
        seq = kseq_load_genome(genome_fp, &seq_n, &seq_m);
        err_gzclose(genome_fp); 
    }
    // 0. parse input bam file names
    bam_aux_t **bam_aux = (sgp->in_list ? sg_par_input_list(sgp, argv[optind+1]) : sg_par_input(sgp, argv[optind+1]));
    sgp->fp_n = sgp->tot_rep_n + 1;
    if (sgp->tot_rep_n <= 0) return cand_asm_usage();

    // 1. set cname --- 1 thread
    chr_name_t *cname = chr_name_init();
    bam_set_cname(bam_aux[0]->h, cname);
    // 2. build splice-graph --- 1 thread (Optional in future, infer exon-intron boundaries by bam records)
    // build from GTF file
    // SG_group *sg_g = construct_SpliceGraph(argv[optind], cname);
    gene_group_t *gg = gene_group_init();
    int g_n = read_gene_group(argv[optind], cname, gg);
    chr_name_free(cname);

    // 3. core process
    sgp->out_fp = iso_output(sgp, out_dir); gen_bit_table16();
    cand_asm_core(gg, g_n, sgp, bam_aux);

    gene_group_free(gg);
    if (seq_n > 0) {
        for (i = 0; i < seq_n; ++i) { free(seq[i].name.s), free(seq[i].seq.s); }
        free(seq);
    }
    for (i = 0; i < sgp->tot_rep_n; ++i) bam_aux_destroy(bam_aux[i]); 
    if (sgp->gtf_fp != NULL) err_fclose(sgp->gtf_fp);
    free(bam_aux); sg_free_para(sgp);
    return 0;
}
