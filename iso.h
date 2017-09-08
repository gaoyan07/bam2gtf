#ifndef _ISO_H
#define _ISO_H
#include "splice_graph.h"
#include "parse_bam.h"

#define ALL_TYPE 0
#define SE_TYPE 1
#define A5SS_TYPE 2
#define A3SS_TYPE 3
#define MXE_TYPE 4
#define RI_TYPE 5

#define _2SE_TYPE 6
#define AIE_TYPE 7
#define ATE_TYPE 8

#define ISO_EDGE_MIN_WEI 0.1
#define ISO_EXON_ALL_CNT 10
#define ISO_EXON_MAX 50
#define ISO_CNT_MAX 20
//#define ISO_READ_CNT_MIN 1
#define NOVEL_JUNC_READ_CNT_MIN 3
#define JUNC_READ_CNT_MIN 1

#define cmptb_map_t uint64_t
#define MAP_STEP 64
#define MAP_STEP_N 6
#define MAP_STEP_M 0x3fULL // 63

typedef struct {
    int weight; // number of read compatible to these exons
    gec_t *exon_id, exon_n, exon_m;
} cmptb_read_exon;

typedef struct {
    int weight;
    int *iso_id, iso_n, iso_m;
} cmptb_read_iso;

/* example:
 * read#1: exon#2, exon#4, start=2, end=4
 * iso#1:  exon#0, exon#2, exon#4, exon#6
 * 
 * read_exon_map: { 0,2 0,4, 00010100 }
 * exon_iso_map:             01010101
 * 
 * read_iso_map =               ...
 * compatibale!
 */
typedef struct {
    int weight; int rlen;
    gec_t map_s, map_si, map_e, map_ei; // map_e - map_s : number of map of each read
    cmptb_map_t *map;
} read_exon_map;

typedef struct {
    int weight, rlen, map_n;
    cmptb_map_t *map;
} read_iso_map;

//cmptb_map_t *iso_exon_map;


void sg_free_asm_group(SGasm_group *asm_g);
SGasm_group *gen_asm(SG_group *sg_g, sg_para *sgp);
int cal_asm_exon_cnt(SG_group *sg_g, samFile *in, bam_hdr_t *h, bam1_t *b);
cmptb_map_t *gen_iso_exon_map(gec_t *node_id, gec_t l, int map_n, int sg_node_n, int *se);
void insert_iso_exon_map(cmptb_map_t **iso_map, int **iso_se, int *iso_i, int map_n, cmptb_map_t *iso_m, int *se);
read_exon_map *bam_sg_cmptb(bam_aux_t *bam_aux, hts_itr_t *itr, double **wei_matrix, int *b_n, SG *sg, sg_para *sgp);
void add_pseu_wei(SG *sg, double **W, uint8_t **con_matrix);
int asm_output(char *fn, char *prefix, SG_group *sg_g, SGasm_group *asm_g, sg_para *sgp);
int bam_infer_exon(bam_aux_t *bam_aux, hts_itr_t *itr, exon_t **bam_e, int *bam_e_n, int *bam_e_m, int **don, int *don_n, int *don_m, sg_para *sgp);
int cand_asm(int argc, char *argv[]);
int full_iso(int argc, char *argv[]);


int dag_longest_path(SG *sg, uint8_t **con_matrix, int src, int sink);

int bias_flow_gen_cand_asm(SG *sg, double **rep_weight, uint8_t **con_matrix, int src, int sink, cmptb_map_t **iso_map, int **iso_se, int map_n, sg_para *sgp);
int enum_gen_cand_asm(SG *sg, uint8_t **con_matrix, int src, int sink, cmptb_map_t **iso_map, int **iso_se, int map_n, sg_para *sgp);
int bias_flow_gen_full_iso(SG *sg, char **cname, double **weight, uint8_t **con_matrix, int src, int sink, sg_para *sgp);

int anno_gen_cand_asm(SG *sg, gene_t *gene, uint8_t **con_matrix, int src, int sink, cmptb_map_t **iso_map, int **iso_se, int map_n, sg_para *sgp);

#endif
