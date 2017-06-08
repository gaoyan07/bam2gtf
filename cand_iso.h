#ifndef _CAND_ISO_H
#define _CAND_ISO_H
#include "splice_graph.h"

#define ISO_EXON_ALL_CNT 10
#define ISO_EXON_MAX 50
#define ISO_CNT_MAX 20
#define ISO_READ_CNT_MIN 1

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
cmptb_map_t *gen_iso_exon_map(gec_t *node_id, gec_t l, int map_n, int sg_node_n);
void insert_iso_exon_map(cmptb_map_t **iso_map, int *iso_i, int map_n, cmptb_map_t *iso_m);
int asm_output(char *fn, char *prefix, SG_group *sg_g, SGasm_group *asm_g, sg_para *sgp);
int cand_iso(int argc, char *argv[]);

#endif
