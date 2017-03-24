#ifndef _GTF_H
#define _GTF_H

#include <stdint.h>
#include <stdio.h>
#include "htslib/htslib/sam.h"

/* TAB-separated standard GTF columns */
/*
 * col-num      content                                     value/format
 * -------------------------------------------------------------------------------------------------------------------------
 *    1     chromosome name                             chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,M}
 *    2     annotation source                           {ENSEMBL,HAVANA,etc.}
 *    3     feature-type                                {gene,transcript,exon,CDS,UTR,start_codon,stop_codon,Selenocysteine}
 *    4     genomic start location                      integer-value (1-based)
 *    5     genomic end location                        integer-value
 *    6     score (not used)                            .
 *    7     genomic strand                              {+,-}
 *    8     genomic phase (for CDS features)            {0,1,2,.}
 *    9     additional information as key-value pairs   see below
 */
/*
 *  key name                value format
 * ------------------------------------------------------------------------------------
 * gene_id              ENSGXXXXXXXXXXX *
 * transcript_id        ENSTXXXXXXXXXXX *
 * gene_type            list of biotypes
 * gene_status          {KNOWN, NOVEL, PUTATIVE}
 * gene_name            string
 * transcript_type      list of biotypes
 * transcript_status    {KNOWN, NOVEL, PUTATIVE}
 * transcript_name      string
 * exon_number          indicates the biological position of the exon in the transcript
 * exon_id              ENSEXXXXXXXXXXX *
 * level                1 (verified loci),
 *                      2 (manually annotated loci),
 *                      3 (automatically annotated loci)
 */
#define MAX_SITE 2147483647
#define DON_SITE_F 0
#define ACC_SITE_F 1

typedef struct {
    int32_t tid; uint8_t is_rev;
    int32_t start, end; //1-based, ref
                        //0: start of init exon
                        //MAX: end of term exon
} exon_t;

typedef struct {
    int32_t tid; uint8_t strand; // 0:undefined, 1:+, 2:-
    int32_t don, acc;
    uint8_t motif, is_anno;
    int32_t uniq_c, multi_c, max_over;
} sj_t;

typedef struct {
    exon_t up, se, down;
} SE_t;   // skipped exon

typedef struct {
    exon_t up, lon, shor;
} A3SS_t; // alternative 3' splice site

typedef struct {
    exon_t lon, shor, down;
} A5SS_t; // alternative 3' splice site

typedef struct {
    exon_t up, fir, sec, down;
} MXE_t; // mutually exclusive exon

typedef struct {
    exon_t up, down;
} RI_t;  // retained intron

typedef struct {
    SE_t *se; int32_t se_n, se_m;
    A5SS_t *a5ss; int32_t a5ss_n, a5ss_m;
    A3SS_t *a3ss; int32_t a3ss_n, a3ss_m;
    MXE_t *mxe; int32_t mxe_n, mxe_m;
    RI_t *ri; int32_t ri_n, ri_m;
} ASE_t;

#define set_l_iden(map) (map |= 0x4)
#define set_r_iden(map) (map |= 0x2)
#define set_b_iden(map) (map |= 0x1)

#define check_l_iden(map) (map & 0x4)
#define check_r_iden(map) (map & 0x2)
#define check_b_iden(map) (map & 0x1)

typedef struct {
    exon_t *exon; int exon_n, exon_m;
    uint8_t *novel_exon_map, *novel_sj_map; // 3-bit map: l-iden | r-iden | both-iden
    int32_t tid; uint32_t is_rev;
    int32_t start, end;
    char tname[100];
    char gname[100];
    int novel_gene_flag, cov;
    uint8_t full, lfull, lnoth, rfull, rnoth; 
    uint8_t novel, all_novel, all_iden;
} trans_t;

typedef struct {
    int32_t tid; uint8_t is_rev;
    int32_t start, end; //1-based, ref
    uint8_t is_canon, is_anno;
    int32_t uniq_c, multi_c;
} intron_t;

typedef struct {
    intron_t *intron; int intron_n, intron_m;
} intron_group_t;

typedef struct {
    trans_t *t; int trans_n, trans_m;
} read_trans_t;

typedef struct {
    trans_t *trans; int trans_n, anno_tran_n, trans_m;
    int32_t tid; uint8_t is_rev;
    int32_t start, end;
    char gname[100];
} gene_t;

typedef struct {
    gene_t *g; int gene_n, gene_m;
    int32_t tid, start, end;
} gene_group_t;

typedef struct {
    char **chr_name;
    int32_t chr_n, chr_m;
} chr_name_t;

exon_t *exon_init(int n);
void exon_free(exon_t *e);

chr_name_t *chr_name_init(void);
void chr_name_free(chr_name_t *cname);
int read_sj_group(FILE *sj_fp, chr_name_t *cname, sj_t **sj_group, int sj_m);
int bam_set_cname(bam_hdr_t *h, chr_name_t *cname);

trans_t *trans_init(int n);
int add_exon(trans_t *t, int32_t tid, int32_t start, int32_t end, uint8_t is_rev);
int set_trans(trans_t *t, char *qname);
trans_t *exon_realloc(trans_t *t);
void trans_free(trans_t *t);

read_trans_t *read_trans_init(void);
void add_read_trans(read_trans_t *r, trans_t t);
read_trans_t *read_trans_realloc(read_trans_t *r);
void novel_read_trans_free(read_trans_t *r);
void read_trans_free(read_trans_t *r);
//int set_read_trans(read_trans_t *r);

intron_t *intron_init(int n);
intron_group_t *intron_group_init(void);
void add_intron(intron_group_t *i, intron_t i1);
void intron_group_free(intron_group_t *i);

gene_t *gene_init(void);
void add_trans(gene_t *g, trans_t t, int novel_gene_flag);
gene_t *trans_realloc(gene_t *g);
void gene_free(gene_t *g);
void gtf_add_info(char add_info[], char tag[], char *info);

gene_group_t *gene_group_init(void);
gene_group_t *gene_group_realloc(gene_group_t *gg);
void add_gene(gene_group_t *gg, gene_t g, int novel_gene_flag);
void set_gene_group(gene_group_t *gg);
void gene_group_free(gene_group_t *gg);
int read_gene_group(FILE *gtf, chr_name_t *cname, gene_group_t *gg);

int print_exon(exon_t e, FILE *out);
int print_trans(trans_t t, bam_hdr_t *h, char *src, FILE *out);
int print_read_trans(read_trans_t r, bam_hdr_t *h, char *src, FILE *out);
int print_gene(gene_t g, FILE *out);
void print_gene_group(gene_group_t gg, bam_hdr_t *h, char *src, FILE *out, char **group_line, int *group_line_n);
void print_gtf_trans(gene_t g, bam_hdr_t *h, char *src, FILE *out);

#define INTRON_MIN_LEN 25
#define INTER_EXON_MIN_LEN 6
#define SPLICE_DISTANCE 0
#define MIN_INTRON_NUM 0

int check_iden(trans_t t1, trans_t t2, int dis);

#endif
