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

typedef struct {
    int32_t tid; uint8_t is_rev;
    int32_t start, end; //1-base, ref
} exon_t;

typedef struct {
    exon_t *exon; int exon_n, exon_m;
    int32_t tid; uint32_t is_rev;
    int32_t start, end;
    char qname[1024];
} trans_t;

typedef struct {
    trans_t *t; int trans_n, trans_m;
} read_trans_t;

typedef struct {
    trans_t *trans; int trans_n, trans_m;
    int32_t tid; uint32_t is_rev;
    int32_t start, end;
} gene_t;

exon_t *exon_init(int n);
void exon_free(exon_t *e);

trans_t *trans_init(int n);
int add_exon(trans_t *t, int32_t tid, int32_t start, int32_t end, uint8_t is_rev);
int set_trans(trans_t *t, char *qname);
trans_t *exon_realloc(trans_t *t);
void trans_free(trans_t *t);

read_trans_t *read_trans_init(void);
void add_read_trans(read_trans_t *r, trans_t t);
read_trans_t *read_trans_realloc(read_trans_t *r);
void read_trans_free(read_trans_t *r);
//int set_read_trans(read_trans_t *r);

gene_t *gene_init(void);
void add_trans(gene_t *g, trans_t t);
gene_t *trans_realloc(gene_t *g);
void gene_free(gene_t *g);

int print_exon(exon_t e, FILE *out);
int print_trans(trans_t t, bam_hdr_t *h, char *src, FILE *out);
int print_read_trans(read_trans_t r, bam_hdr_t *h, char *src, FILE *out);
int print_gene(gene_t g, FILE *out);

#define INTRON_MIN_LEN 50

#endif
