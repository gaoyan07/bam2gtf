/* bam2update_gtf.c
 *   generate junction information based on sam/bam file
 *   then, update existing GTF file
 *   currently, only for single-end long read data
 * 
 * Author:  Yan Gao
 * Contact: yangao07@hit.edu.cn                             */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "htslib/htslib/sam.h"
#include "utils.h"
#include "gtf.h"

#define bam_unmap(b) ((b)->core.flag & BAM_FUNMAP)

extern const char PROG[20];
int bam2update_gtf_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s bam2update_gtf [option] <in.bam> <old.gtf> > new.gtf\n", PROG);
    err_printf("\n");
    err_printf("Notice:  the BAM and GTF files should be sorted in advance.\n");
    err_printf("Options:\n\n");
    err_printf("         -s --source       [STR]    source field in GTF, program, database or project name.\n");
    err_printf("                                    [NONE]\n");
	err_printf("\n");
	return 1;
}

int gen_exon(trans_t *t, bam1_t *b, uint32_t *c, uint32_t n_cigar)
{
    t->exon_n = 0;
    int32_t tid = b->core.tid; int32_t start = b->core.pos+1, end = start-1, qstart = 1, qend = 0;/*1-base*/
    uint8_t is_rev = bam_is_rev(b);
    uint32_t i;
    for (i = 0; i < n_cigar; ++i) {
        int l = bam_cigar_oplen(c[i]);
        switch (bam_cigar_op(c[i])) {
            case BAM_CREF_SKIP:  // N/D(0 1)
            case BAM_CDEL :
                if (l >= INTRON_MIN_LEN) {
                    add_exon(t, tid, start, end, qstart, qend, is_rev);
                    qstart = qend+1;
                    start = end + l + 1;
                }
                end += l;
                break;
            case BAM_CMATCH: // 1 1
            case BAM_CEQUAL:
            case BAM_CDIFF:
                qend += l, end += l;
                break;
            case BAM_CINS: // 1 0
            case BAM_CSOFT_CLIP:
            case BAM_CHARD_CLIP:
                qend += l;
                break;
            case BAM_CPAD: // 0 0
            case BAM_CBACK:
                break;
            default:
                err_printf("Error: unknown cigar type: %d.\n", bam_cigar_op(c[i]));
                break;
        }
    }
    add_exon(t, tid, start, end, qstart, qend, is_rev);
    return 0;
}

int gen_trans(bam1_t *b, trans_t *t)
{
    if (bam_unmap(b)) return 0;

    uint32_t *c = bam_get_cigar(b), n_cigar = b->core.n_cigar;
    gen_exon(t, b, c, n_cigar);
    return 1;
    /*
       if (0) {
       uint8_t *p = bam_aux_get(b, "XA"); // alternative alignment
       if (p) { 
       if (*p == 'Z') {
       char *s = bam_aux2Z(p);
       printf("XA:Z:%s\n", s);
       }
       }
       }
    */
}

char *fgets_gene(FILE *gfp, char *line, int size)
{
    char s[100];
    while (fgets(line, size, gfp) != NULL) {
        if ((sscanf(line, "%*s\t%*s\t%s", s) == 1) && strcmp(s, "gene") == 0) {
            return line;
        }
        else return NULL;
    }
    return NULL;
}

int bam2update_gtf(int argc, char *argv[])
{
    int c;
    char src[1024]="NONE";
	while ((c = getopt(argc, argv, "b:s:")) >= 0)
    {
        switch(c)
        {
           case 's':
                strcpy(src, optarg);
                break;
            default:
                err_printf("Error: unknown option: %s.\n", optarg);
                return bam2update_gtf_usage();
                break;
        }
    }
    if (argc - optind != 2) return bam2update_gtf_usage();

    samFile *in; bam_hdr_t *h; bam1_t *b;

    in = sam_open(argv[optind], "rb");
    if (in == NULL) err_fatal(__func__, "Cannot open \"%s\"\n", argv[optind]);
    h = sam_hdr_read(in);
    if (h == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", argv[optind]);
    b = bam_init1();

    FILE *gfp = fopen(argv[optind+1], "r"); char *gtf_line = (char*)malloc(sizeof(char)*1024);

    read_trans_t *r = read_trans_init();
    trans_t *t = trans_init(1);

    int sam_ret = sam_read1(in, h, b); gtf_line = fgets_gene(gfp, gtf_line, 1024);
    while (sam_ret >= 0 && gtf_line != NULL) {
    }
    while (sam_read1(in, h, b) >= 0) {
        r->trans_n = 0;
        if (gen_trans(b, t)) add_read_trans(r, *t, bam_get_qname(b));
        set_read_trans(r);
        print_read_trans(*r, h, src, stdout);
    }

    read_trans_free(r), trans_free(t); free(gtf_line);
    bam_destroy1(b); bam_hdr_destroy(h); sam_close(in); fclose(gfp);
    return 0;
}
