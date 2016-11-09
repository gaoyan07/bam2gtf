/* bam2gtf.c
 *   generate junction information based on sam/bam file
 *   currently, only for single-end long read data
 * 
 * Author:  Yan Gao
 * Contact: yangao07@hit.edu.cn                             */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include "htslib/htslib/sam.h"
#include "utils.h"
#include "gtf.h"


#define bam_unmap(b) ((b)->core.flag & BAM_FUNMAP)

int bam2gtf_usage(void)
{
    err_printf("\n");
	err_printf("Program: bam2gtf\n");
    err_printf("Usage:   bam2gtf <in.bam> > out.gtf\n");
	//err_printf("Usage:   bam2gtf <command> [options]\n\n");
	//err_printf("Commands: \n");
	//err_printf("         unipath     generate unipath seq from bwt-str\n");
    //err_printf("         index       index unipath's bwt-str\n");
	//err_printf("         query       query the unipath with the bwt index\n");
	err_printf("\n");
	return 1;
}

int gen_exon(trans_t *t, bam1_t *b, uint32_t *c, uint32_t n_cigar)
{
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

    if (0) {
        uint8_t *p = bam_aux_get(b, "XA"); // alternative alignment
        if (p) { 
            if (*p == 'Z') {
                char *s = bam_aux2Z(p);
                //printf("XA:Z:%s\n", s);
            }
        }
    }
    // construct trans based all the exons
    // XXX two split-alignments => two trans
    int res = 0;
    return res;
}

int main(int argc, char *argv[])
{
    samFile *in; bam_hdr_t *h; bam1_t *b; int r;

    if (argc < 2) return bam2gtf_usage();

    in = sam_open(argv[optind], "rb");
    if (in == NULL) err_printf("Error opening \"%s\"\n", argv[optind]);
    h = sam_hdr_read(in);
    if (h == NULL) err_printf("Couldn't read header for \"%s\"\n", argv[optind]);
    b = bam_init1();

    char qname[1024], lqname[1024]="\0"; 
    trans_t *t = trans_init(1);

    while ((r = sam_read1(in, h, b)) >= 0) {
        stdout_printf("%s %d %d\n", bam_get_qname(b), b->core.tid, b->core.pos);
        // gen trans{exon,exon}
        strcpy(qname, bam_get_qname(b));
        if (strcmp(qname, lqname) != 0) {
            set_trans(&t);
            print_trans(*t, stdout);
            strcpy(lqname, qname);
            t->exon_n = 0;
        }
        r = gen_trans(b, t);
    }
    set_trans(&t);
    print_trans(*t, stdout);

    bam_destroy1(b); bam_hdr_destroy(h); sam_close(in);
    return 0;
}
