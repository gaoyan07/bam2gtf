#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include "utils.h"
#include "htslib/htslib/sam.h"

#define bam_unmap(b) ((b)->core.flag & BAM_FUNMAP)
#define COV_RATE 0.9

extern const char PROG[20];
int filter_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s filter [option] <in.bam/sam> | samtools sort > out.sort.bam\n\n", PROG);
    err_printf("Options:\n");
    err_printf("         -b --only-best            only keep the BEST ONE alignment record for each query. [false]\n");
    err_printf("         -v --cov-rate  [FLOAT]    minimum coverage of output alignment record. [%.2f]\n", COV_RATE);
    err_printf("\n");
    return 1;
}

void add_pathid(bam1_t *b, int id)
{
    char path[20]; int l;
    sprintf(path, ".path%d", id);
    l = strlen(path);

    if (b->l_data + l > b->m_data) { //.path1
        b->m_data =  b->l_data;
        kroundup32(b->m_data);
        b->data = (uint8_t*)_err_realloc(b->data, b->m_data*sizeof(uint8_t));
    }
    memmove(b->data+b->core.l_qname+l-1, b->data+b->core.l_qname-1, b->l_data-b->core.l_qname+1);
    memcpy(b->data+b->core.l_qname-1, path, l);
    b->l_data += l;
    b->core.l_qname += l;
}

int gtf_filter(bam1_t *b, int *score, float cov_rate)
{
    if (bam_unmap(b)) return 1;
    uint32_t *c = bam_get_cigar(b), n_c = b->core.n_cigar;
    // cover len/rate
    int cigar_qlen = b->core.l_qseq;
    int op0 = bam_cigar_op(c[0]), op1 = bam_cigar_op(c[n_c-1]);
    if (op0 == BAM_CSOFT_CLIP || op0 == BAM_CHARD_CLIP) cigar_qlen -= bam_cigar_oplen(c[0]);
    if (n_c > 1 && (op1 == BAM_CSOFT_CLIP || op1 == BAM_CHARD_CLIP)) cigar_qlen -= bam_cigar_oplen(c[n_c-1]);
    if ((cigar_qlen+0.0) / b->core.l_qseq < cov_rate) return 1;
    // NM 
    uint8_t *p = bam_aux_get(b, "NM"); // Edit Distance
    int ed; 
    ed = bam_aux2i(p);
    *score = cigar_qlen - ed;
    return 0;
}

const struct option filter_long_opt [] = {
    { "only-best", 0, NULL, 'b' },
    { "cov-rate", 1, NULL, 'v' },

    { 0, 0, 0, 0}
};

int bam_filter(int argc, char *argv[])
{
    int c; int only_best=0; float cov_rate=COV_RATE;
    while ((c = getopt_long(argc, argv, "bv:", filter_long_opt, NULL)) >= 0) {
        switch (c) {
            case 'b': only_best = 1; break;
            case 'v': cov_rate = atof(optarg); break;
            default : return filter_usage();
        }
    }

    if (argc - optind != 1) return filter_usage();

    samFile *in, *out; bam_hdr_t *h; bam1_t *b;
    bam1_t *best_b; int best_score=0, score; // for only-best
    if ((in = sam_open(argv[optind], "rb")) == NULL) err_fatal(__func__, "Cannot open \"%s\"\n", argv[optind]);
    if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", argv[optind]);
    b = bam_init1();  best_b = bam_init1();

    if ((out = sam_open_format("-", "wb", NULL)) == NULL) err_fatal_simple("Cannot open \"-\"\n");
    if (sam_hdr_write(out, h) != 0) err_fatal_simple("Error in writing SAM header\n"); //sam header
    char lqname[1024]="\0"; int id=1;
    while (sam_read1(in, h, b) >= 0) {
        if (gtf_filter(b, &score, cov_rate)) continue;

        if (only_best) {
            if (strcmp(bam_get_qname(b), lqname) == 0) {
                if (score > best_score) {
                    bam_copy1(best_b, b);
                    best_score = score;
                }
            } else { 
                if (strcmp(lqname, "\0") != 0) {
                    add_pathid(best_b, 1);
                    if (sam_write1(out, h, best_b) < 0) err_fatal_simple("Error in writing SAM record\n");
                }
                bam_copy1(best_b, b);
                best_score = score;
                strcpy(lqname, bam_get_qname(b)); 
            }
        } else {
            if (strcmp(bam_get_qname(b), lqname) == 0) id++;
            else { id = 1; strcpy(lqname, bam_get_qname(b)); }
            add_pathid(b, id);
            if (sam_write1(out, h, b) < 0) err_fatal_simple("Error in writing SAM record\n");
        }
    }

    bam_destroy1(b); bam_destroy1(best_b); bam_hdr_destroy(h); sam_close(in); sam_close(out);
    return 0;
}
