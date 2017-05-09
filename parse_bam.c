#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>
#include "htslib/htslib/sam.h"
#include "bam2gtf.h"
#include "parse_bam.h"
#include "utils.h"
#include "gtf.h"
#include "kseq.h"
#include "build_sg.h"

extern const char PROG[20];
const int intron_motif_n = 6;
const char intron_motif[6][10] = {
    "GTAG", "CTAC", 
    "GCAG", "CTGC", 
    "ATAC", "GTAT"
};
const int intron_motif_strand[6] = {
    1, 2, 1, 2, 1, 2
};

int bam2sj_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s bam2sj [option] <in.bam> > out.sj\n\n", PROG);
    err_printf("Note:    in.bam should be sorted in advance\n\n");
    err_printf("Input Options:\n\n");
    err_printf("         -G --gtf-anno    [STR]    GTF annotation file, indicating known splice-junctions. \n");
    err_printf("         -g --genome-file [STR]    genome.fa. Use genome sequence to classify intron-motif. \n");
    err_printf("                                   If no genome file is give, intron-motif will be set as 0\n");
    err_printf("                                   (non-canonical) [None]\n");
    err_printf("\nFilter Options:\n\n");
    err_printf("         -p --prop-pair            set -p to force to filter out reads mapped in improper pair. [False]\n");
    err_printf("         -a --anchor-len  [INT,INT,INT,INT,INT]\n");
    err_printf("                                   minimum anchor length for junction read, [annotated, non-canonical,\n");
    err_printf("                                    GT/AG, GC/AG, AT/AC]. [%d,%d,%d,%d,%d]\n", ANCHOR_MIN_LEN, NON_ANCHOR, ANCHOR1, ANCHOR2, ANCHOR3);
    err_printf("         -U --uniq-map    [INT,INT,INT,INT,INT]\n");
    err_printf("                                   minimum uniq-map read count for junction read, [annotated,\n");
    err_printf("                                   non-canonical, GT/AG, GC/AG, AT/AC]. [%d,%d,%d,%d,%d]\n", UNIQ_MIN, NON_UNIQ_MIN, UNIQ_MIN1, UNIQ_MIN2, UNIQ_MIN3);
    err_printf("         -A --all-map     [INT,INT,INT,INT,INT]\n");
    err_printf("                                   minimum total uniq-map and multi-map read count for junction\n");
    err_printf("                                   read, [annotated, non-canonical, GT/AG, GC/AG, AT/AC].\n");
    err_printf("                                   [%d,%d,%d,%d,%d]\n", ALL_MIN, NON_ALL_MIN, ALL_MIN1, ALL_MIN2, ALL_MIN3);
    err_printf("         -i --intron-len  [INT]    minimum intron length for junction read. [%d]\n", INTRON_MIN_LEN);
	err_printf("\n");
	return 1;
}

const struct option bam2sj_long_opt [] = {
    { "proper-pair", 1, NULL, 'p' },
    { "gtf-anno", 1, NULL, 'G' },
    { "genome-file", 1, NULL, 'g' },
    { "anchor-len", 1, NULL, 'a' },
    { "uniq-map", 1, NULL, 'U' },
    { "all-map", 1, NULL, 'A' },
    { "intron-len", 1, NULL, 'i' },

    { 0, 0, 0, 0}
};

// bam operation
uint8_t bam_is_uniq_NH(bam1_t *b)
{
    uint8_t *p = bam_aux_get(b, "NH");
    if (p == 0) {
        err_printf("No \"NH\" tag.\n");
        return 0;
    }
    return (bam_aux2i(p) == 1);
}

int bam_cigar_opn(int n_cigar, const uint32_t *cigar, uint32_t op, uint32_t len)
{
    int i, j;
    for (i = j = 0; i < n_cigar; ++i)
        if (bam_cigar_op(cigar[i]) == op && bam_cigar_oplen(cigar[i]) >= len) ++j;
    return j;
}

int add_sj(sj_t **sj, int sj_i, int *sj_m, int tid, int don, int acc, uint8_t strand, uint8_t motif_i, uint8_t is_anno, uint8_t is_uniq)
{
    if (sj_i == *sj_m) {
        _realloc(*sj, *sj_m, sj_t)
    }
    (*sj)[sj_i].tid = tid;
    (*sj)[sj_i].don = don;
    (*sj)[sj_i].acc = acc;
    (*sj)[sj_i].strand = strand;
    (*sj)[sj_i].motif = motif_i;
    (*sj)[sj_i].is_anno = is_anno;
    (*sj)[sj_i].uniq_c = is_uniq; 
    (*sj)[sj_i].multi_c = 1-is_uniq;
    return 0;
}

void free_ad_group(ad_t *ad, int ad_n)
{
    int i; for (i = 0; i < ad_n; ++i) {
        free(ad[i].exon_end); 
        if (ad[i].intr_end != NULL) free(ad[i].intr_end);
    }
    free(ad);
}

uint8_t intr_deri_str(kseq_t *seq, int seq_n, int tid, int start, int end, uint8_t *motif_i)
{
    *motif_i = 0;
    if (seq_n == 0) return 0;
    if (tid >= seq_n) err_fatal(__func__, "unknown tid: %d\n", tid); 
    char intron[10]="";
    intron[0] = toupper(seq[tid].seq.s[start-1]);
    intron[1] = toupper(seq[tid].seq.s[start]);
    intron[2] = toupper(seq[tid].seq.s[end-2]);
    intron[3] = toupper(seq[tid].seq.s[end-1]);
    int i;
    for (i = 0; i < intron_motif_n; ++i) {
        if (strcmp(intron, intron_motif[i]) == 0) {
            *motif_i = i+1;
            return intron_motif_strand[i];
        }
    }
    return 0;
}

int sj_sch_group(sj_t *SJ, int SJ_n, sj_t sj, int *hit)
{
    *hit = 0;
    if (SJ_n == 0) return 0;

    int i; int tid, don, acc;
    for (i = SJ_n-1; i >= 0; i--) {
        tid = SJ[i].tid, don = SJ[i].don, acc = SJ[i].acc;
        if (tid == sj.tid &&  don == sj.don && acc == sj.acc) { *hit = 1; return i; }
        else if (tid < sj.tid || don < sj.don || (don == sj.don && acc < sj.acc)) return i+1; // SJ[i] < sj
    }
    return 0;
}

int sj_update_group(sj_t **SJ_group, int *SJ_n, int *SJ_m, sj_t *sj, int sj_n)
{
    if (sj_n + *SJ_n > *SJ_m) _realloc(*SJ_group, *SJ_m, sj_t)
    int i, hit=0;
    for (i = 0; i < sj_n; ++i) {
        int sj_i = sj_sch_group(*SJ_group, *SJ_n, sj[i], &hit);
        if (hit == 0) {
            if ((*SJ_n)++ >= *SJ_m) _realloc(*SJ_group, *SJ_m, sj_t)
            // memmove
            if (sj_i <= *SJ_n-2)
                memmove(*SJ_group+sj_i+1, *SJ_group+sj_i, (*SJ_n-sj_i-1) * sizeof(sj_t));
            // set sj
            (*SJ_group)[sj_i].tid = sj[i].tid;
            (*SJ_group)[sj_i].don = sj[i].don;
            (*SJ_group)[sj_i].acc = sj[i].acc;
            (*SJ_group)[sj_i].strand = sj[i].strand;
            (*SJ_group)[sj_i].motif = sj[i].motif;
            (*SJ_group)[sj_i].is_anno = sj[i].is_anno;
            (*SJ_group)[sj_i].uniq_c = sj[i].uniq_c;
            (*SJ_group)[sj_i].multi_c = sj[i].multi_c;
        } else {
            (*SJ_group)[sj_i].uniq_c += sj[i].uniq_c;
            (*SJ_group)[sj_i].multi_c += sj[i].multi_c;
            if ((*SJ_group)[sj_i].strand != sj[i].strand) (*SJ_group)[sj_i].strand = 0; // undefined
        }
    }
    return 0;
}

kseq_t *kseq_load_genome(gzFile genome_fp, int *_seq_n, int *_seq_m)
{
    int seq_n = 0, seq_m = 30;
    kseq_t *kseq = kseq_init(genome_fp), *seq = (kseq_t*)_err_malloc(30 * sizeof(kseq_t));

    err_func_format_printf(__func__, "loading genome fasta file ...\n");
    while (kseq_read(kseq) >= 0) {
        kseq_copy(seq+seq_n, *kseq);
        seq_n++;
        if (seq_n == seq_m) {
            seq_m <<= 1;
            seq = (kseq_t*)_err_realloc(seq, seq_m * sizeof(kseq_t));
        }
    }
    err_func_format_printf(__func__, "loading genome fasta file done!\n");
    kseq_destroy(kseq);
    *_seq_n = seq_n; *_seq_m = seq_m;
    return seq;
}

int gen_sj(uint8_t is_uniq, int tid, int start, int n_cigar, uint32_t *c, kseq_t *seq, int seq_n, sj_t **sj, int *sj_m, sg_para *sgp)
{
    int end = start - 1; /* 1-base */
    uint8_t strand, motif_i;
    
    int i; int min_intr_len = sgp->intron_len, sj_i = 0;

    for (i = 0; i < n_cigar; ++i) {
        int l = bam_cigar_oplen(c[i]);
        switch (bam_cigar_op(c[i])) {
            case BAM_CREF_SKIP: // N(0 1)
                if (l >= min_intr_len) {
                    strand = intr_deri_str(seq, seq_n, tid, end+1, end+l, &motif_i);
                    // XXX filter with anchor length
                    add_sj(sj, sj_i, sj_m, tid, end+1, end+l, strand, motif_i, 1, is_uniq); ++sj_i;
                    start = end+l+1;
                }
                end += l;
                break;
            case BAM_CMATCH: // 1 1
            case BAM_CEQUAL:
            case BAM_CDIFF:
                end += l;
                break;
            case BAM_CDEL : // D(0 1)
                end += l;
                break; // XXX only M&N allowed
            case BAM_CINS: // 1 0
            case BAM_CSOFT_CLIP:
            case BAM_CHARD_CLIP:
            case BAM_CPAD: // 0 0
            case BAM_CBACK:
                break; // XXX only M&N allowed
            default:
                err_printf("Error: unknown cigar type: %d.\n", bam_cigar_op(c[i]));
                break;
        }
    }

    return sj_i;
}

// only junction-read are kept in AD_T
int parse_bam(int tid, int start, int *_end, int n_cigar, const uint32_t *c, uint8_t is_uniq, kseq_t *seq, int seq_n, ad_t **ad_g, int *ad_n, int *ad_m, sj_t **sj, int *sj_n, int *sj_m, sg_para *sgp)
{
    int i, min_intr_len = sgp->intron_len, SJ_n, sj_i = 0;
#ifdef _RMATS_
    for (i = 0; i  < n_cigar; ++i) {
        uint32_t op = bam_cigar_op(c[i]);
        if (op != BAM_CMATCH && op != BAM_CREF_SKIP) return -1;
        //&& op != BAM_CEQUAL && op != BAM_CDIFF) return 0;
    }
#endif
    SJ_n = bam_cigar_opn(n_cigar, c, BAM_CREF_SKIP, min_intr_len);
    int end = start - 1; /* 1-base */
    uint8_t strand, motif_i;

    ad_t *ad;
    if (SJ_n > 0) {
        if (*ad_n == *ad_m) _realloc(*ad_g, *ad_m, ad_t)
        ad = (*ad_g)+(*ad_n);
        ad->intv_n = 0;
        ad->exon_end = (int*)_err_malloc((SJ_n+1) * sizeof(int));
        ad->intr_end = (int*)_err_malloc(SJ_n * sizeof(int));
        ad->tid = tid; ad->start = start;
        ad->is_uniq = is_uniq; ad->is_splice = 1;
    }
    for (i = 0; i < n_cigar; ++i) {
        int l = bam_cigar_oplen(c[i]);
        switch (bam_cigar_op(c[i])) {
            case BAM_CREF_SKIP: // N(0 1)
                if (l >= min_intr_len) {
                    // sj
                    strand = intr_deri_str(seq, seq_n, tid, end+1, end+l, &motif_i);
                    add_sj(sj, sj_i, sj_m, tid, end+1, end+l, strand, motif_i, 1, is_uniq); ++sj_i;
                    // ad
                    ad->intr_end[ad->intv_n] = end+l;
                    ad->exon_end[(ad->intv_n)++] = end;

                    start = end+l+1;
                }
                end += l;
                break;
            case BAM_CMATCH: // 1 1
            case BAM_CEQUAL:
            case BAM_CDIFF:
                end += l;
                break;
            case BAM_CDEL : // D(0 1)
#ifdef _RMATS_
                return 0; // for rMATS, only M&N allowed
#else
                end += l;
                break;
#endif
            case BAM_CINS: // 1 0
            case BAM_CSOFT_CLIP:
            case BAM_CHARD_CLIP:
            case BAM_CPAD: // 0 0
            case BAM_CBACK:
#ifdef _RMATS_
                return 0; // for rMATS, only M&N allowed
#endif
                break; 
            default:
                err_printf("Error: unknown cigar type: %d.\n", bam_cigar_op(c[i]));
                break;
        }
    }
    *sj_n = SJ_n; *_end = end;

    // ad
    if (SJ_n > 0) {
        ad->exon_end[(ad->intv_n)++] = end;
        ad->end = end;
        (*ad_n)++;
    }
    return SJ_n;
}

// parse alignment details via BAM cigar
// build sg_ad_index to calculate cnt
int parse_bam_record(samFile *in, bam_hdr_t *h, bam1_t *b, kseq_t *seq, int seq_n, SG_group *sg_g, int *sg_ad_idx, ad_t **AD_group, int *AD_n, int AD_m, sj_t **SJ_group, int *SJ_n, int SJ_m, sg_para *sgp)
{
    err_printf(__func__, "parsing bam records...\n");
    int n_cigar; uint32_t *cigar;
    uint8_t is_uniq; int tid, bam_start, bam_end;
    // alignment detail
    *AD_n = 0;
    // junction
    *SJ_n = 0; int sj_n, sj_m = 1; sj_t *sj = (sj_t*)_err_malloc(sizeof(sj_t));
    int last_sg_i = 0, sg_i, idx_sg_i = 0;
    // read bam record
    int ret;
    while (1) {
        ret = sam_read1(in, h, b);
        if (ret == -1) break;
        else if (ret < 0) err_fatal_simple("bam file error!\n");

        if (bam_unmap(b)) continue; // unmap (0)
        is_uniq = bam_is_uniq_NH(b); // uniq-map (1)
#ifdef _RMATS_
        if (is_uniq == 0) continue;
#endif
        if (bam_is_prop(b) != 1 && sgp->read_type == PAIR_T) continue; // prop-pair (2)

        tid = b->core.tid; n_cigar = b->core.n_cigar; cigar = bam_get_cigar(b);
        bam_start = b->core.pos+1;
        // alignment details
        if (parse_bam(tid, bam_start, &bam_end, n_cigar, cigar, is_uniq, seq, seq_n, AD_group, AD_n, &AD_m, &sj, &sj_n, &sj_m, sgp) < 0) continue;
       
        if (sj_n > 0) {
            // junction read
            sj_update_group(SJ_group, SJ_n, &SJ_m, sj, sj_n);
            // set sg_ad_idx
            while (idx_sg_i < sg_g->SG_n) {
                if (tid == sg_g->SG[idx_sg_i]->tid && bam_end >= sg_g->SG[idx_sg_i]->start && bam_start <= sg_g->SG[idx_sg_i]->end) {
                    sg_ad_idx[idx_sg_i] = *AD_n; // 0: no idx
                    idx_sg_i++;
                } else if (tid > sg_g->SG[idx_sg_i]->tid || (tid == sg_g->SG[idx_sg_i]->tid && bam_start > sg_g->SG[idx_sg_i]->end))
                    idx_sg_i++;
                else break; // bam_end < sg.start
            }
        } else {
            // exon-body read
            if (last_sg_i == sg_g->SG_n) continue;
            for (sg_i = last_sg_i; sg_i < sg_g->SG_n; ++sg_i) {
                SG *sg = sg_g->SG[sg_i];
                int start = sg->start, end = sg->end;
                int j;
                if (sg->tid < tid || start == MAX_SITE || end == 0 || (sg->tid == tid && end <= bam_start)) {
                    if (sg_i == last_sg_i) last_sg_i++; continue;
                } else if (sg->tid > tid || (sg->tid == tid && start >= bam_end)) break;
                else {
                    for (j = 0; j < sg->node_n; ++j) {
                        SGnode n = sg->node[j];
                        if (n.start <= bam_start && n.end >= bam_end) {
                            if (is_uniq) sg->node[j].uniq_c++;
                            else sg->node[j].multi_c++;
                        }
                    }
                }
            }
        }
    }
    free(sj);
    err_func_format_printf(__func__, "parsing bam records done!\n");

    return *SJ_n;
}

int bam2cnt_core(samFile *in, bam_hdr_t *h, bam1_t *b, kseq_t *seq, int seq_n, sj_t **SJ_group, int SJ_m, sg_para *sgp) {
    err_func_format_printf(__func__, "calculating junction- and exon-body-read count ...\n");
    int n_cigar; uint32_t *cigar;
    uint8_t is_uniq; int tid, bam_start, bam_end;
    // junction
    int SJ_n = 0, sj_n, sj_m = 1; sj_t *sj = (sj_t*)_err_malloc(sizeof(sj_t));
    // read bam record
    int ret;
    while (1) {
        ret = sam_read1(in, h, b);
        if (ret == -1) break;
        else if (ret < 0) err_fatal_simple("bam file error!\n");

        if (bam_unmap(b)) continue; // unmap (0)
        is_uniq = bam_is_uniq_NH(b); // uniq-map (1)
#ifdef _RMATS_
        if (is_uniq == 0) continue;
#endif
        if (bam_is_prop(b) != 1 && sgp->read_type == PAIR_T) continue; // prop-pair (2)

        tid = b->core.tid; n_cigar = b->core.n_cigar, cigar = bam_get_cigar(b);
        bam_start = b->core.pos+1, bam_end = b->core.pos+bam_cigar2rlen(n_cigar, cigar);
        // junction read
        if ((sj_n = gen_sj(is_uniq, tid, bam_start, n_cigar, cigar, seq, seq_n, &sj, &sj_m, sgp)) > 0) sj_update_group(SJ_group, &SJ_n, &SJ_m, sj, sj_n);
    }
    free(sj);
    err_func_format_printf(__func__, "calculating junction- and exon-body-read count done!\n");

    return SJ_n;
}

int bam2sj_core(samFile *in, bam_hdr_t *h, bam1_t *b, kseq_t *seq, int seq_n, sj_t **SJ_group, int SJ_m, sg_para *sgp)
{
    err_func_format_printf(__func__, "generating splice-junction with BAM file ...\n");
    int n_cigar; uint32_t *cigar;
    uint8_t is_uniq; int tid, bam_start, bam_end;
    int SJ_n = 0, sj_n, sj_m = 1; sj_t *sj = (sj_t*)_err_malloc(sizeof(sj_t));

    int ret;
    while (1) {
        ret = sam_read1(in, h, b);
        if (ret == -1) break;
        else if (ret < 0) err_fatal_simple("bam file error!\n");

        if (bam_unmap(b)) continue; // unmap (0)
        is_uniq = bam_is_uniq_NH(b); // uniq-map (1)
#ifdef _RMATS_
        if (is_uniq == 0) continue;
#endif
        if (bam_is_prop(b) != 1 && sgp->read_type == PAIR_T) continue; // prop-pair (2)

        tid = b->core.tid; n_cigar = b->core.n_cigar, cigar = bam_get_cigar(b);
        bam_start = b->core.pos+1, bam_end = b->core.pos+bam_cigar2rlen(n_cigar, cigar);
        if ((sj_n = gen_sj(is_uniq, tid, bam_start, n_cigar, cigar, seq, seq_n, &sj, &sj_m, sgp)) > 0) sj_update_group(SJ_group, &SJ_n, &SJ_m, sj, sj_n);
    }
    free(sj);
    err_func_format_printf(__func__, "generating splice-junction with BAM file done!\n");

    return SJ_n;
}

/*int sj_filter(sj_t *sj_group, int sj_n, sg_para *sgp, SG_group *sg_g)
{
    err_func_format_printf(__func__, "filtering splice-junctions ...\n");

    int sj_i = 0, last_sg_i = 0, sg_i;
    int known, hit;
    sj_t *sj;
    int i, anc_len;

    while (sj_i < sj_n) {
        sj = sj_group + sj_i;
        known = 0;
        while (last_sg_i < sg_g->SG_n) {
            if (sg_g->SG[last_sg_i]->node_n <= 3 || sg_g->SG[last_sg_i]->start == MAX_SITE || sg_g->SG[last_sg_i]->end == 0) { ++last_sg_i; continue; }
            int comp_res = comp_sj_sg(*sj, *(sg_g->SG[last_sg_i]));
            if (comp_res < 0) goto FILTER;
            else if (comp_res > 0) { ++last_sg_i; continue; }
            for (sg_i = last_sg_i; sg_i < sg_g->SG_n; ++sg_i) {
                if (comp_sj_sg(*sj, *(sg_g->SG[sg_i])) < 0) goto FILTER;
                SG *sg = sg_g->SG[sg_i];
                SGsite *don_site = sg->don_site, *acc_site = sg->acc_site; int32_t acc_n = sg->acc_site_n, don_n = sg->don_site_n;
                // search site/edge: (GTF_don_site_id, GTF_acc_site_id) => GTF_edge_id
                int GTF_don_site_id = sg_bin_sch_site(don_site, don_n, sj->don, &hit); if (hit == 0) continue;
                int GTF_acc_site_id = sg_bin_sch_site(acc_site, acc_n, sj->acc, &hit); if (hit == 0) continue;
                sg_bin_sch_edge(sg, GTF_don_site_id, GTF_acc_site_id, &hit);
                known = hit;
                goto FILTER;
            }
        }
FILTER:
        if (known == 1) { // known
            sj->is_anno = 1;
            anc_len = sgp->anchor_len[0];
        } else {  // novel
            sj->is_anno = 0;
            anc_len = sgp->anchor_len[(sj->motif+3)/2];// 0->1, 1,2->2, 3,4->3, 5,6->4 
        }
        // 1. anchor-len
        int filter_out=0;
        for (i = 0; i < sj->uniq_c; ++i) {
            if (sj->uniq_anc[i].left_anc_len < anc_len || sj->uniq_anc[i].right_anc_len < anc_len)
                filter_out++;
        }
        sj->uniq_c -= filter_out;
        filter_out=0;
        for (i = 0; i < sj->multi_c; ++i) {
            if (sj->multi_anc[i].left_anc_len < anc_len || sj->multi_anc[i].right_anc_len < anc_len)
                filter_out++;
        }
        sj->multi_c -= filter_out;
        // next sj
        sj_i++;
    }

    err_func_format_printf(__func__, "filtering splice-junctions done!\n");
    return sj_n;
}*/

void print_sj(sj_t *sj_group, int sj_n, FILE *out, char **cname)
{
    int i;
    fprintf(out, "###STRAND 0:undefined, 1:+, 2:-\n");
    fprintf(out, "###ANNO 0:novel, 1:annotated\n");
    fprintf(out, "###MOTIF 0:non-canonical, 1:GT/AG, 2:CT/AC, 3:GC/AG, 4:CT/GC, 5:AT/AC, 6:GT/AT\n");
    fprintf(out, "#CHR\tSTART\tEND\tSTRAND\tANNO\tUNIQ_C\tMULTI_C\tMOTIF\n");
    for (i = 0; i < sj_n; ++i) {
        sj_t sj = sj_group[i];
        fprintf(out, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", cname[sj.tid], sj.don, sj.acc, sj.strand, sj.is_anno, sj.uniq_c, sj.multi_c, sj.motif);
    }
}

int bam2sj(int argc, char *argv[])
{
    int c; char *p; char ref_fn[1024]="";
    sg_para *sgp = sg_init_para();
    FILE *gtf_fp=NULL;

    while ((c = getopt_long(argc, argv, "G:g:pa:i:A:U:", bam2sj_long_opt, NULL)) >= 0) {
        switch (c) {
            case 'g': strcpy(ref_fn, optarg); break;
            case 'G': gtf_fp = xopen(optarg, "r"); break;
            case 'p': sgp->read_type = PAIR_T; break;
            case 'a': sgp->anchor_len[0] = strtol(optarg, &p, 10);
                      if (*p != 0) sgp->anchor_len[1] = strtol(p+1, &p, 10); else return bam2sj_usage();
                      if (*p != 0) sgp->anchor_len[2] = strtol(p+1, &p, 10); else return bam2sj_usage();
                      if (*p != 0) sgp->anchor_len[3] = strtol(p+1, &p, 10); else return bam2sj_usage();
                      if (*p != 0) sgp->anchor_len[4] = strtol(p+1, &p, 10); else return bam2sj_usage();
                      break;
            case 'U': sgp->uniq_min[0] = strtol(optarg, &p, 10);
                      if (*p != 0) sgp->uniq_min[1] = strtol(p+1, &p, 10); else return bam2sj_usage();
                      if (*p != 0) sgp->uniq_min[2] = strtol(p+1, &p, 10); else return bam2sj_usage();
                      if (*p != 0) sgp->uniq_min[3] = strtol(p+1, &p, 10); else return bam2sj_usage();
                      if (*p != 0) sgp->uniq_min[4] = strtol(p+1, &p, 10); else return bam2sj_usage();
                      break; 
            case 'A': sgp->all_min[0] = strtol(optarg, &p, 10);
                      if (*p != 0) sgp->all_min[1] = strtol(p+1, &p, 10); else return bam2sj_usage();
                      if (*p != 0) sgp->all_min[2] = strtol(p+1, &p, 10); else return bam2sj_usage();
                      if (*p != 0) sgp->all_min[3] = strtol(p+1, &p, 10); else return bam2sj_usage();
                      if (*p != 0) sgp->all_min[4] = strtol(p+1, &p, 10); else return bam2sj_usage();
                      break;
            case 'i': sgp->intron_len = atoi(optarg); break;

            default: err_printf("Error: unknown option: %s.\n", optarg); return bam2sj_usage();
        }
    }
    if (argc - optind != 1) return bam2sj_usage();

    int seq_n = 0, seq_m; kseq_t *seq;
    if (strlen(ref_fn) != 0) {
        gzFile genome_fp = gzopen(ref_fn, "r");
        if (genome_fp == NULL) { err_fatal(__func__, "Can not open genome file. %s\n", ref_fn); }
        seq = kseq_load_genome(genome_fp, &seq_n, &seq_m);
        err_gzclose(genome_fp); 
    }
    // set cname
    chr_name_t *cname = chr_name_init();

    // open bam and parse bam header
    samFile *in; bam_hdr_t *h; bam1_t *b;
    if ((in = sam_open(argv[optind], "rb")) == NULL) err_fatal_core(__func__, "Cannot open \"%s\"\n", argv[optind]);
    if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", argv[optind]);
    bam_set_cname(h, cname);
    b = bam_init1(); 
    // build splice-graph
    SG_group *sg_g = NULL;
    if (gtf_fp != NULL) {
        sg_g = construct_SpliceGraph(gtf_fp, cname);
        err_fclose(gtf_fp);
    } chr_name_free(cname);

    sj_t *sj_group = (sj_t*)_err_malloc(10000 * sizeof(sj_t)); int sj_m = 10000;
    int sj_n = bam2sj_core(in, h, b, seq, seq_n, &sj_group, sj_m, sgp);
    //sj_filter(sj_group, sj_n, sgp, sg_g);

    print_sj(sj_group, sj_n, stdout, h->target_name);

    bam_destroy1(b); sam_close(in); bam_hdr_destroy(h); 
    sg_free_para(sgp); free(sj_group);
    int i; for (i = 0; i < seq_n; ++i) { free(seq[i].name.s); free(seq[i].seq.s); } free(seq);
    return 0;
}
