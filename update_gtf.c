/* update_gtf.c
 *   generate junction information based on sam/bam file
 *   then, update existing GTF file
 *   currently, only work for single-end long read data
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
extern int read_bam_trans(samFile *in, bam_hdr_t *h, bam1_t *b, int exon_min, read_trans_t *T);
extern int read_intron_group(intron_group_t *I, FILE *fp);
extern int read_anno_trans1(read_trans_t *T, FILE *fp);

int update_gtf_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   %s update-gtf [option] <in.bam> <old.gtf> > new.gtf\n\n", PROG);
    err_printf("Notice:  the BAM and GTF files should be sorted in advance.\n\n");
    err_printf("Options:\n\n");
    err_printf("         -i --intron      [STR]    intron information file output by STAR(*.out.tab). [NONE]\n");
    err_printf("         -e --min-exon    [INT]    minimum length of internal exon. [%d]\n", INTER_EXON_MIN_LEN);
    err_printf("         -d --distance    [INT]    consider same if distance between two splice site is not bigger than d. [%d]\n", SPLICE_DISTANCE);
    err_printf("         -l --full-length [INT]    level of strict criterion for considering full-length transcript. \n");
    err_printf("                                   (1->5, most strict->most relaxed) [%d]\n", 5);
    err_printf("         -u --unclassified         output UNCLASSIFIED novel transcript. [false]\n");
    err_printf("         -s --source      [STR]    source field in GTF, program, database or project name. [NONE]\n");
    err_printf("         -f --full-gtf    [STR]    use this option to output the full GTF information of SAM/BAM to file. [false].\n");
	err_printf("\n");
	return 1;
}

int exon_overlap(exon_t e1, exon_t e2)
{
    if (e1.start > e2.end || e2.start > e1.end) return 0;
    return 1;
}

int check_full(trans_t t, trans_t anno_t, int level)
{
    int i = t.exon_n-1, j = anno_t.exon_n-1;
    if (level == 1) { // identical first and last splice-site
        if (t.is_rev) {
            if (t.exon[i].end != anno_t.exon[j].end) return 0;
            if (t.exon[0].start != anno_t.exon[0].start) return 0;
        } else {
            if (t.exon[0].end != anno_t.exon[0].end) return 0;
            if (t.exon[i].start != anno_t.exon[j].start) return 0;
        }
    } else if (level == 2) { // overlapping first and last exon
        if (!exon_overlap(t.exon[0], anno_t.exon[0])) return 0;
        if (!exon_overlap(t.exon[i], anno_t.exon[j])) return 0;
    } else if (level == 3) { // overlapping first and last exon, or overlapping nothing
        int ii, a1, a2, a3, a4;
        a1 = 0, a3 = 0;
        if (exon_overlap(t.exon[0], anno_t.exon[0])) a1=1;
        if (exon_overlap(t.exon[i], anno_t.exon[j])) a3=1;
        if (a1 == 1 && a3 == 1) return 1;

        a2 = 1, a4 = 1;
        for (ii = 0; ii < anno_t.exon_n; ++ii) {
            if (exon_overlap(t.exon[0], anno_t.exon[ii])) a2 = 0;
            if (exon_overlap(t.exon[i], anno_t.exon[ii])) a4 = 0;
            if (a2 == 0 && a4 == 0) return 0;
        }
        if (a1+a2 >= 1 && a3+a4 >= 1) return 1;
        else return 0;
    } else if (level == 4) { // most 5' exon meets #3, most 3' exon has a polyA+ tail of 15bp or longer XXX
        if (t.is_rev) {
            i = t.exon_n-1; j = anno_t.exon_n-1;
        } else {
            i = 0; j = 0;
        }
        if (exon_overlap(t.exon[i], anno_t.exon[j])) return 1;

        int ii;
        for (ii = 0; ii < anno_t.exon_n; ++ii) {
            if (exon_overlap(t.exon[i], anno_t.exon[ii])) return 0;
        }
        return 1;
    }
    return 1;
}

// check if t is novel and has identical splice site
// if t has all identical splice sites with other novel t, merge two ends
// @return value
//    0: novel, NOT share any identical splice site (unclassified)
//    1: novel, and share identical splice site (gene_id)
//    2: totally identical, can NOT be added to any anno
//    3: other cases that cannot be added to this anno(not full-length to any anno-trans)
int check_novel1(trans_t *bam_t, trans_t anno_t, int dis, int l)
{
    if (bam_t->is_rev != anno_t.is_rev || bam_t->exon_n < 2 || check_full(*bam_t, anno_t, l) == 0) return 3;

    int left=0,right=0, last_j=-1;
    int i, j, iden_n=0, iden_intron_n = 0, not_iden_iden=0;
    if (bam_t->is_rev) { // '-' strand
        for (i = 0; i < bam_t->exon_n-1; ++i) {
            for (j = 0; j < anno_t.exon_n-1; ++j) {
                if (abs(bam_t->exon[i].start - anno_t.exon[j].start) <= dis) left=1;
                if (abs(bam_t->exon[i+1].end - anno_t.exon[j+1].end) <= dis) right=1;
                if (left+right==2) {
                    if (last_j != -1 && j != last_j+1) not_iden_iden = 1;
                    last_j = j;
                    iden_intron_n += 1;
                }
                iden_n += (left+right);
                left=right=0;

                if (anno_t.exon[j+1].end < bam_t->exon[i+1].end) break;
            }
        }
    } else { // '+' strand
        for (i = 0; i < bam_t->exon_n-1; ++i) {
            for (j = 0; j < anno_t.exon_n-1; ++j) {
                if (abs(bam_t->exon[i].end - anno_t.exon[j].end) <= dis) left=1;
                if (abs(anno_t.exon[j+1].start - bam_t->exon[i+1].start) <= dis) right=1;
                if (left+right==2) {
                    if (last_j != -1 && j != last_j+1) not_iden_iden = 1;
                    last_j = j;
                    iden_intron_n += 1;
                }
                iden_n += (left+right);
                left=right=0;

                if (anno_t.exon[j+1].start > bam_t->exon[i+1].start) break;
            }
        }
    }
    if (iden_intron_n == bam_t->exon_n-1 && not_iden_iden == 0) return 2;
    if (iden_n > 0) {
        strcpy(bam_t->gname, anno_t.gname);
        return 1;
    } else {
        return 0;
    }
}

int check_intron1(int tid, int start, int end, uint8_t is_rev, intron_group_t I, int i_start, int dis)
{
    int i = i_start;
    while (i < I.intron_n) {
        if (I.intron[i].tid > tid || I.intron[i].start > start) return 0;
        if (I.intron[i].is_rev != is_rev) { i++; continue; }

        if (abs(I.intron[i].start-start)<=dis && abs(I.intron[i].end-end)<=dis) return 1;
        else i++;
    }
    return 0;
}

int check_intron(trans_t bam_t, int *intron_map, intron_group_t I, int *intron_i, int dis)
{
    int i = *intron_i, j;
    while (i < I.intron_n) {
        if (I.intron[i].tid < bam_t.tid || (I.intron[i].tid == bam_t.tid && I.intron[i].end <= bam_t.start)) {
            i++; *intron_i = i;
        } else if (I.intron[i].tid > bam_t.tid) return 0;
        else {
            if (bam_t.is_rev) { // '-' strand
                for (j = 0; j < bam_t.exon_n-1; ++j) {
                    if (intron_map[j] == 0 && check_intron1(bam_t.tid, bam_t.exon[j+1].end+1, bam_t.exon[j].start-1, bam_t.is_rev, I, i, dis) == 0)
                        return 0;
                }
            } else { // '+' strand
                for (j = 0; j < bam_t.exon_n-1; ++j) {
                    if (intron_map[j] == 0 && check_intron1(bam_t.tid, bam_t.exon[j].end+1, bam_t.exon[j+1].start-1, bam_t.is_rev, I, i, dis) == 0)
                        return 0;
                }
            }
            return 1;
        }
    }
    return 0;
}

int check_novel_intron(trans_t *bam_t, trans_t anno_t, intron_group_t I, int *intron_i, int dis, int l)
{
    if (bam_t->is_rev != anno_t.is_rev || bam_t->exon_n < 2 || check_full(*bam_t, anno_t, l) == 0) return 3;

    int left=0,right=0, last_j=-1;
    int i, j, iden_n=0, iden_intron_n = 0, not_iden_iden=0, ret;
    int *intron_map = (int*)_err_calloc((bam_t->exon_n-1), sizeof(int));
    if (bam_t->is_rev) { // '-' strand
        for (i = 0; i < bam_t->exon_n-1; ++i) {
            for (j = 0; j < anno_t.exon_n-1; ++j) {
                if (abs(bam_t->exon[i].start - anno_t.exon[j].start) <= dis) left=1;
                if (abs(bam_t->exon[i+1].end - anno_t.exon[j+1].end) <= dis) right=1;
                if (left+right==2) {
                    if (last_j != -1 && j != last_j+1) not_iden_iden = 1;
                    last_j = j;
                    iden_intron_n += 1;
                    intron_map[i] = 1;
                }
                iden_n += (left+right);
                left=right=0;

                if (anno_t.exon[j+1].end < bam_t->exon[i+1].end) break;
            }
        }
    } else { // '+' strand
        for (i = 0; i < bam_t->exon_n-1; ++i) {
            for (j = 0; j < anno_t.exon_n-1; ++j) {
                if (abs(bam_t->exon[i].end - anno_t.exon[j].end) <= dis) left=1;
                if (abs(anno_t.exon[j+1].start - bam_t->exon[i+1].start) <= dis) right=1;
                if (left+right==2) {
                    if (last_j != -1 && j != last_j+1) not_iden_iden = 1;
                    last_j = j;
                    iden_intron_n += 1;
                    intron_map[i] = 1;
                }
                iden_n += (left+right);
                left=right=0;

                if (anno_t.exon[j+1].start > bam_t->exon[i+1].start) break;
            }
        }
    }

    if (iden_intron_n == bam_t->exon_n-1 && not_iden_iden == 0) ret=2;
    else {
        if (iden_n > 0) {
            if (check_intron(*bam_t, intron_map, I, intron_i, dis)) {
                strcpy(bam_t->gname, anno_t.gname);
                ret=1;
            } else ret=3;
        } else {
            ret=0;
        }
    }
    free(intron_map);
    return ret;
}

int merge_trans1(trans_t t1, trans_t *t2, int dis)
{
    int i = t1.exon_n-1, j = t2->exon_n-1;
    if (check_iden(t1, *t2, dis)) {
        t2->cov++;
        if (t2->is_rev) { // '-'
            if (t1.exon[0].end > t2->exon[0].end) {
                t2->exon[0].end = t1.exon[0].end;
                t2->end = t1.exon[0].end;
            }
            if (t1.exon[i].start < t2->exon[j].start) {
                t2->exon[j].start = t1.exon[i].start;
                t2->start = t1.exon[i].start;
            }
        } else { // '+'
            if (t1.exon[0].start < t2->exon[0].start)  {
                t2->exon[0].start = t1.exon[0].start;
                t2->start = t1.exon[0].start;
            }
            if (t1.exon[i].end > t2->exon[j].end) {
                t2->exon[j].end = t1.exon[i].end;
                t2->end = t1.exon[i].end;
            }
        }
        return 1;
    } else return 0;
}

int merge_trans(trans_t t, read_trans_t *T, int dis)
{
    int i; 
    for (i = T->trans_n-1; i >= 0; --i) {
        if (merge_trans1(t, T->t+i, dis)) return 1;
        if (t.tid > T->t[i].tid || t.start > T->t[i].end) return 0;
    }
    return 0;
}

int check_novel_trans(read_trans_t bam_T, read_trans_t anno_T, intron_group_t I, 
                      int uncla, int dis, int l, read_trans_t *novel_T)
{
    int i=0, j=0, last_j=0, k=0, ret;
    int all_novel=0, novel=0;
    while (i < bam_T.trans_n && j < anno_T.trans_n) {
        //int x;
        //if (strcmp(bam_T.t[i].tname, "m130614_000849_42175_c100535482550000001823081711101343_s1_p0/87730/ccs.path1")==0)
            //x=i;
        if (merge_trans(bam_T.t[i], novel_T, dis)) { 
            err_printf("merge: %s\n", bam_T.t[i].tname);
            i++; continue; 
        }
        if (bam_T.t[i].tid > anno_T.t[j].tid || (bam_T.t[i].tid == anno_T.t[j].tid && bam_T.t[i].start > anno_T.t[j].end)) {
            j++;
            last_j = j;
        } else if (anno_T.t[j].tid > bam_T.t[i].tid || (anno_T.t[j].tid == bam_T.t[i].tid && anno_T.t[j].start > bam_T.t[i].end)) {
            if (novel != 1) err_printf("unrecover: %s\n", bam_T.t[i].tname);
            if (novel == 1) {
                add_read_trans(novel_T, bam_T.t[i]);
                set_trans(novel_T->t+novel_T->trans_n-1, NULL);
            } else if (uncla == 1 && all_novel == 1) {
                add_read_trans(novel_T, bam_T.t[i]);
                set_trans(novel_T->t+novel_T->trans_n-1, NULL);
            }
            i++;
            novel = 0;
        } else {
            if (I.intron_n > 0) ret = check_novel_intron(bam_T.t+i, anno_T.t[j], I, &k, dis, l);
            else ret = check_novel1(bam_T.t+i, anno_T.t[j], dis, l);
            if (ret == 0) { // all novel
                all_novel = 1;
                j++; continue;
            } else if (ret == 1) { // novel
                novel = 1;
                j++; continue;
            } else if (ret == 2) { // all identical
                novel = 0;
                err_printf("all-iden: %s\n", bam_T.t[i].tname);
                i++;
            } else {
                j++; continue;
            }
        }
        j = last_j;
    }
    return 0;
}

// read
int read_anno_trans(FILE *fp, bam_hdr_t *h, read_trans_t *T)
{
    char line[1024], ref[100]="\0", type[20]="\0"; int start, end; char strand, add_info[1024], gname[100];
    trans_t *t = trans_init(1);
    while (fgets(line, 1024, fp) != NULL) {
        sscanf(line, "%s\t%*s\t%s\t%d\t%d\t%*s\t%c\t%*s\t%[^\n]", ref, type, &start, &end, &strand, add_info);
        uint8_t is_rev = (strand == '-' ? 1 : 0);
        if (strcmp(type, "transcript") == 0) {
            if (t->exon_n > 1) {
                add_read_trans(T, *t);
                set_trans(T->t+T->trans_n-1, NULL);
            }
            t->exon_n = 0;
            char tag[10]="gene_id";
            gtf_add_info(add_info, tag, gname);
            strcpy(t->gname, gname);
        } else if (strcmp(type, "exon") == 0) { // exon
            add_exon(t, bam_name2id(h, ref), start, end, is_rev);
        }
    }
    if (t->exon_n != 0) {
        add_read_trans(T, *t);
        set_trans(T->t+T->trans_n-1, NULL);
    }
    trans_free(t);
    return T->trans_n;
}

const struct option update_long_opt [] = {
    { "intron", 1, NULL, 'i' },
    { "min-exon", 1, NULL, 'e' },
    { "distance", 1, NULL, 'd' },
    { "unclassified", 0, NULL, 'u' },
    { "source", 1, NULL, 's' },
    { "full-gtf", 1, NULL, 'f' },

    { 0, 0, 0, 0}
};

int update_gtf(int argc, char *argv[])
{
    int c; int exon_min = INTER_EXON_MIN_LEN, dis=SPLICE_DISTANCE, l=5, uncla = 0; char src[100]="NONE"; FILE *new_gfp=stdout, *full_gfp=NULL, *intron_fp=NULL;
	while ((c = getopt_long(argc, argv, "i:e:d:l:us:f:", update_long_opt, NULL)) >= 0) {
        switch(c)
        {
            case 'i': if ((intron_fp = fopen(optarg, "r")) == NULL) {
                          err_fatal(__func__, "Can not open intron file \"%s\"\n", optarg);
                          return update_gtf_usage();
                      } 
                      break;
            case 'e': exon_min = atoi(optarg); break;
            case 'd': dis = atoi(optarg); break;
            case 'l': l = atoi(optarg); break;
            case 'u': uncla = 1;
            case 's': strcpy(src, optarg); break;
            case 'f': if ((full_gfp = fopen(optarg, "w")) == NULL) {
                          err_fatal(__func__, "Can not open full-gtf output file \"%s\"\n", optarg);
                          return update_gtf_usage();
                      }
                      break;
            default:
                      err_printf("Error: unknown option: %s.\n", optarg);
                      return update_gtf_usage();
                      break;
        }
    }
    if (argc - optind != 2) return update_gtf_usage();

    samFile *in; bam_hdr_t *h; bam1_t *b; read_trans_t *anno_T, *bam_T, *novel_T;
    if ((in = sam_open(argv[optind], "rb")) == NULL) err_fatal(__func__, "Cannot open \"%s\"\n", argv[optind]);
    if ((h = sam_hdr_read(in)) == NULL) err_fatal(__func__, "Couldn't read header for \"%s\"\n", argv[optind]);
    b = bam_init1(); 
    anno_T = read_trans_init(); bam_T = read_trans_init(); novel_T = read_trans_init();
    intron_group_t *I = intron_group_init();

    FILE *gfp = fopen(argv[optind+1], "r");
    // read all gene
    read_anno_trans(gfp, h, anno_T);
    // read all transcript
    read_bam_trans(in, h, b, exon_min, bam_T);
    // read intron file
    read_intron_group(I, intron_fp);
    
    // merge loop
    check_novel_trans(*bam_T, *anno_T, *I, uncla, dis, l, novel_T);

    // print
    print_read_trans(*novel_T, h, src, new_gfp);

    read_trans_free(anno_T); read_trans_free(bam_T); read_trans_free(novel_T); intron_group_free(I);
    bam_destroy1(b); bam_hdr_destroy(h); sam_close(in); fclose(gfp); if(full_gfp) fclose(full_gfp); if (intron_fp) fclose(intron_fp);
    return 0;
}
