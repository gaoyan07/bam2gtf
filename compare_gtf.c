#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "gtf.h"

int usage(void)
{
    err_printf("\n");
    err_printf("Usage:   comp_gtf 1.gtf 2.gtf short.intron.out\n\n");
	err_printf("         1.gtf is smaller than 2.gtf\n");
	err_printf("\n");
	return 1;
}

int check_iden_intron(trans_t t, intron_group_t I, int start, int dis)
{
    int i = 0, j = start;

    if (t.is_rev == 0) { // '+' strand
        while (i < t.exon_n-1 && j < I.intron_n) {
            if (I.intron[j].start >= t.end) break;
            if (I.intron[j].is_rev != t.is_rev) { j++; continue; }

            if (t.exon[i].end+1 < (I.intron[j].start-dis)) return 0;

            if (t.exon[i].end+1 > (I.intron[j].start+dis)) j++;
            else if (abs(t.exon[i+1].start-1 - I.intron[j].end) <= dis) {
                i++;
                j++;
            } else j++;
        }
        return (i == t.exon_n-1);
    } else { // '-' strand
        i = t.exon_n-1;
        while (i > 0 && j < I.intron_n) {
            if (I.intron[j].start >= t.end) break;
            if (I.intron[j].is_rev != t.is_rev) { j++; continue; }

            if (t.exon[i].end+1 < (I.intron[j].start-dis)) return 0;

            if (t.exon[i].end+1 > (I.intron[j].start+dis)) j++;
            else if (abs(t.exon[i-1].start-1 - I.intron[j].end) <= dis) {
                i--;
                j++;
            } else j++;
        }
        return (i==0);
    }
}

int comp_gtf_intron(read_trans_t T, intron_group_t I)
{
    int i, j, iden, dis;

    //for(i = 0; i < T.trans_n; ++i)
        //printf("all: %s\t%d\n", T.t[i].tname, T.t[i].cov);
    i = 0, j = 0, iden=0, dis=0;
    while (i < T.trans_n && j < I.intron_n) {
        //if (strcmp(T.t[i].tname, "m130614_092349_42175_c100535482550000001823081711101347_s1_p0/96492/ccs.path1") == 0)
        //    printf("ok");
        if (I.intron[j].tid < T.t[i].tid || (I.intron[j].tid == T.t[i].tid && I.intron[j].end <= T.t[i].start)) {
            j++;
        } else if (I.intron[j].tid > T.t[i].tid || (I.intron[j].tid == T.t[i].tid && I.intron[j].start >= T.t[i].end)) {
            // un-recovered
            printf("un-recovered: %s\n", T.t[i].tname);
            i++;
        } else {
            if (check_iden_intron(T.t[i], I, j, dis)) {
                iden++;
                // recovered
                //printf("recover: %s\n", T.t[i].tname);
            } else // un-recovered
                printf("un-recovered: %s\n", T.t[i].tname);
            i++;
        }
    }
    printf("Recovered: %d\n#1 Only: %d\n", iden, T.trans_n-iden);
    return 0;
}

// number of trans in T1 is less than T2
int comp_gtf_core(read_trans_t T1, read_trans_t T2)
{
    int i, j, last_j = 0, iden=0, dis=0;

    for(i = 0; i < T2.trans_n; ++i)
        printf("all: %s\t%d\n", T2.t[i].tname, T2.t[i].cov);
    i = 0, j = 0;
    while (i < T1.trans_n && j < T2.trans_n) {
        if (T1.t[i].tid > T2.t[j].tid || (T1.t[i].tid == T2.t[j].tid && T1.t[i].start > T2.t[j].end)) {
            j++;
            last_j = j;
        } else if (T2.t[j].tid > T1.t[i].tid || (T2.t[j].tid == T1.t[i].tid && T2.t[j].start > T1.t[i].end)) {
            err_printf("%s\t%d\t%d\t%d\n", T1.t[i].tname, T1.t[i].tid, T1.t[i].start, T1.t[i].end);
            i++;
        } else {
            if (check_iden(T1.t[i], T2.t[j], dis)) {
                iden++;
                i++;
                printf("both: %s\t%d\n", T2.t[j].tname, T2.t[j].cov);
                //last_j = j+1;
            } else {
                j++;
                continue;
            }
        }
        j = last_j;
    }
    printf("Overlap: %d\n#1 Only: %d\n#2 Only: %d\n", iden, T1.trans_n-iden, T2.trans_n-iden);
    return 0;
}

int name2id(char ref[])
{
    if (ref[3] == 'X') return 23;
    else if (ref[3] == 'Y') return 24;
    else if (ref[3] == 'M') return 25;
    else return atoi(ref+3);
}

int read_anno_trans1(read_trans_t *T, FILE *fp)
{
    char line[1024], ref[100]="\0", type[20]="\0"; int start, end, cov=1; char strand, add_info[500], tname[100], qual[10];
    trans_t *t = trans_init(1);
    while (fgets(line, 1024, fp) != NULL) {
        sscanf(line, "%*s\t%*s\t%s", type);
        if (strcmp(type, "transcript") == 0) {
            if (t->exon_n != 0) {
                add_read_trans(T, *t);
                set_trans(T->t+T->trans_n-1, tname);
                T->t[T->trans_n-1].cov = cov;
            }
            t->exon_n = 0;
        } else if (strcmp(type, "exon") == 0) { // exon
            sscanf(line, "%s\t%*s\t%s\t%d\t%d\t%s\t%c\t%*s\t%[^\n]", ref, type, &start, &end, qual, &strand, add_info);
            if (qual[0] != '-') cov = (atoi(qual)-450)/50;
            uint8_t is_rev = (strand == '-' ? 1 : 0);
            char tag[20]="transcript_id";
            gtf_add_info(add_info, tag, tname);
            add_exon(t, name2id(ref), start, end, is_rev);
        }
    }
    if (t->exon_n != 0) {
        add_read_trans(T, *t);
        set_trans(T->t+T->trans_n-1, tname);
        T->t[T->trans_n-1].cov = cov;
    }
    trans_free(t);
    return T->trans_n;
}

int read_intron_group(intron_group_t *I, FILE *fp)
{
    char line[1024], ref[100]; int start, end, nstrand, canon, anno, uniq_map, multi_map, overlang;
    intron_t *i = intron_init(1);
    while (fgets(line, 1024, fp) != NULL) {
        sscanf(line, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d", ref, &start, &end, &nstrand, &canon, &anno, &uniq_map, &multi_map, &overlang);
        i->tid = name2id(ref); i->start = start, i->end = end;
        i->is_rev = (nstrand == 1 ? 0 : (nstrand == 2 ? 1 : -1)); i->is_canon = canon;
        add_intron(I, *i);
    }
    free(i);
    return I->intron_n;
}

#ifdef COMP_MAIN
int main(int argc, char *argv[])
{
    if (argc != 4) return usage();
    FILE *fp1 = fopen(argv[1], "r"), *fp2 = fopen(argv[2], "r"), *fp3 = fopen(argv[3], "r");
    read_trans_t *T1, *T2; intron_group_t *I;
    T1 = read_trans_init(), T2 = read_trans_init(); I = intron_group_init();
    read_anno_trans1(T1, fp1); /*read_anno_trans1(T2, fp2);*/ read_intron_group(I, fp3);

    //comp_gtf_core(*T1, *T2);
    comp_gtf_intron(*T1, *I);
    //comp_core(*T1, *T2, *I);

    read_trans_free(T1), read_trans_free(T2), intron_group_free(I);
    fclose(fp1), fclose(fp2), fclose(fp3);
    return 0;
}
#endif
