#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "gtf.h"

int usage(void)
{
    err_printf("\n");
    err_printf("Usage:   comp_gtf 1.gtf 2.gtf\n\n");
	err_printf("         1.gtf is smaller than 2.gtf\n");
	err_printf("\n");
	return 1;
}

// number of trans in T1 is less than T2
int comp_gtf_core(read_trans_t T1, read_trans_t T2)
{
    int i=0, j=0, last_j = 0, iden=0, dis=10;

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
    char line[1024], ref[100]="\0", type[20]="\0"; int start, end; char strand, add_info[500], tname[100];
    trans_t *t = trans_init(1);
    while (fgets(line, 1024, fp) != NULL) {
        sscanf(line, "%*s\t%*s\t%s", type);
        if (strcmp(type, "transcript") == 0) {
            if (t->exon_n != 0) {
                add_read_trans(T, *t);
                set_trans(T->t+T->trans_n-1, tname);
            }
            t->exon_n = 0;
        } else if (strcmp(type, "exon") == 0) { // exon
            sscanf(line, "%s\t%*s\t%s\t%d\t%d\t%*s\t%c\t%*s\t%[^\n]", ref, type, &start, &end, &strand, add_info);
            uint8_t is_rev = (strand == '-' ? 1 : 0);
            char tag[20]="transcript_id";
            gtf_add_info(add_info, tag, tname);
            add_exon(t, name2id(ref), start, end, is_rev);
        }
    }
    if (t->exon_n != 0) {
        add_read_trans(T, *t);
        set_trans(T->t+T->trans_n-1, tname);
    }
    trans_free(t);
    return T->trans_n;
}

#ifdef COMP_MAIN
int main(int argc, char *argv[])
{
    if (argc != 3) return usage();
    FILE *fp1 = fopen(argv[1], "r"), *fp2 = fopen(argv[2], "r");
    read_trans_t *T1, *T2;
    T1 = read_trans_init(), T2 = read_trans_init();
    read_anno_trans1(T1, fp1); read_anno_trans1(T2, fp2);

    comp_gtf_core(*T1, *T2);

    read_trans_free(T1), read_trans_free(T2);
    fclose(fp1), fclose(fp2);
    return 0;
}
#endif
