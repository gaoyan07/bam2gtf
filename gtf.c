#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gtf.h"
#include "utils.h"
#include "htslib/htslib/sam.h"


// exon
exon_t *exon_init(int n) {
    exon_t *e = (exon_t*)_err_malloc(n * sizeof(exon_t));
    return e;
}
void exon_free(exon_t *e) { free(e); } 

//transcript
trans_t *trans_init(int n) { 
    trans_t *t = (trans_t*)_err_malloc(n * sizeof(trans_t));
    strcpy(t->tname, "");
    t->exon_n = 0; t->exon_m = 2;
    t->exon = exon_init(2);
    return t;
}

int add_exon(trans_t *t, int32_t tid, int32_t start, int32_t end, uint8_t is_rev)
{
    if (t->exon_n == t->exon_m) t = exon_realloc(t);
    t->exon[t->exon_n].tid = tid;
    t->exon[t->exon_n].start = start;
    t->exon[t->exon_n].end = end;
    t->exon[t->exon_n].is_rev = is_rev;
    t->exon_n++;
    return 0;
}

// '+': s->e => S->E
// '-': S->E => s->e
int sort_exon(trans_t *t)
{
    int i, j; exon_t e1, e2; int res=0;
    for (i = 0; i < t->exon_n-1; ++i) {
        for (j = i+1; j < t->exon_n; ++j) {
            e1 = t->exon[i], e2 = t->exon[j];
            if ((e1.is_rev == e2.is_rev && ((!e1.is_rev && e1.start > e2.start) || (e1.is_rev && e1.start < e2.start))) //same strand
             || (e1.is_rev && !e2.is_rev)) { // '-' and '+'
                exon_t tmp = e1;
                t->exon[i] = t->exon[j];
                t->exon[j] = tmp;
            }
        }
    }
    return res;
}

int set_trans(trans_t *t, char *tname)
{
    sort_exon(t);
    t->tid = t->exon[0].tid;
    t->is_rev = t->exon[0].is_rev;
    t->start = t->is_rev ? t->exon[t->exon_n-1].start : t->exon[0].start;
    t->end = t->is_rev ? t->exon[0].end: t->exon[t->exon_n-1].end;
    if (tname) strcpy(t->tname, tname);
    return 0;
}

trans_t *exon_realloc(trans_t *t) {
    t->exon_m <<= 1;
    t->exon = (exon_t*)_err_realloc(t->exon, t->exon_m * sizeof(exon_t));
    return t;
}

void trans_free(trans_t *t) { free(t->exon); free(t); }

//for one read: multi-alignments => multi-transcripts
read_trans_t *read_trans_init(void)
{
    read_trans_t *r = (read_trans_t*)_err_malloc(sizeof(read_trans_t));
    r->trans_n = 0, r->trans_m = 1;
    r->t = trans_init(1);
    return r;
}

void add_read_trans(read_trans_t *r, trans_t t)
{
    if (r->trans_n == r->trans_m) r = read_trans_realloc(r);
    int i;
    r->t[r->trans_n].exon_n = 0;
    for (i = 0; i < t.exon_n; ++i)
        add_exon(r->t+r->trans_n, t.exon[i].tid, t.exon[i].start, t.exon[i].end, t.exon[i].is_rev);
    strcpy(r->t[r->trans_n].tname, t.tname);
    r->trans_n++;
}

read_trans_t *read_trans_realloc(read_trans_t *r)
{
    r->trans_m <<= 1;
    r->t = (trans_t*)_err_realloc(r->t, r->trans_m * sizeof(trans_t));
    int i;
    for (i = (r->trans_m >> 1); i < r->trans_m; ++i) {
        r->t[i].exon_n = 0, r->t[i].exon_m = 2;
        r->t[i].exon = exon_init(2);
    }
    return r;
}

void read_trans_free(read_trans_t *r)
{
    int i;
    for (i = 0; i < r->trans_m; ++i) free(r->t[i].exon);
    free(r->t); free(r);
}

//gene
gene_t *gene_init(void) {
    gene_t *g = (gene_t*)_err_malloc(sizeof(gene_t));
    g->trans_n = 0, g->trans_m = 1;
    g->trans = trans_init(1);
    return g; 
}

void add_trans(gene_t *g, trans_t t, int novel_gene_flag)
{
    if (g->trans_n == g->trans_m) g = trans_realloc(g);
    int i;
    g->trans[g->trans_n].tid = t.tid;
    g->trans[g->trans_n].is_rev = t.is_rev;
    g->trans[g->trans_n].start = t.start;
    g->trans[g->trans_n].end = t.end;
    g->trans[g->trans_n].novel_gene_flag = novel_gene_flag;
    strcpy(g->trans[g->trans_n].tname, t.tname);

    g->trans[g->trans_n].exon_n = 0;
    for (i = 0; i < t.exon_n; ++i) {
        add_exon(g->trans+g->trans_n, t.exon[i].tid, t.exon[i].start, t.exon[i].end, t.exon[i].is_rev);
    }
    g->trans_n++;
}

gene_t *trans_realloc(gene_t *g) {
    g->trans_m <<= 1;
    g->trans = (trans_t*)_err_realloc(g->trans, g->trans_m * sizeof(trans_t));
    int i;
    for (i = (g->trans_m >> 1); i < g->trans_m; ++i) {
        g->trans[i].exon_n = 0, g->trans[i].exon_m = 2;
        g->trans[i].exon = exon_init(2);
    }
    return g;
}

void gene_free(gene_t *g) {
    int i;
    for (i = 0; i < g->trans_m; ++i) {
        free(g->trans[i].exon);
    }
    free(g->trans); free(g);
}

// gtf additional information
void gtf_add_info(char add_info[], char tag[], char *info)
{
    for (int i=0; add_info[i] != '\0'; ++i) {
        if (strncmp(add_info+i, tag, strlen(tag)) == 0) {
            sscanf(add_info+i+strlen(tag)+2, "%[^\"]", info);
            return;
        }
    }
}

// gene_group
gene_group_t *gene_group_init(void)
{
    gene_group_t *gg = (gene_group_t*)_err_malloc(sizeof(gene_group_t));
    gg->start = gg->end = gg->gene_n = 0; gg->gene_m = 1;
    gg->g = gene_init();
    return gg;
}

gene_group_t *gene_group_realloc(gene_group_t *gg)
{
    gg->gene_m <<= 1;
    gg->g = (gene_t*)_err_realloc(gg->g, gg->gene_m * sizeof(gene_t));
    for (int i=gg->gene_m>>1; i < gg->gene_m; ++i) {
        gg->g[i].trans_n = 0; gg->g[i].trans_m = 1;
        gg->g[i].trans = trans_init(1);
    }
    return gg;
}

void add_gene(gene_group_t *gg, gene_t g, int novel_gene_flag)
{
    if (gg->gene_n == gg->gene_m) gg = gene_group_realloc(gg);
    for (int i = 0; i < g.trans_n; ++i)
        add_trans(gg->g+gg->gene_n, g.trans[i], novel_gene_flag);
    strcpy(gg->g[gg->gene_n].gname, g.gname);
    gg->g[gg->gene_n].tid = g.tid;
    gg->g[gg->gene_n].start = g.start;
    gg->g[gg->gene_n].end = g.end;
    gg->g[gg->gene_n].is_rev = g.is_rev;
    gg->gene_n++;
}

void set_gene_group(gene_group_t *gg)
{
    gg->tid =  gg->g[0].tid;
    gg->start = gg->g[0].start; 
    gg->end = gg->g[0].end;
    for (int i = 1; i < gg->gene_n; ++i) {
        if (gg->g[i].start < gg->start) gg->start = gg->g[i].start;
        if (gg->g[i].end > gg->end) gg->end = gg->g[i].end;
    }
}

void gene_group_free(gene_group_t *gg)
{
    int i, j;
    for (i = 0; i < gg->gene_m; ++i) {
        for (j = 0; j < gg->g[i].trans_m; ++j)
            free(gg->g[i].trans[j].exon);
        free(gg->g[i].trans);
    }
    free(gg->g); free(gg);
}

//print
int print_exon(exon_t e, FILE *out)
{
    return 0;
}

int print_trans(trans_t t, bam_hdr_t *h, char *src, FILE *out)
{
    int i;
    fprintf(out, "%s\t%s\t%s\t%d\t%d\t.\t%c\t.\ttranscript_id \"%s\";\n", h->target_name[t.tid], src, "transcript", t.start, t.end, "+-"[t.is_rev], t.tname);
    for (i = 0; i < t.exon_n; ++i)
        fprintf(out, "%s\t%s\t%s\t%d\t%d\t.\t%c\t.\ttranscript_id \"%s\";\n", h->target_name[t.tid], src, "exon", t.exon[i].start, t.exon[i].end, "+-"[t.exon[i].is_rev], t.tname);
    return 0;
}

// tid source feature start end score(.) strand phase(.) additional
int print_read_trans(read_trans_t r, bam_hdr_t *h, char *src, FILE *out)
{
    int i, j;
    for (i = 0; i < r.trans_n; ++i) {
        fprintf(out, "%s\t%s\t%s\t%d\t%d\t.\t%c\t.\ttranscript_id \"%s\";\n", h->target_name[r.t[i].tid], src, "transcript", r.t[i].start, r.t[i].end, "+-"[r.t[i].is_rev], r.t[i].tname);
        for (j = 0; j < r.t[i].exon_n; ++j)
            fprintf(out, "%s\t%s\t%s\t%d\t%d\t.\t%c\t.\ttranscript_id \"%s\";\n", h->target_name[r.t[i].exon[j].tid], src, "exon", r.t[i].exon[j].start, r.t[i].exon[j].end, "+-"[r.t[i].exon[j].is_rev], r.t[i].tname);
    }
    return 0;
}

int print_gene(gene_t g, FILE *out)
{
    return 0;
}

void print_gtf_trans(gene_t g, bam_hdr_t *h, char *src, FILE *out)
{
    if (g.trans_n <= g.anno_tran_n) return;
    int i, j; char gene_name[1024];
    for (i = g.anno_tran_n; i < g.trans_n; ++i) {
        if (g.trans[i].novel_gene_flag) strcpy(gene_name, "UNCLASSIFIED");
        else strcpy(gene_name, g.gname);
        trans_t t = g.trans[i];
        fprintf(out, "%s\t%s\t%s\t%d\t%d\t.\t%c\t.\tgene_id \"%s\"; transcript_id \"%s\";\n", h->target_name[t.tid], src, "transcript", t.start, t.end, "+-"[t.is_rev], gene_name, t.tname);
        for (j = 0; j < t.exon_n; ++j)
            fprintf(out, "%s\t%s\t%s\t%d\t%d\t.\t%c\t.\tgene_id \"%s\"; transcript_id \"%s\";\n", h->target_name[t.tid], src, "exon", t.exon[j].start, t.exon[j].end, "+-"[t.exon[j].is_rev], gene_name, t.tname);
    }
}

void print_gene_group(gene_group_t gg, bam_hdr_t *h, char *src, FILE *out, char **group_line, int *group_line_n)
{
    int l_i = 0;
    for (int i = 0; i < gg.gene_n; ++i) {
        // print anno
        for (int j = 0; j < group_line_n[i]; ++i)
            fprintf(out, "%s\n", group_line[l_i++]);
        // print novel trans
        print_gtf_trans(gg.g[i], h, src, out);        
    }
}

