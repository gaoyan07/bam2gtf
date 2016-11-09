#include <stdio.h>
#include <stdlib.h>
#include "gtf.h"
#include "utils.h"


exon_t *exon_init(int n) {
    exon_t *e = (exon_t*)_err_malloc(n * sizeof(exon_t));
    return e;
}
void exon_free(exon_t *e) { free(e); } 

trans_t *trans_init(int n) { 
    trans_t *t = (trans_t*)_err_malloc(n * sizeof(trans_t));
    t->exon_n = 0, t->exon_m = 2;
    t->exon = exon_init(2);
    return t;
}

// '+': s->e => S->E
// '-': S->E => s->e
int sort_exon(trans_t t)
{
    int i, j; exon_t e1, e2; int res=0;
    for (i = 0; i < t.exon_n-1; ++i) {
        for (j = i+1; j < t.exon_n; ++j) {
            e1 = t.exon[i], e2 = t.exon[j];
            if ((e1.is_rev == e2.is_rev && ((!e1.is_rev && e1.start > e2.start) || (e1.is_rev && e1.start < e2.start))) //same strand
             || (e1.is_rev && !e2.is_rev)) { // '-' and '+'
                exon_t tmp = e1;
                t.exon[i] = t.exon[j];
                t.exon[j] = tmp;
            }
        }
    }
    return res;
}

// for split and/or alternative-alignments:
//   construct full transcript with multi exons
//   filter with specific strategy XXX
int set_trans(trans_t **t)
{
    trans_t *_t = *t;

    sort_exon(*_t);
    int i;
    for (i = 0; i < _t->exon_n; ++i) {
        
    }
    return 0;
}

int add_exon(trans_t *t, int32_t tid, int32_t start, int32_t end, int32_t qstart, int32_t qend, uint8_t is_rev)
{
    if (t->exon_n == t->exon_m) t = exon_realloc(t);
    t->exon[t->exon_n].tid = tid;
    t->exon[t->exon_n].start = start;
    t->exon[t->exon_n].end = end;
    t->exon[t->exon_n].qstart = qstart;
    t->exon[t->exon_n].qend = qend;
    t->exon[t->exon_n].is_rev = is_rev;
    t->exon_n++;
    return 0;
}

trans_t *exon_realloc(trans_t *t) {
    t->exon_m <<= 1;
    t->exon = (exon_t*)_err_realloc(t->exon, t->exon_m * sizeof(exon_t));
    return t;
}

void trans_free(trans_t *t) { free(t->exon); free(t); }

gene_t *gene_init(void) { 
    gene_t *g = (gene_t*)_err_malloc(sizeof(gene_t));
    g->trans_n = 0, g->trans_m = 1;
    g->trans = trans_init(1);
    return g; 
}

gene_t *trans_realloc(gene_t *g) {
    g->trans_m <<= 1;
    g->trans = (trans_t*)_err_realloc(g->trans, g->trans_m * sizeof(trans_t));
    int i;
    for (i = g->trans_m >> 1; i < g->trans_m; ++i) {
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

int print_exon(exon_t e, FILE *out)
{
    return 0;
}

int print_trans(trans_t t, FILE *out)
{
    int i;
    for (i = 0; i < t.exon_n; ++i)
        fprintf(out, "\t[exon]\t%d\t%s\t%s\t%d\t%d\t.\t%c\t.\n", t.exon[i].tid+1, "TEST", "exon", t.exon[i].start, t.exon[i].end, "+-"[t.exon[i].is_rev]);
    return 0;
}

int print_gene(gene_t g, FILE *out)
{
    return 0;
}
