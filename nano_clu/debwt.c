#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <math.h>
#include "debwt.h"
#include "rest_index.h"
#include "kmer_hash.h"
#include "kstring.h"
#include "utils.h"

void debwt_index_init(debwt_t *db, hash_idx *h)
{
    // offset 
    //db->n_unipath = h->uni_tol_count;
    //db->n_uni_pos = h->off_tol_count;

    db->uni_pos_c = (uni_sa_t*)_err_malloc((db->n_unipath+1) * sizeof(uni_sa_t));
    db->uni_pos = (ref_off_t*)_err_malloc(db->n_uni_pos * sizeof(ref_off_t));
    db->uni_pos_strand = (uint8_t*)_err_calloc((db->n_uni_pos+7)>>3, sizeof(uint8_t));

    // init debwt_bwt
    db->bwt_l = h->kmer_real_count + h->hp.k * (debwt_count_t)db->n_unipath;
    debwt_count_t n_occ = (db->bwt_l + _OCC_INV - 1) / _OCC_INV; // do NOT store last C[5]
    db->bwt_size = (db->bwt_l + _BWT_INV - 1) / _BWT_INV + n_occ * _OCC_C;
    db->bwt = (debwt_int_t*)_err_calloc(db->bwt_size, sizeof(debwt_int_t));
    db->bwt_unit = 0, db->bwt_i = 0, db->bwt_k = 0;

    // init C[]
    int i; for (i = 0; i < _OCC_C; ++i) db->C[i] = 0;

    // init debwt_sa
    // XXX db->bwt_l - db->n_unipath
    db->n_sa = db->bwt_l / _SA_INV;
    db->n_s_sa = db->n_unipath;

    db->sa_uid = (uni_sa_t*)_err_calloc(db->n_sa, sizeof(uni_sa_t));
    db->sa_u_off = (uni_sa_t*)_err_calloc(db->n_sa, sizeof(uni_sa_t));
    db->s_sa_uid = (uni_sa_t*)_err_calloc(db->n_s_sa, sizeof(uni_sa_t));
    db->uni_len = (ref_off_t*)_err_calloc(db->n_unipath, sizeof(ref_off_t));
}

void debwt_index_init0(debwt_t *db)
{
    // offsets
    //db->n_unipath = h->uni_tol_count;
    //db->n_uni_pos = h->off_tol_count;

    db->uni_pos_c = (uni_sa_t*)_err_calloc((db->n_unipath+1), sizeof(uni_sa_t));
    db->uni_pos = (ref_off_t*)_err_malloc(db->n_uni_pos * sizeof(ref_off_t));
    db->uni_pos_strand = (uint8_t*)_err_calloc((db->n_uni_pos+7)>>3, sizeof(uint8_t));
}

void debwt_index_init1(debwt_t *db, hash_idx *h)
{
    // init debwt_bwt
    db->bwt_l = h->kmer_real_count + h->hp.k * (debwt_count_t)db->n_unipath;
    debwt_count_t n_occ = (db->bwt_l + _OCC_INV - 1) / _OCC_INV; // do NOT store last C[5]
    db->bwt_size = (db->bwt_l + _BWT_INV - 1) / _BWT_INV + n_occ * _OCC_C;
    db->bwt = (debwt_int_t*)_err_calloc(db->bwt_size, sizeof(debwt_int_t));
    db->bwt_unit = 0, db->bwt_i = 0, db->bwt_k = 0;

    // init C[]
    int i; for (i = 0; i < _OCC_C; ++i) db->C[i] = 0;
}

void debwt_index_init2(debwt_t *db)
{
    // init debwt_sa
    // XXX db->bwt_l - db->n_unipath
    db->n_sa = db->bwt_l / _SA_INV;
    //db->n_s_sa = h->uni_tol_count;
    db->n_s_sa = db->n_unipath;

    db->sa_uid = (uni_sa_t*)_err_calloc(db->n_sa, sizeof(uni_sa_t));
    db->sa_u_off = (uni_sa_t*)_err_calloc(db->n_sa, sizeof(uni_sa_t));
    db->s_sa_uid = (uni_sa_t*)_err_calloc(db->n_s_sa, sizeof(uni_sa_t));
    db->uni_len = (ref_off_t*)_err_calloc(db->n_unipath, sizeof(ref_off_t));
}

void debwt_index_free(debwt_t *db_idx)
{
    free(db_idx->bwt);

    free(db_idx->sa_uid); free(db_idx->sa_u_off); free(db_idx->s_sa_uid); free(db_idx->uni_len);

    free(db_idx->uni_pos_c); free(db_idx->uni_pos); free(db_idx->uni_pos_strand);

    free(db_idx->bwt_hash); free(db_idx->bwt_hash_il);
}

// argv:
// sa_uid_i: 0-base sa index
// off: 0-base offset from start of the unipath
uni_sa_t debwt_sa(const debwt_t *db, debwt_count_t sa_uid_i, ref_off_t *off);
void print_debwt(debwt_t *db_idx)
{
    debwt_count_t i; uint32_t j;
    stdout_printf("UID:\n");
    for (i = db_idx->n_unipath; i < db_idx->bwt_l; ++i) {
        stdout_printf("%lld: %lld\t", (long long)i, (long long)debwt_sa(db_idx, i, &j));
        stdout_printf("%lld\n",(long long)j);
    }
    /*
    for (i = 0; i < db_idx->bwt_l-db_idx->n_unipath; ++i) {
        fprintf(stdout, "# %d:\t%d\n", i+1, db_idx->uni_id[i]); 
    }*/
    stdout_printf("Uni_offset_c:\n");
    for (i = 0; i < db_idx->n_unipath; ++i) stdout_printf("U_c %lld:\t%lld\n", (long long)i+1, (long long)db_idx->uni_pos_c[i+1]);
    stdout_printf("Uni_offset:\n");
    for (i = 0; i < db_idx->n_unipath; ++i) {
        stdout_printf("U %lld:\t", (long long)i+1);
        for (j = db_idx->uni_pos_c[i]; j < db_idx->uni_pos_c[i+1]; ++j)
            stdout_printf("%c%lld ","+-"[_debwt_get_strand(db_idx->uni_pos_strand, j)], (long long)db_idx->uni_pos[j]);
        stdout_printf("\n");
    }

    //stdout_printf("BWT-HASH:\n");
    //for (i = 0; i < db_idx->bwt_hash_size; ++i) {
    //    stdout_printf("%lld: %lld %lld\n", (long long)i, (long long)db_idx->bwt_hash[i], (long long)db_idx->bwt_hash_il[i]);
    //}
}

// XXX (sa_i)>>_OCC_INV_B
//#define debwt_bwt_occ_a(i) (((i) << _BWT_OCC_B) + (i))
#define debwt_bwt(db, sa_i) ((db)->bwt[_OCC_C + debwt_bwt_occ_a((sa_i)>>_OCC_INV_B) + (((sa_i)&(_OCC_INV_M))>>_BWT_INV_B)])
#define debwt_bwt_nt(db, sa_i) ((debwt_bwt(db, sa_i)>>((~(sa_i)&_BWT_INV_M)<<_BWT_NT_B))&_BWT_NT_M)
#define debwt_occ_a(db, sa_i) ((db)->bwt + debwt_bwt_occ_a((sa_i)>>_OCC_INV_B))
#define __occ_cnt4(table, b) (table[(b)&0x1111] + table[(b)>>16&0x1111] + table[(b)>>32&0x1111] + table[(b)>>48&0x1111])
#define __occ_cnt2(table, b) (table[(b)&0x11111111] + table[(b)>>32&0x11111111])

// count numbre of '1'
static inline int __occ_aux(debwt_int_t b, const debwt_t *db, uint8_t c)
{
    b = ((c&4)?b:~b)>>2 & ((c&2)?b:~b)>>1 & ((c&1)?b:~b) & 0x1111111111111111ull;
    return __occ_cnt4(db->bit_table16, b);
}

// count number of 'A/C/G/T/N'
static inline debwt_count_t debwt_occ(const debwt_t *db, debwt_count_t sa_uid_i, uint8_t nt)
{
    debwt_int_t *p, *end;
    debwt_count_t occ = (p = debwt_occ_a(db, sa_uid_i))[nt];
    p += _OCC_C;

    // OCC up to sa_uid_i/_BWT_INV
    end = p + ((sa_uid_i & _OCC_INV_M) >> _BWT_INV_B);
    for (; p < end; ++p) occ += __occ_aux(*p, db, nt);

    occ += __occ_aux((*p >> ((~sa_uid_i&_BWT_INV_M)<<_BWT_NT_B))>>(1<<_BWT_NT_B), db, nt);
    if (nt == 0) occ -= ((~sa_uid_i&_BWT_INV_M)+1);
    return occ;
}

// requiring k <= l
void debwt_2occ(const debwt_t *db, debwt_count_t k, debwt_count_t l, uint8_t nt, debwt_count_t *ok, debwt_count_t *ol)
{
    if (k >> _OCC_INV_B != l >> _OCC_INV_B) {
        *ok = debwt_occ(db, k, nt);
        *ol = debwt_occ(db, l, nt);
    } else {
        debwt_count_t occ_k, occ_l, *p0, *p, *end_k, *end_l;
        occ_k = (p = debwt_occ_a(db, k))[nt];
        p += _OCC_C, p0 = p;
        // calculate occ for k
        end_k = p + ((k & _OCC_INV_M) >> _BWT_INV_B);
        for (; p < end_k; ++p) occ_k += __occ_aux(*p, db, nt);
        occ_l = occ_k;
        occ_k += __occ_aux((*p >> ((~k&_BWT_INV_M)<<_BWT_NT_B))>>(1<<_BWT_NT_B), db, nt);
        // calculate occ for l
        end_l = p0 + ((l & _OCC_INV_M) >> _BWT_INV_B);
        for (; p < end_l; ++p) occ_l += __occ_aux(*p, db, nt);
        occ_l += __occ_aux((*p >> ((~l&_BWT_INV_M)<<_BWT_NT_B))>>(1<<_BWT_NT_B), db, nt);

        if (nt == 0) occ_k -= ((~k&_BWT_INV_M)+1), occ_l -= ((~l&_BWT_INV_M)+1);
        *ok = occ_k, *ol = occ_l;
    }
}


/*
#define __occ_aux4(debwt, b)                                           \
   ((debwt)->cnt_table8[(b)&0xff] + (debwt)->cnt_table8[(b)>>8&0xff]      \
  + (debwt)->cnt_table8[(b)>>16&0xff] + (debwt)->cnt_table8[(b)>>24&0xff]  \
  + (debwt)->cnt_table8[(b)>>32&0xff] + (debwt)->cnt_table8[(b)>>40&0xff]  \
  + (debwt)->cnt_table8[(b)>>48&0xff] + (debwt)->cnt_table8[(b)>>56])
#define __occ_aux1(debwt, b, c) ((__occ_aux4(debwt, b) >> ((c)<<3)) & 0x1f )
*/

#define debwt_set_intv(d, nt, k, l) {(k) = (d)->C[(nt)], (l) = (d)->C[(nt)+1]-1;}

debwt_count_t debwt_exact_match(const debwt_t *db, int qlen, const uint8_t *query, debwt_count_t *sa_k, debwt_count_t *sa_l)
{
    if (qlen < 1) return 0;
    debwt_count_t k, l, occ_k, occ_l;

    if (query[qlen-1] >= nt_N) return 0;
    debwt_set_intv(db, query[qlen-1], k, l);

    int i;
    for (i = qlen-2; i >= 0; --i) {
        uint8_t c = query[i];
        if (c >= nt_N) return 0; // 'N' is not in bwt
        debwt_2occ(db, k, l+1, c, &occ_k, &occ_l);
        if (occ_k < occ_l) {
            k = db->C[c] + occ_k;
            l = db->C[c] + occ_l - 1;
        } else return 0;
    }
    *sa_k = k, *sa_l = l;
    return l-k+1;
}

debwt_count_t debwt_exact_match_alt(const debwt_t *db, int qlen, const uint8_t *query, debwt_count_t *sa_k, debwt_count_t *sa_l)
{
    debwt_count_t k=*sa_k, l=*sa_l, occ_k, occ_l;
    int i;
    for (i = qlen-1; i >= 0; --i) {
        uint8_t c = query[i];
        if (c >= nt_N) return 0; // 'N' is not in bwt
        debwt_2occ(db, k, l+1, c, &occ_k, &occ_l);
        if (occ_k < occ_l) {
            k = db->C[c] + occ_k;
            l = db->C[c] + occ_l - 1;
        } else return 0;
    }
    *sa_k = k, *sa_l = l;
    return l-k+1;
}

static inline debwt_int_t debwt_pre_sa_i(const debwt_t *db, debwt_count_t sa_uid_i, uint8_t bwt_nt)
{
    if (bwt_nt == nt_N) return sa_uid_i;
    else return db->C[bwt_nt] + debwt_occ(db, sa_uid_i, bwt_nt);
}

void debwt_gen_bit_table16(debwt_t *db_idx)
{
    int i;
    db_idx->bit_table16[0] = 0;
    for (i = 0; i != 65536; ++i) 
        db_idx->bit_table16[i] = (i&1) + db_idx->bit_table16[i/2];
}

// ACGTN -> 0000|0001|0010|0011|0100 
// cnt_table8[1] = 00000001|00000001|00000001|00000001|00000001
//                        N        T        G        C        A
void debwt_gen_cnt_table8(debwt_t *db_idx)
{
    int i, j;
    for (i = 0; i != 256; ++i) {
        uint64_t x = 0;
        for (j = 0; j != _OCC_C; ++j)
            x |= (((i&7) == j) + ((i>>4&7)==j)) << (j<<3);
        db_idx->cnt_table8[i] = x;
    }
}

void debwt_cal_sa(debwt_t *db_idx)
{
    err_printf("[%s] Calculating SA ...\n", __func__);
    debwt_count_t sa_uid_i; uni_sa_t sa_uid;
    { //update C[5]
        int i, j;
        for (i = 4; i >= 0; --i) {
            db_idx->C[i] = 0;
            for (j = i-1; j >= 0; --j)
                db_idx->C[i] += db_idx->C[j];
            db_idx->C[i] += db_idx->n_unipath;
        }
    }
    // calculate cnt_table
    debwt_gen_bit_table16(db_idx);
    debwt_gen_cnt_table8(db_idx);
    // calculate uni_id
    debwt_count_t i, j;
    uint8_t bwt_nt=0; debwt_count_t off, off_n, off_m, *off_a;
    off_m = 100; off_a = (debwt_count_t*)_err_malloc(off_m * sizeof(debwt_count_t));
    for (i = 0; i < db_idx->n_unipath; ++i) {
        sa_uid_i = i/*0-base*/, sa_uid = i/*0-base*/;
        off = 0; off_n = 0;
        while ((bwt_nt = debwt_bwt_nt(db_idx, sa_uid_i)) < nt_N) {
            if ((sa_uid_i+1) % _SA_INV == 0) {
                db_idx->sa_uid[((sa_uid_i+1)>>_SA_INV_B)-1] = sa_uid;
                db_idx->sa_u_off[((sa_uid_i+1)>>_SA_INV_B)-1] = off;
                if (off_n >= off_m) {
                    off_m <<= 1;
                    off_a = (debwt_count_t*)_err_realloc(off_a, off_m * sizeof(debwt_count_t));
                }
                off_a[off_n++] = ((sa_uid_i+1)>>_SA_INV_B)-1;
            }
            sa_uid_i = debwt_pre_sa_i(db_idx, sa_uid_i, bwt_nt);
            off++;
#ifdef __DEBUG__
            stdout_printf("sa_i: %d\t%c\n", sa_uid_i, "ACGTN"[bwt_nt]);
#endif
        }
        // bwt_nt == nt_N, sa_uid_i
        db_idx->s_sa_uid[debwt_occ(db_idx, sa_uid_i, nt_N)] = sa_uid;
        db_idx->uni_len[sa_uid] = off;
        // reset all the sa_u_off: p -> P-p
        for (j = 0; j < off_n; ++j) 
            db_idx->sa_u_off[off_a[j]] = off-db_idx->sa_u_off[off_a[j]];
    }
    free(off_a);
    err_printf("[%s] Calculating SA done!\n", __func__);
}

uni_sa_t debwt_sa(const debwt_t *db, debwt_count_t sa_uid_i, ref_off_t *off)
{
    *off = 0;
    uint8_t nt = debwt_bwt_nt(db, sa_uid_i);
    while (nt != nt_N && ((sa_uid_i+1) & _SA_INV_M)) {
        (*off)++;
        sa_uid_i = debwt_pre_sa_i(db, sa_uid_i, nt);

        nt = debwt_bwt_nt(db, sa_uid_i);
    }
    if (nt == nt_N) {
        return db->s_sa_uid[debwt_occ(db, sa_uid_i, nt_N)];
    } else { 
        *off += db->sa_u_off[((sa_uid_i+1)>>_SA_INV_B)-1];
        return db->sa_uid[((sa_uid_i+1)>>_SA_INV_B)-1];
    }
}

void push_debwt_bwt(uint8_t bwt_nt, debwt_t *db_idx)
{
    bwt_nt = bwt_nt > nt_N ? nt_N : bwt_nt;
    //stdout_printf("BWT: %c\n", "ACGTN"[bwt_nt]);
    int i;
    if (db_idx->bwt_i % _OCC_INV == 0) {
        for (i = 0; i < _OCC_C; ++i) db_idx->bwt[db_idx->bwt_k++] = db_idx->C[i];
    }
    db_idx->bwt_unit |= (debwt_int_t)(bwt_nt) << ((~db_idx->bwt_i&_BWT_INV_M)<<_BWT_NT_B);
    db_idx->C[bwt_nt]++;
    db_idx->bwt_i++;
                                                // last bwt
    if (db_idx->bwt_i % _BWT_INV == 0 || db_idx->bwt_i == db_idx->bwt_l) {
        db_idx->bwt[db_idx->bwt_k++] = db_idx->bwt_unit;
        db_idx->bwt_unit = 0;
    }
}

// before calculate hash, bwt and sa have been calculated
void debwt_cal_hash(debwt_t *db_idx)
{
    err_printf("[%s] Calculating bwt-hash ...\n", __func__);
    db_idx->bwt_hash_k = _BWT_HASH_K; // 13
    db_idx->bwt_hash_size = (uint32_t)pow(4, _BWT_HASH_K);
    db_idx->bwt_hash = (debwt_count_t*)_err_calloc(db_idx->bwt_hash_size, sizeof(debwt_count_t));
    db_idx->bwt_hash_il = (debwt_count_t*)_err_calloc(db_idx->bwt_hash_size, sizeof(debwt_count_t));

    int i; uint8_t *hash_query = (uint8_t*)_err_malloc(_BWT_HASH_K * sizeof(uint8_t));
    debwt_count_t sa_k=0, sa_l;

    for (i = 0; i < db_idx->bwt_hash_size; ++i) {
        get_hash_nt(hash_query, i, _BWT_HASH_K);
        db_idx->bwt_hash_il[i] = debwt_exact_match(db_idx, _BWT_HASH_K, hash_query, &sa_k, &sa_l);
        db_idx->bwt_hash[i] = sa_k;
        //stdout_printf("%d:\t %d %d\n", i, db_idx->bwt_hash[i], db_idx->bwt_hash_il[i]);
    }
    free(hash_query);
    err_printf("[%s] Calculating bwt-hash done!\n", __func__);
}

void debwt_dump(const char *prefix, debwt_t *db)
{
    char *fn = (char*)_err_calloc(strlen(prefix)+15, sizeof(char));
    // dump bwt
    err_printf("[%s] Writing bwt to file ...\n", __func__);
    strcpy(fn, prefix); strcat(fn, ".rest.debwt");
    FILE *fp = xopen(fn, "wb");
    err_fwrite(&db->bwt_l, sizeof(debwt_count_t), 1, fp);
    err_fwrite(&db->bwt_size, sizeof(debwt_count_t), 1, fp);
    err_fwrite(db->C, sizeof(debwt_count_t), _OCC_C, fp);
    err_fwrite(db->bwt, sizeof(debwt_int_t), db->bwt_size, fp);
    err_fflush(fp);
    err_fclose(fp);

    // dump sa
    err_printf("[%s] Writing sa to file ...\n", __func__);
    strcpy(fn, prefix); strcat(fn, ".rest.usa");
    fp = xopen(fn, "wb");
    err_fwrite(&db->bwt_l, sizeof(debwt_count_t), 1, fp);
    err_fwrite(&db->n_sa, sizeof(debwt_count_t), 1, fp);
    err_fwrite(&db->n_s_sa, sizeof(uni_sa_t), 1, fp); // n_unipath
    err_fwrite(&db->n_uni_pos, sizeof(uni_sa_t), 1, fp);

    err_fwrite(db->sa_uid, sizeof(uni_sa_t), db->n_sa, fp);
    err_fwrite(db->sa_u_off, sizeof(uni_sa_t), db->n_sa, fp);
    err_fwrite(db->s_sa_uid, sizeof(uni_sa_t), db->n_s_sa, fp);
    err_fwrite(db->uni_len, sizeof(ref_off_t), db->n_unipath, fp);

    err_fwrite(db->uni_pos_c+1, sizeof(uni_sa_t), db->n_unipath, fp);
    err_fwrite(db->uni_pos, sizeof(ref_off_t), db->n_uni_pos, fp);
    err_fwrite(db->uni_pos_strand, sizeof(uint8_t), (db->n_uni_pos+7)>>3, fp);
    err_fflush(fp);
    err_fclose(fp);

    // dump bwt-hash
    err_printf("[%s] Writing bwt-hash to file ...\n", __func__);
    strcpy(fn, prefix); strcat(fn, ".rest.dhs");
    fp = xopen(fn, "wb");
    err_fwrite(&db->bwt_hash_size, sizeof(uint32_t), 1, fp);
    err_fwrite(db->bwt_hash, sizeof(debwt_count_t), db->bwt_hash_size, fp);
    err_fwrite(db->bwt_hash_il, sizeof(debwt_count_t), db->bwt_hash_size, fp);
    err_fflush(fp);
    err_fclose(fp);

    free(fn);
    err_printf("[%s] Writing index done.\n", __func__);
}

static void fread_fix(FILE *fp, debwt_count_t size, void *a)
{
	const int bufsize = 0x1000000; // 16M block
	debwt_count_t offset = 0;
	while (size) {
		int x = bufsize < size? bufsize : size;
		if ((x = err_fread_noeof(a + offset, 1, x, fp)) == 0) break;
		size -= x; offset += x;
	}
	//return offset;
}

debwt_t *debwt_restore_index(const char *prefix)
{
    char *fn = (char*)_err_calloc(strlen(prefix)+15, sizeof(char));
    debwt_t *db;
    // restore bwt
    strcpy(fn, prefix); strcat(fn, ".rest.debwt");
    FILE *fp = xopen(fn, "rb");

    db = (debwt_t*)_err_calloc(1, sizeof(debwt_t));
    err_fread_noeof(&db->bwt_l, sizeof(debwt_count_t), 1, fp);
    err_fread_noeof(&db->bwt_size, sizeof(debwt_count_t), 1, fp);
    err_fread_noeof(db->C, sizeof(debwt_count_t), _OCC_C, fp);
    db->bwt = (debwt_int_t*)_err_calloc(db->bwt_size, sizeof(debwt_int_t));
    fread_fix(fp, db->bwt_size * sizeof(debwt_int_t), db->bwt);
    err_fclose(fp); 

    // restore sa
    strcpy(fn, prefix); strcat(fn, ".rest.usa");
    fp = xopen(fn, "rb");
    debwt_count_t bwt_l;
    err_fread_noeof(&bwt_l, sizeof(debwt_count_t), 1, fp);
    xassert(bwt_l == db->bwt_l, "SA-BWT inconsistency: bwt length is not the same.");
    err_fread_noeof(&db->n_sa, sizeof(debwt_count_t), 1, fp);
    err_fread_noeof(&db->n_s_sa, sizeof(uni_sa_t), 1, fp);
    db->n_unipath = db->n_s_sa;
    err_fread_noeof(&db->n_uni_pos, sizeof(uni_sa_t), 1, fp);

    db->sa_uid = (uni_sa_t*)_err_calloc(db->n_sa, sizeof(uni_sa_t));
    fread_fix(fp, db->n_sa * sizeof(uni_sa_t), db->sa_uid);
    db->sa_u_off = (uni_sa_t*)_err_calloc(db->n_sa, sizeof(uni_sa_t));
    fread_fix(fp, db->n_sa * sizeof(uni_sa_t), db->sa_u_off);
    db->s_sa_uid = (uni_sa_t*)_err_calloc(db->n_s_sa, sizeof(uni_sa_t));
    fread_fix(fp, db->n_s_sa * sizeof(uni_sa_t), db->s_sa_uid);
    db->uni_len = (uni_sa_t*)_err_calloc(db->n_unipath, sizeof(uni_sa_t));
    fread_fix(fp, db->n_unipath * sizeof(ref_off_t), db->uni_len);

    db->uni_pos_c = (uni_sa_t*)_err_calloc(db->n_unipath+1, sizeof(uni_sa_t));
    fread_fix(fp, db->n_unipath * sizeof(uni_sa_t), db->uni_pos_c+1);
    db->uni_pos = (ref_off_t*)_err_calloc(db->n_uni_pos, sizeof(ref_off_t));
    fread_fix(fp, db->n_uni_pos * sizeof(ref_off_t), db->uni_pos);
    db->uni_pos_strand = (uint8_t*)_err_calloc((db->n_uni_pos+7)>>3, sizeof(uint8_t));
    fread_fix(fp, (db->n_uni_pos+7)>>3 * sizeof(uint8_t), db->uni_pos_strand);
    err_fclose(fp); 

    // restore bwt-hash
    strcpy(fn, prefix); strcat(fn, ".rest.dhs");
    fp = xopen(fn, "rb");
    err_fread_noeof(&db->bwt_hash_size, sizeof(uint32_t), 1, fp);
    db->bwt_hash = (debwt_count_t*)_err_calloc(db->bwt_hash_size, sizeof(debwt_count_t));
    db->bwt_hash_il = (debwt_count_t*)_err_calloc(db->bwt_hash_size, sizeof(debwt_count_t));
    fread_fix(fp, db->bwt_hash_size * sizeof(debwt_count_t), db->bwt_hash);
    fread_fix(fp, db->bwt_hash_size * sizeof(debwt_count_t), db->bwt_hash_il);
    err_fclose(fp);

    // count table16
    debwt_gen_bit_table16(db);
    debwt_gen_cnt_table8(db);

    free(fn);
    return db;
}

int pac_build_debwt(const char *prefix, debwt_pac_t *db_pac, debwt_count_t l_pac, debwt_count_t f_pac, hash_idx *h_idx, debwt_t *db_idx)
{
    err_printf("[pac_build_debwt] Building de Bruijn-BWT index for genome ...\n");
    debwt_kmer_hash_init_num(h_idx); 

    pac_count_kmer(db_pac, l_pac, h_idx); // count tolal number of kmer
    pac_gen_kmer(db_pac, l_pac, h_idx, db_idx);   // generate normal kmer

    debwt_index_init1(db_idx, h_idx);
    pac_gen_skmer(db_pac, l_pac, f_pac, h_idx, db_idx); // generate special kmer and partial bwt 
    free(db_pac);
    //debwt_kmer_hash_free0(h_idx);

    kmer_merge(h_idx, db_idx); // merge normal kmer and $-contained kmer, generate whole bwt
    debwt_kmer_hash_free(h_idx);

    debwt_index_init2(db_idx);
    debwt_cal_sa(db_idx);   // calculate SA
    debwt_cal_hash(db_idx); // calculate bwt_hash range

    debwt_dump(prefix, db_idx); // dump index files to disk
#ifdef __DEBUG__
    print_debwt(db_idx);
    debwt_t *d = debwt_restore_index(prefix);
    print_debwt(d);
    debwt_index_free(d); free(d);
#endif
    debwt_index_free(db_idx);
    err_printf("[pac_build_debwt] Building de Bruijn-BWT index done!\n");
    return 0;
}
