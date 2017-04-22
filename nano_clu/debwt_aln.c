#include <stdio.h>
#include <stdlib.h>
#include "rest_index.h"
#include "rest_aln.h"
#include "debwt_aln.h"
#include "debwt.h"
#include "kmer_hash.h"
#include "bntseq.h"
#include "utils.h"

#define MEM_LEN 19
#define LOB_LEN 13
#define LOB_DIS 3

void uni_pos_print(uni_sa_t uid, debwt_t *db)
{
    uint32_t m;
    for (m = db->uni_pos_c[uid]; m < db->uni_pos_c[uid+1]; ++m) {
        stdout_printf("\t%c%d\n", "+-"[_debwt_get_strand(db->uni_pos_strand, m)], (int)db->uni_pos[m]);
    }
}

int bi_extend(uint8_t *seq1, uint8_t *seq2, int off1, int off2, int len, int *l1, int *l2)
{
    int i, len1=0, len2=0;

    for (i = off1-1; i >=0; --i) {
        if (seq1[i] != seq2[i]) break;
        ++len1;
    }
    for (i = off2+1; i < len; ++i) {
        if (seq1[i] != seq2[i]) break;
        ++len2;
    }
    *l1 = len1, *l2 = len2;
    return len1+len2+off2-off1+1;
}

/*
 * uni_mem: directly match with pac XXX
 */
/*
 * read         [----- off1 ----- off2 -----]
 * unipath [ ... ----- off1 ----- off2 ---------- ... ]
 */
int uni_mem(uint8_t *read_seq, int read_len, uint32_t read_off1, uint32_t read_off2, ref_off_t uni_off1, ref_off_t uni_off2, debwt_t *db, uni_sa_t uid, bntseq_t *bns, uint8_t *pac, int *l1, int *l2)
{
    uni_sa_t uni_pos_i = db->uni_pos_c[uid];
    ref_off_t uni_pos = db->uni_pos[uni_pos_i], who_uni_len = db->uni_len[uid];
    ref_off_t uni_start, uni_len, uni_len1, uni_len2, off1, off2; // off1/2: for both read and unipath
    ref_off_t pac_coor;
    uint8_t *_read_seq = read_seq;

    if (uni_off1 >= read_off1) {
        off1 = read_off1;
        uni_len1 = read_off1; uni_start = uni_off1-read_off1;
        who_uni_len -= (uni_off1-read_off1), uni_off2 -= (uni_off1-read_off1);
    } else {
        off1 = uni_off1; 
        uni_len1 = uni_off1; uni_start = 0;
        _read_seq = read_seq+(read_off1-uni_off1);
        read_len -= (read_off1-uni_off1), read_off2 -= (read_off1-uni_off1);
    }
    off2 = read_off2; // read_off2 == uni_off2
    if ((who_uni_len-uni_off2) >= (read_len-read_off2)) {
        uni_len2 = read_len-read_off2-1;
    } else {
        uni_len2 = who_uni_len-uni_off2-1;
    }
    uni_len = uni_len1+off2-off1+1+uni_len2;

    uint8_t is_rev = _debwt_get_strand(db->uni_pos_strand, uni_pos_i);
    if (is_rev) pac_coor = uni_pos-uni_start-uni_len;
    else pac_coor = uni_pos+uni_start-1;
    uint8_t *uni_seq = _bns_get_seq(bns->l_pac, pac, pac_coor, uni_len, is_rev);
    int mem_l = bi_extend(_read_seq, uni_seq, off1, off2, uni_len, l1, l2);
    return mem_l;
}

int push_loc(seed_loc_t *loc, uni_loc_t uloc)
{
    if (loc->n == loc->m) realloc_seed_loc(loc);

    loc->loc[loc->n].uid = uloc.uni_id;
    loc->loc[loc->n].uni_off = uloc.uni_off, loc->loc[loc->n].len1 = uloc.uni_loc_len;
    loc->loc[loc->n].read_off = uloc.read_off, loc->loc[loc->n].len2 = uloc.read_loc_len;
    ++loc->n;
    return uloc.read_off;
}

int push_lob(seed_loc_t *loc, lob_t lob)
{
    if (loc->n == loc->m) realloc_seed_loc(loc);

    int l1=lob.lob_flag, l2 = 1-lob.lob_flag;

    loc->loc[loc->n].uid = lob.lob[l1].uni_id;
    loc->loc[loc->n].uni_off = lob.lob[l1].uni_off, loc->loc[loc->n].read_off = lob.lob[l1].read_off;
    // calcu len1,len2 and next_seed_i
    ref_off_t len1 = lob.lob[l2].uni_off - lob.lob[l1].uni_off + lob.lob[l2].uni_loc_len;
    int len2 = lob.lob[l2].read_off - lob.lob[l1].read_off + lob.lob[l2].read_loc_len;
    loc->loc[loc->n].len1 = len1, loc->loc[loc->n].len2 = len2;
    ++loc->n;
    return lob.lob[l1].read_off;
}

int lob_check(uni_loc_t l1, uni_loc_t l2, debwt_t *db)
{
    //int i, j;
    //for (i = 0; i < l1.uni_n; ++i) { // XXX
        //for (j = 0; j < l2.uni_n; ++j) { //XXX
            if (l1.uni_id == l2.uni_id) { // same unipath
                if (abs(l1.read_off-(l2.read_off + l2.read_loc_len)) <= LOB_DIS) {
                    int uni_off2 = l2.uni_off, uni_off1 = l1.uni_off, uni_loc_len2 = l2.uni_loc_len;
                    if (abs((l1.read_off-(l2.read_off+l2.read_loc_len)) - (uni_off1-(uni_off2+uni_loc_len2))) <= LOB_DIS) return 1;
                }
            } else { // different unipath
                //uint32_t m, n;
                //for (m = db->uni_pos_c[l1.uni_id[i]]; m < db->uni_pos_c[l1.uni_id[i]+1]; ++m) {
                //    for (n = db->uni_pos_c[l2.uni_id[j]]; n < db->uni_pos_c[l2.uni_id[j]+1]; ++n) {
                //        if (abs((int)(db->uni_pos[m] + l1.off + l1.len) - (int)(db->uni_pos[n] + l2.off)) <= LOB_DIS) return 1;
                //    }
                //}
            }
        //}
    //}
    return 0;
}

int push_1lob(lob_t *lob, uni_loc_t uni_loc, debwt_t *db)
{
    if (lob->lob_flag == -1) { // NULL
        lob->lob[0] = uni_loc;
        lob->lob_flag = 0;
        return 0;
    } else if (lob->lob_flag == 0 || lob->lob_flag == 1) {
        // check_lob(lob[0/1] and new loc)
        if (lob_check(lob->lob[lob->lob_flag], uni_loc, db)) {
            lob->lob[1-lob->lob_flag] = uni_loc;
            lob->lob_flag = 1-lob->lob_flag;
            return 1;
        } else {
            //lob->lob[1-lob->lob_flag] = uni_loc;
            //lob->lob_flag = 1-lob->lob_flag;
            lob->lob_flag = -1;
            return -1;
        }
    } else {
        err_printf("[push_lob] Error: unknown lob flag: %d.\n", lob->lob_flag);
        exit(-1);
    }
}

void set_uni_loc(uni_loc_t *uni_loc, int read_off, int read_loc_len, uni_sa_t uid, ref_off_t uni_off, ref_off_t uni_loc_len)
{
    uni_loc->read_off = read_off;
    uni_loc->read_loc_len = read_loc_len;
    uni_loc->uni_id = uid;
    uni_loc->uni_off = uni_off;
    uni_loc->uni_loc_len = uni_loc_len;
    uni_loc->uni_n = 1;
    uni_loc->uni_m = 1;
}

int debwt_gen_loc_clu(uint8_t *bseq, int seq_len, debwt_t *db, bntseq_t *bns, uint8_t *pac, rest_aln_para *ap, seed_loc_t *loc_clu)
{
    int cur_i, old_i, old_lob_i;
    debwt_count_t i, uni_occ_thd = ap->debwt_uni_occ_thd, k, l, il;
    lob_t *lob = (lob_t*)_err_malloc(sizeof(lob_t)); lob->lob_flag = -1;
    uni_loc_t uni_loc;

    for (cur_i = seq_len - _BWT_HASH_K; cur_i >= 0; --cur_i) {
        old_i = cur_i;
        // debwt hash
        k = db->bwt_hash[get_hash_value(bseq+cur_i, _BWT_HASH_K)];
        il = db->bwt_hash_il[get_hash_value(bseq+cur_i, _BWT_HASH_K)];
        if (il == 0) continue;
        l = k + il - 1;
        // bwt backtrack
        while (il > uni_occ_thd && cur_i >= 1) {
            il = debwt_exact_match_alt(db, 1, bseq+cur_i-1, &k, &l);
            if (il == 0 || cur_i == 1) break;
            --cur_i;
        }
        if (il == 0) continue;

        // base extend (back/forward) // XXX NOT span unipaths
        int max_len = 0, m_len; ref_off_t uni_off, max_uni_off=0, max_loc_len1=0; uni_sa_t uid, max_uid=-1;
        int l1, l2, max_read_off=0, max_loc_len2=0;
        for (i = 0; i < il; ++i) { // 0 < il <= M
            uid = debwt_sa(db, k+i, &uni_off);
            // XXX
            m_len = uni_mem(bseq, seq_len, cur_i, old_i+_BWT_HASH_K-1, uni_off, uni_off+(old_i-cur_i)+_BWT_HASH_K-1, db, uid, bns, pac, &l1, &l2);
            // for one unipath, only one MEM location
            // for one seed, keep one longest MEM-unipath and (THD-1) secondary MEM-unipath(>l/2)
            // need a heap XXX
            if (m_len > max_len) {
                max_uid = uid;
                max_uni_off = uni_off-l1, max_loc_len1 = m_len;
                max_read_off = cur_i-l1, max_loc_len2 = m_len;
                max_len = m_len;
            }
        }
        // next loop
        cur_i = old_i;
        if (max_len > 0) set_uni_loc(&uni_loc, max_read_off, max_loc_len2, max_uid, max_uni_off, max_loc_len1);
        if (max_len >= MEM_LEN) { // MEM seed
            cur_i = push_loc(loc_clu, uni_loc) - _BWT_HASH_K; // push mem loc
            lob->lob_flag = -1;
            int l_i = loc_clu->n; loc_t l = loc_clu->loc[l_i-1];
            stdout_printf("MEM: id: %d, uni_off: %d, read_off: %d, len: %d\n", l.uid, l.uni_off, l.read_off, l.len1);
            uni_pos_print(l.uid, db);
        } else if (max_len >= LOB_LEN) {
            int res = push_1lob(lob, uni_loc, db);
            if (res == 0) { // lob_flag == -1
                cur_i = uni_loc.read_off - _BWT_HASH_K;
                old_lob_i = old_i;
            } else if (res == 1) { // check == 1
                cur_i = push_lob(loc_clu, *lob) - _BWT_HASH_K; // push lob loc
                old_lob_i = cur_i;
                int lob_i = loc_clu->n; loc_t l = loc_clu->loc[lob_i-1];
                stdout_printf("LOB id: %d, uni_off: %d, read_off: %d, len1: %d, len2: %d\n", l.uid, l.uni_off, l.read_off, l.len1, l.len2);
                uni_pos_print(l.uid, db);
            } else { // check == 0
                cur_i = old_lob_i;
            }
        }
    }
    
    free(lob);
    // return loc
    return 0;
}
