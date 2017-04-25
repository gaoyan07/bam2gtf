#ifndef _BUILD_DE_BWT_
#define _BUILD_DE_BWT_
#include <stdint.h>

#ifdef __LIT__
#define _BWT_HASH_K 2
#else
#define _BWT_HASH_K 10
#endif

#ifdef __LARGE_MEM__
#define uni_sa_t uint64_t
#else
#define uni_sa_t uint32_t
#endif

typedef struct {
    uint8_t k, k_n;        // kmer len
    uint64_t k_m;          // 22: 0xfffffffffff

    uint8_t hash_k, hash_n;// first level's hash-kmer len 
    uint64_t hash_m;
    uint32_t hash_size;    // pow(4, hash_k)

    uint8_t uid_n, uid_ni; // for spe-kmer
    uint64_t uid_m;

    uint8_t remn_n, remn_k, remn_ni; // remain len
    uint64_t remn_m;

    uint8_t in_ni, out_ni, inout_m;  // pos of in/out flag

    uint8_t bwt_char_ni, next_char_ni, char_m; // pos of bwt_char/next_char

    uint8_t sk_ni, sk_n, sk_m; // for special-kmer, len<k

    uint8_t spe_ni, spe_m;     // spe_flag: generate #-kmer

    uint8_t uni_off_flag_ni, uni_off_m;
} hash_para;

#define kmer_int_t uint64_t
#define hash_int_t uint32_t
#define remn_int_t uint32_t
#define _KMER_INT_SIZE 64
#define kmer_num_t uint32_t  // for each hashKey
#define kmer_node_t uint32_t 
#define skmer_node_t uint32_t
#define _KMER_NODE_SIZE 32 

//#define debwt_count_t uint32_t // for forward-only
#define debwt_count_t uint64_t // for forward and reverse
#define ref_off_t uint32_t // for forward-only, < 4G
//#define ref_off_t uint64_t // for forward and reverse, > 4G


#define debwt_pac_t uint8_t

typedef struct {
    hash_para hp;               // hash_para
    
    /* hash index for normal kmer */
    debwt_count_t *kmer_c;       // [0]=0, [i] = cumu_count of kmer
    kmer_num_t    *kmer_num;     // whole count of existing kmer
    debwt_count_t kmer_tol_count;
    debwt_count_t kmer_real_count;
    kmer_node_t   *kmer_node;
                 // same to kmer_node, size: kmer_tol_count
    //XXX debwt_count_t *uni_offset_c; // (cumulative) number of offsets for each unipath
    uni_sa_t *uni_pos_c;
    kmer_int_t *uni_node;
    //XXX debwt_count_t uni_n, uni_m;
    uni_sa_t uni_n, uni_m;
    //XXX debwt_count_t next_offset_c;
    uni_sa_t next_uni_pos_c;
                 // same to kmer_num/skmer_num, size: uni_tol_count
    //debwt_count_t *uni_offset_n;
                 // smae to kmer_node, size: kmer_tol_count
    //XXX debwt_count_t *uni_id;       // UID of normal kmer
    //XXX debwt_count_t last_uid, cur_uid, next_uid;
    uni_sa_t *uni_id;
    uni_sa_t last_uid, cur_uid, next_uid;
                 // smae to skmer_node
    //ref_off_t  *uni_offset;   // 1-base
    //debwt_count_t uni_tol_count; // uni_tol = skmer_tol / (k-1)
    //debwt_count_t off_tol_count;

    /* hash index for special kmer */
    debwt_count_t *skmer_c;      
    kmer_num_t    *skmer_num;   
    debwt_count_t skmer_tol_count;
    skmer_node_t  *skmer_node;   // [UID][NODE]
                                //  1-32 33-64 
} hash_idx;


#define debwt_int_t uint64_t
#define _DEBWT_INT_SIZE  64
#define debwt_hash_t uint32_t

// default schema of bwt
// [5 OCC][4 BWT-INT]                                [5 OCC][4 BWT-INT] [...][...] [5 OCC][last n BWT-INT(n<=4)]
//        [4-4-.-4] [4--4-.-4] [4--4-.-4] [4--4-.-4]                                                ..[4][4]
// index:  0|1| 15| 16|17| 31| 32|33| 47| 48|49| 63|                                          ..|64i+n-1|64i+n|
//
// MACRO of DEBWT_PARA
#define _OCC_C 5         // 5; A/C/G/T/#
#define _BWT_N 4
#define _BWT_NT_K 4      // bits of per bwt_char 
#define _BWT_NT_B 2
#define _BWT_NT_M 0x7    // 0/1/2/3/4
#define _BWT_INV 16      // number of bwt char in one bwt_int; 16=64/4
#define _BWT_INV_B 4  
#define _BWT_INV_M 0xf   // 15

//#define _OCC_INV 64      //store C[5] for every occ_inv bwt_char; 64=16*4
//#define _OCC_INV_B 6
//#define _OCC_INV_M 0x3f  // 63

#define _OCC_INV 128     //store C[5] for every occ_inv bwt_char; 128=16*8
#define _OCC_INV_B 7
#define _OCC_INV_M 0x7f  // 127

//#define _BWT_OCC_B 3     // 9=2^3+1
//#define debwt_bwt_occ_a(i) (((i) << _BWT_OCC_B) + (i))
//#define _BWT_OCC_B       // 8+5=13 = 2^3+2^2+1
#define debwt_bwt_occ_a(i) (((i)<<3) + ((i)<<2) + (i))

#define _SA_INV 32
#define _SA_INV_B 5
#define _SA_INV_M 0x1f

#define _debwt_set_strand(uni_pos_strand, i, s) ((uni_pos_strand)[(i)>>3] |= (s) << (~(i)&7))
#define _debwt_get_strand(uni_pos_strand, i) (((uni_pos_strand)[(i)>>3]>>(~(i)&7))&1) // 0:+, 1:-

typedef struct {
    uint8_t occ_c;              
    uint8_t bwt_nt_k, bwt_nt_b, bwt_nt_m; 
    uint8_t bwt_inv, bwt_inv_b, bwt_inv_m; //  16,4,0xf (16=64/4=2^4)
    uint8_t occ_inv, occ_inv_b, occ_inv_m; //  64,6,0x3f (64=16*4=2^6)
    uint8_t bwt_occ_b;          // 3; 5+4=9=2^3+1

    uint8_t sa_inv;             // 32
    uint8_t sa_inv_b;           // 5
} debwt_para;

typedef struct {
    //debwt_para dp;
    debwt_count_t bwt_l;         // length of bwt_str(seq_len)
    debwt_count_t bwt_size;      // size of *bwt
    debwt_int_t   *bwt;          // bwt_str & OCC
    debwt_int_t   bwt_unit;      // for push_bwt
    debwt_count_t bwt_i, bwt_k;  // index of bwt_str and bwt
    uint8_t bit_table16[65536];  // number of '1'
    uint64_t cnt_table8[256];     // number of 'A/C/G/T/#'

    //XXX//debwt_count_t n_unipath;     // total count of unipath
    //XXX//debwt_count_t n_offset;      // total count of offsets
    uni_sa_t n_unipath, n_uni_pos, n_s_sa;
    debwt_count_t n_sa; 

    //XXX//debwt_count_t *sa_uid, *sa_u_off; 
    uni_sa_t *sa_uid, *sa_u_off; 
    //XXX//debwt_count_t *s_sa_uid;     // normal SA and special SA for all the '#-bwt_char'
    uni_sa_t *s_sa_uid; ref_off_t *uni_len; // length of each unipath

    // uni_offset_c and uni_offset are used in [kmer_hash]
    //XXX//debwt_count_t *uni_offset_c; // cumulative number of offsets for each unipath
    uni_sa_t *uni_pos_c;
    //ref_off_t  *uni_offset;   // offsets of each unipath
    ref_off_t *uni_pos;         // 0-base
    uint8_t *uni_pos_strand;       // 0: forward, 1: reverse

    debwt_count_t C[_OCC_C];     // cumulative count of 'A/C/G/T/#'
                                 // after update: [0]: num of #
                                 //               [1]: num of #, A
                                 //               [2]: num of #, A, C
                                 //               [3]: num of #, A, C, G
                                 //               [4]: num of #, A, C, G, T
    uint8_t bwt_hash_k;
    int bwt_hash_size;
    debwt_count_t *bwt_hash;
    debwt_count_t *bwt_hash_il;  // interval length
} debwt_t;

typedef struct {
    debwt_count_t k, l;
} debwt_intv_t;



int build_debwt(char *prefix, hash_idx *h_idx, debwt_t *de_idx, int rev);
void debwt_index_free(debwt_t *de_idx);
void push_debwt_bwt(uint8_t bwt_nt, debwt_t *db_idx);
int pac_build_debwt(const char *prefix, debwt_pac_t *db_pac, debwt_count_t l_pac, debwt_count_t f_pac, hash_idx *h_idx, debwt_t *db_idx);
void debwt_index_init0(debwt_t *db);

debwt_t *debwt_restore_index(const char *prefix);
debwt_count_t debwt_exact_match(const debwt_t *db, int qlen, const uint8_t *query, debwt_count_t *sa_k, debwt_count_t *sa_l);
debwt_count_t debwt_exact_match_alt(const debwt_t *db, int qlen, const uint8_t *query, debwt_count_t *sa_k, debwt_count_t *sa_l);
uni_sa_t debwt_sa(const debwt_t *db, debwt_count_t sa_uid_i, ref_off_t *off);

#endif
