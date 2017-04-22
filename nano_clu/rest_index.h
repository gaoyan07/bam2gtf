#ifndef _REST_INDEX_H_
#define _REST_INDEX_H_
#include <stdint.h>
#include "kstring.h"

int rest_index(int argc, char *argv[]);

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

#ifdef __LARGE_MEM__
#define uni_sa_t uint64_t
#else
#define uni_sa_t uint32_t
#endif

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

#endif
