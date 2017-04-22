#ifndef _KMER_HASH_H_
#define _KMER_HASH_H_
#include "rest_index.h"
#include "debwt.h"

#define nt_A 0
#define nt_C 1
#define nt_G 2
#define nt_T 3
#define nt_N 4

void debwt_kmer_hash_init_num(hash_idx *h_idx);
void debwt_kmer_hash_init_node(hash_idx *h_idx);
void debwt_kmer_hash_free(hash_idx *h_idx);
void debwt_kmer_hash_free0(hash_idx *h_idx);
void debwt_kmer_hash_free1(hash_idx *h_idx);

int get_hash_nt(uint8_t *hash_nt, int hashKey, uint8_t hash_k);
int get_hash_value(uint8_t *hash_nt, uint8_t hash_k);

int pac_count_kmer(debwt_pac_t *db_pac, debwt_count_t l_pac, hash_idx *h_idx);
int pac_gen_kmer(debwt_pac_t *db_pac, debwt_count_t l_pac, hash_idx *h_idx, debwt_t *db_idx);
int pac_gen_skmer(debwt_pac_t *db_pac, debwt_count_t l_pac, debwt_count_t f_pac, hash_idx *h_idx, debwt_t *db_idx);

int kmer_merge(hash_idx *h_idx, debwt_t *de_idx);

#endif
