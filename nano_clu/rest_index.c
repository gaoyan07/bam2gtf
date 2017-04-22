#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <zlib.h>
#include "rest_index.h"
#include "debwt.h"
#include "bntseq.h"
#include "utils.h"

int rest_index_usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   rest index [option] <ref.fa>\n");
    fprintf(stderr, "                    bulid index for <ref.fa>\n\n");
    fprintf(stderr, "Option:  \n");
    fprintf(stderr, "         -k [INT]     Length of kmer to construct de Bruijn graph. [Def=22]\n");
    //fprintf(stderr, "         -s [INT]     Length of first level's hashed sequence. [Def=12]\n");
    fprintf(stderr, "         -f           Construct index ONLY for the forward strand of reference genome. [Def=false]\n");
    fprintf(stderr, "\n");
    return 0;
}

#ifdef __LIT__
void hash_init_idx32_para(hash_idx *h)
{
    h->hp.k = 4;
    h->hp.k_n = 8;
    h->hp.k_m = 0xff; // 8
    h->hp.hash_k = 2;
    h->hp.hash_n = 4;
    h->hp.hash_m = 0xf; // 4
    h->hp.hash_size = pow(4, 2);
    
    //h->hp.uid_n  = 32;
    //h->hp.uid_ni = 32;
    //h->hp.uid_m = 0xffffffff;

    h->hp.remn_k = 2;
    h->hp.remn_n = 4;
    h->hp.remn_m = 0xf; // 4
    h->hp.remn_ni = 28; // 32-4=28

    //h->hp.in_ni = 27;   // in
    //h->hp.out_ni = 26;  // out
    //h->hp.inout_m = 0x1;    // 1

    h->hp.bwt_char_ni = 5; 
    h->hp.char_m = 0x7;     // 3
    h->hp.next_char_ni = 2;

    h->hp.sk_ni = 0;
    h->hp.sk_n = 5;
    h->hp.sk_m = 0x1f; // 5

    h->hp.spe_ni = 1;  
    h->hp.spe_m = 0x1;

    h->hp.uni_off_flag_ni = 0;
    h->hp.uni_off_m = 0x1;
}

void hash_reset_idx_para(hash_idx *h)
{
    if (h->hp.k != 4) {
        h->hp.k_n = h->hp.k << 1;
        h->hp.k_m = 0;
        uint8_t i=0;
        while (i < h->hp.k << 1) {
            h->hp.k_m <<= 1;
            h->hp.k_m |= 1;
            i++;
        }
        h->hp.hash_k = h->hp.k - 2;
        h->hp.hash_n = h->hp.hash_k << 1;
        h->hp.hash_m = 0;
        i=0;
        while (i < h->hp.hash_n) {
            h->hp.hash_m <<= 1;
            h->hp.hash_m |= 1;
            i++;
        }
        h->hp.hash_size = pow(4, h->hp.hash_k);
    }
}
#else
// NEW for lower memory XXX
// 
// hash-node: (uint32_t)
// ni:  8       5       2   1   0
// [1--24][25--27][28--30][31][32]
//  12-mer  bwt_c  next_c  spe uni_offset_c
//
// for spe-kmer
// hash-node: (uint32_t)
// ni:  8       5       0
// [1--24][25--27][28--32]
//  12-mer  bwt_c   k-len  
//                  k-len: for special-kmer, len<k (2^5=32)
//
//
// hash-node: (uint32_t)
// ni:16         13  12      9       6      1   0
// [1-16][17-18][19--20][21-23][24--26][27-31][32]
//  8-mer none   in/out  bwt_c  next_c  k-len  uni
//                                      k-len: for special-kmer, len<k (2^5=32)
// for spe-kmer
// hash-node: (uint64_t)
// ni:32     16         13  12      9       6      1   0
// [1-32][33-48][49-50][51--52][53-55][56--58][59-63][64]
//  uid   8-mer none   in/out  bwt_c  next_c  k-len  uni
//                                            k-len: for special-kmer, len<k (2^5=32)
void hash_init_idx32_para(hash_idx *h)
{
    h->hp.k = 22;
    h->hp.k_n = 44;
    h->hp.k_m = 0xfffffffffff; // 44
    h->hp.hash_k = 14;
    h->hp.hash_n = 28;
    h->hp.hash_m = 0xfffffff; // 28
    h->hp.hash_size = pow(4, 14);

    h->hp.remn_k = 8;
    h->hp.remn_n = 16;
    h->hp.remn_m = 0xffff; // 16
    h->hp.remn_ni = 16; // 32-16=16

    h->hp.bwt_char_ni = 5; 
    h->hp.char_m = 0x7;     // 3
    h->hp.next_char_ni = 2;

    h->hp.sk_ni = 0;
    h->hp.sk_n = 5;
    h->hp.sk_m = 0x1f; // 5
    
    h->hp.spe_ni = 1;
    h->hp.spe_m = 0x1;

    h->hp.uni_off_flag_ni = 0;
    h->hp.uni_off_m = 0x1;
}

void hash_reset_idx_para(hash_idx *h)
{
    if (h->hp.k != 22) {
        h->hp.k_n = h->hp.k << 1;
        h->hp.k_m = 0;
        uint8_t i=0;
        while (i < h->hp.k_n) {
            h->hp.k_m <<= 1;
            h->hp.k_m |= 1;
            i++;
        }
        h->hp.hash_k = h->hp.k - h->hp.remn_k;
        h->hp.hash_n = h->hp.hash_k << 1;
        h->hp.hash_m = 0;
        i=0;
        while (i < h->hp.hash_n) {
            h->hp.hash_m <<= 1;
            h->hp.hash_m |= 1;
            i++;
        }
        h->hp.hash_size = pow(4, h->hp.hash_k);
    }
}
#endif

int rest_index(int argc, char *argv[])
{
    char *prefix=0; int c, for_only=0;

    hash_idx h_idx; debwt_t db_idx;

    //hash_init_idx_para(&h_idx);
    hash_init_idx32_para(&h_idx);

    while ((c = getopt(argc, argv, "k:f")) >= 0)
    {
        switch (c)
        {
            case 'k': h_idx.hp.k = atoi(optarg); break;
            case 'f': for_only = 1; break;
            default: return rest_index_usage();
        }
    }
    if (optind + 1 > argc) return rest_index_usage();
    hash_reset_idx_para(&h_idx);
    prefix = strdup(argv[optind]);

    { // generate forward.pac
        gzFile fp = xzopen(prefix, "r");
        fprintf(stderr, "[%s] Pack genome FASTA ... ", __func__);
        bns_fasta2bntseq(fp, prefix, 1);
        fprintf(stderr, "done!\n");
        err_gzclose(fp);
    }
    { // generate de Bruijn graph and BWT
        gzFile fp = xzopen(prefix, "r");
        debwt_count_t l_pac, f_pac;
        debwt_pac_t *db_pac = debwt_gen_pac(fp, &l_pac, &f_pac, for_only); 
        err_gzclose(fp);
        pac_build_debwt(prefix, db_pac, l_pac, f_pac, &h_idx, &db_idx);
    }
    free(prefix);
    return 0;
}
