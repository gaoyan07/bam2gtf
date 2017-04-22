#ifndef REST_CLU_H
#define REST_CLU_H
#include <stdint.h>

typedef struct {
    uint64_t ref_off;
    int direction;
    uint64_t ref_start;
    int seed_num;
} s_clu;

typedef struct {
    int n, m, max;
    s_clu *clu;
} seed_clu_t;

#endif
