#ifndef _DEBWT_CLU_
#define _DEBWT_CLU_
#include <stdint.h>
#include "debwt.h"
#include "nano_clu.h"

int debwt_gen_loc_clu(uint8_t *bseq, int seq_len, debwt_t *db, bntseq_t *bns, uint8_t *pac, nano_clu_para *cp, vote_t *v);

#endif
