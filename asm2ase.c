#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "utils.h"
#include "gtf.h"
#include "build_sg.h"


int asm2ase(SG_group *sg_g, SGasm_group *asm_g, ASE_t *ase)
{
    int i, j;
    for (i = 0; i < asm_g->sg_asm_n; ++i) {
        SGasm *sg_asm = asm_g->sg_asm[i];
        int sg_i = sg_asm->SG_id; SG *sg = sg_g->SG[sg_i]; 
        SGnode *node = sg->node; SGsite *acc_site = sg->acc_site; SGsite *don_site = sg->don_site; SGedge *edge = sg->edge;
        int start, end; uint32_t v_s = sg_asm->v_start, v_e = sg_asm->v_end;

        // 

    }
    return ase->se_n+ase->a5ss_n+ase->a3ss_n+ase->mxe_n+ase->ri_n;
}
