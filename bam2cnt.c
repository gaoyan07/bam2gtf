#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "utils.h"
#include "htslib/htslib/sam.h"

/********************************************************
 * bam2cnt
 * @input: chr interval
 * @output: number of reads fully fall into the interval
 ********************************************************/
int bam2cnt_core(samFile *in, bam_hdr_t *h, bam1_t *b, int tid, int start, int end)
{
    int cnt=0;
    while (sam_read1(in, h, b) >= 0) {
         
    }
    return cnt;
}

int bam2cnt(int argc, char *argv[])
{

    return 0;
}
