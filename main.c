#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "bam_filter.h"
#include "update_gtf.h"
#include "bam2gtf.h"
#include "bam2sj.h"
#include "build_sg.h"
#include "pred_sg.h"
#include "pred_asm.h"
#include "asm2ase.h"

const char PROG[20] = "gtools";

static int usage(void)
{
    err_printf("\n");
	err_printf("Program: %s\n", PROG);
    err_printf("Usage:   %s <command> [options]\n\n", PROG);
	err_printf("Commands: \n");
    err_printf("         filter       filter out alignment records with low confidence\n");
	err_printf("         update-gtf   generate new GTF file based on BAM/SAM and existing GTF file\n");
	err_printf("         bam2gtf      generate transcript and exon information based on BAM/SAM file\n");
	err_printf("         bam2sj       generate splice-junction information based on BAM/SAM file\n");
    err_printf("         build-sg     construct splicing graph based on GTF file\n");
    err_printf("         predict-sg   predict splicing graph based on GTF file and short-read splice-junction\n");
    err_printf("         asm          generate alternative splice module from GTF-based splicing graph and\n");
    err_printf("                      short-read splice-junction\n");
    err_printf("         ase          generate 5 types of alternative splice events from GTF-based splicing graph\n");
    err_printf("                      and short-read splice-junction\n");
	err_printf("\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc < 2) return usage();

    if (strcmp(argv[1], "filter") == 0) return bam_filter(argc-1, argv+1);
	else if (strcmp(argv[1], "update-gtf") == 0) return update_gtf(argc-1, argv+1);
	else if (strcmp(argv[1], "bam2gtf") == 0) return bam2gtf(argc-1, argv+1);
    else if (strcmp(argv[1], "bam2sj") == 0) return bam2sj(argc-1, argv+1);
    else if (strcmp(argv[1], "build-sg") == 0) return build_sg(argc-1, argv+1);
    else if (strcmp(argv[1], "predict-sg") == 0) return pred_sg(argc-1, argv+1);
    else if (strcmp(argv[1], "asm") == 0) return pred_asm(argc-1, argv+1);
    else if (strcmp(argv[1], "ase") == 0) return pred_ase(argc-1, argv+1);
	else { fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]); return 1; }
    return 0;
}
