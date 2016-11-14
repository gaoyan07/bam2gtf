#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "bam2gtf.h"
#include "sort_gtf.h"
#include "update_gtf.h"

static int usage(void)
{
    err_printf("\n");
	err_printf("Program: gtools\n");
    err_printf("Usage:   gtools <command> [options]\n\n");
	err_printf("Commands: \n");
	err_printf("         bam2gtf     generate transcript and exon information based on BAM/SAM file\n");
    err_printf("         sort        sort GTF file based on chromosome and coordinate\n");
    err_printf("         update      update existing GTF file with a new GTF file\n");
	err_printf("\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc < 2) return usage();

	if (strcmp(argv[1], "bam2gtf") == 0)      return bam2gtf(argc-1, argv+1);
    else if (strcmp(argv[1], "sort") == 0)    return sort_gtf(argc-1, argv+1);
	else if (strcmp(argv[1], "update") == 0)  return update_gtf(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}
    return 0;
}
