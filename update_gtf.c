#include <stdio.h>
#include <stdlib.h>
#include "utils.h"

int update_gtf_usage(void)
{
    err_printf("\n");
	err_printf("Program: update_gtf\n");
    err_printf("Usage:   update_gtf <old.gtf> <new.gtf> > out.gtf\n");
	err_printf("\n");
	return 1;
}

int update_gtf(int argc, char *argv[])
{
    if (argc != 3) return update_gtf_usage();
    FILE *old_fp, *new_fp;
    old_fp = fopen(argv[1], "r"), new_fp = fopen(argv[2], "r");
    // sort if new is unsorted
    // merge

    fclose(old_fp); fclose(new_fp);
    return 0;
}

