#include <stdio.h>
#include <stdlib.h>
#include "utils.h"

int sort_gtf_usage(void)
{
    err_printf("\n");
    err_printf("Usage:  gtools sort <in.gtf> > out.gtf\n");
    err_printf("\n");
    return 1;
}

int sort_gtf(int argc, char *argv[])
{
    if (argc != 2) return sort_gtf_usage();
    FILE *in_fp = fopen(argv[1], "r");

    fclose(in_fp);
    return 0;
}
