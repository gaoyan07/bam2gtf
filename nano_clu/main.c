#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nano_index.h"
#include "nano_clu.h"

#define VERSION "1.0.0"

static int usage(void)	//main usage
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: NanoClu (Read clustering for Nanopore RNA-seq data)\n");
	fprintf(stderr, "Usage:   nanoclu <command> [options]\n\n");
	fprintf(stderr, "Command: \n");
	fprintf(stderr, "         index      index gene region sequence\n");
	fprintf(stderr, "         clu        cluster RNA-seq based on each gene region\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc < 2) return usage();
	if (strcmp(argv[1], "index") == 0)      return nano_index(argc-1, argv+1);
	else if (strcmp(argv[1], "clu") == 0)   return nano_clu(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}
    return 0;
}
