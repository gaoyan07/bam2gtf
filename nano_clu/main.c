#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rest_index.h"
#include "rest_aln.h"

#define VERSION "1.0.0"

static int usage(void)	//main usage
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: rest (RNA-seq Mapping Tool)\n");
	fprintf(stderr, "Usage:   rest <command> [options]\n\n");
	fprintf(stderr, "Command: \n");
	fprintf(stderr, "         index      index reference sequence\n");
	fprintf(stderr, "         aln        align RNA-seq to reference\n");
	fprintf(stderr, "\n");
	return 1;
}

char rest_pg[1024];
int main(int argc, char *argv[])
{
    int i;
    sprintf(rest_pg, "@PG\tID:rest\tPN:rest\tVN:%s\tCL:%s", VERSION, argv[0]);
	for (i = 1; i < argc; ++i) sprintf(rest_pg+strlen(rest_pg), " %s", argv[i]);
 
	if (argc < 2) return usage();
	if (strcmp(argv[1], "index") == 0)      return rest_index(argc-1, argv+1);
	else if (strcmp(argv[1], "aln") == 0)   return rest_aln(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}

    free(rest_pg);
    return 0;
}
