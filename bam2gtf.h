#ifndef _BAM2GTF_H
#define _BAM2GTF_H

#define bam_unmap(b) ((b)->core.flag & BAM_FUNMAP)

int bam2gtf(int argc, char *argv[]);

#endif
