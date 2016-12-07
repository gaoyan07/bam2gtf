#!/bin/bash
mapL=(0.5 0.6 0.8 0.9)
mapQ=(0.6 0.7 0.8 0.9)
bestR=(0.8 0.9)
full_l=(0 1 2 3 4 5)

for (( i=4; i<=5; ++i ))
do
    echo "./gtools update-gtf -l ${full_l[$i]} ~/DATA/snyder_pnas/GM12878.fa.gmap.sort.bam ~/DATA/snyder_pnas/gencode.v15.annotation.sort.gtf > full$i.gtf"
    ./gtools update-gtf -l ${full_l[$i]} ~/DATA/snyder_pnas/GM12878.fa.gmap.sort.bam ~/DATA/snyder_pnas/gencode.v15.annotation.sort.gtf > full$i.gtf
done
#1. filtered alignment
#   mapL
#   mapQ
#   bestR
#2. Full-length novel
#   full_l
#3. compare

