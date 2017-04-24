#!/bin/bash
# usage: gtf2reg.sh gene.gtf ref.fa output.fa

if [ $# -ne 3 ]; then
    echo "$0 gene.gtf ref.fa output.fa"
    exit
fi

gene_gtf=$1
ref_fa=$2
out_fa=$3

# generate bed file of gene region for each gene
awk 'BEGIN{OFS="\t"} ($3=="gene"){print $1,$4-1,$5}' $gene_gtf > $gene_gtf.bed

bedtools=bedtools
# use bedtools to generate sequence of each gene region
$bedtools getfasta -fi $ref_fa -bed $gene_gtf.bed -fo $out_fa
