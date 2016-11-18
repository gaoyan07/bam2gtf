#!/bin/bash
if [ $# -ne 1 ]; then
    echo "Usage: $0 in.gtf"
    echo "       check if the GTF is sorted"
    exit
fi

awk 'BEGIN{
chr="chr1";start=0;end=0; chr_m=26; unsort=0;
chrom["chr1"]=1; chrom["chr2"]=2; chrom["chr3"]=3; chrom["chr4"]=4; chrom["chr5"]=5;
chrom["chr6"]=6; chrom["chr7"]=7; chrom["chr8"]=8; chrom["chr9"]=9; chrom["chr10"]=10;
chrom["chr11"]=11; chrom["chr12"]=12; chrom["chr13"]=13; chrom["chr14"]=14; chrom["chr15"]=15;
chrom["chr16"]=16; chrom["chr17"]=17; chrom["chr18"]=18; chrom["chr19"]=19; chrom["chr20"]=20;
chrom["chr21"]=21; chrom["chr22"]=22; chrom["chrX"]=23; chrom["chrY"]=24; chrom["chrM"]=25;
}
($3 == "gene") {
if (chr==$1 && ($4 < start || $5 < end)) {
        unsort=1;
        exit;
    } else if ($1 != chr) {
        if (chrom[$1] == "") {
            chrom[$1] = chr_m;
            chr_m++;
        }
        if (chrom[$1] < chrom[chr]) {
            unsort=1;
            exit
        }
    }
    chr=$1; start=$4; end=$5
}
END{
    if (unsort==1) {
        print "Unsorted, please use \"sort_gtf.sh\" to sort the GTF file"
    } else {
        print "Sorted, no need to sort again"
    }
}' $1
