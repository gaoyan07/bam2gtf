#!/bin/bash
if [ $# -ne 2 ]; then
    echo "Usage: $0 in.unsort.gff out.sort.gff"
    echo "       sort GTF file based on 'transcript_id'"
    echo "       only work for \"transcript\" and \"exon\" lines"
    exit
fi

# 1. cat head lines to 
#awk '{if(NR<=2) {print} else {exit}}' $1 > $2
# 2. tag lines based on 'transcript_id'
awk '
BEGIN {
OFS="\t";
id=""; start=3000000000; end=0; chr=0; chr_m=25;
chrom["chr1"]=1; chrom["chr2"]=2; chrom["chr3"]=3; chrom["chr4"]=4; chrom["chr5"]=5;
chrom["chr6"]=6; chrom["chr7"]=7; chrom["chr8"]=8; chrom["chr9"]=9; chrom["chr10"]=10;
chrom["chr11"]=11; chrom["chr12"]=12; chrom["chr13"]=13; chrom["chr14"]=14; chrom["chr15"]=15;
chrom["chr16"]=16; chrom["chr17"]=17; chrom["chr18"]=18; chrom["chr19"]=19; chrom["chr20"]=20;
chrom["chr21"]=21; chrom["chr22"]=22; chrom["chrX"]=23; chrom["chrY"]=24; chrom["chrM"]=25;
} 
{
    if($12 != id) {
        print chr, start, end, id;
        start=3000000000; end=0; id=$12; 
    }
    if (chrom[$1] == "") {
        chr_m += 1;
        chrom[$1] = chr_m;
    }
    chr = chrom[$1]
    if ($4<start) {start = $4}
    if ($5>end) {end = $5}
}
END {
print chr, start, end, id;
} ' $1 > .tmp

# 3. add transcript lines and tag exon lines with transcript start and end
awk '
BEGIN{ OFS="\t"; id=""}
NR==FNR {
    array[$4] = $0
}
NR!=FNR {
{ 
    if ($12 != id) {print array[$12], NR-1, $1, 0, "transcript", "", "", "", "", "", ""; id=$12}
    print array[$12], NR, $0
}
}' .tmp $1 | sort -n -k1 -n -k2 -n -k3 -k4 -n -k5 | awk 'BEGIN{FS="\t"; OFS="\t"} {print $6,$7,$8,$9,$10,$11,$12,$13,$14}' > $2
rm .tmp
