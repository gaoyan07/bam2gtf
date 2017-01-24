#!/bin/bash
OUT_DIR=./data
DIR=/u/home/y/yangao/DATA/snyder_pnas
#INPUT=2013-10
INPUT=2013-03
LOG=$OUT_DIR/${INPUT}.log
SAM=$DIR/GM12878.fa.gmap-${INPUT}.sam
FILT_BAM=$DIR/GM12878.fa.gmap-${INPUT}.filt.bam
rRNA=$DIR/rRNA.gtf
GTF=$DIR/gencode.v15.annotation.gtf
SORT_GTF=$DIR/gencode.v15.annotation.sort.gtf
INTRON=( "-i $DIR/short0SJ.out.tab" "-i $DIR/short1SJ.out.tab" "" )
SNYDER=$DIR/snyder.gtf
SNYDER_SORT=$DIR/snyder.sort.gtf
DIS=0
##0 cpu info
cat /proc/cpuinfo | grep 'model name' | uniq > $LOG
##1 filter
echo "./gtools filter $SAM $rRNA | samtools sort > $FILT_BAM" >> $LOG
{ time ./gtools filter $SAM $rRNA | samtools sort > $FILT_BAM ; } 2>> $LOG
##2 sort gtf
#echo "./sort_gtf.sh $GTF $SORT_GTF" >> $LOG
#{ time ./sort_gtf.sh $GTF $SORT_GTF ; } 2>> $LOG
##
#echo "./sort.sh $SNYDER $SNYDER_SORT" >> $LOG
#{ time ./sort.sh $SNYDER $SNYDER_SORT ; } 2>> $LOG
#3 update
for (( i=1; i<=5; ++i ))
do
    for (( j=0; j<=1; ++j ))
    do
        OUT=$OUT_DIR/${INPUT}-l${i}i${j}
        echo "./gtools update-gtf $FILT_BAM $SORT_GTF ${INTRON[$j]} -l $i > ${OUT}.gtf 2> ${OUT}.eval" >> $LOG
        { time ./gtools update-gtf $FILT_BAM $SORT_GTF ${INTRON[$j]} -l $i > ${OUT}.gtf 2> ${OUT}.eval ; } 2>> $LOG
        echo "./comp-gtf $SNYDER_SORT ${OUT}.gtf $DIS > ${OUT}.comp 2> ${OUT}.un" >> $LOG
        { time ./comp-gtf $SNYDER_SORT ${OUT}.gtf $DIS > ${OUT}.comp 2> ${OUT}.un ; } 2>> $LOG
        python ./extract_name.py ${OUT}.gtf > ${OUT}.name
    done
done
