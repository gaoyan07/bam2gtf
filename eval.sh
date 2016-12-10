#!/bin/bash
#usage

#input
SAM=$1
GTF=$2
SNYDER=$3
DATA=$4

#parameters
mapL=(0.5 0.6 0.67 0.8 0.9)
mapLS=0
mapLN=5
mapQ=(0.6 0.7 0.75 0.8 0.9)
mapQS=0
mapQN=5
bestR=(0.8 0.9 0.98)
bestRS=0
bestRN=3
full_l=(0 1 2 3 4 5)
fullS=1
fullN=5

awk '$3=="transcript"{print $12}' $SNYDER > $SNYDER_T
#1. filter and sort alignment #mapL/mapQ/bestR
for (( i=$mapLS; i<$mapLN; ++i ))
do
    for (( j=$mapQS; j<$mapQN; ++j ))
    do
        for (( k=$bestRS; k<$bestRN; ++k ))
        do
            BAM=$DATA/v${mapL[$i]}_q${mapQ[$j]}_s${bestR[$k]}.sort.bam
            L=$DATA/v${mapL[$i]}_q${mapQ[$j]}_s${bestR[$k]}

            #1. filter
            echo "./gtools filter -v${mapL[$i]} -q${mapQ[$j]} -s${bestR[$k]} $SAM | samtools sort > $BAM"
            ./gtools filter -v${mapL[$i]} -q${mapQ[$j]} -s${bestR[$k]} $SAM | samtools sort > $BAM  2> $BAM.log

            #2 full-length novel
            for (( f=$fullS; f<=$fullN; ++f))
            do
                $FULL=${L}_l$f.gtf
                $FULL_T=${L}_l$f.t
                ./gtools update-gtf -l ${full_l[$f]} $BAM $GTF > $FULL
                #3. compare with snyder
                ./comp-gtf $SNYDER $FULL 2> $FULL.un > $FULL.comp 
                awk -F "[ ;]" '$3=="transcript"{print $5}' $FULL > $FULL_T
                awk '{print("\"%s\"\n", $1)}' $FULL.un > $FULL.un.t
                python ./comp.py t $FULL_T >> $FULL.comp
            done
        done
    done
done
