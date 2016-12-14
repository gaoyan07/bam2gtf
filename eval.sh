#!/bin/bash
#usage
if [ $# -ne 4 ]; then
    echo "Usage:"
    echo "       eval.sh SAM GTF SNYDER DATA_DIR"
    exit
fi

#input
SAM=$1
GTF=$2
SNYDER=$3
DATA=$4

#parameters
mapL=(0.5 0.6 0.67 0.8 0.9)
mapLS=3
mapLN=5
mapQ=(0.6 0.7 0.75 0.8 0.9)
mapQS=0
mapQN=5
bestR=(0.8 0.9 0.98)
bestRS=0
bestRN=3
fullS=1
fullN=5

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
            echo "./gtools filter -v${mapL[$i]} -q${mapQ[$j]} -s${bestR[$k]} $SAM 2> $BAM.log| samtools sort > $BAM"
            ./gtools filter -v${mapL[$i]} -q${mapQ[$j]} -s${bestR[$k]} $SAM 2> $BAM.log | samtools sort > $BAM

            #2 full-length novel
            for (( f=$fullS; f<=$fullN; ++f ))
            do
                FULL=${L}_l$f.gtf
                FULL_T=${L}_l$f.t
                echo "./gtools update-gtf -l $f $BAM $GTF > $FULL"
                ./gtools update-gtf -l $f $BAM $GTF > $FULL
                #3. compare with snyder
                echo "./comp-gtf $SNYDER $FULL 2> $FULL.un > $FULL.comp"
                ./comp-gtf $SNYDER $FULL 2> $FULL.un > $FULL.comp 
                echo "awk -F \"[ ;]\" '{print \$5}' $FULL | uniq > $FULL_T"
                awk -F "[ ;]" '{print $5}' $FULL | uniq > $FULL_T
                echo "awk '{printf("\"%s\"\n", \$1)}' $FULL.un > $FULL.un.t"
                awk '{printf("\"%s\"\n", $1)}' $FULL.un > $FULL.un.t
                echo "python ./comp.py $FULL.un.t $FULL_T >> $FULL.comp"
                python ./comp.py $FULL.un.t $FULL_T >> $FULL.comp
                rm $FULL
                rm $FULL_T
                rm $FULL.un.t
            done
            rm $BAM
        done
    done
done
