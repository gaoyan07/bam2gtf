#!/bin/bash
if [[ $# -ne 1 ]]; then
    echo $0 data_dir
    exit
fi

dir=$1
intron=~/data/bam2gtf/short1SJ.out.tab
bam=~/data/bam2gtf/v0.67_q0.75_s0.98.sort.bam
old_gtf=~/data/bam2gtf/gencode.v15.annotation.sort.gtf
snyder=~/data/bam2gtf/snyder.sort.gtf
dis=0
for ((i=1; i<=5; ++i))
do
    echo "./gtools update-gtf -i $intron $bam $old_gtf -d $dis -l $i > $dir/l$i.gtf 2> $dir/l$i.err"
    ./gtools update-gtf -i $intron $bam $old_gtf -d $dis -l $i > $dir/l$i.gtf 2> $dir/l$i.err
    echo "./comp-gtf $snyder $dir/l$i.gtf $dis > $dir/comp$i.out 2> $dir/comp$i.un"
    ./comp-gtf $snyder $dir/l$i.gtf $dis > $dir/comp$i.out 2> $dir/comp$i.un
done
