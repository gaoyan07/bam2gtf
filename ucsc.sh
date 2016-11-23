#/bin/bash
if [ $# -ne 3 ]; then
    echo "Usage: $0 head.gtf out.gtf ucsc.gtf"
    exit
fi
cat $1 > $3
cat $2 >> $3
