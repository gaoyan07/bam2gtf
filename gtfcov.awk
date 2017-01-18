#! usage: awk -f gtfcov.awk my.gtf
BEGIN{
    cov=0
    num=0
    last=""
}
{
    if ($3 == "exon" && last != $12) {
        cov += ($6-450)/50
        num += 1
        last = $12
    }
}
END{
    print cov, num, cov/num
}
