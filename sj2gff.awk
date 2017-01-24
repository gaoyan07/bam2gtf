#! usage: awk -f sj2gtf.awk in.sj.tab > out.gff
BEGIN{
    OFS="\t"
    l=5
    name="SJ"
}
{
    if ($4 == 1) strand="+"
    else strand="-"
    id="SJ."NR""

    print $1, name, "exon", $2-5, $2-1, ".", strand, ".", id
    print $1, name, "exon", $3+1, $3+5, ".", strand, ".", id
}
