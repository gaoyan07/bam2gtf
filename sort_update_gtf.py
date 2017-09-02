#!/usr/bin/python
import os, sys, re
from collections import OrderedDict

if len(sys.argv) != 2:
    sys.stderr.write('Usage:\n')
    sys.stderr.write('      %s in.update.gtf > out.sort.gtf\n\n' % sys.argv[0])
    sys.exit(1)

old_gtf_fn = sys.argv[1]


old_gtf_fp = open(old_gtf_fn, "r")

gene_group = OrderedDict()
pat = re.compile(r'gene_id \"(.+)\"*')
while True:
    line = old_gtf_fp.readline()
    if not line:
        break
    line = line[:-1]
    try:
        gene = pat.search(line).group(1).split("\"")[0]
    except:
        gene = ''
        continue
    if gene in gene_group:
        gene_group[gene].append(line)
    else:
        gene_group[gene] = [line]

for gene in gene_group:
    for line in gene_group[gene]:
        sys.stdout.writelines(line+'\n')
old_gtf_fp.close()
