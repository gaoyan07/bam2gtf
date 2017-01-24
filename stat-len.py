#!/usr/bin/python
import re
import sys

if len(sys.argv) == 1:
    print "Usage: len.py NUM-OF-COMP full l1.name l2.name l3.name l4.name l5.name\n"
    sys.exit(1)

NUM=int(sys.argv[1])

if len(sys.argv) != NUM+3:
    print "Usage: len.py NUM-OF-COMP full l1.name l2.name l3.name l4.name l5.name\n"
    sys.exit(1)

fp=[]
outname=[]
outlen=[]

name=[]
leng=[]
full=open(sys.argv[2], "r")
line=full.readline()
while line:
    name.append(re.split('\t|\n', line)[0])
    leng.append(re.split('\t|\n', line)[1])
    line=full.readline()
full.close()
#print name, leng

for i in range(NUM):
    fp=open(sys.argv[3+i], "r")
    outname.append([])
    outlen.append([])
    line=fp.readline()
    while line:
        outname[i].append(re.split('\t|\n', line)[0])
        line=fp.readline()
    fp.close()

#for i in range(NUM-2, -1, -1):
#    for j in range(len(outname[i])):
#        if (not (outname[i][j] in outname[i+1])):
#            print i, outname[i][j]

# diff cater
for i in range(NUM-1, 0, -1):
    for j in range(len(outname[i])):
        if (not (outname[i][j] in outname[i-1])):
            if outname[i][j] in name:
                outlen[i].append(leng[name.index(outname[i][j])])
            else :
                print outname[i][j]
                
for i in range(len(outname[0])):
    if outname[0][i] in name:
        outlen[0].append(leng[name.index(outname[0][i])])
    else :
        print outname[0][i]

# output
for i in range(NUM):
    fp=open(sys.argv[3+i]+'.len', "w")
    for j in range(len(outlen[i])):
        fp.write(str(outlen[i][j])+'\n')
    fp.close()
