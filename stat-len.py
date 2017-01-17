#!/usr/bin/python
import re
import sys


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
    name.append(re.split(' |\n', line)[0])
    leng.append(re.split(' |\n', line)[1])
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

# diff cater
for i in range(NUM-1, 0, -1):
    for j in range(len(name[i])):
        if (not (name[i][j] in name[i-1])):
            outlen[i].append(leng[name.index(name[i][j])])
for i in range(len(name[0])):
        outlen[0].append(leng[name.index(name[0][i])])

# output
for i in range(NUM):
    fp=open('l'+str(i+1)+'.len', "w")
    for j in range(len(outlen[i])):
        fp.write(str(outlen[i][j])+'\n')
    fp.close()
