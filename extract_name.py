#!/usr/bin/python
import re
import sys

fp=open(sys.argv[1], "r")
line=fp.readline()
while line:
    if (line.split("\t")[2] == 'transcript'):
        print re.split("\t|\"", line)[11][:-6]

fp.close()
