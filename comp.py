import sys

f1=open(sys.argv[1])
f2=open(sys.argv[2])
a=[]
b=[]

line=f1.readline()
while line:
    a.append(line[:-1])
    line=f1.readline()
line=f2.readline()
while line:
    b.append(line[:-1])
    line=f2.readline()

c=0
for i in range(len(a)):
    if (a[i] in b) == False:
        #print a[i]
        c=c+1

print("pycomp: %d" % c)

f1.close()
f2.close()
