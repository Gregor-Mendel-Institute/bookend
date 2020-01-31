import sys

sj_file = sys.argv[1]
counter=0
for line in open(sj_file):
    counter+=1
    l=line.rstrip().split('\t')
    reads = str(int(l[6])+int(l[7])/2)
    if l[3]=='1':
        print('\t'.join([l[0],str(int(l[1])-1),l[2],'SJ.'+str(counter),reads,'+']))
    elif l[3]=='2':
        print('\t'.join([l[0],str(int(l[1])-1),l[2],'SJ.'+str(counter),reads,'-']))
