import sys
import operator
import os
import csv
import re
from collections import OrderedDict

PARTITIONMEMBER=sys.argv[1]
BLASTRESADDFOLDER=sys.argv[2]
OUTPUTFOLDER=sys.argv[3]
OUTPUTFILE=sys.argv[4]

partition_mem= {} #key is contig name the vlaue is the parition that it is part of
f=open(PARTITIONMEMBER, "rU")
g=csv.reader(f,delimiter=',')
for row in g:#identify major >=50 contigs and minor <50 contigs partitions
	partition=row[1]
	partition_mem[partition]=partition_mem.get(partition,1)+1

dom_partition=OrderedDict((k,v) for k,v in partition_mem.items() if v>=50)
min_partition=OrderedDict((k,v) for k,v in partition_mem.items() if v<50)
dom_partition=OrderedDict(dict.fromkeys(dom_partition, 0))

o=open(OUTPUTFOLDER+"/"+OUTPUTFILE, "w")
#outh=csv.writer(o)
o.write("Partition,"+",".join(dom_partition.keys())+"\n")
for filename in os.listdir(BLASTRESADDFOLDER):
	dom_partition=OrderedDict(dict.fromkeys(dom_partition, 0))
	SAG=re.sub('\.blast.*', '', filename)
	if ".tabout" in filename:
		c=open(BLASTRESADDFOLDER+"/"+filename, "rU")
		d=csv.reader(c,delimiter='\t')
		for row in d:
			partition=row[12]
			if partition in dom_partition:
				dom_partition[partition]+=1
			else: 
				continue
		c.close()
	#write=SAG,",".join(str(x) for x in dom_partition.values())
	#outh.writerow(write)
	o.write(SAG+","+",".join(str(x) for x in dom_partition.values())+"\n")
