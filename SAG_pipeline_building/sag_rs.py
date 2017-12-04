##!/bin/bash

import sys
import operator
import os
import csv
import re
from collections import OrderedDict
from collections import Counter

PARTITIONMEMBER=sys.argv[1]
BLASTRESADDFOLDER=sys.argv[2]
OUTPUTFOLDER=sys.argv[3]
OUTPUTFILE=sys.argv[4]

partition_mem= {} #key is contig name the vlaue is the parition that it is part of
f=open(PARTITIONMEMBER, "rU")
g=csv.reader(f)
for row in g:
	partition=row[1]
	partition_mem[partition]=partition_mem.get(partition,1)+1

dom_partition1=OrderedDict((k,v) for k,v in partition_mem.items() if v>=50)#generate an ordered dictionary where (partition, number of reads)
dom_partition1["unk"]=0	#add unknown source

o=open(OUTPUTFOLDER+"/"+OUTPUTFILE, "w")
outh=csv.writer(o)
for filename in os.listdir(BLASTRESADDFOLDER):#this is my attempt to identify unique reads to each partition and generate a summary csv file
	reads={}
	SAG=re.sub('\.blast.*', "", filename)
	if ".tabout" in filename:
		dom_partition=OrderedDict(dict.fromkeys(dom_partition1,0))
		d=csv.reader(open(BLASTRESADDFOLDER+"/"+filename, "rU"),delimiter="\t")
		for row in d:
			read=row[0]
			partition=row[12]
			if partition in dom_partition:
				readpar=read+partition
				if readpar not in reads:#determine if read+partition combo has been seen before (I might have missed the first read?)
					reads[readpar]=partition
			else: 
				dom_partition["unk"]+=1#increase ordered dictironary value for each new read
		for k,v in reads.iteritems():
			dom_partition[v]+=1
		write=SAG, dom_partition
		outh.writerow(write)

