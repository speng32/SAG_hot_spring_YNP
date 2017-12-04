##!/bin/bash

#Author Jacob Munson-McGee

import operator
import os
import csv
import collections
from collections import OrderedDict

partition_mem= {} #key is contig name the vlaue is the parition that it is part of
f=open("NL10_Archaea_Network_200_1e-30-Partition-Membership.csv", "rU")
g=csv.reader(f)
for row in g:#identify major >=50 contigs and minor <50 contigs partitions
	partition=row[1]
	partition_mem[partition]=partition_mem.get(partition,1)+1

dom_partition=OrderedDict((k,v) for k,v in partition_mem.items() if v>=50)
min_partition=collections.OrderedDict((k,v) for k,v in partition_mem.items() if v<50)

print len(dom_partition)
print len(min_partition)

#print dom_partition

folder_path= "ALL_SAG_2_NL10_vDNA/Contig_partition"
o=open("SAG_2_NL10_vDNA_partition_per_sag.csv", "w")
outh=csv.writer(o)
for filename in os.listdir(folder_path):
	dom_partition=OrderedDict(dict.fromkeys(dom_partition, 0))
	fname=folder_path+"/"+filename
	if "AD" in fname:
		SAG=fname[37:47]
		d=csv.reader(open(fname, "rU"))#add partition to csv file for dominant partitions only
		for row in d:
			partition=row[12]
			if partition in dom_partition:
				dom_partition[partition]=dom_partition.get(partition,0)+1
			else: 
				continue
	write=SAG, dom_partition

	outh.writerow(write)