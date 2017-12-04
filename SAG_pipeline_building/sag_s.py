##!/bin/bash

import sys
import operator
import os
import csv
import re

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

dom_partition=dict((k,v) for k,v in partition_mem.items() if v>=50)
min_partition=dict((k,v) for k,v in partition_mem.items() if v<50)

o=open(OUTPUTFOLDER+"/"+OUTPUTFILE, "w")
outh=csv.writer(o)
header=["SAG", "Total hits", "# Unique reads", "number of partitions", "Dominant Partition hits", "Dominant partition reads", "Minor Partition hits", "Minor partition reads", "Singleton hits"]
outh.writerow(header)

for filename in os.listdir(BLASTRESADDFOLDER): #create variables for every SAG
	Totalhitnum=0 #total hits
	uniquehitnum=0 #unique reads
	numpartitions={} #dict of all partitions
	viralreads={} #dict of reads that are viral
	total_partiton=0 #num of total partitions
	Dom_hits=0 #number of hits to dominant partitions (>=50 contigs)
	unique_dom_reads={} #dict of reads that match to dominant partitions
	dom_read_num=0 #number of reads that hit to dominant partitions
	min_hits=0 #number of hits to minor partitions (<50 contigs)
	unique_min_reads={} #dict of reads to minor partitions
	min_hit_num=0 #number of hits to minor partitions
	single_hits=0 #number of hits to singleton contigs (No blast mach in the database)

	if ".tabout" in filename:
		SAG=re.sub('\.blast.*', "", filename)
		d=csv.reader(open(BLASTRESADDFOLDER+"/"+filename, "rU"),delimiter='\t') #open input csv
		for row in d:
			Totalhitnum+=1
			read=row[0]
			viralreads[read]+=1
			partition=row[12]
			numpartitions[partition]+=1
			if "unk" in partition:
				single_hits+=1
			elif partition in dom_partition:
				unique_dom_reads[read]+=1
				Dom_hits+=1
			elif partition in min_partition:
				min_hits+=1
				unique_min_reads[read]+=1
		
		uniquehitnum=len(viralreads)
		total_partiton=len(numpartitions)
		dom_read_num=len(unique_dom_reads)
		min_hit_num=len(unique_min_reads)


		write=[SAG, Totalhitnum, uniquehitnum, total_partiton, Dom_hits, dom_read_num, min_hits, min_hit_num, single_hits]
		outh.writerow(write)
