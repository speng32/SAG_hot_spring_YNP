#!/usr/bin/python

import sys
import operator
import os
import csv

BLASTFOLDER=sys.argv[1]
PARTITIONFILE=sys.argv[2]
OUTPUTFOLDER=sys.argv[3]
MINLEN=int(sys.argv[4])

#print BLASTFOLDER,PARTITIONFILE,OUTPUTFOLDER,MINLEN

contig_partition= {} #key is contig name the vlaue is the parition that it is part of
f=open(PARTITIONFILE, "rU")
g=csv.reader(f)
for row in g:
	contig=row[0]
	partition=row[1]
	if contig not in contig_partition: #create dict where key is contig name and value is partition num
		contig_partition[contig]=partition
f.close

fnumber=0

if not os.path.isdir(OUTPUTFOLDER):
	os.makedirs(OUTPUTFOLDER)
	

for filename in os.listdir(BLASTFOLDER):
	fname=BLASTFOLDER+"/" +filename
	if ".tabout" in fname:
		#print fname
		d=csv.reader(open(fname, "rB"), delimiter='\t') #open input csv
		outname=OUTPUTFOLDER+"/"+filename[:-6]+"pa."+str(MINLEN)+".tabout"
		o=open(outname, "w")
		outh=csv.writer(o, delimiter="\t")
		for row in d:
			data=row
			length=int(data[3])
			if length >= MINLEN:
				NLcontig=data[1]
				if "NL10_0808" in NLcontig:#there were some contigs in the database that shouldn't have been this ignored them
					continue
				else:
					if NLcontig in contig_partition:
						vpartition=contig_partition[NLcontig]#add partition to row
						data.append(vpartition)
					else:
						data.append("unk")#
					outh.writerow(data)
