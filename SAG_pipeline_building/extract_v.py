##!/bin/bash

import sys
import os
from Bio import SeqIO
import csv

BLASTRESADDFOLDER=sys.argv[1]
FASTAFOLDER=sys.argv[2]
OUTPUTFOLDER=sys.argv[3]

for bf in os.listdir(BLASTRESADDFOLDER):
	if ".tabout" in bf:
		SAG=bf[:13]+".fasta"
		fname=BLASTRESADDFOLDER +"/" + bf
		f=open(fname, "rU")
		g=csv.reader(f, delimiter="\t")
		vreads={}
		
		for row in g:
			read=row[0]
			read.strip('\n')
			if read not in vreads: #generate dict of viral reads
				vreads[read]=0

		if SAG in os.listdir(FASTAFOLDER): #identify fasta file from selected sag
			outputfile_name=OUTPUTFOLDER+"/"+bf[:-6]+"vhits.fasta"
			outputfile=open(outputfile_name, "w")
			for record in SeqIO.parse(FASTAFOLDER+"/"+SAG, "fasta"):
				a=record.id
				a=a.strip("\n")
				if a in vreads: #write viral reads to new fasta file
					print>> outputfile, '>', record.id
					print>> outputfile, record.seq
		else:
			print "The original SAG file: "+SAG+" does not exist!"


		f.close
