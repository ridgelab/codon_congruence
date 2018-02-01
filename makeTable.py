#! /usr/bin/env python

import os
import sys
import argparse
from itertools import product

def parseArgs():
	parser = argparse.ArgumentParser(description='Make TNT file from FASTA Files.')
	parser.add_argument("-i",help="Input Fasta Files",nargs='*',action="store", dest="input", required=False)
	parser.add_argument("-id",help="Input Directory with Fasta Files",action="store", dest="inputDir", required=False)
	parser.add_argument("-o",help="Output for TNT Files",action="store",dest="output", required=True)
	parser.add_argument("-m",help="Minimum Number of Species in Each Group",action="store",dest="min", default=2, type=int, required=False)
	args = parser.parse_args()
	return args

def readInputFiles(args):
	codons = list(product('ACTG',repeat=3))
	codons =["".join(c) for c in codons]
	allChars = dict()
	allInputFiles = []
	if args.input:
		allInputFiles = args.input
	elif args.inputDir:
		allFasta = []
		path = args.inputDir
		allFasta = os.listdir(path)
		if path[-1] != '/':
			path += '/'
		allInputFiles = [path +i for i in allFasta]
	else:
		print "--i or --id is required"
		sys.exit()
	for species in allInputFiles:
		input = ""
		#try:
		if species[-3:] =='.gz':
			import gzip
			input = gzip.open(species,'r')
			species = species[0:-3]
		else:
			input = open(species,'r')
		if '/' in species:
			species = species.split('/')[-1]
		for line in input:
			if line[0] != '>':
				continue
			gene = ""
			if "gene=" in line:
				gene=line.split("gene=")[1].split("]")[0].upper()
			else:
				continue
			#Get the gene name...Note: Important that the gene name is converted to uppercase
			gene = gene.replace("/","_")
			sequence = input.next().strip().upper()
			if len(sequence)%3 != 0:
				continue
			allCodons = set([sequence[i:i+3] for i in xrange(0,len(sequence),3)])
			for code in codons:
				geneCode = gene + "_" + code
				if not geneCode in allChars:
					allChars[geneCode] = {}
				if code in allCodons:
					if not 1 in allChars[geneCode]:
						allChars[geneCode][1] = []
					allChars[geneCode][1].append(species)
				else:
					if not 0 in allChars[geneCode]:
						allChars[geneCode][0] = []
					allChars[geneCode][0].append(species)
		input.close()
		#except Exception:
		#	print "ASDF"
		#	continue
	return allChars

def writeTable(allChars,args):
	output = sys.stdout
	if args.output:
		output = open(args.output,'w')
	output.write("Character	State(0/1)	Species\n")
	for codon in allChars:
		if len(allChars[codon])!=2:
			continue
		if len(allChars[codon][0]) < args.min or len(allChars[codon][1]) < args.min:
			continue
		output.write(codon + "\t" + "0" + "\t" +",".join(allChars[codon][0]) + "\n")
		output.write(codon + "\t" + "1" + "\t" +",".join(allChars[codon][1]) + "\n")
	if args.output:
		output.close()


if __name__ =='__main__':
	args = parseArgs()
	allChars = readInputFiles(args)
	writeTable(allChars,args)


