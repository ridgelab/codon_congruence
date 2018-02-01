#! /usr/bin/env python
import sys
import argparse

def parseArgs():
	'''
	This function parses system arguments.
	A character file and a reference phylogeny are required. 
	An output file path is not required. If it is not supplied, then the output will go to standard out.
	'''
	parser = argparse.ArgumentParser(description='Find Phylogenetic Signal of codon usage and aversion.')
	parser.add_argument("-c",help="Input Character File",action="store", dest="char", required=True)
	parser.add_argument("-r",help="Input Reference File",action="store", dest="ref", required=True)
	parser.add_argument("-o",help="Output File",action="store",dest="output", required=False)
	args = parser.parse_args()
	return args


def getRefTree(refTreePath):
	'''
	Input: Path to a reference phylogeny in Newick format
	Returns: a sorted list of all clades in the phylogeny, 
		where the list is sorted by the number of species in the clade, and the first element in the list has all species in the phylogeny.

	'''
	import pyparsing
	import regex
	import re
	refTreeFile = open(refTreePath,'r')
	clades = []
	thecontent = pyparsing.Word(pyparsing.alphanums) | ' ' | ',' | '_' | '/' | "'" | '[' | ']' 
	parens = pyparsing.nestedExpr('(',')', content= thecontent)
	line = refTreeFile.readline().upper().strip().replace(":1.0","").replace(";","") #single Newick tree
	#Retrieves all clades in a Newick file
	result = regex.search(r'''
	(?<rec> #capturing group rec
	\( #open parenthesis
	(?: #non-capturing group
	[^()]++ #anyting but parenthesis one or more times without backtracking
	| #or
	(?&rec) #recursive substitute of group rec
	 )*
	\) #close parenthesis
	 )
	''',line,flags=regex.VERBOSE)
	
	for clade in (result.captures('rec')):
		cladeSet = set(clade.replace('(','').replace(')','').split(','))
		if len(cladeSet)>1: #Single species clades are not important for this analysis because they are caught in larger clades.
			clades.append(cladeSet)
	c =0
	for x in re.findall(r"[\. \w]+",line):
		c +=1
		clades.append(set([x]))
	#Sorts the number of species in each clade. The smallest clades are first in the list.
	clades.sort(key=len,reverse=False)
	return clades

def makeCharacterDict(characterFilePath):
	'''
	Input: Path to the character matrix file.
	Returns: A dictionary where the key is the gene name, and the value is a dictionary 
		where the key is the species name and the value is the state of that species at that gene (0/1).
	'''
	geneChars = dict()
	characterFile = open(characterFilePath,'r')
	characterFile.readline()
	for line in characterFile:
		info = line.replace("_"," ").upper().strip().split("\t")
		gene = info[0]
		state = info[1] 
		species = []
		if len(info) >2:
			species = info[2].split(",")
		if not gene in geneChars:
			geneChars[gene] = dict()
		for s in species:
			geneChars[gene][s] = state
	characterFile.close()
	return geneChars

def getNumOrigins(clades, geneCharsForGene,gene):
	'''
	Input is a sorted list of reference clades, a dictionary for a gene, 
		where the key is the species name and the value is the state (0/1), and a gene name.
	Returns a formatted string for output to the output file.
	'''
	curState = dict() #clade to ancestral state
	curClade = dict() #species pointing to current clade
	numOriginEvents = 0
	numLossEvents = 0
	speciesWithGene = set(geneCharsForGene.keys())
	refNum0 = 0
	refNum1 = 0
	for species in clades[-1]:
		if not species in speciesWithGene:
			continue
		if geneCharsForGene[species] == '0':
			refNum0 +=1
		else:
			refNum1 +=1
	minRefNum = min(refNum0,refNum1)
	totalSpecies = refNum0 + refNum1
	if totalSpecies == 0 or minRefNum==0:
		return ""
	for clade in clades:
		states = []
		usedSpecies = set()
		for species in clade:
			if not species in speciesWithGene:
				continue
			if species in usedSpecies:
				curClade[species] = tuple(clade)
				continue
			usedSpecies.add(species)
			currentState = geneCharsForGene[species]
			if species in curClade:
				for usedS in curClade[species]:
					usedSpecies.add(usedS)
				currentState = curState[curClade[species]]
			if len(currentState) >1 : 
				for x in currentState:
					states.append(x)
			else:
				states.append(currentState)
			curClade[species] = tuple(clade)
		if len(set(states))==1:
			curState[tuple(clade)] = states[0]
			continue
		zero = states.count("0")
		one = states.count("1")
		if zero==one:
			curState[tuple(clade)] = tuple(states) 
		elif zero < one:
			numLossEvents += states.count("0")
			curState[tuple(clade)] = "1"
		else:
			numOriginEvents += states.count("1")
			curState[tuple(clade)] = "0"
	numUnresolved = len(curState[tuple(clades[-1])]) /2.0
	return (gene.replace(" ","_") + "\t" + str(numOriginEvents) + "\t" +str(numLossEvents) +"\t" + str(int(numUnresolved)) + "\t" + str(minRefNum) + "\t" + str(totalSpecies)+"\t" + str(float(minRefNum)/totalSpecies) + "\t" + str(float(numOriginEvents + numLossEvents + numUnresolved)/minRefNum)  + "\n")#, 1 #1 Unknown 0/1 event at the root node


def writeToFile(geneChars, clades,outputFilePath):
	'''
	Input is dictionary of character states of each gene, list of tuples of species in each clade in the reference phylogeny, and an output file path. 
		In the character states dictionary, the key is the gene name. 
		The value is a dictionary where the key is the species name and the value is the state (0/1).
	This function calls getNumOrigins and writes the number of origin events to the output file.	
	'''
	output = ""
	if outputFilePath:
		output = open(outputFilePath,'w')
		output.write("Gene_And_Codon_Name\tNum_Origin\tNum_Loss\tRoot_Loss(0/1)\tTotal_Species_In_Smaller_Group\tTotal_Species\tPercent_Species_In_Min\tTotal_Origin_And_Loss/Total_Species_In_Smaller_Group\n")
	else:
		sys.stdout.write("Gene_And_Codon_Name\tNum_Origin\tNum_Loss\tRoot_Loss(0/1)\tTotal_Species_In_Smaller_Group\tTotal_Species\tPercent_Species_In_Min\tTotal_Origin_And_Loss/Total_Species_In_Smaller_Group\n")
	for gene in geneChars:
		outString = getNumOrigins(clades, geneChars[gene],gene)
		if outputFilePath:
			output.write(outString)
		else:
			sys.stdout.write(outString)
	if outputFilePath:
		output.close()

if __name__ =='__main__':
	'''
	Main.
	'''
	args = parseArgs()
	clades = getRefTree(args.ref)
	geneChars = makeCharacterDict(args.char)
	writeToFile(geneChars,clades,args.output)	
			
