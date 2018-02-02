#! /usr/bin/env python
import sys
import argparse
import pyparsing
import regex
import re


def parseArgs():
	'''
	This function parses system arguments.
	A character file, a reference phylogeny, and an origin file are required. 
	An output file path is not required. If it is not supplied, then the output will go to standard out.
	'''
	parser = argparse.ArgumentParser(description='Find Phylogenetic Signal of codon usage and aversion.')
	parser.add_argument("-c",help="Input Character File",action="store", dest="char", required=True)
	parser.add_argument("-r",help="Input Reference Phylogeny File",action="store", dest="ref", required=True)
	parser.add_argument("-s",help="Input Origin File From Running scoreCharactersOnTree",action="store", dest="origin", required=True)
	parser.add_argument("-o",help="Output File",action="store",dest="output", required=False)
	args = parser.parse_args()
	return args

def getHomologousClade(refTreePath,cladeToChar):
	revisedClades = dict()
	refTreeFile = open(refTreePath,'r')
	clades = []
	reverse = set() #Where the smaller clade encompasses the root
	thecontent = pyparsing.Word(pyparsing.alphanums) | ' ' | ',' | '_' | '/' | "'" | '[' | ']' 
	parens = pyparsing.nestedExpr('(',')', content= thecontent)
	line = refTreeFile.readline().upper().strip().replace(":1.0","").replace(";","").replace(" ","_") #single Newick tree
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
	refTreeFile.close()	
	for bothClades in cladeToChar:
		clade1,clade2 = bothClades
		setClade1 = set(clade1)
		setClade2 = set(clade2)
		
		for clade in (result.captures('rec')):
			cladeSet = set(clade.replace('(','').replace(')','').split(','))
			if (len(cladeSet & setClade1) == len(setClade1)):
				if not clade1 in revisedClades:
					revisedClades[clade1] = []
				if len(setClade1) > len(setClade2):
					reverse |= set(cladeToChar[bothClades])
				revisedClades[clade1].extend(cladeToChar[bothClades])
				break
			if (len(cladeSet & setClade2) == len(setClade2)):
				if not clade2 in revisedClades:
					revisedClades[clade2] = []
				if len(setClade2) > len(setClade1):
					reverse |= set(cladeToChar[bothClades])
				revisedClades[clade2].extend(cladeToChar[bothClades])
				break


	return revisedClades,reverse

def makeNewRefTree(refTree,refTreePath,numToWrite,curClade):
	'''
	Input: Path to a reference phylogeny in Newick format
	Returns: a sorted list of all clades in the phylogeny, 
		where the list is sorted by the number of species in the clade, and the first element in the list has all species in the phylogeny.

	'''
	refTreeFile = open(refTreePath,'r')
	clades = []
	thecontent = pyparsing.Word(pyparsing.alphanums) | ' ' | ',' | '_' | '/' | "'" | '[' | ']' 
	parens = pyparsing.nestedExpr('(',')', content= thecontent)
	line = refTreeFile.readline().upper().strip().replace(":1.0","").replace(";","").replace(" ","_") #single Newick tree
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
	refTreeFile.close()	
	for clade in (result.captures('rec')):
		cladeSet = set(clade.replace('(','').replace(')','').split(','))
		if (len(cladeSet & curClade) == len(curClade)):
			lastSpecies = clade.split(",")[-1]
			lastSpecies = lastSpecies.replace(")","\)")
			lastSpecies = lastSpecies.replace("\)\)","\)[^)]*\)")
			lastSpecies = lastSpecies.replace("\)\)","\)[^)]*\)")
			match = str(re.findall(lastSpecies,refTree)[0])
			return re.sub(lastSpecies,match + numToWrite,refTree)

def readInputFiles(args,clades,output,outChar):
	charFile = open(args.char) 
	charFile.readline().upper()
	originFile = open(args.origin)
	originFile.readline()
	goodChars = set()
	cladeToChar = {}
	curChar = 1
	unorderedChars = {}
	for line in originFile:
		columns = line.strip().split("\t")
		if float(int(columns[4])) == 1:
			continue
		if int(columns[1]) +int(columns[2]) + int(columns[3]) == 1:
			goodChars.add(columns[0])
	originFile.close()
	line = charFile.readline().upper()
	while line != "":
		columns = line.strip().split("\t")
		if not columns[0] in goodChars:
			line = charFile.readline().upper()
			continue
		otherLine = charFile.readline().upper()
		otherColumn = otherLine.strip().split("\t")
		usedSpecies = set()
		s1 = set(columns[2].split(","))
		s2 = set(otherColumn[2].split(","))
		if len(s1) <len(s2):
			possibleSpecies = s1 & clades[-1] #Allows for species not included in reference phylogeny to be included in the origin file
			for c in clades:
				if len(possibleSpecies & c) == len(possibleSpecies):
					unorderedChars[curChar] = columns[0] + " " +otherColumn[1] + "->" + columns[1] + "\n"
					possibleSpecies = s2 & clades[-1] 
					for d in clades:
						if len(possibleSpecies & d) == len(possibleSpecies):
							tupleClade = tuple([tuple(c),tuple(d)])
							if not tupleClade in cladeToChar:
								cladeToChar[tupleClade] = []
							cladeToChar[tupleClade].append(curChar)
							break
					break
		else:
			possibleSpecies = s2 & clades[-1] 
			for c in clades:
				if len(possibleSpecies & c) == len(possibleSpecies):
					unorderedChars[curChar] =columns[0] + " " +columns[1] + "->" + otherColumn[1] + "\n"
					possibleSpecies = s1 & clades[-1] 
					for d in clades:
						if len(possibleSpecies & d) == len(possibleSpecies):
							tupleClade = tuple([tuple(c),tuple(d)])
							if not tupleClade in cladeToChar:
								cladeToChar[tupleClade] = []
							cladeToChar[tupleClade].append(curChar)
							break
					break
		curChar +=1

		line = charFile.readline().upper()
	num = 1
	refTreeF = open(args.ref,'r')
	refTree = refTreeF.readline().upper().replace(" ","_")
	refTreeF.close()
	
	cladeToChar2,reverse = getHomologousClade(args.ref,cladeToChar)
	for clade in sorted(cladeToChar2.keys(),key=len):
		numToWrite = []
		for x in cladeToChar2[clade]:
			numToWrite.append(num)
			num +=1
			if not x in reverse:
				outChar.write(unorderedChars[x])
			else:
				toWrite = unorderedChars[x]
				if "1->0" in toWrite:
					toWrite = re.sub("1->0","0->1",toWrite)
				else:
					toWrite = re.sub("0->1","1->0",toWrite)
				outChar.write(toWrite)
		if len(numToWrite) >1:
			numToWrite = "\"" + str(numToWrite[0]) + "-" + str(numToWrite[-1]) + "\""
		else:
			numToWrite = "\"" + str(numToWrite[0]) + "\""
		refTree = makeNewRefTree(refTree,args.ref,numToWrite,set(clade))
	#refTree = re.sub(r"([A-Z])([A-Z_]+)",r"\1\L\2",refTree)
	refTree = re.sub(r"([A-Z])([A-Z_]+)",lambda parts: parts.group(1) + parts.group(2).lower().replace("_"," "),refTree)
	output.write(str(refTree))

def getRefTree(refTreePath):
	'''
	Input: Path to a reference phylogeny in Newick format
	Returns: a sorted list of all clades in the phylogeny, 
		where the list is sorted by the number of species in the clade, and the first element in the list has all species in the phylogeny.

	'''
	refTreeFile = open(refTreePath,'r')
	clades = []
	thecontent = pyparsing.Word(pyparsing.alphanums) | ' ' | ',' | '_' | '/' | "'" | '[' | ']' 
	parens = pyparsing.nestedExpr('(',')', content= thecontent)
	line = refTreeFile.readline().upper().strip().replace(":1.0","").replace(";","").replace(" ","_") #single Newick tree
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
	refTreeFile.close()	
	for clade in (result.captures('rec')):
		cladeSet = set(clade.replace('(','').replace(')','').split(','))
		if len(cladeSet)>1: ##Single species clades are not important for this analysis because they are caught in larger clades.
			clades.append(cladeSet)
	c =0
	for x in re.findall(r"[\. \w]+",line):
		c +=1
		clades.append(set([x]))
	#Sorts the number of species in each clade. The smallest clades are first in the list.
	clades.sort(key=len,reverse=False)
	return clades

if __name__ =='__main__':
	'''
	Main.
	'''
	args = parseArgs()
	output = sys.stdout
	outChar = sys.stdout
	if args.output:
		output = open(args.output,'w')
		outChar = open(args.output + "_charactersUsed",'w')
	clades = getRefTree(args.ref)
	readInputFiles(args,clades,output,outChar)
	if args.output:
		output.close()
		outChar.close()
			


