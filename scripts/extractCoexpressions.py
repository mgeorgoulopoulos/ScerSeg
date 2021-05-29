# Copyright 2021 Michael Georgoulopoulos mgeorgoulopoulos@gmail.com

# This code is licensed under the MIT license. See LICENSE for details.

# Extract coexpression scores from zipped archive containing one gene per text file.
# Each text file is in SGDgene format and contains a list of associated genes and their coexpression score.

import os
import zipfile

zipFilename = 'Allgenes_ACS.zip'
aliasesFilename = "aliases.tsv"
coexpressionsFilename = "coexpressions.tsv"

def isInt(s):
	try: 
		int(s)
		return True
	except ValueError:
		return False

def readFile(z, filename):
	data = z.read(filename)
	text = data.decode("utf-8")
	lines = text.splitlines()
	return lines
	
def processFile(z, filename, coexpressionsFile):
	# Figure out gene name by filename
	filenameParts = filename.split(sep='.')
	if len(filenameParts) != 2 or filenameParts[1] != "SGDgene":
		print("Ignoring file: ", filename)
	geneAlias = filenameParts[0]

	lines = readFile(z, filename)
	
	
	# Do a first pass to identify the source gene in long format
	sourceGene = ""
	alias = ""
	for line in lines:
		if "<h1>Application error" in line:
			print("Ignoring file with error: ", filename)
			return
	
		tokens = line.split(sep='\t')
		if len(tokens) != 4: # Keep only lines with 4 tabs
			continue
		# Read alias
		if tokens[0] == 'Q':
			if tokens[2] != geneAlias:
				continue
			sourceGene = tokens[1]
			alias = tokens[2]
			break

	if sourceGene == "":
		print("Failed to find source gene in", filename)
		return
		
	with open(aliasesFilename, 'a') as f:
		f.write(sourceGene + '\t' + alias + '\n')
	
	# Second pass - read scores
	for line in lines:
		tokens = line.split(sep='\t')
		if len(tokens) != 4: # Keep only lines with 4 tabs
			continue
		
		# Read coexpression score
		if not isInt(tokens[0]):
			continue

		destinationGene = tokens[1]
		score = tokens[3]
		coexpressionsFile.write(sourceGene + '\t' + destinationGene + '\t' + score + '\n')	

# Clear output files - we are going to append to them later
try:
	os.remove(aliasesFilename)
	os.remove(coexpressionsFilename)
except:
	pass

# Write column names
with open(aliasesFilename, 'a') as f:
	f.write("Gene\tAlias\n")
with open(coexpressionsFilename, 'a') as f:
	f.write("Gene1\tGene2\tScore\n")

# Read all files and export coexpressions
with open(coexpressionsFilename, 'a') as coexpressionsFile:
	with zipfile.ZipFile(zipFilename) as z:
		for filename in z.namelist():
			if filename.endswith("SGDgene"):
				processFile(z, filename, coexpressionsFile)