#! /usr/bin/env python
# programmer : nshah 
# usage: Translate a fasta file of Repeatmasker based on the parameters listed. This is made to translate all frames in both directions.
#
# There are two main options that can be decided by the user
# (1) Subfamily specific?
# (2) Methionine start specific?
# (3) Protein size minimum length
# (4) Should end/edge cases be counted (Do not have a stop codon within the sequence given)?

import sys
import os
import glob
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import argparse

parser = argparse.ArgumentParser(description='Translate FASTA file in all 6 frames.\nOptimized for use with Repeatmasker TE files, but any FASTA can be tried. The FASTA fieldname orders will be different.')
parser.add_argument('fasta', metavar='<FASTA FILE>', type=str, nargs=1,
                    help='The fasta file with DNA sequences of TEs that will be translated by this program. This script has been optimized to work with Repeatmasker Files.')
parser.add_argument('-m', action='store_true',
                    help='If this flag is placed, then only those peptides starting with a Methionine will be considered')
parser.add_argument('-e', action='store_true',
                    help='If this flag is placed, then peptides with no stop codon will be considered')
parser.add_argument('-l', type=int, nargs=1, metavar='<LENGTH>',
                    help='Minimum peptide length. Default: 1')


#Default arguments
methionineSpecific = False
minimumLength = 1
edgeConsider = False

args = parser.parse_args()

fastafile = args.fasta[0]

if args.m == True:
	methionineSpecific = True

if args.e == True:
	edgeConsider = True

if args.l != None:
	minimumLength = args.l[0]

def getTEProteins(row):
	description = row[0]
	dnasequence = row[1]
	if (len(dnasequence) <= (minimumLength*3)):		
		return('None')
	else:
		frames = [1,2,3,-1,-2,-3]
		pre_coding_dna = Seq(dnasequence, generic_dna)
		pre_coding_dna = pre_coding_dna.complement()
		pre_dnasequence = str(pre_coding_dna)
		for frame in frames:
			if frame > 0:
				coding_dna = Seq(dnasequence[(frame-1):], generic_dna)
				coding_protein = coding_dna.translate()
				stringadd = "F"
			if frame < 0:
				framenew = -1*frame
				coding_dna = Seq(pre_dnasequence[(framenew-1):], generic_dna)
				coding_protein = coding_dna.translate()
				stringadd = "R"
			potentialproteins = str(coding_protein).split('*')
			if edgeConsider == False:
				potentialproteins.pop()
			proteinsfound = []
			proteinstart = []
			proteinend = []
			positioncounter = abs(frame)
			for protein in potentialproteins:
				
				if methionineSpecific == True:
					if 'M' in protein:
						indexstart = protein.find('M')
						positioncounterprot = positioncounter + indexstart*3
						fullprotein = protein[indexstart:]
						if len(fullprotein) >= minimumLength:
							proteinsfound.append(fullprotein)
							proteinstart.append(positioncounterprot)
							proteinend.append(positioncounterprot + (len(fullprotein)*3) - 1)
				else:
					indexstart = 0
					positioncounterprot = positioncounter + indexstart*3
					fullprotein = protein[indexstart:]
					if len(fullprotein) >= minimumLength:
						proteinsfound.append(fullprotein)
						proteinstart.append(positioncounterprot)
						proteinend.append(positioncounterprot + (len(fullprotein)*3) - 1)
					
				positioncounter = positioncounter + ((len(protein) + 1) * 3) #The +1 accounts for the stop codons that are stripped in splitting this
			
			if proteinsfound:
				#Extract relevant fields from description
	
				rawuniqid = description.split('|')[0]
				uniqid = rawuniqid.split('>')[1]
				
				presplit = uniqid.split('rmsk_')[1]
				subfamTE = presplit.split(' ')[0]
				
				presplit = uniqid.split('range=')[1]
				location = presplit.split(' ')[0]
				
				chromosome = location.split(':')[0]
				startends = location.split(':')[1]
				startTE = startends.split('-')[0]
				endTE = startends.split('-')[1]
				
				presplit = uniqid.split('strand=')[1]
				strandTE = presplit.split(' ')[0]
				
				uniqid = "_".join([subfamTE,chromosome,startTE,endTE,strandTE])
				
				for k in range(0,len(proteinsfound)):
					print '\t'.join([uniqid, chromosome, startTE, endTE, strandTE, str(len(dnasequence)), stringadd, str(abs(frame)), str(proteinstart[k]), str(proteinend[k]), str(len(proteinsfound[k])), proteinsfound[k]])
			
	
candidateTable= [i.strip().split('\t') for i in open(fastafile).readlines()]

for candidate in candidateTable:
	getTEProteins(candidate)

	
	