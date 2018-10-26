#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 11:38:05 2018

"""

import argparse
import os
import json
import math
from BWT import BWT
from Seed import Seeds
from readFastQ import readfastq 
from readFasta import readFasta
from dp_alignment import get_alignment
from basicDNAfunctions import reverseComplement
from output_SAM import output_SAM2 



# Parsing command line arguments
parser = argparse.ArgumentParser(
    description="Run alignment on FastQ data")
### basic arguments
parser.add_argument("--fastq", required=True, nargs='+',    #names of FastQ files are gathered to a list
                    help="Path to input FastQ files")
parser.add_argument("--genomeIndexDir", required=True,
                    help="Path do the directory where genome index is/should be stored")
parser.add_argument("--outFile", required=True,
                    help="Path and name of the output SAM file")
parser.add_argument("--runIndexing", action="store_true",
                    help="Run also indexing of the genome; if not specified, the genome index will be looked for in genomeIndexDir")
parser.add_argument("--genomeFasta", 
                    help="Path to the Fasta file with the genome sequence; required if --runIndexing is specified")
### optional arguments, changing the parameters of the algorithm
parser.add_argument("--overhang", default=0.33, type=float,
                    help="Float value specifiying how large overhang over the length of the read will be taken into account when running DP alignment; default 0.33")
parser.add_argument("--seedLength", default=5, type=int,
                    help="Length of the seeds in bp")
parser.add_argument("--numberOfSeeds", default=3, type=int,
                    help="Number of seeds to be tried for each read")

args = parser.parse_args()
#if --runIndexing is present, also --genomeFasta must be present
if args.runIndexing and args.genomeFasta is None:
    parser.error("--runIndexing requires --genomeFasta.")

# Set the paths
fastqs = args.fastq     #list with FastQ files to be processed
indexDir = args.genomeIndexDir
outputFile = args.outFile
genomeFile = args.genomeFasta
overhang = args.overhang
seedLength=args.seedLength
numberOfSeeds=args.numberOfSeeds

# Create the index directory if it does not exist
if not os.path.exists(indexDir):
    os.makedirs(indexDir)

    

################################################################
########### indexing / reading the index
#file names
last_col_file_name = "{}/last_col.txt".format(indexDir)
L_file_name = "{}/L.json".format(indexDir)
F_file_name = "{}/F.json".format(indexDir)
S_file_name = "{}/S.json".format(indexDir)
valid_chars_file_name = "{}/valid_chars.json".format(indexDir)


####### read genome sequence
reference_genome=readFasta(genomeFile,verbose=0)
assert len(reference_genome) == 1 #case with more than one sequence in the genomeFile not yet implemented
reference_genome = next(iter(reference_genome.values()))

                       
if args.runIndexing:
    # run indexing
    print("Running indexing...")
    bwt= BWT(reference_genome+"$")
    [last_col,L,F,S,valid_chars] = bwt.get_transform()
    #saving the output to files
    with open(last_col_file_name, 'w') as outfile:
        outfile.write(last_col)

    with open(L_file_name, 'w') as outfile:
        json.dump(L, outfile)

    with open(F_file_name, 'w') as outfile:
        json.dump(F, outfile)  

    with open(S_file_name, 'w') as outfile:
        json.dump(S, outfile)

    with open(valid_chars_file_name, 'w') as outfile:
        json.dump(valid_chars, outfile)

else:
    # read the required information from file system
    print("Loading the genome index from {}...".format(indexDir))
    with open(last_col_file_name, 'r') as infile:
        last_col = infile.read()

    with open(L_file_name, 'r') as infile:
        L = json.load(infile)

    with open(F_file_name, 'r') as infile:
        F = json.load(infile)

    with open(S_file_name, 'r') as infile:
        S = json.load(infile)

    with open(valid_chars_file_name, 'r') as infile:
        valid_chars = json.load(infile)

  


############## read fastq files
reads={}
for f in fastqs:
    print("Reading file {}...".format(f))
    reads.update(readfastq(f))
    
sequences=list(reads.keys())

############## seeding and alignment
seeding=Seeds(F,L,S,valid_chars)

alignments={}

for i,readName in enumerate(sequences):
    print("Processing read %d...."%i)
    read=reads[readName][0]
    outFw=get_alignment(read,reference_genome,seedLength,numberOfSeeds,overhang,seeding)
    readRev=reverseComplement(read)
    outRev=get_alignment(readRev,reference_genome,seedLength,numberOfSeeds,overhang,seeding)
    if(outFw==None):
        out=outRev
        out['f1'] = 1
    elif(outRev==None):
        out = outFw
        out['f1'] = 0
    elif(outFw['score'] > outRev['score']):
        out = outFw
        out['f1'] = 0
    else:
        out = outRev
        out['f1'] = 1
    out['f2']=0
    out['rname']='ref'
    out['seq']=read
    alignments[readName]=out

################## output in SAM
print("Printing into SAM...")
output_SAM2(outputFile, alignments)

print("JOB SUCCESSFULLY FINISHED.")
