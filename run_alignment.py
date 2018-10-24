#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 11:38:05 2018

"""

import argparse
import os
from ReadFastQ import readfile
from BWT import BWT
from Seed import Seeds
from ReadFastQ import readfile
import random
import math
from dp_alignment import *



# Parsing command line arguments
parser = argparse.ArgumentParser(
    description="Run alignment on FastQ data")
### basic arguments
#parser.add_argument("--fastq", required=True, nargs='+',    #names of FastQ files are gathered to a list
#                    help="Path to input FastQ files")
#parser.add_argument("--genomeIndexDir", required=True,
#                    help="Path do the directory where genome index is/should be stored")
#parser.add_argument("--outdir", required=True,
#                    help="Path to the output directory where resulting .SAM file will be stored")
parser.add_argument("--runIndexing", action="store_true",
                    help="Run also indexing of the genome; if not specified, the genome index will be looked for in genomeIndexDir")
parser.add_argument("--genomeFasta", 
                    help="Path to the Fasta file with the genome sequence; required if --runIndexing is specified")
### optional arguments, changing the parameters of the algorithm
parser.add_argument("--allowForGaps", default=0.33, type=float,
                    help="Float value specifiying how large proportion of the aligned read may be constituted from gaps; default 0.33")
parser.add_argument("--seedLength", default=5, type=int,
                    help="Length of the seeds in bp")
parser.add_argument("--numberOfSeeds", default=3, type=int,
                    help="Number of seeds to be tried for each read")

args = parser.parse_args()
#if --runIndexing is present, also --genomeFasta must be present
if args.runIndexing and args.genomeFasta is None:
    parser.error("--runIndexing requires --genomeFasta.")

# Set the paths
#fastqs = args.fastq     #list with FastQ files to be processed
#indexDir = args.genomeIndexDir
#outputDir = args.outdir
if args.runIndexing:
    genomeFile = args.genomeFasta
allowForGaps = args.allowForGaps
seedLength=args.seedLength
numberOfSeeds=args.numberOfSeeds

# Create the output directory if it does not exist
#if not os.path.exists(outputDir):
#    os.makedirs(outputDir)
    

################################################################
if args.runIndexing:
    # run indexing
    pass

# process each file in fastqs
# output to outputDir
    

reads1=readfile("data_small/output_tiny_30xCov1.fq")
reads2=readfile("data_small/output_tiny_30xCov2.fq")

sequences1=list(reads1.keys())
sequences2=reads2.keys()

reference_genome=(readfile("data_small/genome.chr22.5K.fa",form="fasta")['22_5K'])
#reference_genome="TAGAGC"


bwt= BWT(reference_genome+"$") 
[last_col,L,F,S,valid_chars] = bwt.get_transform()
seeding=Seeds(F,L,S,valid_chars)
alignments=[]

def get_alignment(read,genome,_n_bp=5,overhang=0.5):
    '''
    n_bp: The no. of base pairs to consider for seedong
    read: The read that is required to be aligned
    genome: The reference genome
    overhang: The no. of bp to overhang before and after
    read indexing from 0 treated as str
    genome indexing from 0 treated as str
    '''
    n_bp=_n_bp
    read_len=len(read)
    gen_len=len(reference_genome)
    optimal_alignment={'score':-math.inf}
    for i in range(numberOfSeeds):
        # Generate a random seed of n_bp length 
        idx_seed=random.choice(range(read_len-n_bp+1))
        seed=read[idx_seed:idx_seed+n_bp] 
        locs=seeding.generate_seeds(seed)
        #Find the starting and end positions for the reference genome
        n_before=idx_seed
        n_after=read_len-n_before-n_bp
        print("Seed ",i)
        if(locs==None):
            return None
        for k,pos in enumerate(locs):
            if(pos==None):
                continue
            assert(seed==genome[pos:pos+n_bp])#
            genome_start=pos-n_before-int(n_before*overhang)
            genome_end=pos+n_bp+n_after+int(n_after*overhang)
            if(genome_start)<0:
                genome_start=0
            if(genome_end)>gen_len-1:
                genome_end=gen_len-1
            gen_seq=genome[genome_start:genome_end]
            output=DP_semi_global_alignment(read,gen_seq)
            if(output['score']>optimal_alignment['score']):
                optimal_alignment=output
    return optimal_alignment

for i,read in enumerate(sequences1[:100]):
    print("Processing read %d...."%i)
    out=get_alignment(read,reference_genome,_n_bp=seedLength)
    alignments.append(out)

