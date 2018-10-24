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
from dp_alignment import DP_semi_global_alignment



# Parsing command line arguments
parser = argparse.ArgumentParser(
    description="Run alignment on FastQ data")
### basic arguments
parser.add_argument("--fastq", required=True, nargs='+',    #names of FastQ files are gathered to a list
                    help="Path to input FastQ files")
parser.add_argument("--genomeIndexDir", required=True,
                    help="Path do the directory where genome index is/should be stored")
parser.add_argument("--outdir", required=True,
                    help="Path to the output directory where resulting .SAM file will be stored")
parser.add_argument("--runIndexing", action="store_true",
                    help="Run also indexing of the genome; if not specified, the genome index will be looked for in genomeIndexDir")
parser.add_argument("--genomeFasta", 
                    help="Path to the Fasta file with the genome sequence; required if --runIndexing is specified")
### optional arguments, changing the parameters of the algorithm
parser.add_argument("--allowForGaps", default=0.33, type=float,
                    help="Float value specifiying how large proportion of the aligned read may be constituted from gaps; default 0.33")
parser.add_argument("--seedLength", default=10, type=int,
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
outputDir = args.outdir
if args.runIndexing:
    genomeFile = args.genomeFasta
allowForGaps = args.allowForGaps
seedLength=args.seedLenght
numberOfSeeds=args.numberOfSeeds

# Create the output directory if it does not exist
if not os.path.exists(outputDir):
    os.makedirs(outputDir)
    

################################################################
if args.runIndexing:
    # run indexing

# process each file in fastqs
# output to outputDir