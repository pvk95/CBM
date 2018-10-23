#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 17:07:29 2018

@author: brynja
"""

"""
    Input: 
        Dictionary with
        key: qname = read id
        values: [seq = sequence
                f1 = bitwise flag reporting whether it's coming form the reverse
                    string or the forward string
                f2 = bitwise flag reporting secondary alignments - all false
                rname: reference name (same for everything in our case)
                pos: left most position of the alignment
                CIGAR = cigar string]
    
    Output:
        Sam file with the following header:
            qname seq f1 f2 rname pos CIGAR    
"""
import csv

def output_SAM(outfile, d):
    # USING NO LIBRARIES
    header = ['read id', 'seq', 'reverse', 'secondary_alignment', 'reference',
              'pos', 'cigarstring'] 
    with open(outfile, 'w') as o:
        o.write('\t'.join(header))
        o.write('\n')
        for key, value in d.items():
            o.write(key + '\t')
            for item in value:
                o.write(str(item))
                o.write('\t')
            o.write('\n')
    o.close()
    

def output_SAM2(outfile, d):
    # USING CSV LIBRARY
    header = ['read id', 'seq', 'reverse', 'secondary_alignment', 'reference',
              'pos', 'cigarstring'] 
    with open(outfile, "w") as o:
        writer = csv.writer(o, delimiter='\t')
        writer.writerow(header)
        for key, value in d.items():
            writer.writerow([key, value[0], value[1], value[2], value[3], value[4], value[5]])
            # If all values are strings we can also to '\t'.join(value) and replace       
    o.close()
    