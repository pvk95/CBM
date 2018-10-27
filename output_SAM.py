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
        values: [f1 = bitwise flag reporting whether it's coming form the reverse
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
    
def output_SAM(outfile, header, d):
    with open(outfile, "w") as o:
        writer = csv.writer(o, delimiter='\t')
        o.write(''.join(header))
        for key, value in d.items():
            if value == None:
                writer.writerow([key, 'None'])
            else:
                writer.writerow([key, value[0], value[1], value[2], value[3], value[4]])
    o.close()
    