#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 17:07:29 2018

@author: brynja
"""
import csv  
    
def output_SAM(outfile, header, length, d):
    with open(outfile, "w") as o:
        writer = csv.writer(o, delimiter='\t')
        o.write("@HD VN:1.6\n")
        o.write("@SQ SN:" + header + " LN:" + str(length))
        o.write('\n')
        for key, value in d.items():
            if value == None:
                writer.writerow([key, 4, '*', 0, 0, '*', '*', 0, 0, '*', '*'])
            else:
                if value[1] == 1:
                    value[1] = 16
                writer.writerow([key, value[1], header, value[4], 0, value[5], '*',0, 0, value[0], '*'])
    o.close()
