# -*- coding: utf-8 -*-

import gzip

def readfastq(filename):
    if (filename.endswith(".gz")):
        fq = gzip.open(filename, 'r')
    else:
        fq = open(filename, 'r')
    
    lines = fq.readlines()
    head = [item[1:-1] for item in lines[::4]]
    read = [item[:-1] for item in lines[1::4]]
    qual = [item[:-1] for item in lines[3::4]]

    sequence = dict(zip(head,zip(read, qual)))
    fq.close()
    return sequence;

