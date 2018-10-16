# -*- coding: utf-8 -*-

import numpy as np

def readfastq(filename):
    with open(filename,'r') as fq:
              lines = fq.readlines()
              head = [item[:-1] for item in lines[::4]]
              read = [item[:-1] for item in lines[1::4]]
              qual = [item[:-1] for item in lines[3::4]]
    
              sequence = dict(zip(read, qual))

    return sequence;