"""
Load FastA file, filter the 'N's in the sequence and add '$' at the end of seuqence .
Computes BWT index, obtain the last column of BWT sequence and 2 matrix. 
Returns dictionary with three items:
	L .. The last column of rotation matrix
	C_array .. C: C[k] = total number of occurrences of the characters < c
	matrix .. Occ(c, k) = number of times c occurs in L[1, k]

"""
import numpy as np
import pandas as pd
from Bio import SeqIO
from collections import Counter as count

#Function to load FastA file
def readfasta(filename):
    #Loading FastA data
    fasta_sequences = str(SeqIO.read(filename,'fasta').seq)
    #Delete N and add $
    seq = fasta_sequences.replace("N","") + "$"
    
    return seq

#Function to generate the sorted rotations of genome string
def roatstring(seq):
    #Generate all rotations
    rotate_seq = seq
    str_end = rotate_seq[-1]
    seq_len = len(seq)
    seq_list = []
    seq_list.append(rotate_seq)
    
    for i in range(len(seq)-1):
        rotate_seq = rotate_seq[-1] + rotate_seq[-seq_len:-1]
        seq_list.append(rotate_seq)
    
    #Sort the list
    seq_list_sorted = sorted(seq_list)
    
    #Keep only the last column L
    L = ""
    for sequence in seq_list_sorted:
        L = L + sequence[-1]
    
    return L   

#C[k] = total number of occurrences of the characters < c
def C(seq, L):
    #Get the condons in the sequence 
    char = sorted(count(seq))
    #Create the array with occurrences of the characters < c 
    C_array = pd.DataFrame([np.zeros(len(char))], columns = char)

    for codon in char:
        occ = 0
        for other_chr in L:
            if other_chr < codon:
                occ = occ + 1
                C_array[codon] = occ
    
    return C_array

##########Occ(c, k) = number of times c occurs in L[1, k]
#Count the number of times c occirs in L[1,k]   
def occ(c,k):
    if k == 0:
        return 0
    occ_c = 0
    for i in range(k):
        if L[i] == c:
            occ_c = occ_c + 1
    
    return occ_c

def occ_matrix(seq, L):
    char = sorted(count(seq))
    matr = pd.DataFrame(index = char, columns = range(1,len(L)+1)) 
    for c in char:
        for k in range(1,len(L)+1):
            matr[k][c] = occ(c,k)
    
    return matr

############Produce results:L column, C array and Occ matirx
seq = readfasta("genome.chr22.5K.fa")
L = roatstring(seq) 
C_array = C(seq, L)
matrix = occ_matrix(seq, L)

############Write out the results
f = open("BWT_L_column.txt", "a")
f.write(L)
f.close
C_array.to_csv("C_array.csv",index=False,sep=',')
matrix.to_csv("Occurrences_matirx.csv",index=True,sep=',')