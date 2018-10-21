import numpy as np


def get_total_occurences_list(s, n, list):
    """
    IN: string s, int n=len(s), array list - empty array of length = len(alphabet)
    OUT: occurences list C, where
    C[k] = total number of occurrences of the characters < c
    s='banana$'
    {'a': 0, 'b': 3, 'n': 4}
    """
    for i in range(n):
        if s[i] != '$':
            for char in alphabet:
                if s[i] < char:
                    list[char] = list[char] + 1
    return list


def get_sorted_rotations_and_BWT(s):
    """"
    IN: string s,
    OUT: BWT of s, sorted rotations of s.
    Example input:
    ['banana$', 'anana$b', 'nana$ba', 'ana$ban', 'na$bana', 'a$banan', '$banana']
    ['anana$b', 'ana$ban', 'a$banan', 'banana$', 'nana$ba', 'na$bana', '$banana'] wiki solution
    ['$banana', 'a$banan', 'ana$ban', 'anana$b', 'banana$', 'na$bana', 'nana$ba'] what I do,
    slides too
    """
    n = len(s)
    s_sorted_rotations = sorted([s[i:n]+s[0:i] for i in range(n)],
                                key=None)  # try different key FUN for different sorting order
    s_bwt_index = s_sorted_rotations.index(s)
    s_bwt = ''.join([q[-1] for q in s_sorted_rotations])
    return (s_bwt, s_sorted_rotations)


def get_occurences_mtx(word_list, occurrences_mtx ):
    for i in word_list:
        # now construct the Occ (or C) datastructure
        for char in alphabet:
            if len(occurrences_mtx[char]) == 0:
                prev = 0
            else:
                prev = occurrences_mtx[char][-1]
            if i[-1:] == char:
                occurrences_mtx[char].append(prev + 1)
            else:
                occurrences_mtx[char].append(prev)
    return occurrences_mtx


def get_suffix_array(word_list, array):
    for i in word_list:
        # find a suffix 'xxx$' and append to suffix array
        dollar = i.find('$')
        array.append(i[0:dollar + 1])
    return array


### MAIN ###

# Read reference genome
reference = 'ACAGAATGATGGTGCTAATGCTCAAAGTGGCCCCGCGCAAGGTGGAGAGGGGAGGCCAGG'
reference = 'BANANA'
# Reverse reference genome
# reverse_reference = reference[::-1]#reverse reference


# Get alphabet. For DNA seq: ['A','C','G','T']
alphabet = set()
for char in reference:
    alphabet.add(char)

# Initialize occurence vector C and occurence matrix Occ
C, Occ = [dict() for i in range(2)]
for char in alphabet:
    C[char] = 0
    Occ[char] = list()
    # in Occ, each character has an associated list of integer values
    # (for each index along the reference)

# Append delimiter to end of reference
reference = "%s$" % reference
n = len(reference)

# C = {'A': 0, 'B': 3, 'N': 4}
get_total_occurences_list(reference, n, C)


# Initialize empty string rotations list, suffix array, bwt
# rotation_list =
ref_rotations, suffix_array = [list() for i in range(2)]

# Calculate sorted rotations list, and BWT of reference sequence
# ref_rotations=['$BANANA', 'A$BANAN', 'ANA$BAN', 'ANANA$B', 'BANANA$', 'NA$BANA', 'NANA$BA']
# ref_bwt='ANNB$AA'
ref_bwt, ref_rotations = get_sorted_rotations_and_BWT(reference)


# Calculate suffix array and Occurences matrix
# suffix_array = ['$', 'A$', 'ANA$', 'ANANA$', 'BANANA$', 'NA$', 'NANA$']
get_suffix_array(ref_rotations, suffix_array)

# get_occurences_mtx =
# {'A': [1, 1, 1, 1, 1, 2, 3, 4, 4, 4, 4, 4, 5, 6],
#  'B': [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2],
#  'N': [0, 1, 2, 2, 2, 2, 2, 2, 3, 4, 4, 4, 4, 4]}
get_occurences_mtx(ref_rotations, Occ)


# Save all files
# cannot save nicely for now
