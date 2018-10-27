import gzip

def readFasta(filename,verbose=0):
    """ Parses a classically formatted and possibly 
        compressed FASTA file into a dictionary where the key
        for a sequence is the first part of its header without 
        any white space; if verbose is nonzero then the identifiers 
        together with lengths of the read sequences are printed"""
    if (filename.endswith(".gz")):
        fp = gzip.open(filename, 'r')
    else:
        fp = open(filename, 'r')
    # Get the HD and SQ lines
    header = fp.readlines()[0:2]
    # split at headers
    data = fp.read().split('>')
    fp.close()
    # ignore whatever appears before the 1st header
    data.pop(0)     
    # prepare the dictionary
    D = {}
    for sequence in data:
        lines = sequence.split('\n')
        header = lines.pop(0).split()
        key = header[0]
        D[key] = ''.join(lines)
        if verbose:
            print("Sequence %s of length %d read" % (key,len(D[key])))
    return header, D
