Genome Indexing: BWT  
Seeding: FM-index  
DP alignment: Semi-global alignment  
Score: As in bowtie2: [bowtie manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#end-to-end-alignment-score-example)  


### Done:
- read FASTA
- DP semi-global alignment
    - Computes semi-global alignment of the read against a piece of genome. 
    - Leading and trailing gaps are not penalized. 
    - Scoring according to bowtie2: 0 for match, -6 for mismatch, -5 gap opening, -3 gap extension. 
    - Returns dictionary with three items:
	    - score .. score of the alignment (highest possible: 0, the more negative, the worse)
	    - pos .. position of the left end of the read in given piece of genome, indexing starts from 1
	    - cigar .. CIGAR string representing the alignment
- read FASTQ  
    -read the sequence and filter all 'N' in the sequence.
- BWT index
    - Generate all rotations of genome strings.
    - Sort the matrix lexicographically and only keep the last column L
    - Compute the BWT sequence matrix and obtain the last column 'L' of it.
    - Compute the C array. (C[k] = total number of occurrences of the characters < c)
    - Compute the Occurrence matrix. (Occ(c, k) = number of times c occurs in L[1, k])
    - Ouput three files:
        - TXT file for the last column of Burrows-Wheeler Transformation.
        - CSV file for the C array.
        - CSV file for the Occurrence matrix.