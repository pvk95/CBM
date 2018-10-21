Genome Indexing: BWT  
Seeding: FM-index  
DP alignment: Semi-global alignment  
Score: As in bowtie2: [http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#end-to-end-alignment-score-example](bowtie manual)  


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