import math

'''
Computes semi-global alignment of the read against a piece of genome. Leading and trailing gaps
are not penalized. Scoring according to bowtie2: 0 for match, -6 for mismatch, -5 gap opening,
-3 gap extension. 
Returns dictionary with three items:
	score .. score of the alignment (highest possible: 0, the more negative, the worse)
	pos .. position of the left end of the read in given piece of genome, indexing starts from 1
	cigar .. CIGAR string representing the alignment
'''


### Two helper functions
# match(X,Y), returns 0 if X==Y, 1 otherwise
def match(X,Y):
	if X==Y:
		return 0
	else:
		return 1

# compressCIGAR .. gets a string and compresses repeated chars in it
#	example: MMMDD -> 3M2D
def compressCIGAR(string):
    res = ""
    count = 1
    char = string[0]
   
    #Iterate through loop, skipping last one
    for i in range(len(string)-1):
        if(string[i] == string[i+1]):
            count+=1
        else:
            res += str(count)
            res += string[i]
            count = 1
    #print last one
    if(count == 1):
        res += '1'
        res += string[len(string)-1]
    else:
    	res += str(count)
    	res += string[i]
    return res

######## The main function #########################
def DP_semi_global_alignment(read, genome):
	n = len(read)
	m = len(genome)

	### Scoring: As in Bowtie2. Match 0, mismatch -6, gap opening -5, gap extension -3.
	mismatch = -6
	gapOpen = -5
	gapExtension = -3
	
	### Initialization
	# m .. columns, n .. rows
	# M .. match matrix
	M = [[(-1*math.inf) for x in range(m+1)] for y in range(n+1)]
	M[0][0]=0
	# D .. gap in rows, i.e. deletion in read
	D = [[(-1*math.inf) for x in range(m+1)] for y in range(n+1)]
	for i in range(m+1):
		# leading gaps are free
		D[0][i] = 0
	# I .. gap in columns, i.e. insertion in read
	I = [[(-1*math.inf) for x in range(m+1)] for y in range(n+1)]
	for i in range(1,n+1):
		I[i][0] = gapOpen + (i-1)*gapExtension

	### DP
	for i in range(1,n+1):
		for j in range(1,m+1):
			M[i][j] = match(read[i-1], genome[j-1])*mismatch + max(M[i-1][j-1], D[i-1][j-1], I[i-1][j-1])
			if i < n:
				D[i][j] = max(gapOpen + M[i][j-1], gapExtension + D[i][j-1], gapOpen + I[i][j-1])
			else:
				# free trailing gaps
				D[i][j] = max(M[i][j-1], D[i][j-1], I[i][j-1])
			I[i][j] = max(gapOpen + M[i-1][j], gapExtension + I[i-1][j], gapOpen + D[i-1][j])

	### Final score - at position [n][m], minimum over the three matrices
	score = max(M[n][m], D[n][m], I[n][m])

	### Traceback
	# start traceback at [n][m]
	i = n
	j = m
	if score == M[n][m]:
		matrix = 'M'
	elif score == D[n][m]:
		matrix = 'D'
	else:
		matrix = 'I'

	# most probably there will be several trailing gaps -> we have found the score in D and
	#	now we are moving to the left in the bottom row of D
	while matrix=='D' and D[i][j-1]==score and j>0:
		j = j-1

	# real traceback
	# right end position of the read in the genome
	right = j
	# cigar string (to be built)
	cigar = []
	# while we aren't in the first row, find where did we come from
	while i>0:
		if matrix == 'M':
			if M[i][j] == match(read[i-1], genome[j-1])*mismatch + M[i-1][j-1]:
				matrix='M'
			elif M[i][j] == match(read[i-1], genome[j-1])*mismatch + D[i-1][j-1]:
				matrix = 'D'
			else:
				matrix = 'I'
			cigar.append('M')
			i = i-1
			j = j-1
		elif matrix == 'D':
			if i<n:
				if D[i][j] == gapOpen + M[i][j-1]:
					matrix = 'M'
				elif D[i][j] == gapExtension + D[i][j-1]:
					matrix='D'
				else:
					matrix = 'I'
				cigar.append('D')
			else:
				if D[i][j] == M[i][j-1]:
					matrix = 'M'
				else:	# case D[i][j] = D[i][j-1] should not happen, was already done in the previous while loop 
					matrix = 'I'
			j = j-1
		else:
			if I[i][j] == gapOpen + M[i-1][j]:
				matrix = 'M'
			elif I[i][j] == gapExtension + I[i-1][j]:
				matrix = 'I'
			else:
				matrix = 'D'
			cigar.append('I')
			i = i-1

	# left end position of the read in the genome
	left = j
	# building CIGAR string
	cigarString = compressCIGAR(''.join(cigar)[::-1])

	# return; in SAM, positions are indexed from 1, hence we add 1 to left
	return {'score':score, 'pos':left+1, 'cigar':cigarString}
