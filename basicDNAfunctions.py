complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

def reverseComplement(s):
	s = s.upper()
	revComp = ""
	for c in s:
		revComp += complement[c]
	return revComp[::-1]
