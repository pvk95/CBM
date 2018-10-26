'''
Implements the Burrows Wheeler Transform
'''

from operator import itemgetter

class BWT():
    valid_chars=["$","A","C","G","T"]
    def __init__(self,sequence):
        self.seq=sequence
        self.seq_len=len(sequence)
    
    def rotatestring(self):
        #Generate all rotations
        seq_list = []
        seq_list.append(self.seq)
        
        for i in range(self.seq_len-1):
            rotate_seq = seq_list[i]
            seq_temp=rotate_seq[1:]+rotate_seq[0] 
            seq_list.append(seq_temp)
        
        #Sort the list and get the indices
        indices,seq_list_sorted=zip(*sorted(enumerate(seq_list), key=itemgetter(1)))
        self.indices=list(indices)
        self.seq_list_sorted=list(seq_list_sorted)
    
    def compress(self):
        #Keep only the last column L
        self.last_col = ""
        self.F={'$':0,'A':0,'C':0,'G':0,'T':0}
        for sequence in self.seq_list_sorted:
            self.F[sequence[0]]=self.F[sequence[0]]+1
            self.last_col = self.last_col + sequence[-1]
        
        #Compress the first column
        keys=list(self.F.keys())
        vals=list(self.F.values())
        self.F[keys[0]]=[0,vals[0]]
        for i in range(1,len(keys)):
            temp=self.F[keys[i-1]][1]
            self.F[keys[i]]=[temp,temp+vals[i]]
        
        #Generate the occurence matrix of the last column
        self.last_comp={}
        for char in self.valid_chars:
            self.last_comp[char]=[0]*self.seq_len
        for i,char in enumerate(self.last_col):
            for j in range(i,self.seq_len):
                self.last_comp[char][j]=self.last_comp[char][j]+1
        
    def get_transform(self):
        self.rotatestring()
        self.compress()
        return [self.last_col,self.last_comp,self.F,self.indices,self.valid_chars]
