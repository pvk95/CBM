'''
TO DO: Define exceptions
'''
class Seeds():
    '''
    __first: The first column
    __last: The last column data structure
    __S: The suffix array
    __valid_chars: allowed characters
    '''
    __first={}
    __last={}
    __S={}
    __valid_chars=[]
    def __init__(self,first,last,S,valid_chars):
      self.__first=first
      self.__last=last
      self.__S=S
      self.__valid_chars=valid_chars
    def get_range(self,ch):
        '''
        F is the first column frequences stored as a dict
        '''
        F=self.__first
        start=F[ch][0]
        end=F[ch][1]-1
        return [start,end]
    def get_index(self,ch,r):
        '''
        F is the first column frequences stored as a dict
        ch is the character to search for
        r is the rank
        '''
        F=self.__first
        if(r==0): 
            raise Exception("0th occurence is undefined")
        idx=F[ch][0]+r-1
        return idx
    def get_rank(self,ch,row):
        '''
        L is the rank matrix
        '''
        L=self.__last
        if(row<0):
            return 0
        return L[ch][row]
        
    
    def generate_seeds(self,Q):
        '''
        Q is the query sequence
        S is the Suffix array
        '''
        valid_chars=self.__valid_chars
        S=self.__S
        index=len(Q)-1
        char=Q[index]
        if char not in valid_chars:
            raise Exception("String character mismatch")
        #The range of char in F [start,end] both inclusive
        [start,end]=self.get_range(char)
        locations=[]
        while index>0:
            index=index-1
            char=Q[index]
            if char not in valid_chars:
                raise Exception("Invalid character " + char)
            before_rank=self.get_rank(char,start-1)
            start_rank=self.get_rank(char,start)
            end_rank=self.get_rank(char,end)
            if(end_rank==before_rank):
                # print("Seed not found")
                return None
            if(start_rank==0):
                start_rank=1
            start=self.get_index(char,start_rank)
            end=self.get_index(char,end_rank)
        
        for i in range(start+1,end+1):
            locations.append(S[i])
        assert(locations!=None)
        return locations



'''
F={}
F['A']=[0,2]
F['C']=[2,3]
F['G']=[3,5]
F['T']=[5,6]

L={}
L['A']=[0,0,0,1,2,2,2]
L['C']=[0,0,0,0,0,0,1]
L['G']=[0,1,2,2,2,2,2]
L['T']=[1,1,1,1,1,1,1]

S=[1,3,5,2,4,0,6]

valid_chars=["A","C","G","T","$","-"]
'''

