from ReadFastQ import readfile
from BWT import BWT
from Seed import Seeds

reads1=readfile("data_small/output_tiny_30xCov1.fq")
reads2=readfile("data_small/output_tiny_30xCov2.fq")

sequences1=list(reads1.keys())
sequences2=reads2.keys()

reference_genome=list(readfile("data_small/genome.chr22.5K.fa",form="fasta")['22_5K'])


seq=reference_genome+"$"
bwt= BWT(seq) 
[last_col,L,F,S,valid_chars] = bwt.get_transform()

seeding=Seeds(F,L,S,valid_chars)
reads=["TAG","AGC","TCA"]
locs=[]
for read in sequences1:
    res=seeding.generate_seeds(read)
    locs.append(res)
 


