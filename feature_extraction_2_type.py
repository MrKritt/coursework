import Bio.SeqIO as s
import blosum as bl
from Bio.SeqIO import FastaIO

AA = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']


matrix = bl.BLOSUM(62)

seq = list(s.parse("raw_data\\sequences\\T6SE_Training_Pos_138.fasta","fasta"))

def Calculate_Blosum(sequence):
    res = []
    for i in sequence:
        for j in sequence:
            dipeptide = i + j
            res.append(matrix[dipeptide])
            
    return res

matrr= []
with open('raw_data\\sequences\\T6SE_Training_Pos_138.fasta') as fd:
    for name, sequence in FastaIO.SimpleFastaParser(fd):
        matrr.append(Calculate_Blosum(sequence))

print(matrr)