import Bio.SeqIO as s
import blosum as bl
from Bio.Align import substitution_matrices
from Bio import Align

aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

matrix = bl.BLOSUM(62)

seq = list(s.parse("raw_data\\sequences\\T6SE_Training_Pos_138.fasta","fasta"))

ress = [[] for _ in range(138)]
def Calculate_Blosum_2():
    aa = 0
    for i in seq:
        for j in seq:
            alignments = aligner.align(i.seq, j.seq)
            alignment = alignments[0]
            #print("Score = %.1f" % alignment.score)
            ress[aa].append(alignment.score)
        aa += 1
        print(aa)
    return ress

feture_evolv_type = Calculate_Blosum_2()

pssm_matr = []
with open("raw_data\\all_sorted_pssm_prof.csv") as pssm:
    lis1 = [line.split(",") for line in pssm]
    for elem in lis1:
        elem = [float(x) for x in elem]
        pssm_matr.append(elem)

i = 0
for row in pssm_matr:
    feture_evolv_type[i].extend(row)
    i+=1

