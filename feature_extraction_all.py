import json
import math
import csv

from Bio.SeqIO import FastaIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

import Bio.SeqIO as s
import blosum as bl
from Bio.Align import substitution_matrices
from Bio import Align
import copy

AA = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

massive = {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'E': 0, 'Q': 0, 'G': 0, 'H': 0, 'I': 0, 'L': 0, 'K': 0, 'M': 0,
           'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0}

with open("raw_data\\garntham.json", "r") as f:
    schnaider = json.load(f)

with open("raw_data\\schnider_vreder.json", "r") as f:
    grantham = json.load(f)


def CalculateAAC(sequence_aac):
    X = ProteinAnalysis(sequence_aac)
    for i in range(len(AA)):
        perent_of_amino_acid = X.get_amino_acids_percent()[AA[i]]
        massive[AA[i]] = perent_of_amino_acid
    return massive


def CalculateAAC_ALL(sequence):
    result = []
    X = ProteinAnalysis(sequence)
    for i in range(len(AA)):
        perent_of_amino_acid = X.get_amino_acids_percent()[AA[i]]
        result.append(perent_of_amino_acid)
    return result


def CalculateDPC(sequence_dpc):
    sequence_length = len(sequence_dpc)
    result = []
    for i in AA:
        for j in AA:
            dipeptide = i + j
            result.append((float(sequence_dpc.count(dipeptide)) / (sequence_length - 1) * 100))
    return result


def SequenceOrderNumber(sequence, d: int = 1, distancematrix=schnaider):
    num_protein = len(sequence)
    tau = 0.0
    for i in range(num_protein - d):
        temp1 = sequence[i]
        temp2 = sequence[i + d]
        tau = tau + math.pow(distancematrix[temp1 + temp2], 2)
    return round(tau, 3)


def CalculateQSO1SW(proteinSequence, maxlag=30, weight=0.1, distancematrix=schnaider):
    rightpart = 0.0
    for i in range(maxlag):
        rightpart = rightpart + SequenceOrderNumber(proteinSequence, i + 1, distancematrix)
    AAC = CalculateAAC(proteinSequence)
    result = []
    temp = 1 + weight * rightpart
    for index, aaletter_char in enumerate(AA):
        result.append(AAC[aaletter_char] / temp)

    return result


def CalculateQSO2SW(proteinSequence, maxlag=30, weight=0.1, distancematrix=schnaider):
    rightpart = []
    for i in range(maxlag):
        rightpart.append(SequenceOrderNumber(proteinSequence, i + 1, distancematrix))
    result = []
    temp = 1 + weight * sum(rightpart)
    for index in range(20, 20 + maxlag):
        result.append(weight * rightpart[index - 20] / temp)

    return result


def CalculateQSO1G(proteinSequence, maxlag: int = 30, weight: float = 0.1,
                   distancematrix=grantham, ):
    rightpart = 0.0
    for i in range(maxlag):
        rightpart = rightpart + SequenceOrderNumber(proteinSequence, i + 1, distancematrix)
    AAC = CalculateAAC(proteinSequence)
    result = []
    temp = 1 + weight * rightpart
    for index, aaletter_char in enumerate(AA):
        result.append(AAC[aaletter_char] / temp)

    return result


def CalculateQSO2G(proteinSequence, maxlag=30, weight=0.1,
                   distancematrix=grantham, ):
    rightpart = []
    for i in range(maxlag):
        rightpart.append(SequenceOrderNumber(proteinSequence, i + 1, distancematrix))
    result = []
    temp = 1 + weight * sum(rightpart)
    for index in range(20, 20 + maxlag):
        result.append(weight * rightpart[index - 20] / temp)

    return result


def CalculateQSO(sequence_qso, maxlag=30, weight=0.1):
    result = (CalculateQSO1SW(sequence_qso, maxlag, weight, schnaider))
    result.extend(CalculateQSO2SW(sequence_qso, maxlag, weight, schnaider))
    result.extend(CalculateQSO1G(sequence_qso, maxlag, weight, grantham))
    result.extend(CalculateQSO2G(sequence_qso, maxlag, weight, grantham))
    return result


def Calculate_All(sequence):
    resultt = (CalculateAAC_ALL(sequence))
    resultt.extend(CalculateDPC(sequence))
    resultt.extend(CalculateQSO(sequence))
    return resultt


fetute_seq_type = []
with open('raw_data\\sequences\\T6SE_Training_Neg_1112.fasta') as fd:
    for name, sequence in FastaIO.SimpleFastaParser(fd):
        fetute_seq_type.append(Calculate_All(sequence))

# evolve-information type

aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

matrix = bl.BLOSUM(62)

seq = list(s.parse("raw_data\\sequences\\T6SE_Training_Neg_1112.fasta","fasta"))


def Calculate_Blosum_2():
    ress = [[] for _ in range(1112)]
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
with open("raw_data\\neg_sorted_pssm_prof.csv") as pssm:
    lis1 = [line.split(",") for line in pssm]
    for elem in lis1:
        elem = [float(x) for x in elem]
        pssm_matr.append(elem)

i = 0
for row in pssm_matr:
    feture_evolv_type[i].extend(row)
    i+=1
    
#Physico-chemical type

Hydrophobicity = {"1": "RKEDQN", "2": "GASTPHY", "3": "CLVIMFW"}
NormalizedVDWV = {"1": "GASTPDC", "2": "NVEQIL", "3": "MHKFRYW"}
Polarity = {"1": "LIFWCMVY", "2": "PATGS", "3": "HQRKNED"}
Polarizability = {"1": "GASDT", "2": "CPNVEQIL", "3": "KMHFRYW"}
Charge = {"1": "KR", "2": "ANCQGHILMFPSTWYV", "3": "DE"}
SecondaryStr = {"1": "EALMQKRH", "2": "VIYCWFT", "3": "GNPSD"}
SolventAccessibility = {"1": "ALFCGIVW", "2": "RKQEND", "3": "MPSTHY"}

aa_class = [
    Hydrophobicity,
    NormalizedVDWV,
    Polarity,
    Charge,
    SecondaryStr,
    SolventAccessibility,
    Polarizability,
]

def StringtoNum(ProteinSequence, aa_class):
    hard_protein_sequence = copy.deepcopy(ProteinSequence)
    for k, m in list(aa_class.items()):
        for index in m:
            hard_protein_sequence = hard_protein_sequence.replace(index, k)
    return hard_protein_sequence

def CalculateComposition(ProteinSequence, aa_class):
    tprotein_sequence = StringtoNum(ProteinSequence, aa_class)
    result = []
    num = len(tprotein_sequence)
    result.append(round(float(tprotein_sequence.count("1")) / num, 3))
    result.append(round(float(tprotein_sequence.count("2")) / num, 3))
    result.append(round(float(tprotein_sequence.count("3")) / num, 3))
    return result

def CalculateTransition(ProteinSequence, aa_class):
    tprotein_sequence = StringtoNum(ProteinSequence, aa_class)
    result = []
    num = len(tprotein_sequence)
    ctd = tprotein_sequence
    result.append(round(float(ctd.count("12") + ctd.count("21")) / (num - 1), 3))
    result.append(round(float(ctd.count("13") + ctd.count("31")) / (num - 1), 3))
    result.append(round(float(ctd.count("23") + ctd.count("32")) / (num - 1), 3))
    return result

def CalculateC(ProteinSequence):
    result = []
    result.extend(CalculateComposition(sequence,Hydrophobicity))
    result.extend(CalculateComposition(sequence,NormalizedVDWV))
    result.extend(CalculateComposition(sequence,Polarity))
    result.extend(CalculateComposition(sequence,Polarizability))
    result.extend(CalculateComposition(sequence,Charge))
    result.extend(CalculateComposition(sequence,SecondaryStr))
    result.extend(CalculateComposition(sequence,SolventAccessibility))
    return result

def CalculateT(ProteinSequence):
    result = []
    result.extend(CalculateTransition(sequence,Hydrophobicity))
    result.extend(CalculateTransition(sequence,NormalizedVDWV))
    result.extend(CalculateTransition(sequence,Polarity))
    result.extend(CalculateTransition(sequence,Polarizability))
    result.extend(CalculateTransition(sequence,Charge))
    result.extend(CalculateTransition(sequence,SecondaryStr))
    result.extend(CalculateTransition(sequence,SolventAccessibility))
    return result


def Calculate_CTD(sequence):
    res = []
    res.extend(CalculateC(sequence))
    res.extend(CalculateT(sequence))
    return res

physico_chemical_type = []

with open('raw_data\\sequences\\T6SE_Training_Neg_1112.fasta') as fd:
    for name, sequence in FastaIO.SimpleFastaParser(fd):
        physico_chemical_type.append(Calculate_CTD(sequence))

prepared_features = fetute_seq_type
i = 0
for row in feture_evolv_type:
    prepared_features[i].extend(row)
    i+=1
    
i = 0
for row in physico_chemical_type:
    prepared_features[i].extend(row)
    i += 1
    
with open("prepared_data\\train_data_neg.csv","w") as td:
    csv_writer = csv.writer(td,delimiter = ",", lineterminator="\r")
    for line in prepared_features:
        csv_writer.writerow(line)
        