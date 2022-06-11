import json
import math

from Bio.SeqIO import FastaIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

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
    result = ''
    X = ProteinAnalysis(sequence)
    for i in range(len(AA)):
        perent_of_amino_acid = X.get_amino_acids_percent()[AA[i]]
        result += str(perent_of_amino_acid) + ','
    return result


def CalculateDPC(sequence_dpc):
    sequence_length = len(sequence_dpc)
    result = ''
    for i in AA:
        for j in AA:
            dipeptide = i + j
            result += str(float(sequence_dpc.count(dipeptide)) / (sequence_length - 1) * 100) + ','
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
    result = ''
    temp = 1 + weight * rightpart
    for index, aaletter_char in enumerate(AA):
        result += str(AAC[aaletter_char] / temp) + ','

    return result


def CalculateQSO2SW(proteinSequence, maxlag=30, weight=0.1, distancematrix=schnaider):
    rightpart = []
    for i in range(maxlag):
        rightpart.append(SequenceOrderNumber(proteinSequence, i + 1, distancematrix))
    result = ''
    temp = 1 + weight * sum(rightpart)
    for index in range(20, 20 + maxlag):
        result += str(weight * rightpart[index - 20] / temp) + ','

    return result


def CalculateQSO1G(proteinSequence, maxlag: int = 30, weight: float = 0.1,
                   distancematrix=grantham, ):
    rightpart = 0.0
    for i in range(maxlag):
        rightpart = rightpart + SequenceOrderNumber(proteinSequence, i + 1, distancematrix)
    AAC = CalculateAAC(proteinSequence)
    result = ''
    temp = 1 + weight * rightpart
    for index, aaletter_char in enumerate(AA):
        result += str(AAC[aaletter_char] / temp) + ','

    return result


def CalculateQSO2G(proteinSequence, maxlag=30, weight=0.1,
                   distancematrix=grantham, ):
    rightpart = []
    for i in range(maxlag):
        rightpart.append(SequenceOrderNumber(proteinSequence, i + 1, distancematrix))
    result = ''
    temp = 1 + weight * sum(rightpart)
    for index in range(20, 20 + maxlag):
        result += str(weight * rightpart[index - 20] / temp) + ','

    return result


def CalculateQSO(sequence_qso, maxlag=30, weight=0.1):
    result = ''
    result += str(CalculateQSO1SW(sequence_qso, maxlag, weight, schnaider))
    result += str(CalculateQSO2SW(sequence_qso, maxlag, weight, schnaider))
    result += str(CalculateQSO1G(sequence_qso, maxlag, weight, grantham))
    result += str(CalculateQSO2G(sequence_qso, maxlag, weight, grantham))
    return result


def Calculate_All(sequence):
    resultt = ''
    resultt += str(CalculateAAC_ALL(sequence))
    resultt += str(CalculateDPC(sequence))
    resultt += str(CalculateQSO(sequence))
    return resultt


matrr = []
with open('raw_data\\sequences\\T6SE_Training_Pos_138.fasta') as fd:
    for name, sequence in FastaIO.SimpleFastaParser(fd):
        matrr.append(Calculate_All(sequence))