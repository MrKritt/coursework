import copy
from Bio.SeqIO import FastaIO



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

with open('raw_data\\sequences\\T6SE_Training_Pos_138.fasta') as fd:
    for name, sequence in FastaIO.SimpleFastaParser(fd):
        physico_chemical_type.append(Calculate_CTD(sequence))

