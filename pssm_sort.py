import pandas as pd
import csv

right = "raw_data\\prepare_for_pssm.csv"
need_to_fix = ["raw_data\\finished_profiles\\dpc_pssm.csv", "raw_data\\finished_profiles\\pse_pssm.csv", "raw_data\\finished_profiles\\s_fpssm.csv"]

f3 = open("raw_data\\all_pssm_sort.csv","w")

with open("raw_data\\prepare_for_pssm.csv") as f1:
    lis = [line.split(",") for line in f1]
    with open("raw_data\\finished_profiles\\dpc_pssm.csv") as f2:
        lis2 = [line.split(",") for line in f2]
        for row in lis:
            for row1 in lis2:
                if row[0] == row1[0]:
                    f3.write(str(row1[1:])+"\n")
f3.close()