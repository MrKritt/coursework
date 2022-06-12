right = "raw_data\\prepare_for_pssm.csv"
need_to_fix = ["raw_data\\finished_profiles\\dpc_pssm.csv", "raw_data\\finished_profiles\\pse_pssm.csv", "raw_data\\finished_profiles\\s_fpssm.csv"]

f3 = open("raw_data\\all_pssm_sort.txt","w")

with open("raw_data\\prepare_for_pssm.csv") as f1:
    lis = [line.split(",") for line in f1]
    with open("raw_data\\finished_profiles\\pse_pssm.csv") as f2:
        lis2 = [line.split(",") for line in f2]
        for row in lis:
            for row1 in lis2:
                if row[0] == row1[0]:
                    f3.write(str(row1[1:])+"\n")
f3.close()
   
with open("raw_data\\all_pssm_sort.txt","r") as f3, open("raw_data\\pssm_pse_sort.csv","w") as f4:
    for ii in f3:
        if ii:
            ii = ii.replace("[", "")
            ii = ii.replace("]", "")
            ii = ii.replace("'", "")
            ii = ii.replace("\\n", "")
        f4.write(ii)


final = open("raw_data\\all_sorted_pssm_prof.csv","w")

with open("raw_data\\pssm_dps_sort.csv") as dpc, open("raw_data\\pssm_pse_sort.csv") as pse, open("raw_data\\pssm_f_sort.csv") as s_f:
    for row in dpc:
        for row1 in pse:
            for row2 in s_f:
                final.write(row.replace("\n",",") + row1.replace("\n",",") + row2)
                break
            break

final.close()