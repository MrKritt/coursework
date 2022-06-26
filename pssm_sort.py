right = "raw_data\\prepare_for_pssm_all.csv"
need_to_fix = ["raw_data\\finished_profiles_neg\\dpc_pssm.csv", "raw_data\\finished_profiles_neg\\pse_pssm.csv", "raw_data\\finished_profiles_neg\\s_fpssm.csv"]

def pssm_sorted():
    
    with open("raw_data\\all_pssm_sort.txt","w") as f3:
        with open("raw_data\\prepare_for_pssm_neg.csv") as f1:
            lis = [line.split(",") for line in f1]
            with open("raw_data\\finished_profiles_neg\\s_fpssm.csv") as f2:
                lis2 = [line.split(",") for line in f2]
                for row in lis:
                    for row1 in lis2:
                        if row[0] == row1[0]:
                            f3.write(str(row1[1:])+"\n")

    with open("raw_data\\all_pssm_sort.txt","r") as f3, open("raw_data\\pssm_sort\\pssm_f_sort_neg.csv","w") as f4:
        for ii in f3:
            if ii:
                ii = ii.replace("[", "")
                ii = ii.replace("]", "")
                ii = ii.replace("'", "")
                ii = ii.replace("\\n", "")
            f4.write(ii)

def make_combine_pssm():
    with open("raw_data\\neg_sorted_pssm_prof.csv","w") as final:
        with open("raw_data\\pssm_sort\\pssm_dpc_sort_neg.csv") as dpc, open("raw_data\\pssm_sort\\pssm_pse_sort_neg.csv") as pse, open("raw_data\\pssm_sort\\pssm_f_sort_neg.csv") as s_f:
            for row in dpc:
                for row1 in pse:
                    for row2 in s_f:
                        final.write(row.replace("\n",",") + row1.replace("\n",",") + row2)
                        break
                    break