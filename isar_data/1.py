
def sort_pssm():

    with open("isar_data\\pos_pssm_sort.txt","w") as f3:
        with open("isar_data\\neg\\isar.pair") as f1:
            lis = [line.replace("{.psiblast | .pssm | .ascii.pssm}\n",".ascii").split("\t") for line in f1]
            with open("isar_data\\pssm_neg\\pse_pssm.csv") as f2:
                lis2 = [line.split(",") for line in f2]
                for row in lis:
                    for row1 in lis2:
                        if row[1] == row1[0]:
                            f3.write(str(row1[1:])+"\n")
                
                
                
    with open("isar_data\\pos_pssm_sort.txt","r") as f3, open("isar_data\\pssm_sort_neg\\pse_sort_pos.csv","w") as f4:
        for ii in f3:
            if ii:
                ii = ii.replace("[", "")
                ii = ii.replace("]", "")
                ii = ii.replace("'", "")
                ii = ii.replace("\\n", "")
            f4.write(ii)
            
            
def make_combine_pssm():
    with open("isar_data\\all_sorted_pssm_prof.csv","w") as final:
        with open("isar_data\\pssm_sort_all\\dpc_sort_pos.csv") as dpc, open("isar_data\\pssm_sort_all\\pse_sort_pos.csv") as pse, open("isar_data\\pssm_sort_all\\s_fsort_pos.csv") as s_f:
            for row in dpc:
                for row1 in pse:
                    for row2 in s_f:
                        final.write(row.replace("\n",",") + row1.replace("\n",",") + row2)
                        break
                    break

make_combine_pssm()