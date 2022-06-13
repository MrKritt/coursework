from pssmpro.features  import get_feature

algo_type = ["pse_pssm","s_fpssm"]
pssm_dir_path = "raw_data\\pssm_profiles_all"
output_dir_path = "raw_data\\finished_profiles_all"
for i in range(len(algo_type)):
    get_feature(pssm_dir_path, algo_type[i], output_dir_path)