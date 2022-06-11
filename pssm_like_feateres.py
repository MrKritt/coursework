from pssmpro.features  import get_feature

algo_type = ["dpc_pssm"]
r = [i for i in range(138)]
pssm_dir_path = "raw_data\\pssm_profiles_pos"
output_dir_path = "raw_data\\finished_profiles"
for i in range(len(algo_type)):
    get_feature(pssm_dir_path, algo_type[i], output_dir_path)