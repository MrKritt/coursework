from pssmpro.features  import get_feature

def get_dpc_pse_s_fpssm(algo_type,pssm_dir_path,output_dir_path):
    
    for i in range(len(algo_type)):
        get_feature(pssm_dir_path, algo_type[i], output_dir_path)

get_dpc_pse_s_fpssm(["dpc_pssm","pse_pssm","s_fpssm"],"raw_data\\pssm_profiles_neg", "raw_data\\finished_profiles_neg")