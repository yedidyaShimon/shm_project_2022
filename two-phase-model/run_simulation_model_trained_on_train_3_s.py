import pandas as pd
import torch
from python_code.model.model import TwoPhaseModel
from python_code.model.simulation import simulation
from tqdm import tqdm
from python_code.data_utils.utils import load_multiple_sets
#from python_code.model.simulation_and_inference_fivemers import simulation_and_inference
from python_code.model.simulation_and_inference import simulation_and_inference
#from python_code.model.simulation_and_inference_extended import simulation_and_inference


# --- Load data --- #
# columns_list = ['sequence_alignment', 'ancestor_alignment', 'mutations_all']
#columns_list = ['sequence_alignment', 'ancestor_alignment', 'num_mutations']
columns_list = ['sequence_alignment', 'ancestor_alignment', 'num_mutations']
#all_sets = pd.read_csv('data/final_sets.csv')
#all_sets = all_sets[all_sets.sample_id != 'P4_I19_S1']  # defect set

#all_sets = all_sets[all_sets.sample_id == 'P3_I2_S2']

# --- Choose portion --- #
#paths = all_sets[all_sets.study == 'influenza'].path      # by study
#paths = paths[:1]                                         # by number of sets

#dataset = load_multiple_sets(paths, columns_list)
dataset = pd.read_csv("/home/bcrlab/giladaviv/data/11/test_data_3_v_gene_final.csv", usecols=columns_list)
dataset.ancestor_alignment = dataset.ancestor_alignment.str.replace('.', 'N')

# --- Init model --- #
tpm = TwoPhaseModel()

# --- Set model params --- #
parameters = tpm.state_dict()
parameters['phase1.motifs_prob'] = torch.tensor([0.07385295,0.41977694,0.00207753,0.07376884,0.40280927,0.02771447])
parameters['phase2.replication_prob'] = torch.tensor([0.290345771215344])
parameters['phase2.short_patch_ber_prob'] = torch.tensor([0.161657818615857])
parameters['phase2.lp_ber.profile'] = torch.tensor([3.12888413e-02,1.00000685e-27,4.57585333e-02,1.00000685e-27,1.00000685e-27,8.19214119e-02,1.42577224e-03,2.46750825e-02,6.23929509e-02,1.00000685e-27,4.54112612e-02,3.12098239e-02,4.45668159e-02,6.90164349e-02,1.51027372e-02,1.00000685e-27,1.00000685e-27,8.43248247e-02,4.30101132e-02,3.44957848e-02,2.85203795e-02,3.75771194e-02,1.56566244e-03,8.15386283e-02,4.43226804e-02,9.06660518e-05,3.97438039e-04,3.36009330e-02,2.19470925e-02,1.10586576e-01,2.52524360e-02])
#parameters['phase2.lp_ber.profile'] = torch.concat([torch.zeros(11), torch.ones(9) / 9, torch.zeros(11)])
parameters['phase2.lp_ber.motifs_prob'] = torch.tensor([0.08557295,0.05702175,0.40505158,0.45235372])
tpm.load_state_dict(parameters)


# --- Simulate mutations --- #
tqdm.pandas()
dataset['simulated_sequence'] = dataset.progress_apply(lambda row: simulation(sequence=row.ancestor_alignment,
                                                                              n_mutations=row.num_mutations,
                                                                              model=tpm), axis=1)
dataset.to_csv("/home/bcrlab/giladaviv/data/11/test_data_3_daniel_s_simulate_final.csv",index = False)

