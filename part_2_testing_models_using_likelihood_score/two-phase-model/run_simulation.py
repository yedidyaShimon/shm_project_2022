import pandas as pd

from python_code.data_utils.utils import load_multiple_sets
#from python_code.model.simulation_and_inference_fivemers import simulation_and_inference
from python_code.model.simulation_and_inference import simulation_and_inference
#from python_code.model.simulation_and_inference_extended import simulation_and_inference


# --- Load data --- #
# columns_list = ['sequence_alignment', 'ancestor_alignment', 'mutations_all']
#columns_list = ['sequence_alignment', 'ancestor_alignment', 'num_mutations']
columns_list = ['sequence_alignment', 'ancestor_alignment', 'num_of_mutation']
#all_sets = pd.read_csv('data/final_sets.csv')
#all_sets = all_sets[all_sets.sample_id != 'P4_I19_S1']  # defect set

#all_sets = all_sets[all_sets.sample_id == 'P3_I2_S2']

# --- Choose portion --- #
#paths = all_sets[all_sets.study == 'influenza'].path      # by study
#paths = paths[:1]                                         # by number of sets

#dataset = load_multiple_sets(paths, columns_list)
dataset = pd.read_csv("/home/bcrlab/giladaviv/two-phase-model/data/test_data_3_v_gene_final.csv", usecols=columns_list)
dataset.ancestor_alignment = dataset.ancestor_alignment.str.replace('.', 'N')
dataset.rename(columns = {'num_of_mutation':'num_mutations'}, inplace = True)

# --- Init model --- #
tpm = TwoPhaseModel()

# --- Set model params --- #
parameters = tpm.state_dict()
parameters['phase1.motifs_prob'] = torch.tensor([0.1, 0.4, 0.05, 0.1, 0.3, 0.05])
parameters['phase2.replication_prob'] = torch.tensor([0.5])
parameters['phase2.short_patch_ber_prob'] = torch.tensor([0.2])
parameters['phase2.lp_ber.profile'] = torch.concat([torch.zeros(11), torch.ones(9) / 9, torch.zeros(11)])
parameters['phase2.lp_ber.motifs_prob'] = torch.tensor([1e-19, 1e-19, 0.45, 0.55])
tpm.load_state_dict(parameters)

# --- Simulate mutations --- #
tqdm.pandas()
dataset['simulated_sequence'] = dataset.progress_apply(lambda row: simulation(sequence=row.ancestor_alignment,
                                                                              n_mutations=row.num_mutations,
                                                                              model=tpm), axis=1)
