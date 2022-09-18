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
dataset = pd.read_csv("/home/bcrlab/giladaviv/two-phase-model/data/train_data_3_with_mut_num.csv", usecols=columns_list)
dataset.ancestor_alignment = dataset.ancestor_alignment.str.replace('.', 'N')
dataset.rename(columns = {'num_of_mutation':'num_mutations'}, inplace = True)

simulation_and_inference(dataset, only_synonymous=True, log_postfix='train3_2')