import pandas as pd

from python_code.data_utils.utils import load_multiple_sets
from python_code.definitions import imgt_regions
from python_code.model.model import TwoPhaseModel
from python_code.model.inference import inference


# --- Load data --- #
#dataset = pd.read_csv('/work/daniel/data/P3/P3_I1_S1_cloned_w_filtered_seqs.tsv', sep='\t')
columns_list = ['sequence_alignment', 'ancestor_alignment', 'num_mutations', 'length_equal']
dataset = pd.read_csv("/home/bcrlab/giladaviv/data/12/data_french_all_with_vcall_with_gap_insteadof_N_final2.csv", usecols=columns_list)
dataset.ancestor_alignment = dataset.ancestor_alignment.str.replace('.', 'N')
dataset.rename(columns = {'num_of_mutation':'num_mutations'}, inplace = True)
dataset = dataset[dataset['length_equal']]
# --- Filter by V gene family --- #
#v_gene_family = 'IGHV' + '4'
#dataset = dataset[dataset.germline_v_call.apply(lambda x: x.split('-')[0] == v_gene_family)]

# --- Filter sequences with too many mutations --- #
# dataset = dataset[dataset.mutations_all.apply(len) < 9]

# --- Look only on V gene --- #

# -- Print final dataset stats --- #
#n_sequences = dataset.shape[0]
#n_mutations = dataset.mutations_synonymous.apply(len).sum()
#print(f'Start inference! total number of sequences: {n_sequences}, total number of mutations: {n_mutations}')

# --- Init model --- #
tpm = TwoPhaseModel()

# --- Infer model params --- #
inference(tpm, 
          dataset, 
          ancestor_column='ancestor_alignment', 
          descendant_column='sequence_alignment', 
          only_synonymous=False,
          log_postfix="inference_french_rs")
          #log_postfix=f'-{v_gene_family}-v_only-synonymous-no_mmr')
