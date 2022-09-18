import pandas as pd

from python_code.data_utils.utils import load_multiple_sets
from python_code.definitions import imgt_regions
from python_code.model.model import TwoPhaseModel
from python_code.model.inference import inference


# --- Load data --- #
#dataset = pd.read_csv('/work/daniel/data/P3/P3_I1_S1_cloned_w_filtered_seqs.tsv', sep='\t')
dataset = pd.read_csv("/home/bcrlab/giladaviv/data/shm_oof_french_research/_alignment_updated.data_all_alignment_updated_2.csv")

# --- Filter by V gene family --- #
#v_gene_family = 'IGHV' + '4'
#dataset = dataset[dataset.germline_v_call.apply(lambda x: x.split('-')[0] == v_gene_family)]

# --- Filter sequences with too many mutations --- #
# dataset = dataset[dataset.mutations_all.apply(len) < 9]

# --- Look only on V gene --- #
v_gene_end = imgt_regions['FR3'][1]
dataset.sequence_alignment = dataset.sequence_alignment.apply(lambda x: x[:v_gene_end])
dataset.ancestor_alignment = dataset.ancestor_alignment.apply(lambda x: x[:v_gene_end])


# -- Print final dataset stats --- #
#n_sequences = dataset.shape[0]
#n_mutations = dataset.mutations_synonymous.apply(len).sum()
#print(f'Start inference! total number of sequences: {n_sequences}, total number of mutations: {n_mutations}')

# --- Init model --- #
tpm = TwoPhaseModel()

# --- Infer model params --- #
inference(tpm, 
          dataset, 
          ancestor_column='ancestorseq',
          descendant_column='ALIGNMENT',
          only_synonymous=True,
          log_postfix="_inference_franch_s")
          #log_postfix=f'-{v_gene_family}-v_only-synonymous-no_mmr')
