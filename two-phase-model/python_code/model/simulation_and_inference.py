import torch
import pandas as pd
from tqdm import tqdm

from python_code.model.model import TwoPhaseModel
from python_code.model.simulation import simulation
from python_code.model.inference import inference
from python_code.model.model_utils import normalize



def simulation_and_inference(dataset, only_synonymous=False, log_postfix=''):
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
    # --- Reset model params --- #
    tpm = TwoPhaseModel()
    
    # --- Infer model params --- #
    inference(tpm, dataset, ancestor_column='ancestor_alignment', descendant_column='simulated_sequence', only_synonymous=only_synonymous, log_postfix=log_postfix)
