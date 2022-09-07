import torch
import pandas as pd
import csv
from tqdm import tqdm

from python_code.model.model import TwoPhaseModel
from python_code.model.simulation import simulation
from python_code.model.inference import inference
from python_code.model.model_utils import normalize


# --- load dataset, take only ancestor column --- #
columns_list = ['ancestor_alignment']
dataset = pd.read_csv("/home/bcrlab/giladaviv/data/11/test_data_3_S5F_s_simulate_final.csv", usecols=columns_list)
#dataset = dataset.iloc[range(100)]
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



list_of_lists_of_prob = []


# --- Calculate mutation prob vec--- #
for s in dataset['ancestor_alignment']:
    targeting_probs, _ = tpm(s)
    list_of_lists_of_prob.append(targeting_probs.tolist())

with open("/home/bcrlab/giladaviv/data/11/daniel_probs/probs_daniel_on_data_3_S5F_s_simulate.csv",'w') as f:
    writer = csv.writer(f)
    writer.writerows(list_of_lists_of_prob)