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
dataset = pd.read_csv("/home/bcrlab/giladaviv/data/12/test_data_3_v_gene_final.csv", usecols=columns_list)
#dataset = dataset.iloc[[12,13]]
#dataset.ancestor_alignment = dataset.ancestor_alignment.str.replace('.', 'N')

# --- Init model --- #
tpm = TwoPhaseModel()

# --- Set model params --- #
parameters = tpm.state_dict()
parameters['phase1.motifs_prob'] = torch.tensor([1.47181172e-01, 3.62559368e-01, 4.45706209e-05, 3.45785841e-02, 4.55636305e-01, 1.00002184e-27])
parameters['phase2.replication_prob'] = torch.tensor([0.255284014056178])
parameters['phase2.short_patch_ber_prob'] = torch.tensor([1e-27])
parameters['phase2.lp_ber.profile'] = torch.tensor([0.05344895, 0.00760633, 0.02895327, 0.03101141, 0.01451822, 0.04031003
 ,0.04528457, 0.03792347, 0.02693105, 0.04197933, 0.0327675, 0.0358634
 ,0.02929906, 0.04990636, 0.02009623, 0.04087273, 0.03234872, 0.03153395
 ,0.0034317, 0.00997112, 0.02293392, 0.04552184, 0.02003084, 0.03830806
 ,0.04506391, 0.03102976, 0.04516974, 0.03654953, 0.03519223, 0.03031428
 ,0.03582852]
)
#parameters['phase2.lp_ber.profile'] = torch.concat([torch.zeros(11), torch.ones(9) / 9, torch.zeros(11)])
parameters['phase2.lp_ber.motifs_prob'] = torch.tensor([0.1382648, 0.20380293, 0.37163629, 0.28629598])
tpm.load_state_dict(parameters)


list_of_lists_of_prob = []


# --- Calculate mutation prob vec--- #
for s in dataset['ancestor_alignment']:
    targeting_probs, _ = tpm(s)
    list_of_lists_of_prob.append(targeting_probs.tolist())

with open("/home/bcrlab/giladaviv/data/12/daniel_probs_for_franch/probs_daniel_rs_train_franch_test_daniel.csv",'w') as f:
    writer = csv.writer(f)
    writer.writerows(list_of_lists_of_prob)