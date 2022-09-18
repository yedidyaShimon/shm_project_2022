import os
from datetime import datetime
import pandas as pd
import torch
import torch.nn.functional as F


from python_code.model.inference_funcs import calc_sequence_likelihood
from python_code.model.model_utils import normalize, probablize

# --- Hyperparameters --- #
schedule_steps = 3
batch_size_schedule = [50, 50, 500]
batches_each_schedule_step = [2000, 2000, 2000]
learning_rate_schedule = [1e-2, 1e-3, 1e-4]

log_path = 'results/model/log'
  def inference(model, data, ancestor_column='ancestor_alignment', descendant_column='sequence_alinment', only_synonymous=False, log_postfix=''):
    now = datetime.now().strftime("%d_%m_%Y-%H:%M:%S")
    log_dir = os.path.join(log_path, now + log_postfix)
    os.makedirs(log_dir)
    log_csv = os.path.join(log_dir, 'log.csv')

    print('-------------------------- Training start! ---------------------------')
    for name, param in model.named_parameters():
        if not name.count('motifs_prob'):  # Not too long for print
            print(f'{name}: {param}')  
    print('----------------------------------------------------------------------')

    for step_counter in range(schedule_steps):

        batch_size = batch_size_schedule[step_counter]
        learning_rate = learning_rate_schedule[step_counter]
        batches_this_step = batches_each_schedule_step[step_counter]
        batch_counter = 0
        sample_counter = 0

        optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
        loss = torch.tensor([0.], requires_grad=True)
        mutations_counter = 0

        while batch_counter < batches_this_step:
            for _, row in data.sample(frac=1, ignore_index=True).iterrows():
                sample_counter += 1

                # Compute likelihood and accumulate loss
                targeting_probs, replication_probs = model(row[ancestor_column])
                negative_log_likelihood, n_mutations = calc_sequence_likelihood(row[ancestor_column], 
                                                                                row[descendant_column],
                                                                                targeting_probs, 
                                                                                replication_probs,
                                                                                only_synonymous=only_synonymous)

                loss = loss + negative_log_likelihood
                mutations_counter += n_mutations
                
                if sample_counter == batch_size:
                    sample_counter = 0

                    # Scale loss magnitude by the number of mutations
                    loss = loss / mutations_counter
                    mutations_counter = 0

                    # Backpropagation
                    optimizer.zero_grad()
                    loss.backward(retain_graph=True)
                    optimizer.step()

                    # Prints
                    print(f'Learning step: {step_counter}, Batch: {batch_counter} / {batches_this_step}, Loss: {loss}')
                    
                    # Casting to parameter space 
                    with torch.no_grad():
                        for name, param in model.named_parameters():
                            if param.requires_grad:
                                param.copy_(probablize(param))

                    print('----------------------------------------------------------------------')
                    for name, param in model.named_parameters():
                        print(f'{name}: {param}')
                    print('----------------------------------------------------------------------')

                    # Log
                    log_line = {'loss': [loss.data.detach().item()]}

                    for name, param in model.named_parameters():
                        value = param.detach().numpy()
                        if len(param) > 1:
                            value = [value]
                        log_line[name] = value

                    pd.DataFrame(log_line).to_csv(log_csv, index=False, mode='a', header=not os.path.isfile(log_csv))

                    # Save parameters
                    if not (batch_counter % 100):
                        save_parame_file = os.path.join(log_dir, 'state_dict_' + str(batch_counter))
                        torch.save(model.state_dict(), save_parame_file)

                    # Promote batch counter
                    batch_counter += 1
                    
                    # Reset loss
                    loss = torch.tensor([0.], requires_grad=True)

                    # Exit loop condition
                    if batch_counter == batches_this_step:
                        break





