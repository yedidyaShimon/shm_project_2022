# SHM_project_2022
Code files used for creating our B.Sc. electrical engineering final project. Work done under the supervision of Mr. Daniel Eini & Prof. Gur Yaari.

## Part 1: correlation between s (silent) and rs (replacement + silent) models

Here we added the codes used to produce the correlation coefficient, together with the two files including the model parameters of the s and rs S5F model.
Filtering model parameters included in the correlation calculation, as described in the project book, is done through setting the 'minNum_S' and 'minRatio_RS_to_S' parameters at the beginning of the file

## Part 2: testing models using likelihood score

This part is consisted of four main parts:
 - codes: The main codes used for calculating the likelihood score for the SHM models on test data.
 - datasets: (links to) datasets used for training and evaluating the models.
 - model_inferred_parameters: the parameters of the model we evaluated, in the format used in 'codes' section.
-  In addition, we added the "two-phase-model" (D. Eini) code we used in an additional folder. 
 The S5F model was run using the "createMutabilityMatrix" command in the SHazaM R package. Documentation at: 
 https://shazam.readthedocs.io/en/latest/topics/createMutabilityMatrix/
 The model of N. Spisak et al. was run using codes sent to us by N. Spisak. 
 
 ### Codes
 
 - calc_score.R: calculates the likelihood score, according to model and test data. The model type and the parameters file address are defined within "if" statements along the beginning of the codes. 
 Test data address is to be defined too along the beginning of the code.
 
 - simulate with S5F.R: code used for simulating mutations in V-genes, according to S5F model, which is run on train data.
 
 - calc_mut_vec_for_score.py: a code used in addition to the "calc_score.R" file when one evaluates the 2PM model. Since the 2PM model is more complex, we needed to add 
another file, to calculate the expected mutation rates along the test data. The file creates a .csv file, which is subsequently used in "calc_score.R" file to calculate the likelihood score. 

### Datasets

- a URL link to OneDrive folder, with train data, test data, and simulated data.

### model_inferred_parameters

the folder includes model parameters after training on functional data (accounting either all mutation or only silent mutations), and on OOF data.
for S5F model we have the final parameters (mutability per 5-mer), for the 2PM we have the log file with the parameters along the run of the model. we used the "bottom line" parameters, which were inserted manually into the "calc_mut_vec_for_score.py" file.





