# -*- coding: utf-8 -*-
"""
Created on Fri May 24 12:51:24 2024

@author: b1095820
"""

# %% Set-up: Import necessary libraries and functions
# Standard libraries
# os -  Working with file or folder directories
import os
# random - Pick random values between the lower and the upper bound
import random
# warnings - Ignore warnings about infeasible solutions
# import warnings
# datetime - Display the time it took for the analysis to run
from datetime import datetime
# pandas - Work with dataframes
import pandas as pd
# numpy - Working with arrays
import numpy as np
# Matplotlib - For plotting figures
import matplotlib.pyplot as plt

# Third-party library imports
# cobrapy- Run FBA
import cobra

# Local application/library-specific imports
# escher plots
# from escher.plots import Builder

# %% Set directory

os.system("pwd")
# current_wd = os.getcwd()
# print(current_wd) 
os.chdir('c:\\Users\\b1095820\\Documents\\DGTX\\BR24_E13-20')
current_wd = os.getcwd()
# os.chdir(current_wd)

# %% load cho model

models = {
    # "iCHO1766": "iCHOv1_final.xml",
    "iCHO2441": "iCHO2441.xml",
    # "CHO-K1": "iCHOv1_K1_final.xml",
    # "CHOmpact": "CHOsmallmodel.json"
    # "CHOmpact_small": "CHOsmallmodel_activity4.json"
    # "K1par-0mMCD": "iCHO_K1par-0mMCD.xml"
    }

for model_name, model_file in models.items():
    print(model_name)
    MODEL_PATH = "cho_gems/" + model_file
    _, file_extension = os.path.splitext(MODEL_PATH)

    if file_extension == '.xml':
        model_orig = cobra.io.read_sbml_model(MODEL_PATH)
    elif file_extension == '.json':
        model_orig = cobra.io.load_json_model(MODEL_PATH)
    else:
        # Handle unsupported file types or other extensions
        print(f"Unsupported file extension: {file_extension}")
        continue
    print(f"Loaded {model_name} model successfully!")

# %% create copy of model
model = model_orig.copy()

# %% update igg production reaction equations
# hc equations
# print("Original reaction in the model:", model.reactions.igg_hc.reaction)

# igg_hc_new_equation = '20.0 ala_L_c + 11.0 arg_L_c + 19.0 asn_L_c + 21.0 asp_L_c + 460.0 atp_c + 11.0 cys_L_c + 16.0 gln_L_c + 20.0 glu_L_c + 25.0 gly_c + 916.0 gtp_c + 917.0 h2o_c + 10.0 his_L_c + 9.0 ile_L_c + 35.0 leu_L_c + 35.0 lys_L_c + 6.0 met_L_c + 15.0 phe_L_c + 37.0 pro_L_c + 49.0 ser_L_c + 41.0 thr_L_c + 10.0 trp_L_c + 16.0 tyr_L_c + 44.0 val_L_c --> adp_c + 459.0 amp_c + 916.0 gdp_c + 917.0 h_c + igg_hc_r + 917.0 pi_c + 459.0 ppi_c'  # Replace with the new equation you want to set

# # Access the existing reaction in the model
# existing_reaction = model.reactions.get_by_id('igg_hc')
# # Update the reaction equation
# existing_reaction.reaction = igg_hc_new_equation

# print("Updated reaction in the model:", model.reactions.igg_hc.reaction)


# #lc equations
# print("Original reaction in the model:", model.reactions.igg_lc.reaction)

# igg_lc_new_equation = '12.0 ala_L_c + 6.0 arg_L_c + 5.0 asn_L_c + 10.0 asp_L_c + 217.0 atp_c + 5.0 cys_L_c + 12.0 gln_L_c + 9.0 glu_L_c + 16.0 gly_c + 430.0 gtp_c + 431.0 h2o_c + 3.0 his_L_c + 6.0 ile_L_c + 14.0 leu_L_c + 14.0 lys_L_c + 2.0 met_L_c + 10.0 phe_L_c + 1.0 pro_L_c + 32.0 ser_L_c + 19.0 thr_L_c + 2.0 trp_L_c + 10.0 tyr_L_c + 15.0 val_L_c --> adp_c + 216.0 amp_c + 430.0 gdp_c + 431.0 h_c + igg_lc_r + 431.0 pi_c + 216.0 ppi_c'  # Replace with the new equation you want to set

# # Access the existing reaction in the model
# existing_reaction = model.reactions.get_by_id('igg_lc')
# # Update the reaction equation
# existing_reaction.reaction = igg_lc_new_equation

# print("Updated reaction in the model:", model.reactions.igg_lc.reaction)
# %% load aa and metabolite data
#rates = pd.read_csv("data/aa_rates_reordered_nottshifted_metabolites.csv")
# rates = pd.read_csv("data/aa_rates_reordered_data2_nottshifted.csv")
rates = pd.read_csv("data/aa_rates_reordered_data2_nottshifted.csv")

# Define the specific values to be removed
# values_to_remove = ['Growth_rate', 'Titer']  # Replace these with the actual values you want to remove
values_to_remove = ['Growth_rate'] 
# Remove rows where 'aa_rates' column contains any of the specified values
# rates_filtered = rates[~rates['AA_meta'].isin(values_to_remove)]

# Alternatively, you can modify the existing DataFrame in place
rates.drop(rates[rates['AA_meta'].isin(values_to_remove)].index, inplace=True)

# %% check all names of reactions for which you want to set the bounds
#iCHO1766
# uptake_names = {
#     "Ala": "EX_ala_L_e_",
#     "NH3": "EX_nh4_e_",
#     "Arg": "EX_arg_L_e_",
#     "Asn": "EX_asn_L_e_",
#     "Asp": "EX_asp_L_e_",
#     #"Cysteine": "EX_cys_L_e_",
#     "GLC": "EX_glc_e_",
#     "Glu": "EX_glu_L_e_",
#     "Gln": "EX_gln_L_e_",
#     "Gly": "EX_gly_e_",
#     "His": "EX_his_L_e_",
#     "Ile": "EX_ile_L_e_",
#     "LAC": "EX_lac_L_e_",
#     "Leu": "EX_leu_L_e_",
#     "Lys": "EX_lys_L_e_",
#     "Met": "EX_met_L_e_",
#     "Phe": "EX_phe_L_e_",
#     "Pro": "EX_pro_L_e_",
#     "Ser": "EX_ser_L_e_",
#     "Thr": "EX_thr_L_e_",
#     "Trp": "EX_trp_L_e_",
#     "Tyr": "EX_tyr_L_e_",
#     "Val": "EX_val_L_e_",
#     "Titer": "DM_igg_g_"
# }
# #CHOmpact
# uptake_names = {
#     "Ala": "F110",
#     "NH3": "F108",
#     "Arg": "F111",
#     "Asn": "F112",
#     "Asp": "F113",
#     #"Cysteine": "EX_cys_L_e_",
#     "GLC": "F105",
#     "Glu": "F115",
#     "Gln": "F114",
#     "Gly": "F116",
#     "His": "F117",
#     "Ile": "F118",
#     "LAC": "F107",
#     "Leu": "F119",
#     "Lys": "F120",
#     "Met": "F121",
#     "Phe": "F122",
#     "Pro": "F123",
#     "Ser": "F124",
#     "Thr": "F125",
#     "Trp": "F126",
#     "Tyr": "F127",
#     "Val": "F128",
#     "Titer": "F143"
# }

# for iCHO2441
uptake_names = {
    "Ala": "EX_ala_L(e)",
    "NH3": "EX_nh4(e)",
    "Arg": "EX_arg_L(e)",
    "Asn": "EX_asn_L(e)",
    "Asp": "EX_asp_L(e)",
    #"Cysteine": "EX_cys_L_e_",
    "GLC": "EX_glc(e)",
    "Glu": "EX_glu_L(e)",
    "Gln": "EX_gln_L(e)",
    "Gly": "EX_gly(e)",
    "His": "EX_his_L(e)",
    "Ile": "EX_ile_L(e)",
    "LAC": "EX_lac_L(e)",
    "Leu": "EX_leu_L(e)",
    "Lys": "EX_lys_L(e)",
    "Met": "EX_met_L(e)",
    "Phe": "EX_phe_L(e)",
    "Pro": "EX_pro_L(e)",
    "Ser": "EX_ser_L(e)",
    "Thr": "EX_thr_L(e)",
    "Trp": "EX_trp_L(e)",
    "Tyr": "EX_tyr_L(e)",
    "Val": "EX_val_L(e)",
    "Titer": "DM_igg[g]"
}

# taken_up = [
#     "EX_gln_L_e_",
#     "EX_cys_L_e_",
#     "EX_arg_L_e_",
#     "EX_asn_L_e_",
#     "EX_asp_L_e_",
#     "EX_glc_e_",
#     "EX_glu_L_e_",
#     "EX_h_e_",
#     "EX_h2o_e_",
#     "EX_his_L_e_",
#     "EX_ile_L_e_",
#     "EX_leu_L_e_",
#     "EX_lys_L_e_",
#     "EX_met_L_e_",
#     "EX_o2_e_",
#     "EX_phe_L_e_",
#     "EX_pi_e_",
#     "EX_pro_L_e_",
#     "EX_ser_L_e_",
#     "EX_thr_L_e_",
#     "EX_trp_L_e_",
#     "EX_tyr_L_e_",
#     "EX_val_L_e_",
#     "EX_lnlc_e_",
#     "EX_lnlnca_e_",
#     "EX_Tyr_ggn_e_"
# ]
# for ex in model.reactions:
#     if ex.reversibility or ex.id in taken_up:
#         print(ex.bounds)

# # Turn off igg and epo production
# reaction = model.reactions.get_by_id('DM_igg[g]')

# # Set the lower and upper bounds to 0
# reaction.lower_bound = 0
# reaction.upper_bound = 0
# print(f"Reaction {reaction.id} bounds set to: [{reaction.lower_bound}, {reaction.upper_bound}]")


reaction = model.reactions.get_by_id('DM_epo[g]')

# Set the lower and upper bounds to 0
reaction.lower_bound = 0
reaction.upper_bound = 0
print(f"Reaction {reaction.id} bounds set to: [{reaction.lower_bound}, {reaction.upper_bound}]")

# model.reactions.DM_igg_g_.lower_bound = 0
# model.reactions.DM_igg_g_.upper_bound = 0
# model.reactions.DM_epo_g_.lower_bound = 0
# model.reactions.DM_epo_g_.upper_bound = 0


# Set the objective function
#if strain in producers:
# model.objective = "biomass_cho_producing" #index 6618


# model.reactions.biomass_cho.upper_bound = 0
# model.reactions.biomass_cho.lower_bound = 0
#else:
# model.objective = "biomass_cho" #index 6627
# model.reactions.biomass_cho_producing.upper_bound = 0
# model.reactions.biomass_cho_producing.lower_bound = 0

model.objective = "biomass_cho_prod" 
model.reactions.biomass_cho.upper_bound = 0
model.reactions.biomass_cho.lower_bound = 0

# # Specific for CHOmpact
# model.objective = "F90" #index 6618

# model.reactions.F129.lower_bound = 0.6*0.9
# model.reactions.F129.upper_bound = 0.6*1.1
# model.reactions.F130.lower_bound = 0.05*0.9
# model.reactions.F130.upper_bound = 0.05*1.1
# model.reactions.F131.lower_bound = 0.11*0.9
# model.reactions.F131.upper_bound = 0.11*1.1
# model.reactions.F132.lower_bound = 0.45*0.9
# model.reactions.F132.upper_bound = 0.45*1.1
# model.reactions.F133.lower_bound = 0.46*0.9
# model.reactions.F133.upper_bound = 0.46*1.1
# model.reactions.F134.lower_bound = 1.28*0.9
# model.reactions.F134.upper_bound = 1.28*1.1
# model.reactions.F135.lower_bound = 0.02*0.9
# model.reactions.F135.upper_bound = 0.02*1.1
# model.reactions.F136.lower_bound = 3.41*0.9
# model.reactions.F136.upper_bound = 3.41*1.1
# model.reactions.F137.lower_bound = 0.84*0.9
# model.reactions.F137.upper_bound = 0.84*1.1
# model.reactions.F138.lower_bound = 4.1*0.9
# model.reactions.F138.upper_bound = 4.1*1.1
# model.reactions.F139.lower_bound = 0.19*0.9
# model.reactions.F139.upper_bound = 0.19*1.1
# model.reactions.F140.lower_bound = 0.09*0.9
# model.reactions.F140.upper_bound = 0.09*1.1
# # data from the compact paper SUpplementary Data for additional reaction constrains
# # F_129	0.6	0.6	0.6	0.46	0.46
# # F_130	0.05	0.05	0.05	0.02	0.02
# # F_131	0.11	0.11	0.11	0.04	0.04
# # F_132	0.45	0.45	0.45	0.24	0.24
# # F_133	0.46	0.46	0.46	0.46	0.46
# # F_134	1.28	1.28	1.28	1.56	1.56
# # F_135	0.02	0.02	0.02	2.08	2.08
# # F_136	3.41	3.41	3.41	1.84	1.84
# # F_137	0.84	0.84	0.84	0.41	0.41
# # F_138	4.1	4.1	4.1	5.14	5.14
# # F_139	0.19	0.19	0.19	0.22	0.22
# # F_140	0.09	0.09	0.09	0.08	0.08


# %%
# startTime = datetime.now()
# N = 10

# mus = {}
# sets = set(rates.Window)

# for s in sets:
#     mus[s] = np.zeros(N)


# for s in sets:
#     # reset default bounds on all reactions
#     for ex in model.reactions:
#         ex.upper_bound = 1000
#         if ex.reversibility or ex.id in taken_up:
#             ex.lower_bound = -1000
#         else:
#             ex.lower_bound = 0
#     # uptake and secretion rates for one strain
#     print(s)
#     one_set = rates[rates.Window == s]
    
#     n = 0
#     while n < N:
#         print(n)
#         with model:
#             for idx, row in one_set.iterrows():
#                 uptake = row.AA_meta
#                 qp = row.Rate
#                 err = row.SD
#                 ID = uptake_names[uptake]#[2:-1]
#                 r = model.reactions.get_by_id(ID)
                
#                 # sample LB and UB
#                 picked1 = random.uniform(qp - err, qp + err)
#                 picked2 = random.uniform(qp - err, qp + err)
#                 # picked1 = qp + err
#                 # picked2 = qp - err
#                 picked = sorted([picked1, picked2])
#                 # Set bounds
#                 r.bounds = (picked[0], picked[1])
#                 print(r.bounds)
                
#             FBA = model.optimize()
    
#             mus[s][n] = FBA.objective_value
#             n += 1
#             print(FBA.status)

# # Print how long the script ran
# print(datetime.now() - startTime)

# %%
startTime = datetime.now()
N = 10

# Initialize dictionaries
mus = {}
reaction_data = {}  # To store LB, UB, and flux for all reactions

# Get unique sets of Experiments and Windows
experiments = set(rates.Experiment)
windows = set(rates.Window)

# Initialize the mus and reaction_data dictionaries for experiments and windows
for ex in experiments:
    mus[ex] = {}
    reaction_data[ex] = {}
    for w in windows:
        mus[ex][w] = np.zeros(N)
        reaction_data[ex][w] = [{} for _ in range(N)]  # Initialize a list of empty dictionaries for each iteration

# Iterate over each experiment
for ex in experiments:
    # Iterate over each window within the current experiment
    for w in windows:
        # # Reset default bounds on all reactions for each experiment and window
        # for reaction in model.reactions:
        #     reaction.upper_bound = 1000
        #     if reaction.reversibility or reaction.id in taken_up:
        #         reaction.lower_bound = -1000
        #     else:
        #         reaction.lower_bound = 0

        # Select data for the current experiment and window
        one_set = rates[(rates.Experiment == ex) & (rates.Window == w)]

        print(f"Experiment: {ex}, Window: {w}")

        # Iterate N times for this experiment and window combination
        n = 0
        while n < N:
            print(f"Iteration {n} for Experiment {ex} and Window {w}")
            with model:
                # Initialize dictionary for storing bounds and fluxes for all reactions
                for reaction in model.reactions:
                    reaction_data[ex][w][n][reaction.id] = {
                        'LB': reaction.lower_bound,
                        'UB': reaction.upper_bound,
                        'flux': None  # Placeholder for flux after optimization
                    }

                # Apply uptake and secretion rates for the current strain
                for idx, row in one_set.iterrows():
                    uptake = row.AA_meta
                    qp = row.Rate
                    err = row.SD
                    ID = uptake_names[uptake]
                    r = model.reactions.get_by_id(ID)
                    
                    # Sample lower bound (LB) and upper bound (UB)
                    picked1 = random.uniform(qp - err, qp + err)
                    picked2 = random.uniform(qp - err, qp + err)
                    picked = sorted([picked1, picked2])

                    # Set bounds for the reaction
                    r.bounds = (picked[0], picked[1])
                    print(r.bounds)

                    # Store the updated LB and UB in the reaction_data dictionary
                    reaction_data[ex][w][n][ID]['LB'] = picked[0]
                    reaction_data[ex][w][n][ID]['UB'] = picked[1]

                # Perform Flux Balance Analysis (FBA)
                FBA = model.optimize()

                # Store the FBA objective value in the mus dictionary
                mus[ex][w][n] = FBA.objective_value

                # After optimization, collect all fluxes for each reaction
                for reaction in model.reactions:
                    reaction_data[ex][w][n][reaction.id]['flux'] = reaction.flux

                n += 1
                print(FBA.status)

# Print how long the script ran
print(datetime.now() - startTime)

# %%
print("\nResults (mus):")
for ex in mus:
    for w in mus[ex]:
        print(f"Experiment: {ex}, Window: {w}, mu values: {mus[ex][w]}")

# %%
# # Display reaction data (bounds and fluxes) in the command line
# print("\nReaction Data (LB, UB, and Fluxes):")
# for ex in reaction_data:
#     for w in reaction_data[ex]:
#         for n in range(N):
#             print(f"Experiment: {ex}, Window: {w}, Iteration: {n}")
#             for reaction_id, data in reaction_data[ex][w][n].items():
#                 print(f"  Reaction: {reaction_id}, LB: {data['LB']}, UB: {data['UB']}, Flux: {data['flux']}")
# %%
# Convert mus dictionary to a pandas DataFrame
mus_list = []
for ex in mus:
    for w in mus[ex]:
        for n, mu_value in enumerate(mus[ex][w]):
            mus_list.append({
                'Experiment': ex,
                'Window': w,
                'Iteration': n,
                'mu': mu_value
            })
mus_df = pd.DataFrame(mus_list)

# Convert reaction_data dictionary to a pandas DataFrame
reaction_data_list = []
for ex in reaction_data:
    for w in reaction_data[ex]:
        for n in range(N):
            for reaction_id, data in reaction_data[ex][w][n].items():
                reaction_data_list.append({
                    'Experiment': ex,
                    'Window': w,
                    'Iteration': n,
                    'Reaction': reaction_id,
                    'LB': data['LB'],
                    'UB': data['UB'],
                    'Flux': data['flux']
                })
reaction_data_df = pd.DataFrame(reaction_data_list)

# Save DataFrames to CSV files
mus_df.to_csv('fba_results/mus_results_icho2441.csv', index=False)
reaction_data_df.to_csv('fba_results/reaction_data_results_icho2441.csv', index=False)

#%%
# Calculate mean and standard deviation of mus for each phase
mean_mus = {}
std_mus = {}
for window, value in mus.items():
    mean_mu = np.mean(value)
    std_mu = np.std(value)
    mean_mus[window] = mean_mu
    std_mus[window] = std_mu

# Convert dictionaries to DataFrame
df = pd.DataFrame({
    'Window': mean_mus.keys(),
    'Mean': mean_mus.values(),
    'Std': std_mus.values()
})

# Sort the DataFrame by the 'Window' column
df = df.sort_values(by='Window')

# Plot the bar plot with error bars
plt.figure(figsize=(10, 6))
plt.bar(df['Window'], df['Mean'], yerr=df['Std'], capsize=5, color='green', ecolor='black')
plt.xlabel('Window')
plt.ylabel('Mean $mu$')
plt.title('Mean $mu$ for Each Window with Standard Deviation')
plt.xticks(rotation=45)
plt.tight_layout()

# Show the plot
plt.show()
