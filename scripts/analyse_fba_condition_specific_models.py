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
# import random
# warnings - Ignore warnings about infeasible solutions
# import warnings
# datetime - Display the time it took for the analysis to run
from datetime import datetime
# pandas - Work with dataframes
import pandas as pd
# numpy - Working with arrays
import numpy as np
# Matplotlib - For plotting figures
# import matplotlib.pyplot as plt

# Third-party library imports
# cobrapy- Run FBA
import cobra

# Local application/library-specific imports
# escher plots
from escher.plots import Builder
# %% Set directory

os.system("pwd")
os.chdir('c:\\Users\\b1095820\\Documents\\DGTX\\BR24_E13-20')
current_wd = os.getcwd()

# %% load cho model

models = {
    "iCHO1766": "iCHOv1_final.xml",
    # "iCHO2441": "iCHO2441.xml",
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
model = model_orig

# %% load aa and metabolite data
rates = pd.read_csv("data/aa_rates_reordered_data2.csv")

# Define the specific values to be removed
# values_to_remove = ['Growth_rate', 'Titer']  # Replace these with the actual values you want to remove
values_to_remove = ['Titer']
# Modify the existing DataFrame in place
rates.drop(rates[rates['AA_meta'].isin(values_to_remove)].index, inplace=True)

# # Filter rows where Experiment == 'E13' and Window == 1
# filtered_rates = rates.loc[(rates['Experiment'] == 'E13') & (rates['Window'] == 1)]
#
# # Display the filtered DataFrame
# print(filtered_rates)
#
# rates = filtered_rates

# %% check all names of reactions for which you want to set the bounds
taken_up = [
    "EX_gln_L_e_",
    "EX_cys_L_e_",
    "EX_arg_L_e_",
    "EX_asn_L_e_",
    "EX_asp_L_e_",
    "EX_glc_e_",
    "EX_glu_L_e_",
    "EX_h_e_",
    "EX_h2o_e_",
    "EX_his_L_e_",
    "EX_ile_L_e_",
    "EX_leu_L_e_",
    "EX_lys_L_e_",
    "EX_met_L_e_",
    "EX_o2_e_",
    "EX_phe_L_e_",
    "EX_pi_e_",
    "EX_pro_L_e_",
    "EX_ser_L_e_",
    "EX_thr_L_e_",
    "EX_trp_L_e_",
    "EX_tyr_L_e_",
    "EX_val_L_e_",
    "EX_lnlc_e_",
    "EX_lnlnca_e_",
    "EX_Tyr_ggn_e_"
]
#iCHO1766
uptake_names = dict(Ala="EX_ala_L_e_", NH3="EX_nh4_e_", Arg="EX_arg_L_e_", Asn="EX_asn_L_e_", Asp="EX_asp_L_e_",
                    GLC="EX_glc_e_", Glu="EX_glu_L_e_", Gln="EX_gln_L_e_", Gly="EX_gly_e_", His="EX_his_L_e_",
                    Ile="EX_ile_L_e_", LAC="EX_lac_L_e_", Leu="EX_leu_L_e_", Lys="EX_lys_L_e_", Met="EX_met_L_e_",
                    Phe="EX_phe_L_e_", Pro="EX_pro_L_e_", Ser="EX_ser_L_e_", Thr="EX_thr_L_e_", Trp="EX_trp_L_e_",
                    Tyr="EX_tyr_L_e_", Val="EX_val_L_e_", Growth_rate="biomass_cho_producing" )

# Growth_rate="biomass_cho_producing"
# Titer="DM_igg_g_"
# print(f"Reaction {reaction.id} bounds set to: [{reaction.lower_bound}, {reaction.upper_bound}]")

# Set the objective function
# model.objective = "biomass_cho_producing" #index 6618
model.objective = "DM_igg_g_" #index 6618

# model.reactions.biomass_cho.upper_bound = 0
# model.reactions.biomass_cho.lower_bound = 0

# %%
startTime = datetime.now()
# N = 1 # only one repetition

# Initialize dictionaries
mus_fba = {}
mus_pfba = {}
reaction_data_fba = {}  # To store LB, UB, and flux for all reactions
reaction_data_pfba = {}  # To store LB, UB, and flux for all reactions


# Get unique sets of Experiments and Windows
experiments = set(rates.Experiment)
windows = set(rates.Window)

# Initialize the mus and reaction_data dictionaries for experiments and windows
for ex in experiments:
    mus_fba[ex] = {}
    mus_pfba[ex] = {}
    reaction_data_fba[ex] = {}
    reaction_data_pfba[ex] = {}

    for w in windows:
        # Placeholders for single iteration data

        mus_fba[ex][w] = None
        mus_pfba[ex][w] = None
        reaction_data_fba[ex][w] = {} # A single dictionary for reactions
        reaction_data_pfba[ex][w] = {} # A single dictionary for reactions

# # Iterate over each experiment
# for ex in experiments:
#     # Iterate over each window within the current experiment
#     for w in windows:
        # Reset default bounds on all reactions for each experiment and window
        for reaction in model.reactions:
            reaction.upper_bound = 1000
            if reaction.reversibility or reaction.id in taken_up:
                reaction.lower_bound = -1000
            else:
                reaction.lower_bound = 0

        # Select data for the current experiment and window
        one_set = rates[(rates.Experiment == ex) & (rates.Window == w)]

        # print(f"Experiment: {ex}, Window: {w}")

        # switching off epo production and cho_biomass (for non-producers)
        model.reactions.DM_epo_g_.lower_bound = 0
        model.reactions.DM_epo_g_.upper_bound = 0
        model.reactions.biomass_cho.upper_bound = 0
        model.reactions.biomass_cho.lower_bound = 0

        print(f"Processing Experiment {ex}, Window {w}")

        # # Iterate N times for this experiment and window combination
        # n = 0
        # while n < N:
        #     print(f"Iteration {n} for Experiment {ex} and Window {w}")

        # Perform FBA and pFBA within a context (to avoid modifying the model permanently)
        with model:
                # Initialize dictionary for storing bounds and fluxes for all reactions
                for reaction in model.reactions:
                    reaction_data_fba[ex][w][reaction.id] = {
                        'LB': reaction.lower_bound,
                        'UB': reaction.upper_bound,
                        'flux': None  # Placeholder for flux after optimization
                    }
                    reaction_data_pfba[ex][w][reaction.id] = {
                        'LB': reaction.lower_bound,
                        'UB': reaction.upper_bound,
                        'flux': None  # Placeholder for flux after optimization
                    }

                # Apply uptake and secretion rates for the current strain
                for _, row in one_set.iterrows():
                    uptake = row.AA_meta
                    qp = row.Rate
                    err = row.SD
                    ID = uptake_names[uptake]
                    r = model.reactions.get_by_id(ID)

                    picked1 = qp + err
                    picked2 = qp - err
                    picked = sorted([picked1, picked2])

                    # Set bounds for the reaction
                    r.bounds = (picked[0], picked[1])

                    # Store the updated LB and UB in the reaction_data dictionary
                    reaction_data_fba[ex][w][ID]['LB'] = picked[0]
                    reaction_data_fba[ex][w][ID]['UB'] = picked[1]

                try:
                    # Perform Flux Balance Analysis (FBA)
                    FBA = model.optimize()
                    print("FBA Results:", FBA)

                    # Store the FBA objective value and fluxes
                    mus_fba[ex][w] = FBA.objective_value
                    for reaction in model.reactions:
                        reaction_data_fba[ex][w][reaction.id]['flux'] = reaction.flux

                    # Perform parsimonious FBA (pFBA)
                    pfba_solution = cobra.flux_analysis.pfba(model, fraction_of_optimum=1.0)
                    print("pFBA Results:", pfba_solution)

                    # Store the pFBA objective value in the mus dictionary
                    mus_pfba[ex][w] = pfba_solution.objective_value
                    for reaction in model.reactions:
                        reaction_data_pfba[ex][w][reaction.id]['flux'] = reaction.flux

                except cobra.exceptions.Infeasible:
                    print(f"Infeasible solution encountered in iteration {n} for Experiment {ex} and Window {w}. Skipping to next iteration.")
                    # n += 1
                    continue

                # n += 1

# Print script runtime
print("Runtime:", datetime.now() - startTime)


# %% Escher
# save sbml model to json format
cobra.io.save_json_model(model, "cho_gems/iCHO1766_E13_w1.json")

# Escher
builder = Builder(
    map_name='iCHO1766_E13_w1.map',
    model_name="cho_gems/iCHO1766_E13_w1.json"
)
# builder = Builder(model_json="cho_gems/iCHO1766_E13_w1.json")
builder.reaction_data = FBA.fluxes

# Alternatively, save the map to an HTML file
builder.save_html('iCHO1766_E13_w1_map_with_fluxes.html')

# %%
# Convert mus dictionary (FBA and pFBA objective values) to a pandas DataFrame
mus_list = []
for ex in mus_fba:
    for w in mus_fba[ex]:
        mus_list.append({
            'Experiment': ex,
            'Window': w,
            'FBA_mu': mus_fba[ex][w],  # FBA objective value
            'pFBA_mu': mus_pfba[ex][w]  # pFBA objective value
        })

mus_df = pd.DataFrame(mus_list)  # Create DataFrame for objective values
mus_df.to_csv('fba_results/mus_results_icho1766_FBA_pFBA.csv', index=False)  # Save to CSV
print("Saved FBA and pFBA objective data to CSV.")

# Convert reaction_data dictionaries (FBA and pFBA flux data) to a pandas DataFrame
reaction_data_list = []
for ex in reaction_data_fba:
    for w in reaction_data_fba[ex]:
        for reaction_id, fba_data in reaction_data_fba[ex][w].items():
            # Gather data for both FBA and pFBA
            pfba_data = reaction_data_pfba[ex][w].get(reaction_id, {})
            reaction_data_list.append({
                'Experiment': ex,
                'Window': w,
                'Reaction': reaction_id,
                'FBA_LB': fba_data['LB'],  # FBA lower bound
                'FBA_UB': fba_data['UB'],  # FBA upper bound
                'FBA_Flux': fba_data['flux'],  # FBA flux value
                'pFBA_LB': pfba_data.get('LB'),  # pFBA lower bound (if available)
                'pFBA_UB': pfba_data.get('UB'),  # pFBA upper bound (if available)
                'pFBA_Flux': pfba_data.get('flux')  # pFBA flux value (if available)
            })

reaction_data_df = pd.DataFrame(reaction_data_list)  # Create DataFrame for reaction data
reaction_data_df.to_csv('fba_results/reaction_data_results_icho1766_FBA_pFBA.csv', index=False)  # Save to CSV
print("Saved FBA and pFBA reaction data to CSV.")


# %%
# FVA analysis
startTime = datetime.now()
N = 1 # only one repetition

# Initialize dictionaries
mus = {}
fva_results_dict = {}  # To store LB, UB, and flux for all reactions

# Get unique sets of Experiments and Windows
experiments = set(rates.Experiment)
windows = set(rates.Window)

# Initialize the mus and reaction_data dictionaries for experiments and windows
for ex in experiments:
    mus[ex] = {}
    fva_results_dict[ex] = {}
    for w in windows:
        mus[ex][w] = np.zeros(N)
        fva_results_dict[ex][w] = [None] * N # Initialize a list to store FVA results for each iteration
        # reaction_data[ex][w] = [{} for _ in range(N)]  # Initialize a list of empty dictionaries for each iteration

# Iterate over each experiment
for ex in experiments:
    # Iterate over each window within the current experiment
    for w in windows:
        # Reset default bounds on all reactions for each experiment and window
        for reaction in model.reactions:
            reaction.upper_bound = 1000
            if reaction.reversibility or reaction.id in taken_up:
                reaction.lower_bound = -1000
            else:
                reaction.lower_bound = 0

        # Select data for the current experiment and window
        one_set = rates[(rates.Experiment == ex) & (rates.Window == w)]

        print(f"Experiment: {ex}, Window: {w}")

        # switching off epo production and cho_biomass (for non-producers)
        model.reactions.DM_epo_g_.lower_bound = 0
        model.reactions.DM_epo_g_.upper_bound = 0
        model.reactions.biomass_cho.upper_bound = 0
        model.reactions.biomass_cho.lower_bound = 0

        # Iterate N times for this experiment and window combination
        n = 0
        while n < N:
            print(f"Iteration {n} for Experiment {ex} and Window {w}")
            with model:
                # # Initialize dictionary for storing bounds and fluxes for all reactions
                # for reaction in model.reactions:
                #     reaction_data[ex][w][n][reaction.id] = {
                #         'LB': reaction.lower_bound,
                #         'UB': reaction.upper_bound,
                #         'flux': None  # Placeholder for flux after optimization
                #     }

                # Apply uptake and secretion rates for the current strain
                for idx, row in one_set.iterrows():
                    uptake = row.AA_meta
                    qp = row.Rate
                    err = row.SD
                    ID = uptake_names[uptake]
                    r = model.reactions.get_by_id(ID)

                    picked1 = qp + err
                    picked2 = qp - err
                    picked = sorted([picked1, picked2])

                    # Set bounds for the reaction
                    r.bounds = (picked[0], picked[1])

                    # # Store the updated LB and UB in the reaction_data dictionary
                    # reaction_data[ex][w][n][ID]['LB'] = picked[0]
                    # reaction_data[ex][w][n][ID]['UB'] = picked[1]

                # try:
                #     # Perform Flux Balance Analysis (FBA)
                #     FBA = model.optimize()
                #
                #     # Perform parsimonious FBA
                #     pfba_solution = cobra.flux_analysis.pfba(model, fraction_of_optimum=1.0)
                #
                #     # Store the FBA objective value in the mus dictionary
                #     mus[ex][w][n] = pfba_solution.objective_value
                #
                #     # After optimization, collect all fluxes for each reaction
                #     for reaction in model.reactions:
                #         reaction_data[ex][w][n][reaction.id]['flux'] = reaction.flux
                #
                # except cobra.exceptions.Infeasible:
                #     print(f"Infeasible solution encountered in iteration {n} for Experiment {ex} and Window {w}. Skipping to next iteration.")
                #     n += 1
                #     continue
                # Perform FVA
                try:
                    fva_results = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=1.0)
                    fva_results_dict[ex][w][n] = fva_results  # Store FVA results
                    print(f"FVA completed for Experiment {ex}, Window {w}, Iteration {n}")
                except Exception as e:
                    print(f"Error during FVA for Experiment {ex}, Window {w}, Iteration {n}: {e}")
            n += 1

# Print how long the script ran
print(datetime.now() - startTime)
# %%

# Save FVA results to Excel files
for ex in experiments:
    for w in windows:
        for n in range(N):
            fva_results = fva_results_dict[ex][w][n]
            if fva_results is not None:
                df = pd.DataFrame({
                    'Reaction': [rxn.id for rxn in model.reactions],  # Reaction IDs
                    'Min. Flux': fva_results['minimum'],  # Min. flux values
                    'Max. Flux': fva_results['maximum']  # Max. flux values
                })
                filename = f"FVA_Experiment_{ex}_Window_{w}_Iteration_{n}.xlsx"
                filepath = os.path.join("fba_results", filename)
                df.to_excel(filepath, index=False)
                print(f"Saved FVA results to {filepath}")

# df = pd.DataFrame({
#     'Reaction': [rxn.id for rxn in model.reactions],  # reaction ID
#     'Min. Flux': fva_results['minimum'],  # min. flux values
#     'Max. Flux': fva_results['maximum']  # max. flux values
# })
#
# filename = "CHO_solution_FVA_fluxes.xlsx"
# filepath = os.path.join("fba_results/", filename)
# df.to_excel(filepath, index=False)