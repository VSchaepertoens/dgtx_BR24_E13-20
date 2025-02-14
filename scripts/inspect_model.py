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
import warnings
# datetime - Display the time it took for the analysis to run
from datetime import datetime
# pandas - Work with dataframes
import pandas as pd
# numpy - Working with arrays
import numpy as np

# Third-party library imports
# cobrapy- Run FBA
import cobra

# Local application/library-specific imports
# escher plots
from escher.plots import Builder

# %% Set directory

os.system("pwd")
current_wd = os.getcwd()
os.chdir(current_wd)

# %% load cho model

models = {
    # "iCHO1766": "iCHOv1_final.xml",
    # "iCHO2441": "iCHO2441.xml",
    # "CHO-K1": "iCHOv1_K1_final.xml",
    # "CHOmpact": "CHOmpact_generic_producing_CK.json",
    "CHOmpact": "CHOmpact_generic_producing_Published_Fixed.json",
    # "CHOmpact": "CHOsmallmodel.json"
    # "CHOmpact_small": "CHOsmallmodel_activity4.json"
    # "K1par-0mMCD": "iCHO_K1par-0mMCD.xml"
    }

for model_name, model_file in models.items():
    print(model_name)
    MODEL_PATH = "cho_gems/" + model_file
    _, file_extension = os.path.splitext(MODEL_PATH)

    if file_extension == '.xml':
        model = cobra.io.read_sbml_model(MODEL_PATH)
    elif file_extension == '.json':
        model = cobra.io.load_json_model(MODEL_PATH)
    else:
        # Handle unsupported file types or other extensions
        print(f"Unsupported file extension: {file_extension}")
        continue
    print(f"Loaded {model_name} model successfully!")


# %%  Explore the model by exporting it to an Excel file

# Create a list of reactions and their coefficients
reaction_list = []
for reaction in model.reactions:
    reaction_data = {
        "Reaction": reaction.id,
        "Equation": reaction.reaction,
        "Objective Coefficient": reaction.objective_coefficient,
        "Reversibility": reaction.reversibility, #the reversibility attribute of a reaction object represents whether the reaction is reversible. This attribute is a boolean value (True or False), where True indicates that the reaction is reversible, meaning it can proceed in both the forward and reverse directions, while False indicates that the reaction is irreversible, meaning it can only proceed in one direction.
        "Lower bound": reaction.lower_bound,
        "Upper bound": reaction.upper_bound
    }
    reaction_list.append(reaction_data)
# Create a DataFrame from the list
reaction_df = pd.DataFrame(reaction_list)

# Export the model to excel format
DIRECTORYPATH = "C:/Users/b1095820/Documents/DGTX/BR24_E13-20/cho_gems"
FILENAME = "CHOmpact_generic_producing_Published_Fixed.xlsx"
FILEPATH = os.path.join(DIRECTORYPATH, FILENAME)
reaction_df.to_excel(FILEPATH, index=True)

# %% check all names of reactions for which you want to set the bounds

uptake_names = np.load("cho_gems/CHO_uptake_names.npy", allow_pickle=True).item()
uptake_names = {
    "Alanine": "EX_ala_L_e_",
    "Ammonia": "EX_nh4_e_",
    "Arginine": "EX_arg_L_e_",
    "Asparagine": "EX_asn_L_e_",
    "Aspartic_acid": "EX_asp_L_e_",
    "Cysteine": "EX_cys_L_e_",
    "Glucose": "EX_glc_e_",
    "Glutamic_acid": "EX_glu_L_e_",
    "Glutamine": "EX_gln_L_e_",
    "Glycine": "EX_gly_L_e_",
    "Histidine": "EX_his_L_e_",
    "Isoleucine": "EX_ile_L_e_",
    "Lactate": "EX_lac_L_e_",
    "Leucine": "EX_leu_L_e_",
    "Lysine": "EX_lys_L_e_",
    "Methionine": "EX_met_L_e_",
    "Phenylalanine": "EX_phe_L_e_",
    "Proline": "EX_pro_L_e_",
    "Serine": "EX_ser_L_e_",
    "Threonine": "EX_thr_L_e_",
    "Tryptophan": "EX_trp_L_e_",
    "Tyrosine": "EX_tyr_L_e_",
    "Valine": "EX_val_L_e_"
}

uptake_names = {
    "Ala": "EX_ala_L_e_",
    #"Ammonia": "EX_nh4_e_",
    "Arg": "EX_arg_L_e_",
    "Asn": "EX_asn_L_e_",
    "Asp": "EX_asp_L_e_",
    #"Cysteine": "EX_cys_L_e_",
    #"Glucose": "EX_glc_e_",
    "Glu": "EX_glu_L_e_",
    "Gln": "EX_gln_L_e_",
    "Gly": "EX_gly_L_e_",
    "His": "EX_his_L_e_",
    "Ile": "EX_ile_L_e_",
    #"Lactate": "EX_lac_L_e_",
    "Leu": "EX_leu_L_e_",
    "Lys": "EX_lys_L_e_",
    "Met": "EX_met_L_e_",
    "Phe": "EX_phe_L_e_",
    "Pro": "EX_pro_L_e_",
    "Ser": "EX_ser_L_e_",
    "Thr": "EX_thr_L_e_",
    "Trp": "EX_trp_L_e_",
    "Tyr": "EX_tyr_L_e_",
    "Val": "EX_val_L_e_"
}

# Turn off igg and epo production
model.reactions.DM_igg_g_.lower_bound = 0
model.reactions.DM_igg_g_.upper_bound = 0
model.reactions.DM_epo_g_.lower_bound = 0
model.reactions.DM_epo_g_.upper_bound = 0

# Set the objective function
#if strain in producers:
model.objective = "biomass_cho_producing" #index 6618
model.reactions.biomass_cho.upper_bound = 0
model.reactions.biomass_cho.lower_bound = 0
#else:
model.objective = "biomass_cho" #index 6627
model.reactions.biomass_cho_producing.upper_bound = 0
model.reactions.biomass_cho_producing.lower_bound = 0











