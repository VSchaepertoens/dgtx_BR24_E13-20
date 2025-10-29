# %% Set-up: Import necessary libraries and functions

# Standard libraries
# os -  Working with file or folder directories
import os
# pandas - Work with dataframes
import pandas as pd
# numpy - Working with arrays
import numpy as np
# Matplotlib - For plotting figures
# import matplotlib.pyplot as plt

# Third-party library imports
# cobrapy- Run FBA
import cobra
# escher plots
import escher
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

# %% save sbml model to json format
cobra.io.save_json_model(model_orig, "cho_gems/iCHO1766.json")

# %% Escher
builder = Builder(model_json="cho_gems/iCHO1766.json")
builder.reaction_data = global_result.fluxes

# %%
escher.list_available_maps()