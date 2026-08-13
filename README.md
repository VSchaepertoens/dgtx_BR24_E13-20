# DGTX Project NISTCHO Characterization 
# Analysis of mAb subunits and N-glycans on intact mAb level
## 📄 Overview
(Abstract of the paper)

## 🗂 Repository Structure

- `/scripts`: This directory contains all R scripts used for the analysis of CHO fed-batch mass spectrometry data, including N-glycans quantification, correction for hexosylation bias and the quantification of subunits from the UV chromatograms. 
Below is an overview of the purpose and outputs of each script:
  - mAb subunits analysis:
    - `/00_subunits_plot_heatmaps.R` : Collects all data for subunits (LC, LC2, intact mAb) and plots lineplot. Saves data as charrun_E13-E20_subunit_V0x_YYYYMMDD_VS.RData. 
    - `/01_analyse_linear_models_subunit.R`: Uses linear models to test for differences in subunit abundances between conditions and timepoints. 
  
  - N-glycans abundance analysis:
    - `/01_analyse_all_files.R`: Using the package fragquaxi [fragquaxi](https://github.com/cdl-biosimilars/fragquaxi), quantifies the abundance of N-glycans in the input mzml files. For quantification of glycation in PNGaseF-digested mzml files, change line 29 to "pngase". 
    - `/02_plot_abundances.R`: Collects N-glycan abundances from all files & plots as barplots for first glimpse of the data. Change line 7 to "pngase" to visualise & assemble glycation data from PNGaseF-digested samples. Produces intermediary files abundance_data_none.RData and abundance_data_pngase.RData 
    - `/03_prepare_data_cafog.R`: Assembles all data required for the CAFOG analysis in the folder analysis/cafog. 
      - `/subprocess_cafog.ipynb`: Uses the hexose bias correction algorithm [cafog](https://github.com/cdl-biosimilars/cafog) to correct N-glycan abundances for hexosylation bias.
    - `/04_plot_abundance_cafog_corrected.R`: Plots the corrected N-glycan abundances for all experiments and individual biological replicates and saves all corrected N-glycan abundances as corr_abundance_data.RData and charrun_E13-E20_CQA_V0x_YYYYMMDD_VS.RData.
    - `/05_analyse_splineomics.R`: Uses linear models to test for differences in N-glycan abundances between conditions and timepoints. 

    - Calculate galactosylation, fucosylation and glycation indices:
      - `/galactosylation_index.R`:
      - `/fucosylation_index.R`:
      - `/glycation_index.R`:
    - `/06_analyse_linear_models_index.R`: Uses linear models to test for differences in indices between conditions and timepoints. 



**Note: Copy [subprocess_cafog.ipynb](subprocess_cafog.ipynb) to cafog folder to run directly from the source & base_folder in .ipynb must be changed to match the directory of `analysis/cafog`.

## 📦 Data Access (Zenodo)

The full raw dataset for this study is archived and publicly available on Zenodo:
  
  Zenodo DOI:[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21643681.svg)](https://doi.org/10.5281/zenodo.21643681  [Add to Citavi project by DOI] )

The processed data provided in this repository were generated entirely from the Zenodo data using the scripts in /scripts.

## 📝 Citation
If you use the code or data, please cite:
  
  ## 🤝 Contact
  
  #### Mass spectrometry analysis & Computational analysis
  **Veronika Schäpertöns**  
  University of Salzburg
📧 **veronika.schaepertoens@plus.ac.at**  
  🔗 **GitHub:** [@VSchaepertoens](https://github.com/VSchaepertoens) 


