library(tidyverse)

measured_fluxes <- read_csv("data/aa_rates_reordered_data2.csv")

color_mapping_experiment <- c(
  "E13" = "#FD8D3C",
  "E14" = "#9E9AC8",
  "E15" = "#F16913",
  "E16" = "#807DBA",
  "E17" = "#D94801",
  "E18" = "#6A51A3",
  "E19" = "#A63603",
  "E20" = "#54278F"
)


# Alanine -----------------------------------------------------------------

measured_fluxes_ala <- measured_fluxes %>%
  filter(AA_meta %in% "Ala") %>%
  mutate(Type = "measured",
         Reaction = "EX_ala_L_e_") %>%
  select(!c(Condition, AA_meta))

predicted_fluxes_ala <- read_csv("fba_results/condition_specific/iCHO1766_biomass_producing/aa_validation/reaction_data_results_icho1766_FBA_pFBA_Ala.csv") %>%
  filter(Reaction %in% "EX_ala_L_e_") %>%
  pivot_longer(cols = c("FBA_Flux", "pFBA_Flux"), 
               names_to = "Type",
               values_to = "Rate") %>%
  select(!c(FBA_LB, FBA_UB, pFBA_LB, pFBA_UB)) %>%
  mutate(SD = NA)

# merge measured and predicted
# Bind predicted and measured by rows for plotting bar plots
mp_fluxes_ala <- measured_fluxes_ala %>%
  rbind(predicted_fluxes_ala) 

# plot


ggplot(mp_fluxes_ala, aes(x = Window, y = Rate, color = Experiment, group = interaction(Window, Experiment))) +
  geom_point(
    aes(shape = Type),
    position = position_dodge(width = 0.2),  
    size = 2
  ) +
  geom_errorbar(
    aes(
      ymin = Rate - SD,  
      ymax = Rate + SD   
    ),
    position = position_dodge(width = 0.2),  
    width = 0.2,  
    linewidth = 0.5 
  ) + 
  scale_color_manual(values = color_mapping_experiment) +
  ggtitle("EX_ala_L_e_") +
  theme_bw() +
  theme(
    legend.position = "top"  
  )

# NH3 -----------------------------------------------------------------

measured_fluxes_ala <- measured_fluxes %>%
  filter(AA_meta %in% "NH3") %>%
  mutate(Type = "measured",
         Reaction = "EX_nh4_e_") %>%
  select(!c(Condition, AA_meta))

predicted_fluxes_ala <- read_csv("fba_results/condition_specific/iCHO1766_biomass_producing/aa_validation/reaction_data_results_icho1766_FBA_pFBA_NH3.csv") %>%
  filter(Reaction %in% "EX_nh4_e_") %>%
  pivot_longer(cols = c("FBA_Flux", "pFBA_Flux"), 
               names_to = "Type",
               values_to = "Rate") %>%
  select(!c(FBA_LB, FBA_UB, pFBA_LB, pFBA_UB)) %>%
  mutate(SD = NA)

# merge measured and predicted
# Bind predicted and measured by rows for plotting bar plots
mp_fluxes_ala <- measured_fluxes_ala %>%
  rbind(predicted_fluxes_ala) 

# plot 

ggplot(mp_fluxes_ala, aes(x = Window, y = Rate, color = Experiment, group = interaction(Window, Experiment))) +
  geom_point(
    aes(shape = Type),
    position = position_dodge(width = 0.2),  
    size = 2
  ) +
  geom_errorbar(
    aes(
      ymin = Rate - SD,  
      ymax = Rate + SD  
    ),
    position = position_dodge(width = 0.2),  
    width = 0.2,  
    linewidth = 0.5  
  ) + 
  scale_color_manual(values = color_mapping_experiment) +
  ggtitle("EX_nh4_e_") +
  theme_bw() +
  theme(
    legend.position = "top" 
  )


# generic -----------------------------------------------------------------
# Create a named list 
metabolite_reaction_pairs <- list(
  Ala = "EX_ala_L_e_",
  NH3 = "EX_nh4_e_",
  Arg = "EX_arg_L_e_",
  Asn = "EX_asn_L_e_",
  Asp = "EX_asp_L_e_",
  GLC = "EX_glc_e_",
  Glu = "EX_glu_L_e_",
  Gln = "EX_gln_L_e_",
  Gly = "EX_gly_e_",
  His = "EX_his_L_e_",
  Ile = "EX_ile_L_e_",
  LAC = "EX_lac_L_e_",
  Leu = "EX_leu_L_e_",
  Lys = "EX_lys_L_e_",
  Met = "EX_met_L_e_",
  Phe = "EX_phe_L_e_",
  Pro = "EX_pro_L_e_",
  Ser = "EX_ser_L_e_",
  Thr = "EX_thr_L_e_",
  Trp = "EX_trp_L_e_",
  Tyr = "EX_tyr_L_e_",
  Val = "EX_val_L_e_"
)

# Loop over the named list
for (metabolite in names(metabolite_reaction_pairs)) {
  
  # Get the corresponding reaction for the metabolite
  reaction <- metabolite_reaction_pairs[[metabolite]]
  
  # Process and plot 
  measured_fluxes_metabolite <- measured_fluxes %>%
    filter(AA_meta == metabolite) %>%
    mutate(Type = "measured",
           Reaction = reaction) %>%
    select(!c(Condition, AA_meta))
  
  predicted_fluxes_metabolite <- read_csv(paste0("fba_results/condition_specific/iCHO1766_biomass_producing/aa_validation/reaction_data_results_icho1766_FBA_pFBA_",
                                                 metabolite,
                                                 ".csv")) %>%
    filter(Reaction == reaction) %>%
    pivot_longer(cols = c("FBA_Flux", "pFBA_Flux"), 
                 names_to = "Type",
                 values_to = "Rate") %>%
    select(!c(FBA_LB, FBA_UB, pFBA_LB, pFBA_UB)) %>%
    mutate(SD = NA)
  
  mp_fluxes_metabolite <- bind_rows(measured_fluxes_metabolite, predicted_fluxes_metabolite)
  
  p <- ggplot(mp_fluxes_metabolite, aes(x = Window, y = Rate, color = Experiment, group = interaction(Window, Experiment))) +
    geom_point(aes(shape = Type), position = position_dodge(width = 0.4), size = 2) +
    geom_errorbar(aes(ymin = Rate - SD, ymax = Rate + SD), position = position_dodge(width = 0.4), width = 0.2, linewidth = 0.5) + 
    scale_color_manual(values = color_mapping_experiment) +
    ggtitle(paste0("Reaction ",reaction, "\nMetabolite ",metabolite)) +
    theme_bw() +
    theme(legend.position = "top")
  
  # plot(p)
  
  ggsave(filename = paste0("figures/fba/pfba_fba_fva/validation_fluxes/",metabolite,".png"),
         plot = p,
         height = 120,
         width = 200,
         units = "mm",
         dpi = 600,
         bg = "white") 
}


