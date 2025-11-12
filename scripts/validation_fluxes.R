library(tidyverse)
library(ggpubr)

measured_fluxes <- read_csv("data/rates_05032025/rates_mM_gDCW_h_new.csv")

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

color_mapping_condition <- c(
  "Constant" = "#A63603",
  "Temp. shifted" = "#54278F"
)

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


# plot 1 ------------------------------------------------------------------

# Initialize an empty list to store data frames
mp_fluxes_list <- list()

# Loop over the named list
for (metabolite in names(metabolite_reaction_pairs)) {
  
  # Get the corresponding reaction for the metabolite
  reaction <- metabolite_reaction_pairs[[metabolite]]
  
  # Process and plot 
  measured_fluxes_metabolite <- measured_fluxes %>%
    filter(Metabolite == metabolite) %>%
    mutate(Type = "measured",
           Reaction = reaction,
           Rate = rate_mM_gDCW_h,
           SD = SD_mM_gDCW_h) %>%
    select(!c(Metabolite))
  
  predicted_fluxes_metabolite <- read_csv(paste0("fba_results/window_specific_5/iCHO1766_biomass_cho_producing/validation_uptake_rates/reaction_data_results_icho1766_FBA_pFBA_",
                                                 metabolite,
                                                 ".csv")) %>%
    filter(Reaction == reaction) %>%
    pivot_longer(cols = c("FBA_Flux", "pFBA_Flux"), 
                 names_to = "Type",
                 values_to = "Rate") %>%
    select(!c(FBA_LB, FBA_UB, pFBA_LB, pFBA_UB)) %>%
    mutate(SD = NA)
  
  mp_fluxes_metabolite <- bind_rows(measured_fluxes_metabolite, predicted_fluxes_metabolite)
  
  p <- ggplot(mp_fluxes_metabolite, aes(x = Window, y = Rate, color = Condition, group = interaction(Window, Condition))) +
    geom_point(aes(shape = Type), position = position_dodge(width = 0.4), size = 2) +
    geom_errorbar(aes(ymin = Rate - SD, ymax = Rate + SD), position = position_dodge(width = 0.4), width = 0.2, linewidth = 0.5) + 
    scale_color_manual(values = color_mapping_condition) +
    ggtitle(paste0("Reaction ",reaction, "\nMetabolite ",metabolite)) +
    theme_bw() +
    theme(legend.position = "top")
  
  plot(p)
  
  # ggsave(filename = paste0("figures/fba/pfba_fba_fva/validation_fluxes/",metabolite,".png"),
  #        plot = p,
  #        height = 120,
  #        width = 200,
  #        units = "mm",
  #        dpi = 600,
  #        bg = "white") 
  
  # Store in list
  mp_fluxes_list[[metabolite]] <- mp_fluxes_metabolite
}


# Combine all data frames into a single one
mp_fluxes_all <- bind_rows(mp_fluxes_list, .id = "Metabolite")

mp_fluxes_all_wider <- mp_fluxes_all %>%
  select(Metabolite, Condition, Window, Rate, SD, Type, Reaction) %>%
  pivot_wider(names_from = Type,
              values_from = c("Rate", "SD")) %>%
  mutate(Window = as.character(Window))

# plot 2: dotplot per reaction:metabolite --------------------------------
# mp_fluxes_metabolite_wider <- mp_fluxes_metabolite %>%
#   pivot_wider(names_from = Type,
#               values_from = c("Rate", "SD")) %>%
#   mutate(Window = as.character(Window))


# # Calculate the overall range for both axes, considering SD and adding a buffer
# limits_x <- range(
#   c(mp_fluxes_all_wider$Rate_measured - mp_fluxes_all_wider$SD_measured,
#     mp_fluxes_all_wider$Rate_measured + mp_fluxes_all_wider$SD_measured,
#     mp_fluxes_all_wider$Rate_pFBA_Flux), na.rm = TRUE
# )


# Add a buffer (e.g., 5%) to the range
# buffer_x <- 0.05 * (limits_x[2] - limits_x[1])

# Use the same limits for y-axis (since both axes should have the same range)
# limits_y <- limits_x  # Set y-axis limits same as x-axis limits
# Define abline parameters
abline_params <- list(
  geom_abline(intercept = 0, slope = 1, linewidth = 0.5),
  geom_abline(intercept = 0, slope = 0.75, linetype = 2, color = 'grey40'),
  geom_abline(intercept = 0, slope = 1.25, linetype = 2, color = 'grey40')
)

fba_comp <- ggplot(mp_fluxes_all_wider, aes(x = Rate_measured, y = Rate_FBA_Flux, color = Condition)) +
  abline_params +
  geom_point(aes(shape = Window)) +
  geom_errorbar(
    aes(xmin = Rate_measured - SD_measured,
        xmax = Rate_measured + SD_measured,
        group = Condition),
    position = position_dodge(.9),
    width = .00025,
    linewidth = .25
  ) +
  scale_color_manual(values = color_mapping_condition) +
  # Set the same limits for both axes with buffer
  # scale_x_continuous(limits = c(limits_x[1] - buffer_x, limits_x[2] + buffer_x)) +
  # scale_y_continuous(limits = c(limits_y[1] - buffer_x, limits_y[2] + buffer_x)) +
  labs(title = "iCHO1766 Predicted vs Measured Rates") +
  theme_minimal() +
  # Enhance readability and aesthetics 
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "top"
  ) +
  # Fix legend title for 'Window'
  guides(shape = guide_legend(title = "Window"))

plot(fba_comp)

pfba_comp <- ggplot(mp_fluxes_all_wider, aes(x = Rate_measured, y = Rate_pFBA_Flux, color = Condition)) +
  abline_params +
  geom_point(aes(shape = Window)) +
  geom_errorbar(
    aes(xmin = Rate_measured - SD_measured,
        xmax = Rate_measured + SD_measured,
        group = Condition),
    position = position_dodge(.9),
    width = .00025,
    linewidth = .25
  ) +
  scale_color_manual(values = color_mapping_condition) +
  # Set the same limits for both axes with buffer
  # scale_x_continuous(limits = c(limits_x[1] - buffer_x, limits_x[2] + buffer_x)) +
  # scale_y_continuous(limits = c(limits_y[1] - buffer_x, limits_y[2] + buffer_x)) +
  labs(title = "iCHO1766 Predicted vs Measured Rates") +
  theme_minimal() +
  # Enhance readability and aesthetics 
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "top"
  ) +
  # Fix legend title for 'Window'
  guides(shape = guide_legend(title = "Window"))


ggarrange(fba_comp, 
          pfba_comp,
          ncol = 2,
          common.legend = TRUE)



# plot per metabolite and calculate R2 ------------------------------------

mp_fluxes_metabolite_wider <- mp_fluxes_all_wider %>%
  filter(Metabolite == "Ala")


fba_comp <- ggplot(mp_fluxes_metabolite_wider, aes(x = Rate_measured, y = Rate_FBA_Flux, color = Condition)) +
  abline_params +
  geom_point(aes(shape = Window)) +
  geom_errorbar(
    aes(xmin = Rate_measured - SD_measured,
        xmax = Rate_measured + SD_measured,
        group = Condition),
    position = position_dodge(.9),
    width = .00025,
    linewidth = .25
  ) +
  scale_color_manual(values = color_mapping_condition) +
  # Set the same limits for both axes with buffer
  # scale_x_continuous(limits = c(limits_x[1] - buffer_x, limits_x[2] + buffer_x)) +
  # scale_y_continuous(limits = c(limits_y[1] - buffer_x, limits_y[2] + buffer_x)) +
  labs(title = "iCHO1766 Predicted vs Measured Rates") +
  theme_minimal() +
  # Enhance readability and aesthetics 
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "top"
  ) +
  # Fix legend title for 'Window'
  guides(shape = guide_legend(title = "Window"))

plot(fba_comp)

pfba_comp <- ggplot(mp_fluxes_metabolite_wider, aes(x = Rate_measured, y = Rate_pFBA_Flux, color = Condition)) +
  abline_params +
  geom_point(aes(shape = Window)) +
  geom_errorbar(
    aes(xmin = Rate_measured - SD_measured,
        xmax = Rate_measured + SD_measured,
        group = Condition),
    position = position_dodge(.9),
    width = .00025,
    linewidth = .25
  ) +
  scale_color_manual(values = color_mapping_condition) +
  # Set the same limits for both axes with buffer
  # scale_x_continuous(limits = c(-0.03, 0.04)) +
  # scale_y_continuous(limits = c(-0.025, 0.025)) +
  labs(title = "iCHO1766 Predicted vs Measured Rates") +
  theme_minimal() +
  # Enhance readability and aesthetics 
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "top"
  ) +
  # Fix legend title for 'Window'
  guides(shape = guide_legend(title = "Window"))


plot(pfba_comp)










