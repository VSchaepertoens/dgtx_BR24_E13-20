library(tidyverse)

measured_fluxes <- read_csv("data/aa_rates_reordered_data2.csv")

# GROWTH rates ------------------------------------------
## MEASURED ##
m_growth_rates <- measured_fluxes %>%
  filter(AA_meta == "Growth_rate") %>%
  mutate(Type = "measured") %>%
  select(!Condition)


p_growth_rates <- read_csv("fba_results/condition_specific/iCHO1766_biomass_producing/mus_results_icho1766_FBA_pFBA.csv") %>%
  mutate(AA_meta = "Growth_rate",
         Type = "iCHO1766") %>%
  pivot_longer(cols = c("FBA_mu", "pFBA_mu"), names_to = "method_type") %>%
  mutate(Type = paste(Type, method_type, sep = "_"),
         SD = 0) %>%
  mutate(Type = gsub("_mu", "", Type)) %>%
  rename(Rate = value) %>%
  select(!method_type)

# # merge predicted with measured dfs
mp_growth_rates <- m_growth_rates %>%
  rbind(p_growth_rates) 
  
# 
# 
# save(mp_growth_rates, file = "fba_results/both_conditions/case_2_1_mus.RData")
# load("fba_results/both_conditions/case_2_1_mus.RData")

# Ensure that rows align by joining based on common columns for dotplot
mp_growth_rates_joined <- m_growth_rates  %>%
  full_join(p_growth_rates %>% 
              filter(Type != "iCHO1766_pFBA"), 
            by = c("Experiment", "Window", "AA_meta")) %>%
  # full_join(p_growth_rates_icho2441, by = c("Experiment", "Window", "Condition", "AA_meta")) %>%
  rename(measured = Rate.x,
         measured_error = SD.x,
         icho1766 = Rate.y,
         # icho1766_error = SD.y,
         # icho2441 = Rate,
         # icho2441_error = SD
         ) %>%
  select(!c("Type.x", "Type.y"))


# TITER rates ----------------------------
## MEASURED ##
  m_titer_rates <- measured_fluxes %>%
    filter(AA_meta == "Titer") %>%
    mutate(Type = "measured")%>%
    select(!Condition)

## PREDICTED iCHO1766 ##
p_titer_rates <- read_csv("fba_results/condition_specific/iCHO1766_igg/mus_results_icho1766_FBA_pFBA.csv") %>%
  mutate(AA_meta = "Titer",
         Type = "iCHO1766") %>%
  pivot_longer(cols = c("FBA_mu", "pFBA_mu"), names_to = "method_type") %>%
  mutate(Type = paste(Type, method_type, sep = "_"),
         SD = 0) %>%
  mutate(Type = gsub("_mu", "", Type)) %>%
  rename(Rate = value) %>%
  select(!method_type)

mp_titer_rates <- m_titer_rates %>%
  rbind(p_titer_rates) 

  # Ensure that rows align by joining based on common columns for dotplot
  mp_titer_rates_joined <- m_titer_rates %>%
    full_join(p_titer_rates %>% 
                filter(Type != "iCHO1766_FBA"), 
              by = c("Experiment", "Window", "AA_meta")) %>%
    # full_join(p_titer_rates_icho2441, by = c("Experiment", "Window", "Condition", "AA_meta")) %>%
    rename(measured = Rate.x,
           measured_error = SD.x,
           icho1766 = Rate.y,
           # icho1766_error = SD.y,
           # icho2441 = Rate,
           # icho2441_error = SD
           ) %>%
    select(!c("Type.x", "Type.y"))
  
# function to plot data barplots ---------------------------------------------------
plot_rates <- function(data_to_plot,
                       y_axis_label,
                       y_axis_limits,
                       plot_growth_rate = FALSE){
   ggplot(data_to_plot, aes(x = Experiment)) +
     geom_col(
       aes(y = Rate,
           fill = Type),
       position = position_dodge(width = 0.9)
       ) +
     geom_errorbar(
       aes(ymin = Rate - SD,
           ymax = Rate + SD,
           group = Type),
       position = position_dodge(.9),
       width = .5,
       linewidth = .25
       ) +
     facet_wrap(~Window, ncol =  5) +
    theme_light() +
    theme(text = element_text(size = 8,
                              face = "plain",
                              family = "sans")
          ) +
     if (plot_growth_rate) {
       scale_y_continuous(name = y_axis_label,
                          limits = y_axis_limits)
     }else{
       scale_y_continuous(name = y_axis_label,
                          labels = function(x) x * 1e5,
                          limits = y_axis_limits)
     }
}
  
# growth rate
y_label_growth_rate <- expression(paste("growth rate (",h^-1,")"))
growth_rate_limits = c(-0.015, 0.045)

plot_rates(data_to_plot = mp_growth_rates %>% filter(Type != "iCHO1766_pFBA"),
           y_axis_label = y_label_growth_rate,
           y_axis_limits = growth_rate_limits,
           plot_growth_rate = TRUE)

ggsave(filename = "figures/growth_rates_case2_1.png",
       height = 70,
       width = 250,
       units = "mm",
       dpi = 600)


# titer
y_label_titer_rate <- expression(paste("specific antibody productivity (mM g DC",W^-1,h^-1,") x",~ 10^-5))
titer_rate_limits = c(-1e-5, 8e-5)

plot_rates(data_to_plot = mp_titer_rates %>% filter(Type != "iCHO1766_pFBA"),
           y_axis_label = y_label_titer_rate,
           y_axis_limits = titer_rate_limits,
           plot_growth_rate = FALSE)

ggsave(filename = "figures/titer_rates_case2_2.png",
       height = 70,
       width = 250,
       units = "mm",
       dpi = 600)


# dotplots exp_growth vs pred_growth plots -------------------------------------------
# Plotting measured vs predicted rates
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

#growth rate
ggplot(mp_growth_rates_joined %>% filter(Window %in% c("1", "2", "3")), aes(x = measured, y = icho1766, color = Experiment)) +
  geom_abline(intercept = 0, slope = 1, linewidth = 1) +
  geom_abline(intercept = 0, slope = 0.75, linetype = 2, color = 'grey40') +
  geom_abline(intercept = 0, slope = 1.25, linetype = 2, color = 'grey40') +
  geom_point(alpha = 0.7) +
  geom_errorbar(
    aes(xmin = measured - measured_error,
        xmax = measured + measured_error,
        group = Experiment),
    position = position_dodge(.9),
    width = .001,
    linewidth = .25
  ) +
  # geom_errorbar(
  #   aes(ymin = icho1766 - icho1766_error,
  #       ymax = icho1766 + icho1766_error,
  #       group = Experiment),
  #   position = position_dodge(.9),
  #   width = .001,
  #   linewidth = .25
  # ) +
  scale_color_manual(values = color_mapping_experiment) +
  scale_y_continuous(limits = c(-0.005, 0.04)) +
  scale_x_continuous(limits = c(-0.005, 0.04)) +
  coord_fixed() +
  labs(x = "Measured (Experimental) Rates", 
       y = "Predicted Rates", 
       title = "iCHO1766 Predicted vs Measured Rates") +
  facet_wrap(~Window, nrow = 1) 

 
ggsave(filename = "figures/fba/dotplots_growth_rates_case2_1_icho1766_window1_2_3_axes_fixed.png",
       height = 100,
       width = 200,
       units = "mm",
       dpi = 600) 

#titer
ggplot(mp_titer_rates_joined %>% filter(Window %in% c("3", "4", "5")), aes(x = measured, y = icho1766, color = Experiment)) +
  geom_abline(intercept = 0, slope = 1, linewidth = 1) +
  geom_abline(intercept = 0, slope = 0.75, linetype = 2, color = 'grey40') +
  geom_abline(intercept = 0, slope = 1.25, linetype = 2, color = 'grey40') +
  geom_point(alpha = 0.7) +
  geom_errorbar(
    aes(xmin = measured - measured_error,
        xmax = measured + measured_error,
        group = Experiment),
    position = position_dodge(.9),
    width = .001,
    linewidth = .25
  ) +
  # geom_errorbar(
  #   aes(ymin = icho1766 - icho1766_error,
  #       ymax = icho1766 + icho1766_error,
  #       group = Experiment),
  #   position = position_dodge(.9),
  #   width = .001,
  #   linewidth = .25
  # ) +
  scale_color_manual(values = color_mapping_experiment) +
  scale_y_continuous(limits = c(-0.2e-5, 3.5e-5), labels = function(x) x * 1e5) +
  scale_x_continuous(limits = c(-0.2e-5, 3.5e-5), labels = function(x) x * 1e5) +
  coord_fixed() +
  labs(x = "Measured (Experimental) Rates",
       y = "Predicted Rates",
       title = "iCHO1766 Predicted vs Measured Rates") +
  facet_wrap(~Window, nrow = 1)


ggsave(filename = "figures/fba/dotplots_titer_rates_case2_1_icho1766_window3_4_5_coord_fixed.png",
       height = 100,
       width = 200,
       units = "mm",
       dpi = 600)


  
  
  
  
  
  
  
  
  

