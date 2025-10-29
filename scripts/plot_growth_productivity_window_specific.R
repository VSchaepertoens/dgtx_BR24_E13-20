library(tidyverse)

# measured_fluxes <- read_csv("data/aa_rates_reordered_data2.csv")

# GROWTH rates ------------------------------------------
## MEASURED ##
m_growth_rates <- read_csv('data/rates_05032025/growth_rate.csv') %>%
  # filter(AA_meta == "Growth_rate") %>%
  mutate(Type = "measured",
         Metabolite = "Growth_rate")
  # select(!Condition)


p_growth_rates <- read_csv("fba_results/window_specific/iCHO1766_biomass_cho_producing_notbiomasscorrected/mus_results_icho1766_FBA_pFBA.csv") %>%
  mutate(Metabolite = "Growth_rate",
         Type = "iCHO1766"
         ) %>%
  pivot_longer(cols = c("FBA_mu", "pFBA_mu"), 
               names_to = "method_type"
               ) %>%
  mutate(Type = paste(Type, method_type, sep = "_"),
         SD = 0
         ) %>%
  mutate(Type = gsub("_mu", "", Type)) %>%
  rename(Rate = value) %>%
  select(!method_type)

# Bind predicted and measured by rows for plotting bar plots
mp_growth_rates <- m_growth_rates %>%
  rbind(p_growth_rates) 

# Join predicted and measured for plotting dot plots
mp_growth_rates_joined <- m_growth_rates  %>%
  full_join(p_growth_rates %>% filter(Type != "iCHO1766_FBA"), 
            by = c("Condition", "Window", "Metabolite")
            ) %>%
  rename(measured = Rate.x,
         measured_error = SD.x,
         icho1766 = Rate.y
         ) %>%
  select(!c("Type.x", "Type.y")) %>%
  mutate(Window = as.character(Window))


# TITER rates ----------------------------
## MEASURED ##
  m_titer_rates <- measured_fluxes %>%
    filter(AA_meta == "Titer") %>%
    mutate(Type = "measured")%>%
    select(!Condition)

## PREDICTED iCHO1766 ##
p_titer_rates <- read_csv("fba_results/condition_specific/iCHO1766_igg/mus_results_icho1766_FBA_pFBA.csv") %>%
  mutate(AA_meta = "Titer",
         Type = "iCHO1766"
         ) %>%
  pivot_longer(cols = c("FBA_mu", "pFBA_mu"), 
               names_to = "method_type"
               ) %>%
  mutate(Type = paste(Type, method_type, sep = "_"),
         SD = 0
         ) %>%
  mutate(Type = gsub("_mu", "", Type)
         ) %>%
  rename(Rate = value) %>%
  select(!method_type)

# Bind predicted and measured by rows for plotting bar plots
mp_titer_rates <- m_titer_rates %>%
  rbind(p_titer_rates) 

# Join predicted and measured for plotting dot plots
mp_titer_rates_joined <- m_titer_rates %>%
  full_join(p_titer_rates %>% filter(Type != "iCHO1766_pFBA"), 
            by = c("Experiment", "Window", "AA_meta")
            ) %>%
  rename(measured = Rate.x,
         measured_error = SD.x,
         icho1766 = Rate.y
         ) %>%
  select(!c("Type.x", "Type.y")) %>%
  mutate(Window = as.character(Window))
  
# function to plot barplots ---------------------------------------------------
plot_barplot <- function(data_to_plot,
                       y_axis_label,
                       y_axis_limits,
                       plot_growth_rate = FALSE){
   ggplot(data_to_plot, aes(x = Condition)) +
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
  
## GROWTH RATE
y_label_growth_rate <- expression(paste("growth rate (",h^-1,")"))
growth_rate_limits <-  c(round(min(mp_growth_rates$Rate - mp_growth_rates$SD, na.rm = TRUE), digits = 7), 
                         round(max(mp_growth_rates$Rate + mp_growth_rates$SD, na.rm = TRUE), digits = 7))

plot_barplot(data_to_plot = mp_growth_rates %>% filter(Type != "iCHO1766_pFBA"),
             y_axis_label = y_label_growth_rate,
             y_axis_limits = c(-0.00153, 0.038),
             plot_growth_rate = TRUE)

# ggsave(filename = "figures/fba/pfba_fba_fva/condition_specific_biomass_FBA.png",
#        height = 70,
#        width = 250,
#        units = "mm",
#        dpi = 600)


# # TITER RATE
# y_label_titer_rate <- expression(paste("specific antibody productivity (mM g DC",W^-1,h^-1,") x",~ 10^-5))
# titer_rate_limits <-  c(round(min(mp_titer_rates$Rate - mp_titer_rates$SD, na.rm = TRUE), digits = 7), 
#                       round(max(mp_titer_rates$Rate + mp_titer_rates$SD, na.rm = TRUE), digits = 7))
# 
# plot_barplot(data_to_plot = mp_titer_rates %>% filter(Type != "iCHO1766_FBA"),
#            y_axis_label = y_label_titer_rate,
#            y_axis_limits = titer_rate_limits,
#            plot_growth_rate = FALSE)
# 
# ggsave(filename = "figures/fba/pfba_fba_fva/condition_specific_igg_pFBA.png",
#        height = 70,
#        width = 250,
#        units = "mm",
#        dpi = 600)


# dotplots exp_growth vs pred_growth plots -------------------------------------------
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
  # "E13" = "#FD8D3C",
  # "E14" = "#9E9AC8",
  # "E15" = "#F16913",
  # "E16" = "#807DBA",
  # "E17" = "#D94801",
  # "E18" = "#6A51A3",
  "Constant" = "#A63603",
  "Temp. shifted" = "#54278F"
)

# function to plots dot plots of measured vs predicted rates
plot_dotplot <- function(in_data = mp_growth_rates_joined %>% filter(Window %in% c("1", "2", "3")),
                         x_y_lim = c(-0.005, 0.04),
                         facet_window = TRUE,
                         xlab,
                         ylab,
                         plot_titer = FALSE
                         ){
  
  # Define point layer
  point_layer <- geom_point(aes(shape = if (facet_window) NULL else Window))
  
  # Define facet layer
  facet_layer <- if (facet_window) facet_wrap(~Window, nrow = 1) else NULL

  # Define axis labels for Titer
  axis_labels <- if (plot_titer) function(x) x * 1e5 else identity
  
  # Define abline parameters
  abline_params <- list(
    geom_abline(intercept = 0, slope = 1, linewidth = 0.5),
    geom_abline(intercept = 0, slope = 0.75, linetype = 2, color = 'grey40'),
    geom_abline(intercept = 0, slope = 1.25, linetype = 2, color = 'grey40')
  )
  
  ggplot(in_data, aes(x = measured, y = icho1766, color = Condition)) +
    abline_params +
    point_layer +
    geom_errorbar(
      aes(xmin = measured - measured_error,
          xmax = measured + measured_error,
          group = Condition),
      position = position_dodge(.9),
      width = .00025,
      linewidth = .25
    ) +
    scale_color_manual(values = color_mapping_condition) +
    scale_y_continuous(name = ylab,
                       limits = x_y_lim,
                       labels = axis_labels) +
    scale_x_continuous(name = xlab,
                       limits = x_y_lim,
                       labels = axis_labels) +
    coord_fixed() +
    labs(title = "iCHO1766 Predicted vs Measured Rates") +
    theme_minimal() +
    facet_layer +
    # Enhance readability and aesthetics 
    theme(
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "top"
    ) +
    # Fix legend title for 'Window'
    guides(shape = guide_legend(title = "Window"))
}


## GROWTH rates
# faceted by window
growth_rate_limits <-  c(round(min(mp_growth_rates$Rate - mp_growth_rates$SD, na.rm = TRUE), digits = 7), 
                         round(max(mp_growth_rates$Rate + mp_growth_rates$SD, na.rm = TRUE), digits = 7))

plot_dotplot(in_data = mp_growth_rates_joined,
             x_y_lim = c(-0.00153, 0.038),
             facet_window = TRUE,
             xlab = expression(paste("Measured (Experimental) Rates [",h^-1 ,"]")), 
             ylab = expression(paste("Predicted Rates [",h^-1 ,"]")),
             plot_titer = FALSE
             )
ggsave(filename = "figures/fba/pfba_fba_fva/window_specific_notbiomasscorrected/dotplots_biomass_FBA_facet_window.png",
       height = 100,
       width = 200,
       units = "mm",
       dpi = 600,
       bg = "white") 


mp_growth_rates_w123 <- mp_growth_rates %>% filter(Window %in% c("1", "2", "3"))
growth_rate_limits_w123 <-  c(round(min(mp_growth_rates_w123$Rate - mp_growth_rates_w123$SD, na.rm = TRUE), digits = 7), 
                              round(max(mp_growth_rates_w123$Rate + mp_growth_rates_w123$SD, na.rm = TRUE), digits = 7))

plot_dotplot(in_data = mp_growth_rates_joined %>% filter(Window %in% c("1", "2", "3")),
             x_y_lim = c(-0.00153, 0.038),
             facet_window = TRUE,
             xlab = expression(paste("Measured (Experimental) Rates [",h^-1 ,"]")), 
             ylab = expression(paste("Predicted Rates [",h^-1 ,"]")),
             plot_titer = FALSE)

ggsave(filename = "figures/fba/pfba_fba_fva/window_specific_notbiomasscorrected/dotplots_biomass_FBA_facet_w123.png",
       height = 100,
       width = 200,
       units = "mm",
       dpi = 600,
       bg = "white") 

# not faceted by window
plot_dotplot(in_data = mp_growth_rates_joined ,
             x_y_lim = c(-0.00153, 0.038),
             facet_window = FALSE,
             xlab = expression(paste("Measured (Experimental) Rates [",h^-1 ,"]")), 
             ylab = expression(paste("Predicted Rates [",h^-1 ,"]")),
             plot_titer = FALSE)

ggsave(filename = "figures/fba/pfba_fba_fva/window_specific_notbiomasscorrected/dotplots_biomass_FBA_NOfacet_window.png",
       height = 120,
       width = 170,
       units = "mm",
       dpi = 600,
       bg = "white") 

plot_dotplot(in_data = mp_growth_rates_joined %>% filter(Window %in% c("1", "2", "3")),
             x_y_lim = c(-0.00153, 0.038),
             facet_window = FALSE,
             xlab = expression(paste("Measured (Experimental) Rates [",h^-1 ,"]")), 
             ylab = expression(paste("Predicted Rates [",h^-1 ,"]")),
             plot_titer = FALSE)

ggsave(filename = "figures/fba/pfba_fba_fva/window_specific_notbiomasscorrected/dotplots_biomass_FBA_NOfacet_w123.png",
       height = 120,
       width = 150,
       units = "mm",
       dpi = 600,
       bg = "white")



## TITER rates
# faceted by window
titer_rate_limits <-  c(round(min(mp_titer_rates$Rate - mp_titer_rates$SD, na.rm = TRUE), digits = 7), 
                         round(max(mp_titer_rates$Rate + mp_titer_rates$SD, na.rm = TRUE), digits = 7))

plot_dotplot(in_data = mp_titer_rates_joined,
             x_y_lim = titer_rate_limits,
             facet_window = TRUE,
             xlab = expression(paste("Measured (Experimental) Rates [mM g DC",W^-1,h^-1,"] x",~ 10^-5)), 
             ylab = expression(paste("Predicted Rates  [mM g DC",W^-1,h^-1,"] x",~ 10^-5)),
             plot_titer = TRUE
             )

ggsave(filename = "figures/fba/pfba_fba_fva/dotplots_igg_FBA_facet_window.png",
       height = 55,
       width = 200,
       units = "mm",
       dpi = 600,
       bg = "white") 




  
  
  
  
  
  
  
  
  

