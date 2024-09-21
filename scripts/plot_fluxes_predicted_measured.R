library(tidyverse)

# prepare growth rates ------------------------------------------
# measured
measured_fluxes <- read_csv("data/aa_rates_reordered_data2.csv")

m_growth_rates <- measured_fluxes %>%
  filter(AA_meta == "Growth_rate") %>%
  mutate(Type = "measured")

# m_growth_rates_nottshift <- m_growth_rates %>%
#   filter(Condition == "Constant")

m_growth_rates_tshift <- m_growth_rates %>%
  filter(Condition == "Temp. shifted")

#predicted
p_growth_rates_iter <- read_csv("fba_results/tshifted/case_1/mus_results_icho2441.csv")
p_growth_rates_iter <- read_csv("fba_results/tshifted/case_2.1/mus_results_icho2441.csv")
p_growth_rates_iter <- read_csv("fba_results/nottshifted/case_3.1/mus_results.csv")

# p_growth_rates_iter <- read_csv("fba_results/tshifted/case_1/mus_results.csv")
# p_growth_rates_iter <- read_csv("fba_results/tshifted/case_2.1/mus_results.csv")
# p_growth_rates_iter <- read_csv("fba_results/tshifted/case_3.1/mus_results.csv")

p_growth_rates <- p_growth_rates_iter %>%
  group_by(Experiment, Window) %>%
  summarise(Rate = mean(mu),
            SD = sd(mu)) %>%
  mutate(Condition = "Temp. shifted",
         AA_meta = "Growth_rate") %>%
  mutate(Type = "predicted") %>%
  ungroup()

p_growth_rates <- p_growth_rates_iter %>%
  group_by(Experiment, Window) %>%
  summarise(Rate = mean(mu),
            SD = sd(mu)) %>%
  mutate(Condition = "Constant",
         AA_meta = "Growth_rate") %>%
  mutate(Type = "predicted") %>%
  ungroup()

#merge dfs
mp_growth_rates <- m_growth_rates_tshift %>%
  rbind(p_growth_rates)

# Transformation to wider table for dot plots
# Rename the Rate column in each dataframe to distinguish between them
p_growth_rates <- p_growth_rates %>% 
  rename(Predicted = Rate, Predicted_sd = SD) %>%
  select(Experiment, Predicted, Predicted_sd, Window, Condition, AA_meta)
# m_growth_rates_tshift <- m_growth_rates_tshift %>% 
#   rename(Experimental = Rate, Experimental_sd = SD) %>%
#   select(Experiment, Window, Experimental, Experimental_sd)

m_growth_rates_nottshift <- m_growth_rates_nottshift %>% 
  rename(Experimental = Rate, Experimental_sd = SD) %>%
  select(Experiment, Window, Experimental, Experimental_sd)

# Combine the two dataframes by the Group column
# df_wide <- left_join(p_growth_rates, m_growth_rates_tshift, 
#                      by = c("Experiment", "Window"),
#                      keep = FALSE)
 
df_wide <- left_join(p_growth_rates, m_growth_rates_nottshift, 
                     by = c("Experiment", "Window"),
                     keep = FALSE) 

# prepare cell specific productivity ----------------------------
# measured
  m_titer_rates <- measured_fluxes %>%
    filter(AA_meta == "Titer") %>%
    mutate(Type = "measured")
  
  # m_titer_rates_nottshift <- m_titer_rates %>%
  #   filter(Condition == "Constant")
  # 
  m_titer_rates_tshift <- m_titer_rates %>%
    filter(Condition == "Temp. shifted")
  
#predicted
  p_titer_rates_iter <- read_csv("fba_results/tshifted/case_2.2/mus_results_icho2441.csv")
  # p_titer_rates_iter <- read_csv("fba_results/nottshifted/case_3.2/mus_results.csv")
  
  # p_titer_rates_iter <- read_csv("fba_results/tshifted/case_2.2/mus_results.csv")
  # p_titer_rates_iter <- read_csv("fba_results/tshifted/case_3.2/mus_results.csv")
  
  # p_titer_rates <- p_titer_rates_iter %>%
  #   group_by(Experiment, Window) %>%
  #   summarise(Rate = mean(mu),
  #             SD = sd(mu)) %>%
  #   mutate(Condition = "Constant",
  #          AA_meta = "Titer") %>%
  #   mutate(Type = "predicted") %>%
  #   ungroup()
  
  p_titer_rates <- p_titer_rates_iter %>%
    group_by(Experiment, Window) %>%
    summarise(Rate = mean(mu),
              SD = sd(mu)) %>%
    mutate(Condition = "Temp. shifted",
           AA_meta = "Titer") %>%
    mutate(Type = "predicted") %>%
    ungroup()
  
# merge dfs
  # mp_titer_rates <- m_titer_rates_nottshift %>%
  #   rbind(p_titer_rates)

  mp_titer_rates <- m_titer_rates_tshift %>%
    rbind(p_titer_rates)

  
  # Transformation to wider table for dot plots
  # Rename the Rate column in each dataframe to distinguish between them
  p_titer_rates <- p_titer_rates %>% 
    rename(Predicted = Rate, Predicted_sd = SD) %>%
    select(Experiment, Predicted, Predicted_sd, Window, Condition, AA_meta)

  
  m_titer_rates_tshift <- m_titer_rates_tshift %>% 
    rename(Experimental = Rate, Experimental_sd = SD) %>%
    select(Experiment, Window, Experimental, Experimental_sd)
  
  m_titer_rates_nottshift <- m_titer_rates_nottshift %>% 
    rename(Experimental = Rate, Experimental_sd = SD) %>%
    select(Experiment, Window, Experimental, Experimental_sd)
  # Combine the two dataframes by the Group column
  # df_wide <- left_join(p_growth_rates, m_growth_rates_tshift, 
  #                      by = c("Experiment", "Window"),
  #                      keep = FALSE)
  
  df_wide <- left_join(p_titer_rates, m_titer_rates_nottshift, 
                       by = c("Experiment", "Window"),
                       keep = FALSE) 

# function to plot data ---------------------------------------------------
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
growth_rate_limits = c(-0.02, 0.06)

plot_rates(data_to_plot = mp_growth_rates,
           y_axis_label = y_label_growth_rate,
           y_axis_limits = growth_rate_limits,
           plot_growth_rate = TRUE) 

ggsave(filename = "figures/growth_rates_tshifted_case2_1_icho2441.png",
       height = 70,
       width = 200,
       units = "mm",
       dpi = 600)


# titer
y_label_titer_rate <- expression(paste("specific antibody productivity (mM g DC",W^-1,h^-1,") x",~ 10^-5))
titer_rate_limits = c(-8e-5, 8e-5)

plot_rates(data_to_plot = mp_titer_rates,
           y_axis_label = y_label_titer_rate,
           y_axis_limits = titer_rate_limits,
           plot_growth_rate = FALSE)

ggsave(filename = "figures/titer_rates_tshifted_case2_2_icho2441.png",
       height = 70,
       width = 200,
       units = "mm",
       dpi = 600)

 

# exp_growth vs pred_growth plots -------------------------------------------


# Plotting measured vs predicted rates
ggplot(df_wide, aes(x = Experimental, y = Predicted, color = Experiment)) +
  geom_abline(intercept = 0, slope = 1, size = 1) +
  geom_abline(intercept = 0, slope = 0.75, linetype = 2, color = 'grey40') +
  geom_abline(intercept = 0, slope = 1.25, linetype = 2, color = 'grey40') +
  geom_point(size = 3) +
  geom_errorbar(
    aes(xmin = Experimental - Experimental_sd,
        xmax = Experimental + Experimental_sd,
        group = Experiment),
    position = position_dodge(.9),
    width = .001,
    linewidth = .25
  ) +
  geom_errorbar(
    aes(ymin = Predicted - Predicted_sd,
        ymax = Predicted + Predicted_sd,
        group = Experiment),
    position = position_dodge(.9),
    width = .001,
    linewidth = .25
  ) +
  scale_y_continuous(limits = c(-8e-5, 8e-5), labels = function(x) x * 1e5) +
  scale_x_continuous(limits = c(-8e-5, 8e-5), labels = function(x) x * 1e5) +
  labs(x = "Measured (Experimental) Rates", 
       y = "Predicted Rates", 
       title = "Predicted vs Measured Rates") +
  facet_wrap(~Window, nrow = 1)

 
ggsave(filename = "figures/fba/dotplots_titer_rates_nottshifted_case3_1.png",
       height = 70,
       width = 200,
       units = "mm",
       dpi = 600) 

  
  
  
  
  
  
  
  
  
  
  

