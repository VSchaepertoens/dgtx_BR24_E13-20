library(tidyverse)
library("ggradar")

aa_rates <- read_csv("data/AA_rates.csv") %>%
  mutate(Window = str_replace_all(Window, c("76 - 120 h" = "window_1",
                                            "124 - 168 h" = "window_2",
                                            "172 - 216 h" = "window_3",
                                            "220 - 264 h" = "window_4",
                                            "268 - 312 h" = "window_5"))
         )

aa_rates_reordered <- aa_rates %>%
  # mutate(across(Asn_SE, as.numeric)) %>%
  pivot_longer(cols = -c(Condition, Window), 
               names_to = "Variable", 
               values_to = "Value") %>%
  separate(Variable, into = c("amino_acid", "Type"), sep = "_", extra = "merge") %>%
  pivot_wider(names_from = Type, values_from = Value) %>%
  rename("qp" = `NA`, "se" = "SE", "condition" = "Condition", "window" = "Window")


aa_rates_reordered_tshifted <- aa_rates_reordered %>%
  filter(condition == "Temperature shifted")

aa_rates_reordered_nottshifted <- aa_rates_reordered %>%
  filter(condition == "Not temperature shifted")


write_csv(aa_rates_reordered, "data/aa_rates_reordered.csv")
write_csv(aa_rates_reordered_nottshifted, "data/aa_rates_reordered_nottshifted.csv")
write_csv(aa_rates_reordered_tshifted, "data/aa_rates_reordered_tshifted.csv")
  