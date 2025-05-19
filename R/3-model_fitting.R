
# GitHub repo: github.com/dario-ssm/phenodev_predict

# 0. Load ----
library(terra)
library(tidyverse)
#library(lubridate)
library(here)
library(rTPC)
library(nls.multstart)
library(ggthemes)
library(chillR)
library(sf)
library(viridis)
source(here("R/1-functions_phenodev_predict.R"))
library(mappestRisk)

# 1. Model fitting ------------------------------------------------------

#####  a) von Schmalensee 2023 ---------------------------------------------------

schmalensee_pieris_data <- read_delim("~/Data and scripts von Schmalensee et al. 2023/Data and scripts von Schmalensee et al. 2023/Data/full_data.txt", 
                                      delim = "\t") |> 
  filter(species == "rapae",
         life.stage == "pupa") |> 
  select(temp, dev.rate) |> 
  drop_na() 

fit_models_pieris_rapae <- mappestRisk::fit_devmodels(temp = schmalensee_pieris_data$temp,
                                                      dev_rate = schmalensee_pieris_data$dev.rate,
                                                      model_name = "all")
my_species_name <- expression(~italic("Pieris rapae"))

lm_schmalensee <- lm(dev.rate~temp,
                     data = schmalensee_pieris_data)
params_lm_schmalensee <- coef(lm_schmalensee)
  
tbl_preds_lm_schmalensee <- tibble(
  temp_preds = seq(min(schmalensee_pieris_data$temp) -15, 
                  max(schmalensee_pieris_data$temp) + 15, 0.01)) |> 
  mutate(preds = map_dbl(.x = temp_preds,
                         .f = ~params_lm_schmalensee[1] + params_lm_schmalensee[2]*.x)) |> 
  filter(preds >= 0)

plot_preds_lm_schmalensee <- ggplot(data = tbl_preds_lm_schmalensee,
                                    aes(x = temp_preds,
                                        y = preds))+
  geom_point(data = schmalensee_pieris_data,
             aes(x = temp,
                 y = dev.rate),
             color = "darkslategray",
             alpha = 0.8, size = 1.5)+
  geom_line(color = "lightcoral", linewidth = 1.3)+
  theme_bw()+
  labs(x = "Temperature",
       y = expression(italic(R(T)) ~ 
                        (d^-1)))
ggsave(here("figures/schmalensee_lm.svg"),
       height = 800,
       width = 800,
       units = "px")

ggsave(here("figures/schmalensee_lm.png"),
       height = 800,
       width = 800,
       units = "px")

plot_devmodels(temp = schmalensee_pieris_data$temp,
               dev_rate = schmalensee_pieris_data$dev.rate,
               fitted_parameters = fit_models_pieris_rapae,
               species = "Pieris rapae",
               life_stage = "Pupa")+
  labs(title = my_species_name)+
  khroma::scale_color_batlow(discrete = TRUE)+
  khroma::scale_fill_batlow(discrete = TRUE)+
  theme_few()+
  theme(legend.position = "none")

ggsave(here("figures/schmalensee_tpcs.png"),
       height = 2600,
       width = 2600,
       units = "px")
  

boots_pieris_rapae <- mappestRisk::predict_curves(temp = schmalensee_pieris_data$temp,
                                                  dev_rate = schmalensee_pieris_data$dev.rate,
                                                  fitted_parameters = fit_models_pieris_rapae,
                                                  model_name_2boot = c("schoolfield", "wang", "mod_polynomial",
                                                                       "ratkowsky", "thomas","lrf", "boatman",
                                                                       "lactin2", "joehnk", "briere2", "beta", 
                                                                       "mod_weibull",                                                                       "oneill", "kamykowski", "lactin1", "pawar"),
                                                  propagate_uncertainty = TRUE,
                                                  n_boots_samples = 100) 
                                                  
mappestRisk::plot_uncertainties(bootstrap_uncertainties_tpcs = boots_pieris_rapae,
                                temp = schmalensee_pieris_data$temp,
                                dev_rate = schmalensee_pieris_data$dev.rate,
                                species = "Pieris rapae",
                                life_stage = "Pupa")
# exclude unrealistic TPC shapes or unreasonable uncertainties 
selected_models_pieris_rapae <- fit_models_pieris_rapae #|> 
 # filter(!model_name %in% c("mod_polynomial", "oneill", "mod_weibull",
  #                          "ratkowsky", "schoolfield", "wang"))
plot_devmodels(temp = schmalensee_pieris_data$temp,
               dev_rate = schmalensee_pieris_data$dev.rate,
               fitted_parameters = selected_models_pieris_rapae,
               species = "Pieris rapae",
               life_stage = "Pupa")+
  labs(title = my_species_name)+
  khroma::scale_color_batlow(discrete = TRUE)+
  khroma::scale_fill_batlow(discrete = TRUE)

write_rds(selected_models_pieris_rapae,
          here("data/selected_models_schmalensee.rds"))

#####  b) Gilbert, 1987 ---------------------------------------------------

gilbert_pieris_data <- read_delim(here("data/gilbert_pupa.csv"),
                                  delim = ";") |> 
  rename(temp = temperature,
         dev.rate = devrate,
         life_stage = ...3,
         notes = ...4,
         reference = interact) |> 
  filter(life_stage == "pupa") |> 
  dplyr::select(temp, dev.rate)
  
  
fit_models_pieris_rapae <- mappestRisk::fit_devmodels(temp = gilbert_pieris_data$temp,
                                                      dev_rate = gilbert_pieris_data$dev.rate,
                                                      model_name = "all")
my_species_name <- expression(~italic("Pieris rapae"))

plot_devmodels(temp = gilbert_pieris_data$temp,
               dev_rate = gilbert_pieris_data$dev.rate,
               fitted_parameters = fit_models_pieris_rapae,
               species = "Pieris rapae",
               life_stage = "Pupa")+
  labs(title = my_species_name)+
  khroma::scale_color_batlow(discrete = TRUE)+
  khroma::scale_fill_batlow(discrete = TRUE)+
  theme_few()+
  theme(legend.position = "none")

ggsave(here("figures/gilbert_tpcs.png"),
       width = 2600,
       height = 2600,
       units = "px")
lm_gilbert <- lm(dev.rate~temp,
                     data = gilbert_pieris_data)
params_lm_gilbert <- coef(lm_gilbert)

tbl_preds_lm_gilbert <- tibble(
  temp_preds = seq(min(gilbert_pieris_data$temp) -15, 
                   max(gilbert_pieris_data$temp) + 15, 0.01)) |> 
  mutate(preds = map_dbl(.x = temp_preds,
                         .f = ~params_lm_gilbert[1] + params_lm_gilbert[2]*.x)) |> 
  filter(preds >= 0)

plot_preds_lm_gilbert <- ggplot(data = tbl_preds_lm_gilbert,
                                    aes(x = temp_preds,
                                        y = preds))+
  geom_point(data = gilbert_pieris_data,
             aes(x = temp,
                 y = dev.rate),
             color = "darkslategray",
             alpha = 0.8, size = 1.5)+
  geom_line(color = "lightcoral", linewidth = 1.3)+
  theme_bw()+
  labs(x = "Temperature",
       y = expression(italic(R(T)) ~ 
                        (d^-1)))
plot_preds_lm_gilbert

ggsave(here("figures/gilbert_lm.svg"),
       height = 800,
       width = 800,
       units = "px")

ggsave(here("figures/gilbert_lm.png"),
       height = 800,
       width = 800,
       units = "px")


boots_pieris_rapae <- mappestRisk::predict_curves(temp = gilbert_pieris_data$temp,
                                                  dev_rate = gilbert_pieris_data$dev.rate,
                                                  fitted_parameters = fit_models_pieris_rapae,
                                                  model_name_2boot = c("oneill", "mod_weibull", "ratkowsky",
                                                                       "lrf", "thomas", "briere2", "beta", "boatman",
                                                                       "wang", "joehnk", "lactin1", "kamykowski", "flextpc"),
                                                  propagate_uncertainty = TRUE,
                                                  n_boots_samples = 100) 

mappestRisk::plot_uncertainties(bootstrap_uncertainties_tpcs = boots_pieris_rapae,
                                temp = gilbert_pieris_data$temp,
                                dev_rate = gilbert_pieris_data$dev.rate,
                                species = "Pieris rapae",
                                life_stage = "Pupa")
# exclude unrealistic TPC shapes or unreasonable uncertainties 
selected_models_pieris_rapae_gilbert <- fit_models_pieris_rapae |> 
  filter(!model_name %in% c("mod_polynomial"))
plot_devmodels(temp = gilbert_pieris_data$temp,
               dev_rate = gilbert_pieris_data$dev.rate,
               fitted_parameters = selected_models_pieris_rapae_gilbert,
               species = "Pieris rapae",
               life_stage = "Pupa")+
  labs(title = my_species_name)+
  khroma::scale_color_batlow(discrete = TRUE)+
  khroma::scale_fill_batlow(discrete = TRUE)
ggsave(here("figures/gilbert_tpcs.png"),
       width = 2600,
       height = 2600,
       units = "px")

write_rds(selected_models_pieris_rapae_gilbert,
           here("data/selected_models_gilbert.rds"))

