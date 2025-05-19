library(terra)
library(tidyverse)
#library(lubridate)
library(here)
library(rTPC)
library(nls.multstart)
library(ggthemes)
library(chillR)
#library(sf)
#library(viridis)
source(here("Scripts/1-functions_phenodev_predict.R"))
library(mappestRisk)

## 1. Temperatures mainland Spain ----

load(here("Data/daily_tmax_df.RData"))
print(daily_tmax_df)
load(here("Data/daily_tmin_df.RData"))
print(daily_tmin_df)
daily_temperatures <- inner_join(daily_tmin_df, daily_tmax_df) |> 
  mutate(year = year(date),
         daily_tavg = map2_dbl(.x = daily_tmin,
                               .y = daily_tmax,
                               .f = ~mean(c(.x, .y))
         )
  )

daily_temperatures_spaincase_vert <- daily_temperatures |> 
  rename(Tmin = daily_tmin,
         Tmax = daily_tmax) |> 
  pivot_longer(cols = c("Tmin", "Tmax"),
               names_to = "var",
               values_to = "temperature") |> 
  rename()

spaincase_daily_temperatures_rainplot <- ggplot(data = daily_temperatures_spaincase_vert,
                                                aes(x = temperature,
                                                    fill = var,
                                                    color = var))+
  ggdist::stat_dist_halfeye(alpha = .5)+
  theme_few()+
  labs(x = "Temperature (ºC)",
       y = NULL,
       fill = "Variable",
       color = "Variable")+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

spaincase_daily_temperatures_rainplot
ggsave(plot = spaincase_daily_temperatures_rainplot,
       filename = here("figures/spaincase_daily_temperatures_rainplot.png"),
       width = 2600, height = 2600, units = "px")




# 2. Gilbert's data ---------------------------------------------------------
selected_models_pieris_rapae <- readRDS(here("data/selected_models_gilbert.rds"))
gilbert_pieris_data <- read_delim(here("data/gilbert_pupa.csv"),
                                  delim = ";") |> 
  rename(temp = temperature,
         dev.rate = devrate,
         life_stage = ...3,
         notes = ...4,
         reference = interact) |> 
  filter(life_stage == "pupa") |> 
  dplyr::select(temp, dev.rate)

##### 2.1. Daily resolution ----

###### a) nonlinear rate summation ----

daily_rate_preds <- tibble()

for(model_i in unique(selected_models_pieris_rapae$model_name)) {
  rate_daily_i <- rate_pred(model2ratesum = model_i,
                            temperature = daily_temperatures,
                            fitted_params = selected_models_pieris_rapae,
                            res = "daily")
  daily_rate_preds <- bind_rows(daily_rate_preds, rate_daily_i)
}

rate_sum_preds_daily <- daily_rate_preds |> 
  mutate(rate_pred = case_when(rate_pred < 0 ~ 0, # <- avoid negative summation
                               is.na(rate_pred) ~0, # <- zeros where model yields NA
                               .default = rate_pred), 
         rate_pred = 100*rate_pred) |>  # as rate_dev equals 1/days_dev
  group_by(year, model_name) |> 
  mutate(rate_sum = cumsum(rate_pred), #heat summation
         emergence_lgl = logic_ratesum(rate_sum), #threshold point
         doy = yday(date)) |> 
  summarise(day_of_emergence = sum(emergence_lgl))  #number of days until threshold is reach

rate_sum_preds_daily

###### b) linear DDs ----
lm_ratedev <- lm(dev.rate ~ temp,
                 data = gilbert_pieris_data) #or gilbert pieris data


dd_preds_daily <- daily_temperatures |> 
  mutate(model_name = "linear") |> 
  group_by(year, model_name) |> 
  mutate(rate_pred = map_dbl(.x = daily_tavg,
                             .f = ~lm_pred(lm_fit = lm_ratedev,
                                           .x))) |> 
  mutate(rate_pred = if_else(rate_pred <0,
                             0,
                             rate_pred)) |>
  mutate(rate_sum = map_dbl(.x = daily_tavg,
                            .f = ~calc_dd(lm_fit = lm_ratedev,
                                          tavg = .x)), #heat summation
         rate_sum = cumsum(rate_sum),
         emergence_lgl = map_dbl(.x = rate_sum,
                                 .f = ~logic_dd(heat_units = 1/coef(lm_ratedev)[2],
                                                rate_cumsum = .x)), #threshold point
         doy = yday(date)) |> 
  summarise(day_of_emergence = sum(emergence_lgl))

###### c) combined  ----
doe_preds_daily <- rate_sum_preds_daily |> 
  full_join(dd_preds_daily) |> 
  mutate(time_slot = "daily")

##### 2.2. Hourly res. ----------------------------------------------------------
##  we use chillR package to simulate temperature variation across hours within days.
###### a) nonlinear rate summation ----

daily_temps_doy <- daily_temperatures  |> 
  mutate(JDay = yday(date)) |>  
  rename(Tmin = daily_tmin,
         Tmax = daily_tmax) 

hourly_temp_doy <- chillR::make_hourly_temps(latitude = 40.5,
                                             year_file = daily_temps_doy,
                                             keep_sunrise_sunset = TRUE )  |> 
  pivot_longer(cols = 10:33,
               names_to = "hour_of_day",
               values_to = "temperature") |>  
  rename(tmax = Tmax,
         tmin = Tmin,
         doy = JDay,
         sunrise = Sunrise,
         sunset = Sunset,
         daylength = Daylength) 


hourly_rate_preds <- tibble()

for(model_i in unique(selected_models_pieris_rapae$model_name)) {
  rate_hourly_i <- rate_pred(model2ratesum = model_i,
                             temperature = hourly_temp_doy,
                             fitted_params = selected_models_pieris_rapae,
                             res = "hourly")
  hourly_rate_preds <- bind_rows(hourly_rate_preds, rate_hourly_i)
}
hourly_rate_preds

rate_sum_preds_hourly <- hourly_rate_preds |> 
  mutate(rate_pred = case_when(rate_pred < 0 ~ 0, # <- avoid negative summation
                               is.na(rate_pred) ~0, # <- zeros where model yields NA
                               .default = rate_pred), 
         rate_pred = rate_pred*100/24) |>  # 24h a day
  group_by(year, model_name) |> 
  mutate(rate_sum = cumsum(rate_pred), #heat summation
         emergence_lgl = logic_ratesum(rate_sum), #threshold point
         doy = yday(date)) |> 
  summarise(day_of_emergence = round(sum(emergence_lgl)/24))  #number of days until threshold is reach

rate_sum_preds_hourly

###### b) linear DDs ----
dd_preds_hourly <- hourly_temp_doy |> 
  mutate(model_name = "linear") |> 
  group_by(year, model_name) |> 
  mutate(rate_pred = map_dbl(.x = temperature,
                             .f = ~lm_pred(lm_fit = lm_ratedev,
                                           .x))) |> 
  mutate(rate_pred = if_else(rate_pred <0,
                             0,
                             rate_pred)) |> #fraction of 1/24 hours
  mutate(rate_sum = map_dbl(.x = temperature,
                            .f = ~calc_dd(lm_fit = lm_ratedev,
                                          tavg = .x)), #heat summation
         rate_sum = cumsum(rate_sum),
         emergence_lgl = map_dbl(.x = rate_sum,
                                 .f = ~logic_dd(heat_units = 24*1/coef(lm_ratedev)[2],
                                                rate_cumsum = .x)), #threshold point
         doy = yday(date)) |> 
  summarise(day_of_emergence = round(sum(emergence_lgl)/24))

###### c) combined  ----
doe_preds_hourly <- rate_sum_preds_hourly |> 
  full_join(dd_preds_hourly) |> 
  mutate(time_slot = "hourly")

### 2.3. Validation ---------------------------------------
###### a) Visual comparison ----

# incorporate observational data from Gordo & Sanz (2006)
pieris_obs_trends <- read_csv(here("data/gordo_sanz_2006_doe.csv")) |> 
  rename(day_of_emergence = doe) |> 
  relocate(year, day_of_emergence) |> 
  group_by(year) |> 
  slice(1)

rate_sum_validation <- doe_preds_daily |> 
  bind_rows(doe_preds_hourly)

###### b) RMSE ----
preds_vs_obs_trends <- rate_sum_validation  |>
  filter(year > 1951 & year < 2005) |> 
  group_by(model_name, time_slot) |> 
  mutate(observed_emergence = pieris_obs_trends$day_of_emergence) |> 
  summarise(rmsep = chillR::RMSEP(predicted = day_of_emergence,
                                  observed = observed_emergence)) |>  
  arrange(rmsep) %>% 
  print() 

## see which temperature temporal resolution predicts better
transfer_rank_resolution <- preds_vs_obs_trends %>% 
  group_by(time_slot) %>% 
  summarise(rmsep = mean(rmsep)) %>% 
  arrange(rmsep) %>% 
  print()

## see which models predict better
transfer_rank_model <- preds_vs_obs_trends %>% 
  group_by(model_name) %>% 
  summarise(rmsep = mean(rmsep)) %>% 
  arrange(rmsep) %>% 
  print() # 

# Make Figure 1

preds_obs_compare_rmsep <- rate_sum_validation |> 
  full_join(preds_vs_obs_trends) |> 
  group_by(year, model_name, time_slot) |> 
  summarise(rmsep = min(rmsep)) |>  
  arrange(rmsep) |> 
  print()
rmsep_text <-  preds_obs_compare_rmsep |>  
  group_by(model_name, time_slot) |>  
  summarise(rmsep = mean(rmsep)) |> 
  arrange(rmsep)
rmsep_order <- rmsep_text |> 
  pull(model_name)
rmsep_order_min <- rmsep_text |> 
  slice_min(rmsep) |>
  arrange(rmsep) |> 
  pull(model_name)
rmsep_values <- rmsep_text |> 
  mutate(rmsep =  case_when(time_slot == "daily" ~paste("d-RMSEP =",round(rmsep,2)),
                            time_slot == "hourly" ~paste("h-RMSEP =",round(rmsep,2))),
         year = 2004,
         day_of_emergence = case_when(time_slot == "daily" ~165,
                                      time_slot == "hourly" ~145))
preds_obs_compare_plot <- ggplot()+
  geom_point(data = rate_sum_validation |> 
               filter(year %in% pieris_obs_trends$year),
             aes(x = year, y = day_of_emergence, color = time_slot),
             alpha = .66)+
  geom_line(data = rate_sum_validation |> 
              filter(year %in% pieris_obs_trends$year),
            aes(x = year, y = day_of_emergence, color = time_slot),
            linetype = "longdash",
            alpha = .66)+
  geom_smooth(data =  rate_sum_validation |> 
                filter(year %in% pieris_obs_trends$year),
              aes(x = year, y = day_of_emergence, 
                  color = time_slot, 
                  fill = time_slot))+
  scale_color_manual(values = c("#3A848B", "#E85038"))+
  scale_fill_manual(values = c("#3A848B", "#E85038"))+
  geom_point(data = pieris_obs_trends, 
             aes(x = year, y = day_of_emergence),
             color = "#E9C86B",
             fill = "#E9C86B",
             alpha = .66)+
  geom_line(data = pieris_obs_trends, 
            aes(x = year, y = day_of_emergence),
            color = "#E9C86B",
            linetype = "longdash",
            alpha = .66)+
  geom_smooth(data = pieris_obs_trends, 
              aes(x = year, y = day_of_emergence),
              color = "#E9C86B",
              fill = "#E9C86B")+
  labs(x = "Year",
       y = "Day of first adult emergence",
       color = "Resolution (model)",
       fill = "Resolution (model)")+
  facet_wrap(~factor(model_name, levels = rmsep_order_min))+
  theme_few()+
  theme(strip.text.x = element_text(size = 10, 
                                    color = "grey42", 
                                    face = "bold"))+
  geom_text(data = rmsep_values, 
            aes(label = rmsep,
                x = year,
                y = day_of_emergence), 
            size = 3,
            fontface = "italic",hjust = "right")
preds_obs_compare_plot  

## save for gilbert
ggsave(filename = here("figures/gilbert_modelpreds_validate_rmsep.png"),
       width = 25,
       height = 25,
       units = "cm")
ggsave(filename = here("figures/gilbert_modelpreds_validate_rmsep.svg"),
       width = 25,
       height = 25,
       units = "cm")

# 3. Schmalensee's data ---------------------------------------------------------
selected_models_pieris_rapae <- readRDS(here("data/selected_models_schmalensee.rds"))

schmalensee_pieris_data <- read_delim(here("data/schmalensee2023.txt"), 
                                      delim = "\t") |> 
  filter(species == "rapae",
         life.stage == "pupa") |> 
  select(temp, dev.rate) |> 
  drop_na() 

##### 3.1. Daily resolution ----

###### a) nonlinear rate summation ----

daily_rate_preds <- tibble()

for(model_i in unique(selected_models_pieris_rapae$model_name)) {
  rate_daily_i <- rate_pred(model2ratesum = model_i,
                            temperature = daily_temperatures,
                            fitted_params = selected_models_pieris_rapae,
                            res = "daily")
  daily_rate_preds <- bind_rows(daily_rate_preds, rate_daily_i)
}

rate_sum_preds_daily <- daily_rate_preds |> 
  mutate(rate_pred = case_when(rate_pred < 0 ~ 0, # <- avoid negative summation
                               is.na(rate_pred) ~0, # <- zeros where model yields NA
                               .default = rate_pred), 
         rate_pred = 100*rate_pred) |>  # as rate_dev equals 1/days_dev
  group_by(year, model_name) |> 
  mutate(rate_sum = cumsum(rate_pred), #heat summation
         emergence_lgl = logic_ratesum(rate_sum), #threshold point
         doy = yday(date)) |> 
  summarise(day_of_emergence = sum(emergence_lgl))  #number of days until threshold is reach

rate_sum_preds_daily

###### b) linear DDs ----
lm_ratedev <- lm(dev.rate ~ temp,
                 data = gilbert_pieris_data) #or gilbert pieris data


dd_preds_daily <- daily_temperatures |> 
  mutate(model_name = "linear") |> 
  group_by(year, model_name) |> 
  mutate(rate_pred = map_dbl(.x = daily_tavg,
                             .f = ~lm_pred(lm_fit = lm_ratedev,
                                           .x))) |> 
  mutate(rate_pred = if_else(rate_pred <0,
                             0,
                             rate_pred)) |>
  mutate(rate_sum = map_dbl(.x = daily_tavg,
                            .f = ~calc_dd(lm_fit = lm_ratedev,
                                          tavg = .x)), #heat summation
         rate_sum = cumsum(rate_sum),
         emergence_lgl = map_dbl(.x = rate_sum,
                                 .f = ~logic_dd(heat_units = 1/coef(lm_ratedev)[2],
                                                rate_cumsum = .x)), #threshold point
         doy = yday(date)) |> 
  summarise(day_of_emergence = sum(emergence_lgl))

###### c) combined  ----
doe_preds_daily <- rate_sum_preds_daily |> 
  full_join(dd_preds_daily) |> 
  mutate(time_slot = "daily")

##### 3.2. Hourly res. ----------------------------------------------------------
##  we use chillR package to simulate temperature variation across hours within days.
###### a) nonlinear rate summation ----

daily_temps_doy <- daily_temperatures  |> 
  mutate(JDay = yday(date)) |>  
  rename(Tmin = daily_tmin,
         Tmax = daily_tmax) 

hourly_temp_doy <- chillR::make_hourly_temps(latitude = 40.5,
                                             year_file = daily_temps_doy,
                                             keep_sunrise_sunset = TRUE )  |> 
  pivot_longer(cols = 10:33,
               names_to = "hour_of_day",
               values_to = "temperature") |>  
  rename(tmax = Tmax,
         tmin = Tmin,
         doy = JDay,
         sunrise = Sunrise,
         sunset = Sunset,
         daylength = Daylength) 


hourly_rate_preds <- tibble()

for(model_i in unique(selected_models_pieris_rapae$model_name)) {
  rate_hourly_i <- rate_pred(model2ratesum = model_i,
                             temperature = hourly_temp_doy,
                             fitted_params = selected_models_pieris_rapae,
                             res = "hourly")
  hourly_rate_preds <- bind_rows(hourly_rate_preds, rate_hourly_i)
}
hourly_rate_preds

rate_sum_preds_hourly <- hourly_rate_preds |> 
  mutate(rate_pred = case_when(rate_pred < 0 ~ 0, # <- avoid negative summation
                               is.na(rate_pred) ~0, # <- zeros where model yields NA
                               .default = rate_pred), 
         rate_pred = rate_pred*100/24) |>  # 24h a day
  group_by(year, model_name) |> 
  mutate(rate_sum = cumsum(rate_pred), #heat summation
         emergence_lgl = logic_ratesum(rate_sum), #threshold point
         doy = yday(date)) |> 
  summarise(day_of_emergence = round(sum(emergence_lgl)/24))  #number of days until threshold is reach

rate_sum_preds_hourly

###### b) linear DDs ----
dd_preds_hourly <- hourly_temp_doy |> 
  mutate(model_name = "linear") |> 
  group_by(year, model_name) |> 
  mutate(rate_pred = map_dbl(.x = temperature,
                             .f = ~lm_pred(lm_fit = lm_ratedev,
                                           .x))) |> 
  mutate(rate_pred = if_else(rate_pred <0,
                             0,
                             rate_pred)) |> #fraction of 1/24 hours
  mutate(rate_sum = map_dbl(.x = temperature,
                            .f = ~calc_dd(lm_fit = lm_ratedev,
                                          tavg = .x)), #heat summation
         rate_sum = cumsum(rate_sum),
         emergence_lgl = map_dbl(.x = rate_sum,
                                 .f = ~logic_dd(heat_units = 24*1/coef(lm_ratedev)[2],
                                                rate_cumsum = .x)), #threshold point
         doy = yday(date)) |> 
  summarise(day_of_emergence = round(sum(emergence_lgl)/24))

###### c) combined  ----
doe_preds_hourly <- rate_sum_preds_hourly |> 
  full_join(dd_preds_hourly) |> 
  mutate(time_slot = "hourly")

### 3.3. Validation ---------------------------------------
###### a) Visual comparison ----

# incorporate observational data from Gordo & Sanz (2006)
pieris_obs_trends <- read_csv(here("data/gordo_sanz_2006_doe.csv")) |> 
  rename(day_of_emergence = doe) |> 
  relocate(year, day_of_emergence) |> 
  group_by(year) |> 
  slice(1)

rate_sum_validation <- doe_preds_daily |> 
  bind_rows(doe_preds_hourly)

###### b) RMSE ----
preds_vs_obs_trends <- rate_sum_validation  |>
  filter(year > 1951 & year < 2005) |> 
  group_by(model_name, time_slot) |> 
  mutate(observed_emergence = pieris_obs_trends$day_of_emergence) |> 
  summarise(rmsep = chillR::RMSEP(predicted = day_of_emergence,
                                  observed = observed_emergence)) |>  
  arrange(rmsep) %>% 
  print() 

## see which temperature temporal resolution predicts better
transfer_rank_resolution <- preds_vs_obs_trends %>% 
  group_by(time_slot) %>% 
  summarise(rmsep = mean(rmsep)) %>% 
  arrange(rmsep) %>% 
  print()

## see which models predict better
transfer_rank_model <- preds_vs_obs_trends %>% 
  group_by(model_name) %>% 
  summarise(rmsep = mean(rmsep)) %>% 
  arrange(rmsep) %>% 
  print() # 

# Make Figure 1

preds_obs_compare_rmsep <- rate_sum_validation |> 
  full_join(preds_vs_obs_trends) |> 
  group_by(year, model_name, time_slot) |> 
  summarise(rmsep = min(rmsep)) |>  
  arrange(rmsep) |> 
  print()
rmsep_text <-  preds_obs_compare_rmsep |>  
  group_by(model_name, time_slot) |>  
  summarise(rmsep = mean(rmsep)) |> 
  arrange(rmsep)
rmsep_order <- rmsep_text |> 
  pull(model_name)
rmsep_order_min <- rmsep_text |> 
  slice_min(rmsep) |>
  arrange(rmsep) |> 
  pull(model_name)
rmsep_values <- rmsep_text |> 
  mutate(rmsep =  case_when(time_slot == "daily" ~paste("d-RMSEP =",round(rmsep,2)),
                            time_slot == "hourly" ~paste("h-RMSEP =",round(rmsep,2))),
         year = 2004,
         day_of_emergence = case_when(time_slot == "daily" ~165,
                                      time_slot == "hourly" ~145))
preds_obs_compare_plot <- ggplot()+
  geom_point(data = rate_sum_validation |> 
               filter(year %in% pieris_obs_trends$year),
             aes(x = year, y = day_of_emergence, color = time_slot),
             alpha = .66)+
  geom_line(data = rate_sum_validation |> 
              filter(year %in% pieris_obs_trends$year),
            aes(x = year, y = day_of_emergence, color = time_slot),
            linetype = "longdash",
            alpha = .66)+
  geom_smooth(data =  rate_sum_validation |> 
                filter(year %in% pieris_obs_trends$year),
              aes(x = year, y = day_of_emergence, 
                  color = time_slot, 
                  fill = time_slot))+
  scale_color_manual(values = c("#3A848B", "#E85038"))+
  scale_fill_manual(values = c("#3A848B", "#E85038"))+
  geom_point(data = pieris_obs_trends, 
             aes(x = year, y = day_of_emergence),
             color = "#E9C86B",
             fill = "#E9C86B",
             alpha = .66)+
  geom_line(data = pieris_obs_trends, 
            aes(x = year, y = day_of_emergence),
            color = "#E9C86B",
            linetype = "longdash",
            alpha = .66)+
  geom_smooth(data = pieris_obs_trends, 
              aes(x = year, y = day_of_emergence),
              color = "#E9C86B",
              fill = "#E9C86B")+
  labs(x = "Year",
       y = "Day of first adult emergence",
       color = "Resolution (model)",
       fill = "Resolution (model)")+
  facet_wrap(~factor(model_name, levels = rmsep_order_min))+
  theme_few()+
  theme(strip.text.x = element_text(size = 10, 
                                    color = "grey42", 
                                    face = "bold"))+
  geom_text(data = rmsep_values, 
            aes(label = rmsep,
                x = year,
                y = day_of_emergence), 
            size = 3,
            fontface = "italic",hjust = "right")
preds_obs_compare_plot  

## save for schmalensee
ggsave(filename = here("figures/schmalensee_modelpreds_validate_rmsep.png"),
       width = 25,
       height = 25,
       units = "cm")
ggsave(filename = here("figures/schmalensee_modelpreds_validate_rmsep.svg"),
       width = 25,
       height = 25,
       units = "cm")

# 4. Thermal regimes ---------------------------------------

##only for Gilbert's data
daily_temperatures_doy <- daily_temperatures |> 
  filter(year %in% 1952:2004) |> 
  mutate(doy = yday(date)) |> 
  rename(tmin = daily_tmin,
         tmax = daily_tmax,
         tavg = daily_tavg) |> 
  pivot_longer(cols = c(2, 3, 5),
               names_to = "temp_var",
               values_to = "daily_temperature")

ggplot(daily_temperatures_doy, aes(x = doy, y = daily_temperature)) +
  geom_point(aes(color = temp_var), alpha = .1)+
  scale_color_manual(values = c("#F1CB56", "#E36944", "#7FD3C3"))+
  ggdark::dark_theme_minimal()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")+
  labs(x = "Day of Year",
       y = "Daily temperature",
       color = "Temperature variables")

ggsave(here("figures/daily_temperature_regimes.png"),
       width = 2400,
       height = 2400,
       units = "px")
ggsave(here("figures/daily_temperature_regimes.svg"),
       width = 2400,
       height = 2400,
       units = "px") 

daily_temperatures_doy_decades <- daily_temperatures_doy |> 
  mutate(decade = case_when(year %in% c(1950:1959) ~ "1950's",
                            year %in% c(1960:1969) ~ "1960's",
                            year %in% c(1970:1979) ~ "1970's",
                            year %in% c(1980:1989) ~ "1980's",
                            year %in% c(1990:1999) ~ "1990's",
                            year %in% c(2000:2009) ~ "2000's",
                            year %in% c(2010:2015) ~ "2010's"
  )
  ) |> 
  group_by(decade) |> 
  mutate(decade_count_day = rep(1:(n()/3), each = 3)) |>  # Assign the same value to all rows for the same doy
  ungroup()

ggplot(daily_temperatures_doy_decades, aes(x = decade_count_day, y = daily_temperature)) +
  #geom_point(aes(color = temp_var), alpha =.3)+
  geom_line(aes(color = temp_var),
            linewidth = .5,
            alpha = 1)+
  scale_color_manual(values = c("#F1CB56", "#E36944", "#7FD3C3"))+
  scale_x_reverse()+
  ggdark::dark_theme_minimal()+
  labs(x = "Days",
       y = "Daily temperature",
       color = "Temperature variables")+
  coord_flip()+
  facet_wrap(~decade,
             scales = "free_x",
             nrow = 1
  )+
  theme(legend.position = "none")

ggsave(here("figures/daily_temperature_regimes_decades.png"),
       width = 2500,
       height = 4000,
       units = "px")

ggsave(here("figures/daily_temperature_regimes_decades.svg"),
       width = 2500,
       height = 4000,
       units = "px")

## and plot oneill

plot_devmodels(temp = gilbert_pieris_data$temp,
               dev_rate = gilbert_pieris_data$dev.rate,
               fitted_parameters = selected_models_pieris_rapae |> 
                 filter(model_name == "oneill"),
               species = "Pieris rapae",
               life_stage = "Pupa")+
  labs(title = NULL,
       subtitle = NULL,
       y = NULL,
       x = NULL)+
  scale_fill_manual(values = "#f39e71ff")+
  scale_color_manual(values = "#f39e71ff")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_blank())+
  facet_null()

ggsave(here("figures/oneill.svg"),
       width = 600,
       height = 600,
       units = "px")



