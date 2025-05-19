library(tidyverse)
library(here)
library(lubridate)
library(leaflet)
library(leaflet.extras2)
library(RColorBrewer)
library(terra)
library(sf)
library(easyclimate)
source(here("R/1-functions_phenobraspests.R"))

# 1. CBMS data wrangling ----

## load pieris rapae only data counting
piepra_raw_count_data_csv <- read_delim("~/dev_rapae/data_source/peticio_sansegundo.csv") |> 
  mutate(Nindiv = if_else(Nindiv == 0.5, # <- misread with R
                          1,
                          Nindiv))

## Regular data set

piepra_raw_count_data <- piepra_raw_count_data_csv |> 
  rename(id_transect = IDitin,
         id_species = IDesp,
         year = Any,
         date = Date,
         n_indiv = Nindiv,
         length = longitud,
         lat = LAT,
         lon = LNG) |> 
  select(-year) |> 
  mutate(date = lubridate::dmy(date),
         year = lubridate::year(date),
         id_transect = as_factor(id_transect)) |> 
  relocate(id_transect, id_species, year, date, n_indiv, length, lat, lon) |> 
  glimpse()

## Combined data set
cbms_all_transects <- read_delim("~/dev_rapae/data_source/m_visit_sub.csv") |> 
 rename(id_transect = SITE_ID,
        date = DATE)|> 
  select(id_transect, date) |> 
  mutate(id_transect = as_factor(id_transect),
         year = year(date))

pieris_pres_abs <- full_join(piepra_raw_count_data,
                             cbms_all_transects) |>
  group_by(id_transect) |> 
  mutate(pres_abs = map_chr(.x = n_indiv,
                            ~if_else(is.na(.x),
                                     "absent",
                                     "present"))) |> 
  mutate(n_indiv = map2_dbl(.x = pres_abs,
                            .y = n_indiv,
                            .f = ~if_else(.x == "absent",
                                          0,
                                          .y))) |> 
  mutate(doy = yday(date))

## c) Data wrangling -----------------------------------------------------
## examine trends

## now we will apply exclusion criteria:
##  
## (1) we'll exclude complete transects with >50% data absences than presences.
## since populations at these places may be transient or population dynamics may
## reflect other processes rather than heat accumulation.

## (2) we'll exclude anomalous years for each transect, i.e., those
##  with late first appearances (i.e., > doy150) with no previous or very early absences
##  followed by no data for months that may indicate lack of sampling rather than lack of presence.
##  Similarly, we'll exclude years within transects having exceptionally low abundances 
##  when more  than ten year data have been collected (i.e., years with extraordinary sitations that may 
##  indicate noise, rather than a phenological response to climate (see eg. transect 68 year 2012)

## filter by presence/absence ratio
pieris_cbms_pres <- pieris_pres_abs |> 
  group_by(id_transect) |> 
  count(pres_abs) |>
  mutate(total_counts_transect = sum(n),
         prop_pres = n/total_counts_transect) |> 
  filter(pres_abs == "present" & prop_pres > 0.5) 


sites_selection <- unique(pieris_cbms_pres$id_transect)

#individual plots 
for(transect_i in sites_selection) {
  pieris_selected_sites_i <- pieris_pres_abs |> 
    filter(id_transect == transect_i)
  ggplot_pieris_i <- ggplot(data = pieris_selected_sites_i, 
                            aes(x = doy,
                                y = n_indiv,
                                color = pres_abs)
  )+
    geom_point()+
    geom_line()+
    labs(x = "Day of Year",
         y = "N individuals observed",
         color = NULL,
         title = paste("Transect", transect_i))+
    facet_wrap(.~year)+
    theme_bw()
  ggsave(plot = ggplot_pieris_i, 
         filename = here(paste0("figures/cbms_transects_pres_abs/cbms_pieris_count_", transect_i,".png")),
         width = 2600,
         height = 2600,
         units = "px")
  
}

## filtering the data set:
pieris_cbms_selection <- pieris_pres_abs |> 
  filter(id_transect %in% sites_selection) |> 
  filter(!(id_transect == 17 & year == 2000),
         !(id_transect == 29 & year == 1997),    
         !(id_transect == 47 & year == 2020),
         !(id_transect == 68 & year == 2012),
         !(id_transect == 89 & year == 2020), #likely COVID-19 restrictions
         !(id_transect == 90 & year %in% c(2007, 2021)),
         !(id_transect == 94 & year == 2020), #likely COVID-19 restrictions
         !(id_transect == 106 & year == 2020), #likely COVID-19 restrictions
         !(id_transect == 107 & year == 2020), #likely COVID-19 restrictions
         !(id_transect == 110 & year == 2020), #likely COVID-19 restrictions
         !(id_transect == 114 & year == 2020), #likely COVID-19 restrictions
         !(id_transect == 116 & year == 2020), #likely COVID-19 restrictions
         !(id_transect == 117 & year == 2020)) #likely COVID-19 restrictions

ids_cbms_selection <- unique(pieris_cbms_selection$id_transect)

sites <- pieris_cbms_selection |> 
  rename(ID_coords = id_transect)


first_emergence_adult_sites <- sites |> 
  filter(pres_abs == "present") |> 
  group_by(ID_coords, year) |> 
  slice_min(date) |> 
  mutate(doy = yday(date)) 

years_sites <- sort(unique(first_emergence_adult_sites$year))

sites_coords <- first_emergence_adult_sites |> 
  group_by(ID_coords) |> 
  slice(1) |> 
  ungroup() |> 
  dplyr::select(lat, lon)

easyclimate_temp_transects <- easyclimate::get_daily_climate(coords = sites_coords,
                                                             climatic_var = c("Tmin", "Tmax"),
                                                             period = 1988:2022,
                                                             output = "df",
                                                             version = 4) |> 
  as_tibble() |> 
  mutate(date = ymd(date))


easyclimate_temp_transects

daily_temperatures <- easyclimate_temp_transects |> 
  mutate(daily_tavg = map2_dbl(.x = Tmin,
                               .y = Tmax,
                               .f = ~mean(c(.x, .y), na.rm = TRUE))) |> 
  mutate(year = year(date)) 
  

saveRDS(daily_temperatures, file = here::here("data/easyclimate_data_extracted.rds"))

##re-arrange pieris_cbms_selection to ensure correct joining

pieris_cbms_uniques <- pieris_cbms_selection |> 
  group_by(id_transect) |> 
  mutate(ID_coords = cur_group_id()) |> 
  ungroup() |> 
  select(-id_transect)


pieris_cbms_selection_arranged <- pieris_cbms_selection |> 
  mutate(id_transect = as.numeric(id_transect)) |> 
  arrange(id_transect) |> 
  group_by(id_transect) |> 
  mutate(ID_coords = cur_group_id()) |> 
  ungroup() |> 
  select(-id_transect)
  
daily_temperatures <- read_rds(here::here("data/easyclimate_data_extracted.rds")) |> 
  filter(!ID_coords %in% c(13, 132, 154)) |>  # <- no climate data for these points with easyclimate
  #remove also those transects with unreliable data
  filter(!(ID_coords == 17 & year == 2000),
         !(ID_coords == 29 & year == 1997),    
         !(ID_coords == 47 & year == 2020),
         !(ID_coords == 68 & year == 2012),
         !(ID_coords == 89 & year == 2020), #likely COVID-19 restrictions
         !(ID_coords == 90 & year %in% c(2007, 2021)),
         !(ID_coords == 94 & year == 2020), #likely COVID-19 restrictions
         !(ID_coords == 106 & year == 2020), #likely COVID-19 restrictions
         !(ID_coords == 107 & year == 2020), #likely COVID-19 restrictions
         !(ID_coords == 110 & year == 2020), #likely COVID-19 restrictions
         !(ID_coords == 114 & year == 2020), #likely COVID-19 restrictions
         !(ID_coords == 116 & year == 2020), #likely COVID-19 restrictions
         !(ID_coords == 117 & year == 2020)) #likely COVID-19 restrictions

daily_temperatures_vert <- pivot_longer(daily_temperatures,
                                        cols = c("Tmin", "Tmax"),
                                        names_to = "var",
                                        values_to = "temperature")

cbms_daily_temperatures_rainplot <- ggplot(data = daily_temperatures_vert,
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

cbms_daily_temperatures_rainplot
ggsave(plot = cbms_daily_temperatures_rainplot,
       filename = here("figures/cbms_daily_temperatures_rainplot.png"),
       width = 2600, height = 2600, units = "px")

#2. Gilbert Predictions -----------------------------------------------------
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

###### a) daily  ----
daily_rate_preds_cbms <- tibble()

  for(model_i in unique(selected_models_pieris_rapae$model_name)) {
    print(paste("Fitting model", model_i))
  
  rate_daily_i <- rate_pred(model2ratesum = model_i,
                            temperature = daily_temperatures,
                            fitted_params = selected_models_pieris_rapae,
                            res = "daily")
  daily_rate_preds_cbms <- bind_rows(daily_rate_preds_cbms, rate_daily_i)
}

rate_sum_preds_daily <- daily_rate_preds_cbms |> 
  mutate(rate_pred = case_when(rate_pred < 0 ~ 0, # <- avoid negative summation
                               is.na(rate_pred) ~0, # <- zeros where model yields NA
                               .default = rate_pred), 
         rate_pred = 100*rate_pred) |>  # as rate_dev equals 1/days_dev
  group_by(year, model_name, ID_coords) |> 
  mutate(rate_sum = cumsum(rate_pred), #heat summation
         emergence_lgl = logic_ratesum(rate_sum), #threshold point
         doy = yday(date)) |> 
  summarise(day_of_emergence = sum(emergence_lgl))  #number of days until threshold is reach

lm_ratedev <- lm(dev.rate ~ temp,
                 data = gilbert_pieris_data) 

dd_preds_daily <- daily_temperatures |> 
  mutate(model_name = "linear") |> 
  group_by(year, model_name, ID_coords) |> 
  mutate(rate_pred = lm_pred(lm_fit = lm_ratedev,
                             daily_tavg)) |> 
  mutate(rate_pred = if_else(rate_pred <0,
                             0,
                             rate_pred)) |>
  mutate(rate_sum = calc_dd(lm_fit = lm_ratedev,
                            daily_tavg), #heat summation
         rate_sum = cumsum(rate_sum),
         emergence_lgl = logic_dd(heat_units = 1/coef(lm_ratedev)[2],
                                  rate_cumsum = rate_sum), #threshold point
         doy = yday(date)) |> 
  summarise(day_of_emergence = sum(emergence_lgl))

daily_preds_ratesum_cbms <- rate_sum_preds_daily |> 
  full_join(dd_preds_daily) |> 
  mutate(time_slot = "daily")


###### b) hourly  ----
daily_temps_doy <- daily_temperatures  |> 
  mutate(JDay = yday(date))

hourly_temps_cbms <- tibble()

for(transect_i in unique(daily_temps_doy$ID_coords)) {
  print(paste0(transect_i, "/", length(unique(daily_temps_doy$ID_coords))))
  daily_temps_doy_i <- daily_temps_doy |> 
    filter(ID_coords == transect_i)
  lat_i <- unique(daily_temps_doy_i$lat)
  hourly_temps_i <- chillR::make_hourly_temps(latitude = lat_i,
                                              year_file = daily_temps_doy_i,
                                              keep_sunrise_sunset = TRUE) |> 
    pivot_longer(cols = 13:36,
                 names_to = "hour_of_day",
                 values_to = "temperature") |>  
    rename(tmax = Tmax,
           tmin = Tmin,
           sunrise = Sunrise,
           sunset = Sunset,
           daylength = Daylength,
           doy = JDay
           ) 
  
  hourly_temps_cbms <- bind_rows(hourly_temps_cbms, hourly_temps_i) 
}
  
hourly_temps_cbms

hourly_rate_preds <- tibble()

for(model_i in unique(selected_models_pieris_rapae$model_name)) {
  print(paste0("Predicting phenology from TPC model ", model_i))
  rate_hourly_i <- rate_pred(model2ratesum = model_i,
                             temperature = hourly_temps_cbms,
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
  group_by(year, model_name, ID_coords) |> 
  mutate(rate_sum = cumsum(rate_pred), #heat summation
         emergence_lgl = logic_ratesum(rate_sum), #threshold point
         doy = yday(date)) |> 
  summarise(day_of_emergence = round(sum(emergence_lgl)/24))  #number of days until threshold is reach

dd_preds_hourly <- hourly_temps_cbms |> 
  mutate(model_name = "linear") |> 
  group_by(year, model_name, ID_coords) |>
  mutate(rate_sum = calc_dd(lm_fit = lm_ratedev,
                            tavg = temperature), #heat summation
         rate_sum = cumsum(rate_sum),
         emergence_lgl = logic_dd(heat_units = 24*1/coef(lm_ratedev)[2],
                                  rate_cumsum = rate_sum), #threshold point
         doy = yday(date)) |> 
  summarise(day_of_emergence = round(sum(emergence_lgl)/24))


hourly_preds_ratesum_cbms <- rate_sum_preds_hourly |> 
  bind_rows(dd_preds_hourly) |> 
  mutate(time_slot = "hourly")

###### c) combine resolutions  ----
preds_ratesum_cbms <- bind_rows(daily_preds_ratesum_cbms,
                                hourly_preds_ratesum_cbms)

write_rds(preds_ratesum_cbms, file = here("data/preds_ratesum_cbms.rds"))
#gc()

preds_ratesum_cbms <- readRDS(here("data/preds_ratesum_cbms.rds"))
## and observations:

first_emergence_adult_sites <- pieris_cbms_selection_arranged |> 
  filter(!ID_coords %in% c(13, 132, 154)) |> 
  filter(pres_abs == "present") |> 
  group_by(ID_coords, year) |> 
  slice_min(date) |> 
  mutate(doy = yday(date)) |>
  ungroup() |> 
  group_by(ID_coords) |> 
  mutate(n_years = n()) |> 
  filter(n_years >= 10)

joined_cbms_preds_gilbert <- inner_join(first_emergence_adult_sites,
                                        preds_ratesum_cbms)

#ggplot(first_emergence_adult_sites, aes(x = year, y = doy))+
 # geom_point()+
  # geom_line()

### d) plot no-pooling ----
joined_cbms_preds <- joined_cbms_preds_gilbert |> 
  mutate(rmsep = map2_dbl(.x = doy, .y = day_of_emergence,
                          .f = ~chillR::RMSEP(predicted = .y,
                                              observed = .x))) |> 
  ungroup() 
  
transfer_rank_resolution <- joined_cbms_preds |> 
  group_by(time_slot) |> 
  summarise(rmsep = mean(rmsep)) |>  
  arrange(rmsep) |> 
  print()# overall hourly rmsep is lower

transfer_rank_model <- joined_cbms_preds |> 
  group_by(model_name) |>  
  summarise(rmsep = mean(rmsep)) |> 
  arrange(rmsep) |> 
  print() # overall mod_weibull, beta, flextpc, schoolfield and oneill

#re-plot the graph (manually) annotating rmsep

preds_obs_compare_rmsep <- joined_cbms_preds |> 
  group_by(model_name, time_slot) |> 
  summarise(rmsep = mean(rmsep)) |>  
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
         year = 2022,
         day_of_emergence = case_when(time_slot == "daily" ~170,
                                      time_slot == "hourly" ~155))

preds_obs_compare_plot <- ggplot(joined_cbms_preds,
                                 aes(x = year, 
                                     y = day_of_emergence,
                                 ))+
  geom_point(aes(color = time_slot), 
             alpha = .1)+
  geom_smooth(aes(color = time_slot,
                  fill = time_slot))+
  scale_fill_manual(values = c("#3A848B", "#E85038"))+
  scale_color_manual(values = c("#3A848B", "#E85038"))+
  ggthemes::theme_clean()+
  geom_point(data = first_emergence_adult_sites,
             aes(x = year,
                 y = doy),
             color = "#E9C86B",
             alpha = .1)+
  geom_smooth(data = first_emergence_adult_sites,
              aes(x = year, y  = doy),
              color = "#E9C86B",
              fill = "#E9C86B")+
  labs(x = "Year",
       y = "Day of first adult emergence")+
  facet_wrap(~factor(model_name, levels = rmsep_order_min))+
  ggthemes::theme_few()+
  theme(strip.text.x = element_text(size = 10, 
                                    color = "grey42", 
                                    face = "bold"))+
  geom_text(data = rmsep_values, 
            aes(label = rmsep,
                x = year,
                y = day_of_emergence), 
            size = 3,
            fontface = "italic",
            hjust = "right")
preds_obs_compare_plot  

ggsave(filename = here("figures/preds_obs_compare_plot_rmsep_all.png"),
       width = 2600, height = 2600, units = "px")

ggsave(filename = here("figures/preds_obs_compare_plot_rmsep_all.svg"),
       units = "cm",
       width = 25,
       height = 25)


##### e) plot with pooling -----------------------------------------------

transects_years <- first_emergence_adult_sites |> 
  group_by(year, ID_coords) |> 
  summarise(doy = mean(doy))


###  pooling (avg) across coordinates and then taking RMSEPs

cbms_preds_coords <- preds_ratesum_cbms |> 
  group_by(ID_coords, year, model_name, time_slot) |> 
  summarise(day_of_emergence = mean(day_of_emergence)) |> 
  ungroup() 

cbms_preds_obs_pooled <- transects_years |> 
  left_join(cbms_preds_coords) |> 
  group_by(year, model_name, time_slot) |> 
  summarise(day_of_emergence = mean(day_of_emergence),
            doy = mean(doy)) |> 
  tidyr::drop_na() |> 
  mutate(rmsep = map2_dbl(.x = day_of_emergence,
                          .y = doy,
                          .f = ~chillR::RMSEP(predicted = .x,
                                              observed = .y)))   

transfer_rank_resolution <- cbms_preds_obs_pooled |> 
  group_by(time_slot) |> 
  summarise(rmsep = mean(rmsep)) |>  
  arrange(rmsep) |> 
  print()# overall hourly rmsep is lower

transfer_rank_model <- cbms_preds_obs_pooled |> 
  group_by(model_name) |>  
  summarise(rmsep = mean(rmsep)) |> 
  arrange(rmsep) |> 
  print() # overall mod_weibull, beta, flextpc, schoolfield and oneill

#re-plot the graph (manually) annotating rmsep

preds_obs_compare_rmsep_pooled <- cbms_preds_obs_pooled |> 
  group_by(model_name, time_slot) |> 
  summarise(rmsep = mean(rmsep)) |>  
  arrange(rmsep) |> 
  print()

rmsep_text <-  preds_obs_compare_rmsep_pooled |>  
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
         year = 2022,
         day_of_emergence = case_when(time_slot == "daily" ~130,
                                      time_slot == "hourly" ~120))
preds_obs_compare_plot_pooled <- ggplot(cbms_preds_obs_pooled,
                                 aes(x = year, 
                                     y = day_of_emergence,
                                 ))+
  geom_point(aes(color = time_slot), alpha = .66)+
  geom_line(aes(color = time_slot),
            linetype = "longdash")+
  geom_smooth(aes(color = time_slot,
                  fill = time_slot))+
  scale_fill_manual(values = c("#3A848B", "#E85038"))+
  scale_color_manual(values = c("#3A848B", "#E85038"))+
  ggthemes::theme_clean()+
  geom_point(aes(x = year, y  = doy),
             color = "#E9C86B",
             alpha = .66)+
  geom_line(aes(x = year, y  = doy),
            color = "#E9C86B")+
  geom_smooth(aes(x = year, y  = doy),
              color = "#E9C86B",
              fill = "#E9C86B")+
  labs(x = "Year",
       y = "Day of first adult emergence")+
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
            fontface = "italic",
            hjust = "right")
preds_obs_compare_plot_pooled  

ggsave(filename = here("figures/preds_obs_compare_plot_rmsep_group.png"),
       width = 2600, height = 2600, units = "px")

ggsave(filename = here("figures/preds_obs_compare_plot_rmsep_group.svg"),
       width = 25,
       height = 25,
       units = "cm")



# 3. Schmalensee Predictions -----------------------------------------------------
selected_models_pieris_rapae <- readRDS(here("data/selected_models_schmalensee.rds"))
schmalensee_pieris_data <- read_delim("~/Data and scripts von Schmalensee et al. 2023/Data and scripts von Schmalensee et al. 2023/Data/full_data.txt", 
                                      delim = "\t") |> 
  filter(species == "rapae",
         life.stage == "pupa") |> 
  select(temp, dev.rate) |> 
  drop_na() 
###### a) daily  ----
daily_rate_preds_cbms <- tibble()

for(model_i in unique(selected_models_pieris_rapae$model_name)) {
  print(paste("Fitting model", model_i))
  
  rate_daily_i <- rate_pred(model2ratesum = model_i,
                            temperature = daily_temperatures,
                            fitted_params = selected_models_pieris_rapae,
                            res = "daily")
  daily_rate_preds_cbms <- bind_rows(daily_rate_preds_cbms, rate_daily_i)
}

rate_sum_preds_daily <- daily_rate_preds_cbms |> 
  mutate(rate_pred = case_when(rate_pred < 0 ~ 0, # <- avoid negative summation
                               is.na(rate_pred) ~0, # <- zeros where model yields NA
                               .default = rate_pred), 
         rate_pred = 100*rate_pred) |>  # as rate_dev equals 1/days_dev
  group_by(year, model_name, ID_coords) |> 
  mutate(rate_sum = cumsum(rate_pred), #heat summation
         emergence_lgl = logic_ratesum(rate_sum), #threshold point
         doy = yday(date)) |> 
  summarise(day_of_emergence = sum(emergence_lgl))  #number of days until threshold is reach

lm_ratedev <- lm(dev.rate ~ temp,
                 data = schmalensee_pieris_data) 

dd_preds_daily <- daily_temperatures |> 
  mutate(model_name = "linear") |> 
  group_by(year, model_name, ID_coords) |> 
  mutate(rate_pred = lm_pred(lm_fit = lm_ratedev,
                             daily_tavg)) |> 
  mutate(rate_pred = if_else(rate_pred <0,
                             0,
                             rate_pred)) |>
  mutate(rate_sum = calc_dd(lm_fit = lm_ratedev,
                            daily_tavg), #heat summation
         rate_sum = cumsum(rate_sum),
         emergence_lgl = logic_dd(heat_units = 1/coef(lm_ratedev)[2],
                                  rate_cumsum = rate_sum), #threshold point
         doy = yday(date)) |> 
  summarise(day_of_emergence = sum(emergence_lgl))

daily_preds_ratesum_cbms <- rate_sum_preds_daily |> 
  full_join(dd_preds_daily) |> 
  mutate(time_slot = "daily")

###### b) hourly  ----
daily_temps_doy <- daily_temperatures  |> 
  mutate(JDay = yday(date))

hourly_temps_cbms <- tibble()

for(transect_i in unique(daily_temps_doy$ID_coords)) {
  print(paste0(transect_i, "/", length(unique(daily_temps_doy$ID_coords))))
  daily_temps_doy_i <- daily_temps_doy |> 
    filter(ID_coords == transect_i)
  lat_i <- unique(daily_temps_doy_i$lat)
  hourly_temps_i <- chillR::make_hourly_temps(latitude = lat_i,
                                              year_file = daily_temps_doy_i,
                                              keep_sunrise_sunset = TRUE) |> 
    pivot_longer(cols = 13:36,
                 names_to = "hour_of_day",
                 values_to = "temperature") |>  
    rename(tmax = Tmax,
           tmin = Tmin,
           sunrise = Sunrise,
           sunset = Sunset,
           daylength = Daylength,
           doy = JDay
    ) 
  
  hourly_temps_cbms <- bind_rows(hourly_temps_cbms, hourly_temps_i) 
}

hourly_temps_cbms

hourly_rate_preds <- tibble()

for(model_i in unique(selected_models_pieris_rapae$model_name)) {
  print(paste0("Predicting phenology from TPC model ", model_i))
  rate_hourly_i <- rate_pred(model2ratesum = model_i,
                             temperature = hourly_temps_cbms,
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
  group_by(year, model_name, ID_coords) |> 
  mutate(rate_sum = cumsum(rate_pred), #heat summation
         emergence_lgl = logic_ratesum(rate_sum), #threshold point
         doy = yday(date)) |> 
  summarise(day_of_emergence = round(sum(emergence_lgl)/24))  #number of days until threshold is reach

dd_preds_hourly <- hourly_temps_cbms |> 
  mutate(model_name = "linear") |> 
  group_by(year, model_name, ID_coords) |>
  mutate(rate_sum = calc_dd(lm_fit = lm_ratedev,
                            tavg = temperature), #heat summation
         rate_sum = cumsum(rate_sum),
         emergence_lgl = logic_dd(heat_units = 24*1/coef(lm_ratedev)[2],
                                  rate_cumsum = rate_sum), #threshold point
         doy = yday(date)) |> 
  summarise(day_of_emergence = round(sum(emergence_lgl)/24))


hourly_preds_ratesum_cbms <- rate_sum_preds_hourly |> 
  bind_rows(dd_preds_hourly) |> 
  mutate(time_slot = "hourly")

###### c) combine resolutions  ----
preds_ratesum_cbms <- bind_rows(daily_preds_ratesum_cbms,
                                hourly_preds_ratesum_cbms)

write_rds(preds_ratesum_cbms, file = here("data/preds_ratesum_cbms_schmalensee.rds"))
 #gc()

preds_ratesum_cbms <- readRDS(here("data/preds_ratesum_cbms_schmalensee.rds"))
## and observations:

first_emergence_adult_sites <- pieris_cbms_selection_arranged |> 
  filter(!ID_coords %in% c(13, 132, 154)) |> 
  filter(pres_abs == "present") |> 
  group_by(ID_coords, year) |> 
  slice_min(date) |> 
  mutate(doy = yday(date)) |>
  ungroup() |> 
  group_by(ID_coords) |> 
  mutate(n_years = n()) |> 
  filter(n_years >= 10)


joined_cbms_preds_schmalensee <- inner_join(first_emergence_adult_sites,
                                        preds_ratesum_cbms)

#ggplot(first_emergence_adult_sites, aes(x = year, y = doy))+
# geom_point()+
# geom_line()

### d) plot no-pooling ----

joined_cbms_preds <- joined_cbms_preds_schmalensee |> 
  mutate(rmsep = map2_dbl(.x = doy, .y = day_of_emergence,
                          .f = ~chillR::RMSEP(predicte = .y,
                                              observed = .x))) |> 
  drop_na() |> 
  ungroup() 

transfer_rank_resolution <- joined_cbms_preds |> 
  group_by(time_slot) |> 
  summarise(rmsep = mean(rmsep)) |>  
  arrange(rmsep) |> 
  print()# overall hourly rmsep is lower

transfer_rank_model <- joined_cbms_preds |> 
  group_by(model_name) |>  
  summarise(rmsep = mean(rmsep)) |> 
  arrange(rmsep) |> 
  print() # overall mod_weibull, beta, flextpc, schoolfield and oneill

#re-plot the graph (manually) annotating rmsep

preds_obs_compare_rmsep <- joined_cbms_preds |> 
  group_by(model_name, time_slot) |> 
  summarise(rmsep = mean(rmsep)) |>  
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
         year = 2022,
         day_of_emergence = case_when(time_slot == "daily" ~130,
                                      time_slot == "hourly" ~120))

preds_obs_compare_plot <- ggplot(joined_cbms_preds,
                                 aes(x = year, 
                                     y = day_of_emergence,
                                 ))+
  
  geom_smooth(aes(color = time_slot,
                  fill = time_slot))+
  scale_fill_manual(values = c("#3A848B", "#E85038"))+
  scale_color_manual(values = c("#3A848B", "#E85038"))+
  ggthemes::theme_clean()+
  geom_smooth(aes(x = year, y  = doy),
              color = "#E9C86B",
              fill = "#E9C86B")+
  labs(x = "Year",
       y = "Day of first adult emergence")+
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
            fontface = "italic",
            hjust = "right")+
  theme(legend.position = "none")
preds_obs_compare_plot  

ggsave(filename = here("figures/preds_obs_compare_plot_rmsep_all_schmalensee.png"),
       width = 2600, height = 2600, units = "px")

ggsave(filename = here("figures/preds_obs_compare_plot_rmsep_all_schmalensee.svg"),
       units = "cm",
       width = 25,
       height = 25)


##### e) plot with pooling -----------------------------------------------

## and now rmsep 
### refactor first_emergence_adult_sites (CAUTION!)
transects_years <- first_emergence_adult_sites |> 
  group_by(year, ID_coords) |> 
  summarise(doy = mean(doy))

###  pooling (avg) across coordinates and then taking RMSEPs

cbms_preds_coords <- preds_ratesum_cbms |> 
  group_by(ID_coords, year, model_name, time_slot) |> 
  summarise(day_of_emergence = mean(day_of_emergence)) |> 
  ungroup() 

cbms_preds_obs_pooled <- transects_years |> 
  left_join(cbms_preds_coords) |> 
  group_by(year, model_name, time_slot) |> 
  summarise(day_of_emergence = mean(day_of_emergence),
            doy = mean(doy)) |> 
  tidyr::drop_na() |> 
  mutate(rmsep = map2_dbl(.x = day_of_emergence,
                          .y = doy,
                          .f = ~chillR::RMSEP(predicted = .x,
                                              observed = .y)))   

transfer_rank_resolution <- cbms_preds_obs_pooled |> 
  group_by(time_slot) |> 
  summarise(rmsep = mean(rmsep)) |>  
  arrange(rmsep) |> 
  print()# overall hourly rmsep is lower

transfer_rank_model <- cbms_preds_obs_pooled |> 
  group_by(model_name) |>  
  summarise(rmsep = mean(rmsep)) |> 
  arrange(rmsep) |> 
  print() # overall mod_weibull, beta, flextpc, schoolfield and oneill

#re-plot the graph (manually) annotating rmsep

preds_obs_compare_rmsep_pooled <- cbms_preds_obs_pooled |> 
  group_by(model_name, time_slot) |> 
  summarise(rmsep = mean(rmsep)) |>  
  arrange(rmsep) |> 
  print()

rmsep_text <-  preds_obs_compare_rmsep_pooled |>  
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
         year = 2022,
         day_of_emergence = case_when(time_slot == "daily" ~130,
                                      time_slot == "hourly" ~120))
preds_obs_compare_plot_pooled <- ggplot(cbms_preds_obs_pooled,
                                        aes(x = year, 
                                            y = day_of_emergence,
                                        ))+
  geom_point(aes(color = time_slot), alpha = .66)+
  geom_line(aes(color = time_slot),
            linetype = "longdash")+
  geom_smooth(aes(color = time_slot,
                  fill = time_slot))+
  scale_fill_manual(values = c("#3A848B", "#E85038"))+
  scale_color_manual(values = c("#3A848B", "#E85038"))+
  ggthemes::theme_clean()+
  geom_point(aes(x = year, y  = doy),
             color = "#E9C86B",
             alpha = .66)+
  geom_line(aes(x = year, y  = doy),
            color = "#E9C86B")+
  geom_smooth(aes(x = year, y  = doy),
              color = "#E9C86B",
              fill = "#E9C86B")+
  labs(x = "Year",
       y = "Day of first adult emergence")+
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
            fontface = "italic",
            hjust = "right")+
  theme(legend.position = "none")
preds_obs_compare_plot_pooled  

ggsave(filename = here("figures/preds_obs_compare_plot_rmsep_group_schmalensee.png"),
       width = 25,
       height = 25,
       units = "cm")
ggsave(filename = here("figures/preds_obs_compare_plot_rmsep_group_schmalensee.svg"),
       width = 25,
       height = 25,
       units = "cm")

# 4. Thermal regimes ---------------------------------------
selected_models_pieris_rapae <- readRDS(here("data/selected_models_gilbert.rds"))

daily_temperatures_doy <- daily_temperatures |> 
  mutate(doy = yday(date)) |> 
  rename(tmin = Tmin,
         tmax = Tmax,
         tavg = daily_tavg) |> 
  group_by(year, doy) |> 
  summarise(tmin = mean(tmin),
            tmax = mean(tmax),
            daily_tavg = mean(tavg)) |> 
   pivot_longer(cols = c(3:5),
               names_to = "temp_var",
               values_to = "daily_temperature")


daily_temperatures_doy_decades <- daily_temperatures_doy |> 
  mutate(decade = case_when(year %in% c(1988:2004) ~ "first_half",
                            year %in% c(2005:2022) ~ "second_half",
                            )) |> 
  group_by(decade) |> 
  mutate(decade_count_day = rep(1:(n()/3), each = 3)) |>  # Assign the same value to all rows for the same doy
  ungroup()

ggplot(daily_temperatures_doy_decades, aes(x = decade_count_day, 
                                           y = daily_temperature)) +
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
       width = 1500,
       height = 4000,
       units = "px")

ggsave(here("figures/daily_temperature_regimes_decades.svg"),
       width = 1500,
       height = 4000,
       units = "px")

## and plot lrf

plot_devmodels(temp = gilbert_pieris_data$temp,
               dev_rate = gilbert_pieris_data$dev.rate,
               fitted_parameters = selected_models_pieris_rapae |> 
                 filter(model_name == "lrf"),
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

ggsave(here("figures/lrf.svg"),
       width = 600,
       height = 600,
       units = "px")


