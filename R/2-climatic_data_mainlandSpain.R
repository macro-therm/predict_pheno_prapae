#### Script Info ####

# GitHub repo: github.com/dario-ssm/PhenoBrassicaPests

#previous: load packages
library(raster)
library(here)
library(tidyverse)
library(lubridate)
source(here("Scripts/1-functions_phenodev_predict.R"))

#### 1. Data extraction ####
##### a) Tmax ####
## downloaded data from http://www.meteo.unican.es/datasets/spain02 on Jan-24th-2022
#  (now in 2023 it is not available anymore; might be available upon request to the authors)

Tmax_brick <- brick(x = here("/Spain02_v5.0_010reg_aa3d/Spain02_v5.0_DD_010reg_aa3d_tasmax.nc"))
Tmax_avg <- as.data.frame(na.omit(Tmax_brick[,])) #convert brick to data frame
t_Tmax_avg <- as.data.frame(t(Tmax_avg)) #transpose
Tmax_avg_date <- t_Tmax_avg |> #extract dates
  rownames() |>
  str_sub(2,11)|>
  ymd()
#add date to tiblle
tmax_df <- t_Tmax_avg |>
  bind_cols(date=Tmax_avg_date)|>
  mutate(year = year(date))

#make tibble more easy to use
tmax_daily_avg_spain <- tmax_df |>
  pivot_longer(cols = starts_with("V")) #expand vertically
# summarise mean for each day across all cells
daily_tmax_df <- tmax_daily_avg_spain |>
  group_by(date)|> #grouping variable is "date"
  summarise(daily_tmax=mean(value)) |> 
  mutate(date=ymd(date))
save(daily_tmax_df, file = here("Data/daily_tmax_df.RData"))

##### b) Tmin ####
## downloaded data from http://www.meteo.unican.es/datasets/spain02 on Jan-24th-2022
#  (now in 2023 it is not available anymore)
Tmin_brick <-brick(x = here("/Spain02_v5.0_010reg_aa3d/Spain02_v5.0_DD_010reg_aa3d_tasmin.nc"))
Tmin_avg <- as.data.frame(na.omit(Tmin_brick[,])) #convert brick to data frame
t_Tmin_avg <- as.data.frame(t(Tmin_avg)) #transpose
Tmin_avg_date <- t_Tmin_avg |> #extract dates
  rownames() |>
  str_sub(2,11)|>
  ymd()
#add date to tiblle
tmin_df <- t_Tmin_avg |>
  bind_cols(date=Tmin_avg_date)|>
  mutate(year = year(date))

#make tibble more easy to use
tmin_daily_avg_spain <- tmin_df |>
  pivot_longer(cols = starts_with("V")) #expand vertically

# summarise mean for each day across all cells
daily_tmin_df <- tmin_daily_avg_spain |>
  group_by(date)|> #grouping variable is "date"
  summarise(daily_tmin=mean(value)) |> 
  mutate(date=ymd(date))
save(daily_tmin_df, file = here("Data/daily_tmin_df.RData"))

