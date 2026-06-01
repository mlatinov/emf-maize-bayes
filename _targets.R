
## Workflow Orchestration 
library(targets)
library(tidyverse)
library(cmdstanr)

## Soruce Functions ##
tar_source("R/clean_data_.R")
tar_source("R/maize_height/")

list(
  # Load the Table 1 
  tar_target(
    name = maize_heigh_raw,
    command = readxl::read_excel("data/Table S1 - 2026 Paunov et al. - 868 MHz EMF Maize - Plant Height.xlsx",sheet = 1)
  ),
  # Clean Table 1
  tar_target(
    name = model_data_maize_heigh,
    command = clean_meize_raw(maize_heigh_raw)
  ),
  # Fit the Stan model Plant Height ~ Treatment Model M1
  tar_target(
    name = m1_maize_model,
    command = m1_stan_model(model_data_maize_heigh, prior_only = 0)
  ),
  tar_target(
    name = m1_maize_model_pp,
    command = m1_stan_model(model_data_maize_heigh, prior_only = 1)
  ),
  # Fit the Stan model Plant Height ~ Treatment x Day Model M2
  tar_target(
    name = m2_maize_model,
    command = m2_stan_model(model_data_maize_heigh, prior_only = 0)
  ),
  tar_target(
    name = m2_maize_model_pp,
    command = m2_stan_model(model_data_maize_heigh, prior_only = 1)
  ),
  # Fit the Stan model Plant Height ~ Pot x Day Model M3
  tar_target(
    name = m3_maize_model,
    command = m3_stan_model(model_data_maize_heigh, prior_only = 0)
  ),
  tar_target(
    name = m3_maize_model_pp,
    command = m3_stan_model(model_data_maize_heigh, prior_only = 1)
  ),
  ## M4 Gompertz per treatment ##
  tar_target(
    name = m4_maize_model,
    command = m4_stan_model(model_data_maize_heigh, prior_only = 0)
  ),
  tar_target(
    name = m4_maize_model_pp,
    command = m4_stan_model(model_data_maize_heigh, prior_only = 1)
  )
)