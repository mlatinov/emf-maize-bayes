
## Workflow Orchestration 
library(targets)
library(tidyverse)
library(tarchetypes)
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
  ## Hierarchical cell-means, pots pooled within treatment ##
  tar_target(
    name = m5_recovery_data,
    command = maize_h_m5_dgp(n = 700)
  ),
  tar_target(
    name = m5_recovery_check,
    command = m5_stan_model(m5_recovery_data$data, prior_only = 0)
  ),
  tar_target(
    name = m5_maize_model,
    command = m5_stan_model(model_data_maize_heigh, prior_only = 0)
  ),
  tar_target(
    name = m5_maize_model_pp,
    command = m5_stan_model(model_data_maize_heigh, prior_only = 1)
  ),
  ## Hierarchical Gompertz, per-pot curves pooled within treatment ##
  tar_target(
    name = m6_recovery_data,
    command = maize_h_m6_dgp(n = 700)
  ),
  tar_target(
    name = m6_recovery_check,
    command = m6_stan_model(m6_recovery_data$data, prior_only = 0)
  ),
  tar_target(
    name = m6_maize_model,
    command = m6_stan_model(model_data_maize_heigh, prior_only = 0)
  ),
  tar_target(
    name = m6_maize_model_pp,
    command = m6_stan_model(model_data_maize_heigh, prior_only = 1)
  ), 
  ## Quatro Report 
  tar_quarto(
    name = maize_height_quatro_report,
    path = "analysis/maize_height.qmd",
    quiet = FALSE 
  )
)