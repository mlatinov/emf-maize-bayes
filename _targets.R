
## Workflow Orchestration 
library(targets)
library(tidyverse)
library(cmdstanr)

## Soruce Functions ##
tar_source("R/clean_data_.R")

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
  )
)