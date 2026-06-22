
## Workflow Orchestration 
library(targets)
library(tidyverse)
library(tarchetypes)
library(cmdstanr)

## Soruce Functions ##
tar_source("R/clean_data_.R")
tar_source("R/maize_height/")
tar_source("R/maize_jip/")
tar_source("tar_factory.R")

list(
  #### Maize Height Estimand Pipeline ####
  maize_height_pipeline(),

  # Load the JIP data #
  tar_target(
    name = data_maize_jip,
    command = readxl::read_excel(
      path = "data/Table S5 - 2026 Paunov et al. - 868 MHz EMF Maize - Photosynthetic Performance (JIP-test).xlsx",
      sheet = 4
    )
  ),
  # CLean JIP data #
  tar_target(
    name = model_data_jip,
    command = clean_jip_data(data = data_maize_jip)
  ),
  ## Quatro Report 
  tar_quarto(
    name = maize_height_quatro_report,
    path = "analysis/maize_height.qmd",
    quiet = FALSE 
  )
)