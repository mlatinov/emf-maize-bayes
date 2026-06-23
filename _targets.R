
## Workflow Orchestration 
library(targets)
library(tidyverse)
library(tarchetypes)

## Soruce Functions ##
tar_source("R/clean_data_.R")
tar_source("R/maize_height/")
tar_source("R/maize_jip/")
tar_source("tar_factory.R")

#### Main Pipeline ### ========================================================
list(
  #### Maize Height Estimand Pipeline ####
  maize_height_pipeline(),

  #### Maize JIP GammaRC Estimand Pipeline ####
  maize_jip_gammaRC_pipilene(),

  ## Quatro Report 
  tar_quarto(
    name = maize_height_quatro_report,
    path = "analysis/maize_height.qmd",
    quiet = FALSE 
  )
)