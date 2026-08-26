
## Workflow Orchestration ##
library(targets)
library(tidyverse)
library(tarchetypes)

## Source Functions ##
tar_source("R/maize_bb/")
tar_source("R/maize_height/")
tar_source("R/maize_jip/")
tar_source("tar_factory.R")

#### Main Pipeline ### ========================================================

## INFORMATION: You can Skip a Pipeline with run = FALSE (eg maize_height_pipeline = FALSE). FALSE is the default behavior
# Its recommended because of the size of the project and the computing time it will take to run the entire pipeline
# You will get an Warining message for every skipped pipeline all pipelines produce a Quatro report in end where you can
# see the summary of the results, but i strongly recommend examining it. You can run a specific target ignoring the rest of the pipeline
# with tar_make(name), you can find the name of the specific target you want to run in tar_factory.R file, keep in mind that for targets
# that depend on other targets they will be ran too

## IMPORTANT NOTE !! The Quatro documents are mostly checking steps to summarize and make diagnostics and leave dev notes. At this point they are not intended for scientific communication.
#  The final document per project will be shared personally via email metodilatinov@abv.bg on request or will be in the docs file in this project. 
#  If you dont find any documents there they are in the making process or 
#  I dont have permission to share them openly. The data also will not be shared publicly which can cause problems for the pipeline execution. The solution is to point the pipeline to the 
#  directory where the data is if you have or to point to simulated data directory in this project.
#  For the restoration of the externally used libraries you can use renv::restore() (can possibly take some time to restore all the libraries)

## NOTE : All the model use 4 chains without parallel chains or multithreading you can change this behavers in the call_stan.R file where you can find the functions and the MCMC settings 
list(
  #### Maize Height Estimands Pipeline ####
  maize_height_pipeline(run = FALSE),

  #### Maize JIP Estimands Pipeline ####
  maize_jip_pipeline(run = FALSE),

  ### Maize Biomass/Biochem Pipeline ####
  biochem_biomass_pipeline(run = TRUE)

)