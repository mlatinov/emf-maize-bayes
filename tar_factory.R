#### Maize Height Pipeline Function ####
## This pipeline takes the data from data directory in the root and Table 1 and runs sequentially :
## 1 Cleaning of the data : clean_maize_raw()
## 2 Data Simulation for Cell means model and Gompertz model 
## 3 Test the simulation recovery to both of them 
## 4 Fit the Stan models to both of them 
## 5 Produces a Small Quatro Document for the model summary.
maize_height_pipeline <- function(run = FALSE){
  ## Condition for running the pipeline 
  if(run){
    list(
      # Load the Table 1 
      tar_target(
        name = maize_heigh_raw,
        command = readxl::read_excel("data/Table S1 - 2026 Paunov et al. - 868 MHz EMF Maize - Plant Height.xlsx",sheet = 1)
      ),
      # Clean Table 1
      tar_target(
        name = model_data_maize_heigh,
        command = clean_maize_raw(maize_heigh_raw)
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
  }else{
    warning("Target Pipeline for Maize Height will be skipped")
  }
}

#### Maize JIP Pipeline ####
## Maize JIP pipeline have the goal to work with multiple outcomes and model them together. It runs sequentially
## 1 Cleaning of the data 
## 2 Simulation of the data 
## 3 Recovery Checks 
## 4 Real Data Modeling 
## 5 Small Quatro Document 
maize_jip_pipeline <- function(run = FALSE){
  ## Condition to run the pipeline or get warnings if skipped
  if(run){
    list(
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
      ## Simulation Datasets  
      tar_target(
        name = jip_sim_data,
        command = dgp_combined_jip(n = 188)
      ),
      #### Joined Model ####
      tar_target(
        name = jip_model_pp,
        command = combined_jip(model_data_jip, prior_only = 1)
      ),
      tar_target(
        name = jip_model_check,
        command = combined_jip(jip_sim_data$data, prior_only = 0)
      ),
      tar_target(
        name = jip_model,
        command = combined_jip(model_data_jip, prior_only = 0)
      ),
      ## Quatro Report 
      tar_quarto(
        name = jip_quatro_report,
        path = "analysis/maize_jip.qmd",
        quiet = FALSE 
      )
    )
  }else{
    warning("Targets Pipeline Maize JIT will be skipped")
  }
}
#### Maize Biochem/Biomass Pipeline ####
maize_biochem_biomass_pipeline <- function(run = FALSE){
   ## Condition to run the pipeline or get warnings if skipped
  if(run){
    list(
    ## Load join and clean All the datasets needed for the pipeline 
    tar_target(
      name = data_bb,            # bb - Biochem-Biomass
      command = clean_bb_data()
    ),
    ## Simulation dataset + truth and primaries Only Simulates DW Covariates 
    tar_target(
      name = sim_dw,
      command = simulate_biomass()
    ),
    ## Prior Predictive Checks
    tar_target(
      name = model_bb_pp,
      command = bb_model(data_bb, prior_only = 1)
    ),
    ## Simulation Recovery Check
    tar_target(
      name = model_bb_sim,
      command = bb_model(sim_dw)
    ),
    ## Model Fitting
    tar_target(
      name = model_bb,
      command = bb_model(data_bb)
    ),
    ## Quatro Report
    tar_quarto(
      name = bb_quatro_report,
      path = "analysis/maize_bb.qmd",
      quiet = FALSE 
    )
  )
  }else{
    warning("Target Pipeline Maize Biochem/Biomass will be skipped")
  }
}