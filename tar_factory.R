#### Maize Height Pipeline Function ####
maize_height_pipeline <- function(){
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
    )
  )
}
#### Maize JIP GammaRC Pipeline ####
maize_jip_gammaRC_pipilene <- function(){
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
      name = sim_gammaRC_jip_1,
      command = dgp_gammaRC_jip_1(n = 188)
    ),
    tar_target(
      name = sim_gammaRC_jip_2,
      command = dgp_gammaRC_jip_2(n = 188)
    ),
    tar_target(
      name = jip_sim_data_ncor,
      command = dgp_combined_jip(n = 188)
    ),
    #### Cell Means GammaRC Model 1 Without Heirarcle or Pooling ####
    ## Model Checks 
    tar_target(
      name = gammarc_model_1_pp,
      command = m1_gamma_rc_model(model_data_jip, prior_only = 1)
    ),
    tar_target(
      name = gammarc_model_1_recovery_check,
      command = m1_gamma_rc_model(sim_gammaRC_jip_1$data, prior_only = 0)
    ),
    ## Fit model to the data 
    tar_target(
      name = gammarc_model_1,
      command = m1_gamma_rc_model(model_data_jip, prior_only = 0)
    ),
    #### Cell-Means GammRc Hierarcle Model with partial pooling ####
    ## Model Checks 
    tar_target(
      name = gammarc_model_2_pp,
      command = m2_gamma_rc_model(model_data_jip, prior_only = 1)
    ),
    tar_target(
      name = gammarc_model_2_recovery_check,
      command = m2_gamma_rc_model(sim_gammaRC_jip_2$data, prior_only = 0)
    ),
    ## Fit model to the data 
    tar_target(
      name = gammarc_model_2,
      command = m2_gamma_rc_model(model_data_jip, prior_only = 0)
    )
  )
}