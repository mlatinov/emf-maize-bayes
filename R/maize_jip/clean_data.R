
#### Function to Clean JIP data ####
clean_jip_data <- function(data_maize_jip){
  data_maize_jip %>%
    # Rename The Variables 
    rename(
      treatment = Treatment,
      day      = Day,
      leaf     = Leaf,
      pi_abs   = PIabs,
      pi_total = PItotal,     
      plant    = Plant,
      pot      = Pot,
      rc_cso   = `RC/CSo`,
      tech_rep = Tech_Replicate,
      id       = Sample
    ) %>%
    # Average Over the Tech replicate where they are present 
    group_by(day, treatment, pot, plant, leaf) %>%
    summarise(
      rc_cso   = mean(rc_cso),
      gammaRC  = mean(gammaRC),
      phi_Po   = mean(phi_Po),
      psi_Eo   = mean(psi_Eo),
      delta_Ro = mean(delta_Ro),
      n_reps   = n(),          
      .groups  = "drop"
    ) %>%
    mutate(
      # Add Unique Pot identifier
      pot_idx  = as.integer(factor(paste(treatment, pot))),
      day_idx = as.integer(factor(day,levels = c(13,20,27))),
      treatment = as.integer(factor(treatment,levels = c("Control","Sham","EMF")))
    ) %>%
    # Order by Pot and treatment 
    arrange(treatment,pot)
}