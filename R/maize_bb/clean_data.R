#### Function to Clean SOT data ####
clean_sot <- function(sot_data){
  sot_data %>%
    # Rename the variables with R friendly names  
    rename(
      sod_dw = `SOD, U/mg DW`,
      sod_fw = `SOD, U/mg FW`,
      sod_protein = `SOD, U/mg Protein`,
    ) %>%
    # Average over the replicates
    group_by(Day,Treatment,Leaf, Pot) %>%
    summarise(
        sod_dw = mean(sod_dw),
        sod_fw = mean(sod_fw),
        sod_protein = mean(sod_protein),
        .groups = "drop"
    ) %>%
    mutate(
      id = paste(Day,Treatment,Leaf, Pot) # Recover the Sample Id later for join
    ) %>%
    # Remove repeating covariates across datasets
    select(id, sod_dw, sod_fw, sod_protein)
}

#### Function to clean the TEAC data ####
clean_teac <- function(teac_data){
  teac_data %>%
    # Rename the variables with R friendly names 
    rename(
      trolox_umol_fw = `Trolox Equivalent, µmol/g FW`,
      trolox_umol_dw = `Trolox Equivalent, µmol/g DW`
    ) %>%
    # Average over the replicates
    group_by(Day,Treatment,Leaf, Pot) %>%
    summarise(
      trolox_umol_fw = mean(trolox_umol_fw),
      trolox_umol_dw = mean(trolox_umol_dw),
      .groups = "drop"
    ) %>%
    mutate(
      id = paste(Day,Treatment,Leaf, Pot) # Recover the Sample Id later for join
    ) %>%
    # Remove repeating covariates across datasets
    select(id, trolox_umol_fw, trolox_umol_dw)
}

#### Function to clean the CAT data ####
clean_cat <- function(cat_data){ 
  cat_data %>%
    # Rename the variables with R friendly names 
    rename(
      cat_fw    = `CAT, µmol H2O2/min/mg FW`,
      cat_dw    = `CAT, µmol H2O2/min/mg DW`,
      cat_protein = `CAT, µmol H2O2/min/mg Protein`
    ) %>%
    # Average over the replicates
    group_by(Day,Treatment,Leaf, Pot) %>%
    summarise(
      cat_fw      = mean(cat_fw),
      cat_dw      = mean(cat_dw),
      cat_protein = mean(cat_protein),
      .groups = "drop"
    ) %>%
    mutate(
      id = paste(Day,Treatment,Leaf, Pot) # Recover the Sample Id later for join
    ) %>%
    # Remove repeating covariates across datasets
    select(id, cat_fw, cat_dw, cat_protein)
}

#### Function to clean the H202 data ####
clean_h2o2 <- function(h2_o2_data){
  h2_o2_data %>%
    # Rename the variables with R friendly names 
    rename(
      h2_o2_fw  = `H2O2, nmol/g FW`,
      h2_o2_dw  = `H2O2, nmol/g DW`
    ) %>%
    # Average over the replicates
    group_by(Day, Treatment, Leaf, Pot) %>%
    summarise(
      h2_o2_fw = mean(h2_o2_fw),
      h2_o2_dw = mean(h2_o2_dw),
      .groups = "drop"
  ) %>%
    mutate(
      id = paste(Day,Treatment,Leaf, Pot) # Recover the Sample Id later for join
    ) %>%
    # Remove repeating covariates across datasets
    select(id, h2_o2_fw, h2_o2_dw)
}

#### Function to clean MDA data ####
clean_mda <- function(mda_data){ 
  mda_data %>%
    # Rename the variables with R friendly names 
    rename(
      mda_fw    = `MDA, nmol/g FW`,
      mda_dw    = `MDA, nmol/g DW`
    ) %>%
    mutate(
      # For the names and construction of the Id to match the Rest of the datasets
      id = paste(Day,Treatment,Leaf, Pot)
    ) %>%
    # Remove repeating covariates across datasets
    select(id, mda_fw, mda_dw)
}

#### Function to clean the Sugar data ####
clean_sugar <- function(sugar_data){  
  sugar_data %>%
    # Rename the variables with R friendly names 
    rename(
      reduced_sugar_fw = `Reducing Sugars, mg Glucose/g FW`,
      reduced_sugar_dw = `Reducing Sugars, mg Glucose/g DW`
    ) %>%
    # Average over the replicates
    group_by(Day, Treatment, Leaf, Pot) %>%
    summarise(
      reduced_sugar_fw = mean(reduced_sugar_fw),
      reduced_sugar_dw = mean(reduced_sugar_dw),
      .groups = "drop"
  ) %>%
    mutate(
      id = paste(Day, Treatment, Leaf, Pot) # Recover the Sample Id later for join
    ) %>%
    # Remove repeating covariates across datasets
    select(id, reduced_sugar_fw, reduced_sugar_dw)
}

#### Function to clean the Biomass data ####
clean_biomass <- function(biomass_data){
  biomass_data %>%
    # Rename the variables with R friendly names 
    rename(
      sample = Sample,
      day    = Day,
      treatment = Treatment,
      pot       = Pot,
      leaf      = Leaf,
      total_fresh_mass = `Total Fresh Leaf Biomass, g (1 leaf per plant, 4 plants per pot, combined)`,
      dw_fw            = `DW/FW, %`,
      water_content_dw = `Leaf Water Content Relative to DW, rel. u.`,
      water_percent    = `Leaf Water Content, %`
    ) %>%
    mutate(
      # For the names and construction of the Id to match the Rest of the datasets
      id = paste(day,treatment,leaf, pot) 
  ) 
}

#### Function to Load, Join and clean all the dataset needed for the Biochem-Biomass Pipeline ####
clean_bb_data <- function(){

  ## Load the datasets needed 
  sot_data_raw  <- readxl::read_excel("data/Table S9 - 2026 Paunov et al. - 868 MHz EMF Maize - SOD Activity.xlsx", sheet = 2)
  teac_data_raw <- readxl::read_excel("data/Table S8 - 2026 Paunov et al. - 868 MHz EMF Maize - TEAC.xlsx",         sheet = 2)
  cat_data_raw  <- readxl::read_excel("data/Table S10 - 2026 Paunov et al. - 868 MHz EMF Maize - CAT Activity.xlsx",sheet = 2)  
  h2_o2_data_raw   <- readxl::read_excel("data/Table S7 - 2026 Paunov et al. - 868 MHz EMF Maize - Hidrogen Peroxide.xlsx",sheet = 2) 
  mda_data_raw     <- readxl::read_excel("data/Table S6 - 2026 Paunov et al. - 868 MHz EMF Maize - Malondialdehyde.xlsx",  sheet = 2) 
  biomass_data_raw <- readxl::read_excel("data/Table S3 - 2026 Paunov et al. - 868 MHz EMF Maize - Leaf Biomass.xlsx",     sheet = 2)
  sugars_data_raw  <- readxl::read_excel("data/Table S4 - 2026 Paunov et al. - 868 MHz EMF Maize - Reducing Sugars.xlsx",  sheet = 2) 

  ### Average Over the Tech Replicate and clean the names ###  
  sot_data   <- clean_sot(sot_data_raw)
  cat_data   <- clean_cat(cat_data_raw)
  mda_data   <- clean_mda(mda_data_raw)
  teac_data  <- clean_teac(teac_data_raw)
  h2_o2_data <- clean_h2o2(h2_o2_data_raw)
  sugar_data <- clean_sugar(sugars_data_raw)
  biomass_data <- clean_biomass(biomass_data_raw)

  ## Join the dataset together via day treatment pot combination ##
  combined_data <- biomass_data %>% 
    inner_join(sot_data,   by = "id") %>%
    inner_join(cat_data,   by = "id") %>%
    inner_join(mda_data,   by = "id") %>%
    inner_join(teac_data,  by = "id") %>%
    inner_join(h2_o2_data, by = "id") %>%
    inner_join(sugar_data, by = "id") %>%
    select(-id,-sample) %>%
    # Get Unique identifiers later for modeling in Stan
    mutate(
      treatment_id     = as.integer(factor(treatment, levels = c("Control","Sham","EMF"))),
      pot_treatment_id = as.integer(as.factor(paste(pot, treatment,sep = "_"))),
      cell_id          = as.integer(as.factor(paste(treatment, day,sep = "_"))),
      day_id           = as.integer(day)
    )

  # Return Combined data 
  return(combined_data)
}