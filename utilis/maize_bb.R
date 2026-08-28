library(tidyverse)
library(patchwork)

tar_load(data_bb)

#### Univariate Function ####
get_unif <- function(
  data,
  x,
  title,
  x_title
){
  ## Get density split by day and treatment 
  dens <- ggplot(data, aes(x = .data[[x]], fill = treatment))+
  geom_density() +
  theme_minimal() +
  scale_fill_viridis_d(option = "A", begin = 0.2, end = 0.6) +
  facet_wrap(~treatment + day) +
  labs(
    title = title,
    x     = x_title,
    y     = "Density", 
    fill  = "Treatment"
  ) +
  theme(
    plot.title = element_text(face = "italic", size = 20, hjust = 0.5),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
  ## Get a Boxplot Comparison day and treatment 
  box_plot <- ggplot(data, aes(x = .data[[x]], fill = treatment))+
  geom_boxplot() +
  theme_minimal() +
  scale_fill_viridis_d(option = "C", begin = 0.2, end = 0.6) +
  facet_wrap(~treatment + day) +
  labs(
    title = title,
    x     = x_title,
    y     = "Density", 
    fill  = "Treatment"
  ) +
  theme(
    plot.title = element_text(face = "italic", size = 20, hjust = 0.5),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
  ## Get a regular histogram 
  hist <- ggplot(data, aes(x = .data[[x]])) +
    geom_histogram(color = "black",fill = "lightblue") +
    theme_minimal() +
    labs(
    title = title,
    x     = x_title,
    y     = "Count", 
  ) +
  theme(
    plot.title = element_text(face = "italic", size = 20, hjust = 0.5),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
  ## Return list with the viz
  list(
    density_day_x_treatment = dens,
    box_plot_comparison     = box_plot,
    histogram               = hist
  )
}

#### Function to get the relationship between x and y ####
get_relan <- function(data, x, y, title, x_title, y_title){
  ## Overview 
  plot <- ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm") +
  theme_minimal() +
  labs(
    title = title,
    x     = x_title,
    y     = y_title
  ) +
  theme(
    plot.title   = element_text(face = "italic", size = 20, hjust = 0.5),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
  # By Treatment 
  plot_2 <- ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm") +
  facet_wrap(~treatment) +
  theme_minimal() +
  labs(
    title = title,
    x     = x_title,
    y     = y_title
  ) +
  theme(
    plot.title   = element_text(face = "italic", size = 20, hjust = 0.5),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
  # Return a list 
  return(list(
    overview     = plot,
    by_treatment = plot_2 
  ))
}

#### SOT ####
sot_viz <- get_unif(
  data     = data_bb, 
  x        = "sod_dw", 
  title    = "SOD Dry Weight Primary",
  x_title = "SOD, U/mg DW"
)
compare <- rsims::check_fidelity(real = data_bb, sim = sod_dw$data, vars = "sod_dw")

#### CAT ###
cat_viz <- get_unif(
  data   = data_bb,
  x      = "cat_dw",
  title  = "CAT Dry Weight Primary",
  x_title = "CAT, µmol H2O2/min/mg DW" 
)
compare <- rsims::check_fidelity(real = data_bb, sim = cat_dw$data, vars = "cat_dw")

#### TEAC ####
teac_viz <- get_unif(
  data     = data_bb, 
  x        = "trolox_umol_dw", 
  title    = "Trolox Dry Weight Primary",
  x_title = "Trolox Equivalent, µmol/g DW"
)
compare <- rsims::check_fidelity(real = data_bb, sim = trolox_umol_dw$data, vars = "trolox_umol_dw")

#### Sugar ###
sugar_viz <- get_unif(
  data   = data_bb,
  x      = "reduced_sugar_dw",
  title  = "Reducing Sugars Dry Weight Primary",
  x_title = "Reducing Sugars, mg Glucose/g DW" 
)

#### H202 ####
h2_o2_viz <- get_unif(
  data   = data_bb,
  x      = "h2_o2_dw",
  title  = "H2O2 Dry Weight Secondary",
  x_title = "H2O2, nmol/g DW" 
)
h2_o2_x_sod <- get_relan(
  data = data_bb, 
  x = "sod_dw", 
  y = "h2_o2_dw",
  title = "Relationship between H2O2 and SOD",
  y_title = "H2O2, nmol/g DW",
  x_title = "SOD, U/mg DW"
)
h2_o2_x_cat <- get_relan(
  data = data_bb, 
  x = "cat_dw", 
  y = "h2_o2_dw",
  title = "Relationship between H2O2 and CAT",
  y_title = "H2O2, nmol/g DW",
  x_title = "CAT, µmol H2O2/min/mg DW"
)
h2_o2_x_teac <- get_relan(
  data = data_bb, 
  x = "trolox_umol_dw", 
  y = "h2_o2_dw",
  title = "Relationship between H2O2 and TEAC",
  y_title = "H2O2, nmol/g DW",
  x_title = "Trolox Equivalent, µmol/g DW"
)
h2_o2 <- simulate_h2_02(seed = 32)$sim
compare <- rsims::check_fidelity(
  real = data_bb, 
  sim = h2_o2, 
  vars = c("cat_dw", "trolox_umol_dw", "sod_dw", "h2_o2_dw")
)

#### MAD ###
mad <- get_unif(
  data = data_bb, 
  x       = "mda_dw",
  title   = "MDA Dry Weight", 
  x_title = "MDA, nmol/g DW"
)
mad_x_h2o2 <- get_relan(
  data = data_bb,
  x = "h2_o2_dw",
  y = "mda_dw",
  title = "Relationship MAD vs H202",
  x_title = "H2O2, nmol/g DW",
  y_title = "MDA, nmol/g DW"
)
mda_data <- simulate_mda()$sim_data
compare <- rsims::check_fidelity(
  real = data_bb, 
  sim = mda_data, 
  vars = c("cat_dw", "trolox_umol_dw", "sod_dw", "h2_o2_dw", "mda_dw")
)

#### Sugars ####
sugar_viz <- get_unif(
  data    = data_bb, 
  x       = "reduced_sugar_dw",
  title   = "Reducing Sugars : Glucose ",
  x_title = "Reducing Sugars, mg Glucose/g DW "
)
sugar_sim <- simulate_primaries(
    sdlog_obs   = 0.2,
    cell_mean   = 80,
    cell_sd_log = 0.02,
    zeta_rate   = 15,
    primary_name = "reduced_sugar_dw"
) 
compare <- rsims::check_fidelity(
  real = data_bb,
  sim  = sugar_sim$data,
  vars = "reduced_sugar_dw"
 )

#### Biomass ####
biomass_viz <- get_unif(
  data = data_bb,
  x    = "water_content_dw",
  title   = "Leaf Water Content Relative to DW",
  x_title = "Leaf Water Content"
)

biomass_relan <- get_relan(
  data = data_bb,
  x    = "mda_dw",
  y    = "water_content_dw",
  title = "Relationship Between Water Content Dry Weight and MDA",
  x_title = "MDA, nmol/g DW",
  y_title = "Leaf Water Content Relative to DW"
)

biomass_relan <- get_relan(
  data = data_bb,
  x    = "reduced_sugar_dw",
  y    = "water_content_dw",
  title = "Relationship Between Water Content Dry Weight and Sugars",
  x_title = "Reducing Sugars, mg Glucose/g DW ",
  y_title = "Leaf Water Content Relative to DW"
)
biomass_data <- simulate_biomass()
compare <- rsims::check_fidelity(
  real = data_bb,
  sim  = biomass_data$sim,
  vars = c("cat_dw", "trolox_umol_dw", "sod_dw", "h2_o2_dw", "mda_dw","water_content_dw","reduced_sugar_dw")  
)