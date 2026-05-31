#### Function to clean the Raw data ####
clean_meize_raw <- function(maize_heigh_raw){

  maize_heigh_raw %>%
    # Rename the columns 
    rename(
      day          = Day,
      plant_height = `Plant_Height_(cm)`,
      plant_id     = Plant_Number ,
      pot          = Pot,
      treatment    = Treatment 
    ) %>%
    # Fix the Types
    mutate(
      day = as.integer(factor(day, levels = c(7, 14, 20, 27))),
      pot = as.integer(pot),
      plant_id  = as.integer(plant_id),
      treatment = as.integer(factor(treatment,levels = c("Control","Sham","EMF")))
    ) %>%
    # Remove Plants where Height = 0
    filter(plant_height > 0) %>%
    # Build the pot index 
    arrange(treatment, pot) %>%
    mutate(
      pot_id = as.integer(factor(paste(treatment, pot)))
    )
}