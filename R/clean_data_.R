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
      day = as.integer(day),
      pot = as.integer(pot),
      plant_id = as.integer(plant_id) 
    )
}