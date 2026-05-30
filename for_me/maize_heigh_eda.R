
#### Libraries ####
library(ggstatsplot)
library(patchwork)

#### Basic Checks ####

## Summary Stats by day and treatment ##
model_data_maize_heigh %>%
  group_by(day, treatment) %>%
  summarise(
    n = n(),
    mean_height = mean(plant_height),
    sd_height   = sd(plant_height),
    min_height  = min(plant_height),
    max_height  = max(plant_height)
  )

## Summary Stats by treatment and pot ##
model_data_maize_heigh %>%
  group_by(treatment, pot) %>%
  summarise(
    n = n(),
    mean_height = mean(plant_height),
    sd_height   = sd(plant_height),
    min_height  = min(plant_height),
    max_height  = max(plant_height)
  ) %>%
  arrange(pot)

## Plants By pot by week  ##
model_data_maize_heigh %>%
  group_by(pot,day, treatment) %>%
  summarise(
    n = n()
  )

### Devide by Week ###
## Week 1
week_1 <- model_data_maize_heigh %>% filter(day == 7) %>%
  mutate(
    treatment = as.factor(treatment)
  )
## Week 2
week_2 <- model_data_maize_heigh %>% filter(day == 14) %>%
  mutate(
    treatment = as.factor(treatment)
  )
## Week 3
week_3 <- model_data_maize_heigh %>% filter(day == 20) %>%
  mutate(
    treatment = as.factor(treatment)
  )
## Week 4
week_4 <- model_data_maize_heigh %>% filter(day == 27) %>%
  mutate(
    treatment = as.factor(treatment)
  )

#### Distributions of Plan Heights for the treatment group and week ####

# Week 1 
week_1_dist <- grouped_gghistostats(
  data = week_1,
  x = plant_height,
  grouping.var = treatment,
  results.subtitle = FALSE
)

week_1_dist +
  plot_annotation(
    title = "Week 1 Plant Height Distribution by Treatment"
  ) &
  theme_bw(base_size = 14) &
  theme(
    plot.title = element_text(
      size = 20,
      face = "bold",
      hjust = 0.5
    ),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(
      size = 14,
      face = "bold"
    ),
    panel.grid.minor = element_blank()
  ) &
  labs(
    x = "Plant Height (cm)",
    y = "Count"
  )

# Week 2
grouped_gghistostats(
  data = week_2,
  x = plant_height,
  grouping.var = treatment
)

# Week 3
grouped_gghistostats(
  data = week_3,
  x = plant_height,
  grouping.var = treatment
)

# Week 4
grouped_gghistostats(
  data = week_4,
  x = plant_height,
  grouping.var = treatment
)

#### Diff in Plant Height between the treatment group per week ####

# Weel 1
ggbetweenstats(
  data = week_1,
  x = treatment,
  y = plant_height
)
# Weel 2
ggbetweenstats(
  data = week_2,
  x = treatment,
  y = plant_height
)
# Weel 3
ggbetweenstats(
  data = week_3,
  x = treatment,
  y = plant_height
)
# Weel 4
ggbetweenstats(
  data = week_4,
  x = treatment,
  y = plant_height
)

#### Pot Level Diff in weeks ####
control_week_1 <- week_1 %>%
  filter(treatment == "Control")
emf_week1 <- week_1 %>%
  filter(treatment == "EMF")
sham_week1 <- week_1 %>%
  filter(treatment == "Sham")

## Control Only Week 1 
ggstatsplot::ggdotplotstats(
  data = control_week_1,
  x = plant_height,
  y = pot
)
## Treatment EMF Only Week 1
ggstatsplot::ggdotplotstats(
  data = emf_week1,
  x = plant_height,
  y = pot
)
## Treatment Sham Only Week 1
ggstatsplot::ggdotplotstats(
  data = emf_week1,
  x = plant_height,
  y = pot
)
