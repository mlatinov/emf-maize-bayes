library(ggplot2)
library(ggstatsplot)
library(patchwork)
library(dplyr)

# === Setup: consistent theme + palette ============================

treatment_colors <- c(
    "Control" = "#2E5984",
    "Sham"    = "#B8860B",
    "EMF"     = "#7B2D8E"
)

plot_data <- model_data_maize_heigh %>%
    mutate(
        treatment = factor(treatment, 
                          levels = c(1, 2, 3),
                          labels = c("Control", "Sham", "EMF")),
        day = factor(day, 
                    levels = c(7, 14, 20, 27),
                    labels = c("Day 7", "Day 14", "Day 20", "Day 27"))
    )

theme_paper <- theme_minimal(base_size = 13) +
    theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0),
        plot.subtitle = element_text(color = "gray40", size = 11, hjust = 0),
        strip.text = element_text(face = "bold", size = 11),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
    )

plot_distributions <- ggplot(plot_data, aes(x = plant_height, fill = treatment)) +
    geom_density(alpha = 0.5, color = NA) +
    geom_rug(aes(color = treatment), alpha = 0.3, sides = "b", length = unit(0.04, "npc")) +
    facet_wrap(~ day, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = treatment_colors) +
    scale_color_manual(values = treatment_colors) +
    labs(
        x = "Plant height (cm)",
        y = "Density",
        fill = "Treatment",
        color = "Treatment",
        title = "Distribution of plant heights by day and treatment",
        subtitle = "Each panel shows one measurement day; rug ticks mark individual plants"
    ) +
    theme_paper +
  theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5)    
  )

make_day_comparison <- function(day_label) {

  day_data <- plot_data %>% filter(day == day_label)
  
  p <- ggbetweenstats(
    data = day_data,
    x = treatment,
    y = plant_height,
    pairwise.display = "none",
    bf.message = FALSE,
    results.subtitle = FALSE,
    caption = NULL,
    centrality.label.args = list(size = 3.5),
    violin.args = list(width = 0.7, alpha = 0.4),
    point.args = list(alpha = 0.4, size = 1.5),
    ggplot.component = list(
        scale_color_manual(values = treatment_colors),
        labs(
            title = day_label,
            x = NULL,
            y = "Plant height (cm)"
        )
    )
)
p +
  labs(subtitle = NULL, caption = NULL) +
    theme_paper +
      theme(
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_blank(),
        plot.caption = element_blank()
    )
}

p_day7  <- make_day_comparison("Day 7")
p_day14 <- make_day_comparison("Day 14")
p_day20 <- make_day_comparison("Day 20")
p_day27 <- make_day_comparison("Day 27")

plot_comparisons <- (p_day7 + p_day14) / (p_day20 + p_day27) +
    plot_annotation(
        title = "Treatment comparisons at each measurement day",
        subtitle = "Violin = distribution; dot in centre = median; points = individual plants",
        theme = theme(
            plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
            plot.subtitle = element_text(color = "gray40", size = 11, hjust = 0)
        )
    )


plot_pot_uniformity <- plot_data %>%
    group_by(pot_id, treatment, day) %>%
    summarize(
        mean_height = mean(plant_height),
        n_plants = n(),
        .groups = "drop"
    ) %>%
    ggplot(aes(x = treatment, y = mean_height, color = treatment)) +
    geom_jitter(aes(size = n_plants), width = 0.15, alpha = 0.7) +
    facet_wrap(~ day, scales = "free_y", ncol = 4) +
    scale_color_manual(values = treatment_colors) +
    scale_size_continuous(range = c(2, 5), name = "Plants per pot") +
    labs(
        x = NULL,
        y = "Pot mean height (cm)",
        color = "Treatment",
        title = "Pot-level variation within treatments",
        subtitle = "Each point is one pot's mean height on that day; point size reflects number of plants in pot"
    ) +
    theme_paper +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(face = "bold", size = 20, hjust = 0.5)
    )

pot_trajectories <- plot_data %>%
    group_by(pot_id, treatment, day) %>%
    summarize(
        mean_height = mean(plant_height),
        .groups = "drop"
    ) %>%
    mutate(day_num = as.numeric(as.character(
        factor(day, labels = c(7, 14, 20, 27))
    )))

# Convert day back to numeric for plotting on continuous axis
pot_trajectories <- plot_data %>%
    group_by(pot_id, treatment) %>%
    mutate(day_num = case_when(
        day == "Day 7"  ~ 7,
        day == "Day 14" ~ 14,
        day == "Day 20" ~ 20,
        day == "Day 27" ~ 27
    )) %>%
    group_by(pot_id, treatment, day_num) %>%
    summarize(mean_height = mean(plant_height), .groups = "drop")

plot_pot_trajectories <- ggplot(pot_trajectories,
                                 aes(x = day_num, y = mean_height,
                                     group = pot_id, color = treatment)) +
    geom_line(alpha = 0.6, linewidth = 0.8) +
    geom_point(alpha = 0.8, size = 2) +
    facet_wrap(~ treatment, ncol = 3) +
    scale_color_manual(values = treatment_colors, guide = "none") +
    scale_x_continuous(breaks = c(7, 14, 20, 27)) +
    labs(
        x = "Day from experiment start",
        y = "Pot mean height (cm)",
        title = "Per-pot growth trajectories",
        subtitle = "Each line is one pot; the same pot is followed across all four measurement days"
    ) +
    theme_paper +
  theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5)
  )


make_pot_plot <- function(trt) {

  ggplot(
    filter(pot_trajectories, treatment == trt),
    aes(
      x = day_num,
      y = mean_height,
      group = pot_id
    )
  ) +
    geom_line(
      color = treatment_colors[[trt]],
      alpha = 0.5,
      linewidth = 0.7
    ) +
    geom_point(
      color = treatment_colors[[trt]],
      alpha = 0.8,
      size = 2
    ) +

    # population trend
    geom_smooth(
      aes(group = 1),
      method = "loess",
      se = FALSE,
      color = treatment_colors[[trt]],
      linewidth = 1
    ) +

    scale_x_continuous(
      breaks = c(7, 14, 20, 27)
    ) +

    labs(
      title = trt,
      x = NULL,
      y = "Mean height (cm)"
    ) +

    theme_paper +
    theme(
      plot.title = element_text(face = "bold",hjust = 0.5)
    )
}

p_control <- make_pot_plot("Control")
p_sham    <- make_pot_plot("Sham")
p_emf     <- make_pot_plot("EMF")

plot_pot_trajectories <- (p_control / p_sham) / p_emf +
  plot_annotation(
    title = "Per-pot growth trajectories",
    subtitle = "Each line represents one pot tracked across measurement days",
    theme = theme(
        plot.title = element_text(size = 20, face = "bold",hjust = 0.5)
    )
  )

plot_pot_trajectories
