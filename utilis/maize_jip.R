#### 
library(ggstatsplot)
library(patchwork)

#### Univariate Analysis ####

## GammaRC 
ggplot(data = model_data_jip, aes(x = gammaRC, fill = treatment))+
  geom_density(alpha = 0.2)+
  facet_wrap(~day)

ggplot(data = model_data_jip, aes(x = gammaRC, fill = treatment))+
  geom_boxplot()+
  facet_wrap(~day)

model_data_jip %>%
  group_by(treatment,pot,day) %>%
  summarise(
    mean_gamma = mean(gammaRC),
    median_gamma = median(gammaRC),
    sd_gamma     = sd(gammaRC),
    n_gamma      = n() 
  ) %>%
  ggplot(aes(x = pot, y = mean_gamma, fill = treatment))+
  geom_col(position = "dodge")+
  facet_grid(treatment ~ day)

grouped_ggbetweenstats(
  data = model_data_jip,
  x = treatment,
  y = gammaRC,
  grouping.var = day
)

model_data_jip %>%
  filter(treatment == "Control") %>%
  grouped_ggbetweenstats(
    x = pot,
    y = gammaRC,
    grouping.var = day
  )

model_data_jip %>%
  filter(treatment == "Sham") %>%
  grouped_ggbetweenstats(
    x = pot,
    y = gammaRC,
    grouping.var = day
  )

model_data_jip %>%
  filter(treatment == "EMF") %>%
  grouped_ggbetweenstats(
    x = pot,
    y = gammaRC,
    grouping.var = day
  )

## PhiPo ##
ggplot(data = model_data_jip, aes(x = phi_Po, fill = factor(day)))+
  geom_density(alpha = 0.3) +
  facet_wrap(~factor(treatment,labels = c("Control","Sham", "EMF")))+
  labs(
    title = "Distribution of Phi By Treatment and Day",
    x = "Phi",
    y = "Density",
    fill = "Day"
  )+
  theme_minimal()+
  theme(plot.title = element_text(family = "bold", size = 20, hjust = 0.5))

grouped_ggbetweenstats(
  data = model_data_jip,
  x = treatment,
  y = phi_Po,
  grouping.var = day
)

model_data_jip %>%
  filter(treatment == "1") %>%
  grouped_ggbetweenstats(
    grouping.var = day,
    x            = pot,
    y = phi_Po
)
model_data_jip %>%
  filter(treatment == "2") %>%
  grouped_ggbetweenstats(
    grouping.var = day,
    x            = pot,
    y            = phi_Po
)
model_data_jip %>%
  filter(treatment == "3") %>%
  grouped_ggbetweenstats(
    grouping.var = day,
    x            = pot,
    y            = phi_Po
)

## Psi ##
ggplot(data = model_data_jip, aes(x = psi_Eo, fill = factor(day)))+
  geom_density(alpha = 0.3)+
  facet_wrap(~factor(treatment, labels = c("Control", "Sham", "EMF")))+
  labs(
    title = "Density Distribution of Psi Eo by Day and Treatment", 
    x = "Psi Eo",
    y = "Density",
    fill = "Day"
  )+
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5))

grouped_ggbetweenstats(
  data = model_data_jip,
  x = treatment,
  y = psi_Eo,
  grouping.var = day
)

model_data_jip %>%
  filter(treatment == "1") %>%
  grouped_ggbetweenstats(
    grouping.var = day,
    x            = pot,
    y            = psi_Eo
)
model_data_jip %>%
  filter(treatment == "2") %>%
  grouped_ggbetweenstats(
    grouping.var = day,
    x            = pot,
    y            = psi_Eo
)
model_data_jip %>%
  filter(treatment == "3") %>%
  grouped_ggbetweenstats(
    grouping.var = day,
    x            = pot,
    y            = psi_Eo
)

## Ro Model ##
ggplot(data = model_data_jip, aes(x = delta_Ro, fill = factor(day)))+
  geom_density(alpha = 0.3)+
  facet_wrap(~factor(treatment, labels = c("Control", "Sham", "EMF")))+
  labs(
    title = "Density Distribution of Delta Ro by Day and Treatment",
    x = "Delta Ro",
    y = "Density",
    fill = "Day"
  )+
  theme_minimal()+
  theme(plot.title = element_text(face = "bold",size = 20, hjust = 0.5))

ggplot(
  data = model_data_jip, 
  aes(x = delta_Ro, fill = factor(treatment,labels = c("Control", "Sham", "EMF"))))+
  geom_density(alpha = 0.3)+
  facet_wrap(~factor(day, labels = c("Day 13", "Day 20", "Day 27")))+
  labs(
    title = "Density Distribution of Delta Ro by Day and Treatment",
    x = "Delta Ro",
    y = "Density",
    fill = "treatment"
  )+
  theme_minimal()+
  theme(plot.title = element_text(face = "bold",size = 20, hjust = 0.5))



#### EFA ####
model_data_jip %>%
  group_by(treatment, day, pot_idx) %>%
  summarise(
    cor_phi_psi   = cor(phi_Po, psi_Eo, use = "complete.obs"),
    cor_psi_delta = cor(psi_Eo, delta_Ro, use = "complete.obs"),
    cor_gamma_phi = cor(gammaRC, phi_Po, use = "complete.obs"),
    cor_gamma_psi = cor(gammaRC, psi_Eo, use = "complete.obs"),
    .groups = "drop"
  )
ggplot(model_data_jip, aes(psi_Eo, delta_Ro, color = factor(day))) +
  geom_point() +
  facet_wrap(~treatment)
summary(lm(delta_Ro ~ psi_Eo * day + treatment, data = model_data_jip))

