## intraspecific variability in growth alone vs. neighbours
## Author: Emily H
## Created: November 28, 2024
## Last edited: January 21, 2025

#install.packages("tidyverse")
#install.packages("vegan")
#install.packages("reshape2")
#install.packages("lme4")
#install.packages("lmerTest")
#install.packages("car")
#devtools::install_github("cttobin/ggthemr")

library(vegan)
library(tidyverse)
library(reshape2)
library(lme4)
library(lmerTest)
library(car)
library(ggthemr)

ggthemr("fresh")

#setwd("C:/Users/emilyh/Documents/GitHub/Mesocosm_networks")
setwd("~/Documents/GitHub/Mesocosm networks/Mesocosm_networks")

#### import data ####
## import mesocosm data
meso1 <- read.csv("data/2022 Mesocosm AB and BG biomasses-2.csv") %>%
  dplyr::select(Pot, Block, Nutrients, Position.in.pot, Plant.identity, Neighbour.identity, Dead.at.harvest.,
                Aboveground.biomass..g., Belowground.biomass..g., Total.biomass..g.) %>%
  dplyr::rename(AB.biomass = Aboveground.biomass..g.,
                BG.biomass = Belowground.biomass..g.,
                total.biomass = Total.biomass..g.,
                pot.position = Position.in.pot) %>%
  subset(Plant.identity != "") %>% #removes cells where the plant was not found at harvest
  mutate(Neighbour.identity = ifelse(Neighbour.identity %in% "", "alone", Neighbour.identity)) %>%
  # clarify which plants where grown alone
  mutate(Dead.alive = if_else(Dead.at.harvest. == "Y",1,0)) %>%
  filter(Plant.identity != "Ast lae" & Neighbour.identity != "Ast lae") %>% #remove species with data only in the unfertilized treatment
  filter(Plant.identity != "Ast fal" & Neighbour.identity != "Ast fal") %>% #remove species with data only in the unfertilized treatment
  filter(Plant.identity != "Fra vir" & Neighbour.identity != "Fra vir") %>% #remove species with data only in the unfertilized treatment
  filter(Neighbour.identity != "Hel hoo") %>% #remove Hel hoo (no alone plant data)
  dplyr::select(!Dead.at.harvest.)
str(meso1) #1190 plants
length(unique(meso$Pot)) #649 pots

##import alone plant mean growth data
alone.df <- read_rds("output/alone_df_new.rds")

## create new df of the variables I want
datdat <- meso1 %>%
  left_join(alone.df, by = c("Plant.identity", "Nutrients")) %>% #add alone data to df
  dplyr::select(Pot, Block, pot.position, Plant.identity, Neighbour.identity, Nutrients, alone.replicates, total.biomass, alone.mean.total.biomass) %>% #select columns I want
  mutate(wneighbour = if_else(Neighbour.identity == "alone", "Alone", "With neighbour")) #differentiate plants grown with and without neighbours
str(datdat)

#### calculate CV ####
cvs1 <- datdat %>% 
  group_by(Plant.identity, wneighbour, Nutrients) %>%
  summarise(replicates = length(total.biomass),
            mean.bio = mean(total.biomass, na.rm = TRUE),
            sd.bio = sd(total.biomass, na.rm = TRUE),
            cv.bio = sd.bio/mean.bio)
str(cvs1)

# Calculate the mean cv for alone vs wneighbour plants
mean_cv_bio_by_wneighbour <- cvs1 %>%
  group_by(wneighbour) %>%
  summarize(mean_cv_bio = mean(cv.bio, na.rm = TRUE))

#### make table S8 ####
tables8 <- cvs1 %>%
  dplyr::select(Plant.identity, Nutrients, wneighbour, cv.bio) %>%
  mutate(cv.bio = round(cv.bio, 3)) %>%
  pivot_wider(names_from = wneighbour, values_from = cv.bio) %>%
  rename('Species code' = Plant.identity,
         'CV when grown alone' = Alone,
         'CV when grown with neighbours' = 'With neighbour')

##export
write_csv(tables8, "output/tableS8.csv")

#### model ####
cv.mod <- lmer(cv.bio~wneighbour*Nutrients + (1|Plant.identity), data = cvs1)
summary(cv.mod)
Anova(cv.mod, type = "III")
qqnorm(resid(cv.mod))
qqline(resid(cv.mod))
shapiro.test(resid(cv.mod))

#### plot ####
cv.plots <- ggplot(cvs1, aes(x = Plant.identity, y = cv.bio, fill = Nutrients)) +
  geom_col(position = "dodge") +
  facet_wrap(~wneighbour, scales = "free") +
  #geom_errorbar(position = "dodge", aes(ymin = mean.percent.final.size-se.percent.final.size, ymax = mean.percent.final.size+se.percent.final.size)) +
  ylab("CV of total biomass (g)") +
  xlab("Species identity") +
  theme_classic(base_size = 20) +
  scale_x_discrete(limits = c("Gai ari", "Tra dub", "Hed alp", "Bou gra", "Bro ine", "Ely tra",
                              "Koe mac", "Poa pra", "Sti com", "Sti vir", "Ane pat", "Geu tri", "Pot arg"),
                   guide = guide_axis(angle = 90)) + 
  #order species
  geom_hline(
    data = mean_cv_bio_by_wneighbour,
    aes(yintercept = mean_cv_bio),
    linetype = "dashed", 
    color = "red",
    size = 1)
cv.plots

#export
ggsave(filename = "growth cvs.png", 
       cv.plots,
       path = "figures/",
       width = 8,
       height = 6,
       units = "in"
)

#### bootstrap analysis ####
# Initialize storage for bootstrap results
bootstrap_results <- list()

# Set the number of bootstrap iterations
n_bootstrap <- 1000

# Perform the bootstrap
set.seed(123) # For reproducibility
for (i in 1:n_bootstrap) {
  # Step 1: Sample 3 observations per Plant.identity*Nutrients*wneighbour combination
  boot_sample <- datdat %>%
    group_by(Plant.identity, Nutrients, wneighbour) %>%
    sample_n(size = 3, replace = TRUE) %>%
    ungroup()
  
  # Step 2: Calculate mean, sd, and cv for each group
  boot_cvs <- boot_sample %>%
    group_by(Plant.identity, wneighbour, Nutrients) %>%
    summarise(
      mean.bio = mean(total.biomass, na.rm = TRUE),
      sd.bio = sd(total.biomass, na.rm = TRUE),
      cv.bio = sd.bio / mean.bio,
      .groups = "drop"
    )
  
  # Step 3: Fit the linear mixed-effects model
  model <- tryCatch({
    lmer(cv.bio ~ wneighbour * Nutrients + (1 | Plant.identity), data = boot_cvs)
  }, error = function(e) NULL) # Handle errors gracefully
  
  # Step 4: Store the model results
  if (!is.null(model)) {
    bootstrap_results[[i]] <- summary(model)$coefficients
  } else {
    bootstrap_results[[i]] <- NULL
  }
}

# Analyze the bootstrap results
# Combine results into a single data frame
bootstrap_summary <- do.call(rbind, bootstrap_results)

# Inspect the bootstrap summary
str(bootstrap_summary)

#### compare bootstrap results to initial model ####
# Extract coefficients from initial model
initial_coefficients <- summary(cv.mod)$coefficients

# Combine bootstrap coefficients into a data frame
bootstrap_coefficients <- do.call(rbind, bootstrap_results) %>%
  as.data.frame()
str(bootstrap_coefficients)

#update row names
colnames(bootstrap_coefficients) <- colnames(initial_coefficients)

# Summarize bootstrap distributions
bootstrap_summary <- bootstrap_coefficients %>%
  summarise(across(everything(), list(
    mean = ~mean(.),
    sd = ~sd(.),
    lower = ~quantile(., 0.025),
    upper = ~quantile(., 0.975)
  )))
head(bootstrap_summary)

# Check if initial coefficients fall within bootstrap confidence intervals
comparison <- data.frame(
  Coefficient = rownames(initial_coefficients),
  Initial = initial_coefficients[, "Estimate"],
  Bootstrap_Mean = bootstrap_summary[, "Estimate_mean"],
  Lower_CI = bootstrap_summary[, "Estimate_lower"],
  Upper_CI = bootstrap_summary[, "Estimate_upper"],
  In_CI = initial_coefficients[, "Estimate"] >= bootstrap_summary[, "Estimate_lower"] &
    initial_coefficients[, "Estimate"] <= bootstrap_summary[, "Estimate_upper"]
)

# View comparison results
print(comparison)

##prepare to export
comparison2 <- comparison %>%
  mutate(across(2:5, ~ round(.x, 3)))

##export
write_csv(comparison2, "output/bootstrap estimate comparison.csv")

# Melt the bootstrap coefficients for plotting
library(reshape2)
bootstrap_long <- melt(bootstrap_coefficients, variable.name = "Coefficient", value.name = "Value")

# Plot the bootstrap distributions
ggplot(bootstrap_long, aes(x = Value)) +
  geom_density(fill = "blue", alpha = 0.4) +
  geom_vline(data = data.frame(Coefficient = rownames(initial_coefficients), 
                               Value = initial_coefficients[, "Estimate"]),
             aes(xintercept = Value), color = "red", linetype = "dashed") +
  facet_wrap(~ Coefficient, scales = "free") +
  labs(title = "Bootstrap Distributions vs. Initial Coefficients",
       x = "Coefficient Value", y = "Density") +
  theme_minimal()

#compare distributions
comparison$Z_Score <- (comparison$Initial - comparison$Bootstrap_Mean) / bootstrap_summary[, "Estimate_sd"]

comparison$P_Value <- 2 * pmin(
  mean(bootstrap_coefficients <= comparison$Initial),
  mean(bootstrap_coefficients >= comparison$Initial)
)
comparison$P_Value
