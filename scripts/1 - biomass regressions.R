## Biomass regressions for focal species plants
## Author: Emily H
## Created: July 16. 2024
## Last edited: November 17, 2024

#install.packages("tidyverse")
#install.packages("vegan")

library(vegan)
library(tidyverse)

setwd("/Users/emilyholden/Documents/GitHub/Analyses/Reframing-plant-strategies")

#### import data ####
## read in biomass data for sacrificed plants
bio <- read.csv("data/Copy of 2022 initial measures for meso biomass regressions.csv") %>%
  dplyr::select(!Date.Planted, Count.of.days) %>%
  rename(height = Plant.Height..cm.,
         width = Plant.longest.width..cm.,
         no.leaves = Plant...of.leaves,
         tot.biomass = Total.biomass..g.) %>%
  dplyr::select(Species, height, width, no.leaves, tot.biomass) %>%
  mutate_at(2:5, as.numeric) %>%
  mutate(Species = ifelse(Species == "Poa Pra", "Poa pra", Species)) #fix poa pra spelling
str(bio)

# read in initial size measures for meso plants
meso.bio <- read.csv("data/2022 Mesocosm Initial Size Measures_Data_for analysis.csv") %>%
  rename(height = Height..cm.,
         width = Longest.width..cm.,
         no.leaves = No.leaves) %>%
  dplyr::select(Pot, Block, Plant.identity, pot.position, height, width, no.leaves) %>%
  mutate(Plant.identity = ifelse(Plant.identity == "Poa Pra", "Poa pra", Plant.identity)) #fix poa pra spelling
str(meso.bio)

#### calculate R2 and regression coefficients for all spp*nutrient combinations ####
# Get unique species
species_list <- unique(bio$Species)

# Create an empty data frame to store the coefficients and R^2 values
coefficients_df <- data.frame(
  Species = character(),
  Intercept_Height = numeric(),
  Height_Coeff = numeric(),
  R2_Height = numeric(),
  Intercept_Width = numeric(),
  Width_Coeff = numeric(),
  R2_Width = numeric(),
  Intercept_No_Leaves = numeric(),
  No_Leaves_Coeff = numeric(),
  R2_No_Leaves = numeric(),
  Intercept_Height_Width = numeric(),
  Height_Width_Coeff = numeric(),
  R2_Height_Width = numeric(),
  Intercept_Height_No_Leaves = numeric(),
  Height_No_Leaves_Coeff = numeric(),
  R2_Height_No_Leaves = numeric(),
  Intercept_Width_No_Leaves = numeric(),
  Width_No_Leaves_Coeff = numeric(),
  R2_Width_No_Leaves = numeric(),
  Intercept_Full = numeric(),
  Full_Coeff = numeric(),
  R2_Full = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each species and perform the regressions
for (species in species_list) {
  # Subset the data for the current species
  species_data <- subset(bio, Species == species)
  
  # Fit the linear models
  model_height <- lm(tot.biomass ~ height, data = species_data)
  model_width <- lm(tot.biomass ~ width, data = species_data)
  model_no_leaves <- lm(tot.biomass ~ no.leaves, data = species_data)
  model_height_width <- lm(tot.biomass ~ height + width, data = species_data)
  model_height_no_leaves <- lm(tot.biomass ~ height + no.leaves, data = species_data)
  model_width_no_leaves <- lm(tot.biomass ~ width + no.leaves, data = species_data)
  model_full <- lm(tot.biomass ~ height + width + no.leaves, data = species_data)
  
  # Extract coefficients and R^2 values
  coefficients_height <- coef(model_height)
  r2_height <- summary(model_height)$r.squared
  
  coefficients_width <- coef(model_width)
  r2_width <- summary(model_width)$r.squared
  
  coefficients_no_leaves <- coef(model_no_leaves)
  r2_no_leaves <- summary(model_no_leaves)$r.squared
  
  coefficients_height_width <- coef(model_height_width)
  r2_height_width <- summary(model_height_width)$r.squared
  
  coefficients_height_no_leaves <- coef(model_height_no_leaves)
  r2_height_no_leaves <- summary(model_height_no_leaves)$r.squared
  
  coefficients_width_no_leaves <- coef(model_width_no_leaves)
  r2_width_no_leaves <- summary(model_width_no_leaves)$r.squared
  
  coefficients_full <- coef(model_full)
  r2_full <- summary(model_full)$r.squared
  
  # Create a data frame for the current species
  species_coefficients <- data.frame(
    Species = species,
    Intercept_Height = coefficients_height[1],
    Height_Coeff = coefficients_height[2],
    R2_Height = r2_height,
    Intercept_Width = coefficients_width[1],
    Width_Coeff = coefficients_width[2],
    R2_Width = r2_width,
    Intercept_No_Leaves = coefficients_no_leaves[1],
    No_Leaves_Coeff = coefficients_no_leaves[2],
    R2_No_Leaves = r2_no_leaves,
    Intercept_Height_Width = coefficients_height_width[1],
    Height_Width_Coeff = coefficients_height_width[2],
    R2_Height_Width = r2_height_width,
    Intercept_Height_No_Leaves = coefficients_height_no_leaves[1],
    Height_No_Leaves_Coeff = coefficients_height_no_leaves[2],
    R2_Height_No_Leaves = r2_height_no_leaves,
    Intercept_Width_No_Leaves = coefficients_width_no_leaves[1],
    Width_No_Leaves_Coeff = coefficients_width_no_leaves[2],
    R2_Width_No_Leaves = r2_width_no_leaves,
    Intercept_Full = coefficients_full[1],
    Full_Coeff = coefficients_full[2],
    R2_Full = r2_full,
    stringsAsFactors = FALSE
  )
  
  # Bind the species_coefficients data frame to the main coefficients_df data frame
  coefficients_df <- rbind(coefficients_df, species_coefficients)
}

str(coefficients_df)

#### Identify most predictive models ####
greatest_r2_df <- data.frame(
  Species = character(),
  Best_Variable = character(),
  Best_R2 = numeric(),
  Best_Coefficient = numeric(),
  Best_Intercept = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each species to find the variable with the greatest R^2
for (i in 1:nrow(coefficients_df)) {
  species <- coefficients_df$Species[i]
  
  # Extract R^2 and coefficients
  r2_values <- c(coefficients_df$R2_Height[i], coefficients_df$R2_Width[i], coefficients_df$R2_No_Leaves[i], coefficients_df$R2_Height_Width[i], coefficients_df$R2_Height_No_Leaves[i], coefficients_df$R2_Width_No_Leaves[i], coefficients_df$R2_Full[i])
  coefficients <- c(coefficients_df$Height_Coeff[i], coefficients_df$Width_Coeff[i], coefficients_df$No_Leaves_Coeff[i], coefficients_df$Height_Width_Coeff[i], coefficients_df$Height_No_Leaves_Coeff[i], coefficients_df$Width_No_Leaves_Coeff[i], coefficients_df$Full_Coeff[i])
  intercepts <- c(coefficients_df$Intercept_Height[i], coefficients_df$Intercept_Width[i], coefficients_df$Intercept_No_Leaves[i], coefficients_df$Intercept_Height_Width[i], coefficients_df$Intercept_Height_No_Leaves[i], coefficients_df$Intercept_Width_No_Leaves[i], coefficients_df$Intercept_Full[i])
  variable_names <- c("Height", "Width", "No_Leaves", "Height_Width", "Height_No_Leaves", "Width_No_Leaves", "Full")
  
  # Ensure r2_values are not NA
  if (all(is.na(r2_values))) {
    next  # Skip to the next iteration if all are NA
  }
  
  # Find the index of the maximum R^2
  max_index <- which.max(r2_values)
  
  # Create a new row with the species, best variable, best R^2, best coefficient, and best intercept
  greatest_r2_df <- rbind(greatest_r2_df, data.frame(
    Species = species,
    Best_Variable = variable_names[max_index],
    Best_R2 = r2_values[max_index],
    Best_Coefficient = coefficients[max_index],
    Best_Intercept = intercepts[max_index],
    stringsAsFactors = FALSE
  ))
}

str(greatest_r2_df)

##reformat for export
greatest_r2_df_for_export <- greatest_r2_df %>%
  mutate(across(3:5, round, 3)) %>% #round values to 3 decimal places
  mutate(Best_Variable = ifelse(Best_Variable == "Full", "Height + Width + Number of leaves", Best_Variable)) %>% #fix spelling
  mutate(Best_Variable = ifelse(Best_Variable == "Height_Width", "Height + Width", Best_Variable)) %>%
  mutate(Best_Variable = ifelse(Best_Variable == "Height_No_Leaves", "Height + Number of leaves", Best_Variable)) %>% #fix spelling
  rename("Best size variable" = Best_Variable,
         "Best R2 value" = Best_R2,
         "Best slope coefficient" = Best_Coefficient,
         "Best y-intercept" = Best_Intercept)

##export table
write_csv(greatest_r2_df_for_export, "output/biomass regression r2 values.csv")

#### calculate initial biomass for meso plants ####
meso.bio2 <- meso.bio %>%
  rename(Species = Plant.identity,
         No_Leaves = no.leaves,
         Height= height,
         Width = width) %>% #rename columns to match r2 df
  full_join(greatest_r2_df, by = "Species") %>%
  filter(Species != "",
         Width != "") %>%
  mutate(Height = as.numeric(Height),
         Width = as.numeric(Width),
         No_Leaves = as.numeric(No_Leaves)) %>% #convert variables to numeric
  mutate(initial_biomass = case_when(
    Best_Variable == "Height" ~ Best_Coefficient * Height + Best_Intercept,
    Best_Variable == "Width" ~ Best_Coefficient * Width + Best_Intercept,
    Best_Variable == "No_Leaves" ~ Best_Coefficient * No_Leaves + Best_Intercept,
    Best_Variable == "Height_Width" ~ Best_Coefficient * (Height + Width) + Best_Intercept,
    Best_Variable == "Height_No_Leaves" ~ Best_Coefficient * (Height + No_Leaves) + Best_Intercept,
    Best_Variable == "Full" ~ Best_Coefficient * (Height + Width + No_Leaves) + Best_Intercept)) #calculate initial biomass
str(meso.bio2)

## export df
write_rds(meso.bio2, "output/meso initial biomasses.rds")
