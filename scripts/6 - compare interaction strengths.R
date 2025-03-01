## Compare facilitation vs. competition frequency and interaction strength
## Author: Emily H
## Created: September 20, 2024
## Last edited: January 30, 2025

#install.packages("tidyverse")
#install.packages("vegan")
#install.packages("reshape2")
#install.packages("lme4")
#install.packages("lmerTest")
#install.packages("car")
#devtools::install_github('Mikata-Project/ggthemr', force = TRUE)

library(vegan)
library(tidyverse)
library(reshape2)
library(ggthemr)

ggthemr("fresh")

setwd("/Users/emilyholden/Documents/GitHub/Analyses/Reframing-plant-strategies")

#### import data ####
## import NR (neighbour response) data
NR.df <- read_rds("output/NNR df.rds") %>%
  dplyr::select(Pot, Block, pot.position, Plant.identity, Neighbour.identity, Nutrients, total.biomass, NNR)
str(NR.df)

##import NE (neighbour effect) data
NE.df <- read_rds("output/CE df.rds") %>%
  dplyr::select(Pot, pot.position, Plant.identity, Neighbour.identity, Nutrients, CE)
str(NE.df)

##join dfs 
competitive.df <- NR.df %>%
  full_join(NE.df, by = c("Plant.identity", "Neighbour.identity", "Nutrients", "Pot", "pot.position"), relationship = "many-to-many") %>%
  filter(!is.na(total.biomass)) %>%
  mutate(NR = if_else(NNR<0,"Competition", "Facilitation"), # create columns to id interaction type
         NE = if_else(CE<0,"Competition", "Facilitation"),
         abs.NE = abs(CE), #create columns of absolute interaction strength
         abs.NR = abs(NNR)) %>%
  pivot_longer(cols = NR:NE, names_to = "metric", values_to = "interaction.type") %>%
  pivot_longer(cols = abs.NE:abs.NR, names_to = "abs.metric", values_to = "interaction.strength") %>%
  mutate(metric = replace(metric, metric == 'NE', 'Neighbour effect')) %>%
  mutate(metric = replace(metric, metric == 'NR', 'Neighbour response'))
str(competitive.df)    

#### boxplot ####
abs.interactions.boxplot <- ggplot(competitive.df, aes(x=interaction.type, y=interaction.strength, fill = Nutrients)) + 
  geom_boxplot() +
  facet_wrap(~metric) +
  ylab("Absolute interaction strength") +
  xlab("Interaction type") + 
  theme_classic(base_size = 18)
abs.interactions.boxplot


##export image
ggsave(filename = "interaction strength boxplot.png", 
       abs.interactions.boxplot,
       path = "figures/",
       width = 8,
       height = 6,
       units = "in"
)

#### table ####
interactions.summary <- competitive.df %>%
  group_by(interaction.type, metric, Nutrients) %>%
  summarise("Count of interactions" = length(interaction.strength),
            mean.abs.interaction.strength = round(mean(interaction.strength), digits = 3),
            se.abs.interaction.strength = round(sd(interaction.strength)/sqrt(length(interaction.strength)), digits = 3),
            "Median absolute interaction strength" = round(median(interaction.strength), digits = 3)) %>%
  mutate("Mean absolute interaction strength" = paste(mean.abs.interaction.strength, "\u00B1", se.abs.interaction.strength)) %>%
  rename(Metric = metric,
         "Interaction type" = interaction.type,
         "Nutrient treatment" = Nutrients) %>%
  dplyr::select(!c(mean.abs.interaction.strength, se.abs.interaction.strength))
## write table
write.csv(interactions.summary, "output/interaction summary table.csv")
