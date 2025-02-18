## Calculate neighbour response for individs, create competitive hierarchies, and count instances of facilitation
## Author: Emily H
## Created: July 16. 2024
## Last edited: February 17, 2025

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
palette_colors <- ggthemr::swatch()

#setwd("C:/Users/emilyh/Documents/GitHub/Mesocosm_networks")
setwd("~/Documents/GitHub/Mesocosm networks/Mesocosm_networks")

#### import data ####
alone.df <- read_rds("output/alone_df_new.rds")
str(alone.df)

neighbour.df <- read_rds("output/meso biomass changes.rds")
str(neighbour.df)

#### calculate neighbour response ####
### NOTE that this code uses total biomass, not biomass change
## create new df of the variables I want
dat <- neighbour.df %>%
  left_join(alone.df, by = c("Plant.identity", "Nutrients")) %>% #add alone data to df
  dplyr::select(Pot, Block, pot.position, Plant.identity, Neighbour.identity, Nutrients, alone.replicates, total.biomass, alone.mean.total.biomass) %>% #select columns I want
  mutate(NNR = log(total.biomass/alone.mean.total.biomass)) %>% #calculate net neighbour (competitive) response
  filter(Neighbour.identity != "alone")
str(dat)

#export file
write_rds(dat, "output/NNR df.rds")

#### find number of unique combinations ####
unique_combinations <- unique(dat[, c("Plant.identity", "Neighbour.identity", "Nutrients")]) # Create a data frame with unique combinations of the three columns
num_unique_combinations <- nrow(unique_combinations) # Count the number of unique combinations
num_unique_combinations

#### prevalence of positive vs. negative interactions ####
interaction.counts <- dat %>%
  group_by(Nutrients) %>%
  summarise(facilitation = sum(NNR>0, na.rm = TRUE),
            competition = sum(NNR<0, na.rm = TRUE),
            total.interactions = length(NNR)) %>%
  mutate(percent.pos = facilitation/total.interactions*100,
         percent.neg = competition/total.interactions*100) %>%
  pivot_longer(cols = c(facilitation, competition), names_to = "Interaction.type", values_to = "Interaction.count")
print(interaction.counts)

#find min and max NNR values
min(dat$NNR, na.rm = TRUE)
## as proportional reduction
exp(min(dat$NNR, na.rm = TRUE))*100

max(dat$NNR, na.rm = TRUE)
## as proportional reduction
exp(max(dat$NNR, na.rm = TRUE))*100

#### plot ####
CR.interaction.counts.plot <-ggplot(dat = interaction.counts, aes(x=Interaction.type, y = Interaction.count, fill = Nutrients)) + 
  geom_col(position = "stack") + 
  xlab("Interaction type") +
  ylab("Count of interactions") +
  #ylim(0,800) +
  labs(title = "Neighbour response") +
  theme_classic(base_size = 18)
CR.interaction.counts.plot

## export plot
ggsave(filename = "CR interaction counts.png", 
       CR.interaction.counts.plot,
       path = "figures/",
       width = 8,
       height = 6,
       units = "in"
)

#### find mean facilitative and competitive NR values ####
mean.comp.NR <- dat %>%
  filter(NNR < 0) %>%
  group_by(Nutrients) %>%
  summarise(mean = mean(NNR, na.rm = TRUE),
            se = sd(NNR)/sqrt(length(NNR)),
            median = median(NNR, na.rm = TRUE))
print(mean.comp.NR)

mean.facil.NR <- dat %>%
  filter(NNR > 0) %>%
  group_by(Nutrients) %>%
  summarise(mean = mean(NNR, na.rm = TRUE),
            se = sd(NNR)/sqrt(length(NNR)),
            median = median(NNR, na.rm = TRUE))
print(mean.facil.NR)