## Calculate neighbour effect
## Author: Emily H
## Created: November 15, 2024
## Last edited: February 17, 2025

#install.packages("tidyverse")
#install.packages("vegan")
#install.packages("reshape2")
#install.packages("lme4")
#install.packages("lmerTest")
#install.packages("car")
#install.packages("ggthemr")
#install.packages("cowplot")

library(vegan)
library(tidyverse)
library(reshape2)
library(lme4)
library(lmerTest)
library(car)
library(ggthemr)
library(cowplot)

ggthemr("fresh")
palette_colors <- ggthemr::swatch()

#setwd("C:/Users/emilyh/Documents/GitHub/Mesocosm_networks")
setwd("~/Documents/GitHub/Mesocosm networks/Mesocosm_networks")

#### import data ####
alone.df <- read_rds("output/alone_df_new.rds") %>%
  rename(Neighbour.identity = Plant.identity) #rename to facilitate joining by neighbour identity
str(alone.df)

neighbour.df <- read_rds("output/meso biomass changes.rds")
str(neighbour.df)

# read in initial size measures and biomass estimates for meso plants
initial.bio <- read_rds("output/meso initial biomasses.rds") %>%
  rename(Plant.identity = Species) %>%
  filter(initial_biomass > 0) #remove impossible numbers
str(initial.bio)

#### calculate neighbour effect ####
## create new df of the variables I want
dat.CE <- neighbour.df %>%
  left_join(alone.df, by = c("Neighbour.identity", "Nutrients")) %>% #add alone data to df
  dplyr::select(Pot, Block, pot.position, Plant.identity, Neighbour.identity, Nutrients, total.biomass, alone.mean.total.biomass) %>% #select columns I want
  filter(Neighbour.identity != "alone") %>% #just to be sure
  mutate(CE = log(total.biomass/alone.mean.total.biomass)) # calculate competitive effect for focal plant
str(dat.CE)
hist(dat.CE$CE)

# export
write_rds(dat.CE, "output/CE df.rds")

#### prevalence of positive vs. negative interactions ####
CE.interaction.counts <- dat.CE %>%
  group_by(Nutrients) %>%
  summarise(facilitation = sum(CE>0, na.rm = TRUE),
            competition = sum(CE<0, na.rm = TRUE),
            total.interactions = sum(!is.na(CE))) %>%
  mutate(percent.pos = facilitation/total.interactions*100,
         percent.neg = competition/total.interactions*100) %>%
  pivot_longer(cols = c(facilitation, competition), names_to = "Interaction.type", values_to = "Interaction.count")
print(CE.interaction.counts)

#find min and max CE values
min(dat.CE$CE, na.rm = TRUE)
## as proportional growth reduction
exp(min(dat.CE$CE, na.rm = TRUE))*100

max(dat.CE$CE, na.rm = TRUE)
## as proportional growth enhancement
exp(max(dat.CE$CE, na.rm = TRUE))*100

##plot
CE.interaction.counts.plot <-ggplot(dat = CE.interaction.counts, aes(x=Interaction.type, y = Interaction.count, fill = Nutrients)) + 
  geom_col(position = "stack") + 
  xlab("Interaction type") +
  ylab("Count of interactions") +
  labs(title = "Neighbour effect") +
  ylim(0,800) +
  theme_classic(base_size = 18)
CE.interaction.counts.plot

## export plots
ggsave(filename = "CE interaction counts.png", 
       CE.interaction.counts.plot,
       path = "figures/",
       width = 8,
       height = 6,
       units = "in"
)

#### compare NE to plant initial size ####
#select initial size columns I want
initial.bio.for.CE <- initial.bio %>%
  dplyr::select(Pot, Plant.identity, pot.position, Height, initial_biomass)
str(initial.bio.for.CE)

#join dfs
CE.comparison <- dat.CE %>%
  left_join(initial.bio.for.CE, by = c("Pot","Plant.identity", "pot.position")) %>%
  rename(initial.height = Height)
str(CE.comparison)

#export df
write_rds(CE.comparison, "output/competitive effect df.rds")

#### plots ####
CE.vs.initial.height <- ggplot(data= CE.comparison, aes(x = initial.height, y = CE, colour = Nutrients)) +
  geom_point(position = "jitter") + 
  #geom_smooth(method = lm, se = FALSE, fullrange = FALSE) +
  xlab("Initial height (cm)") + 
  ylab("Neighbour effect") +
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p(height) = 0.162")), hjust = 1, vjust = 1, size = 5) +
  theme_classic(base_size = 18)
CE.vs.initial.height

CE.vs.initial.biomass <- ggplot(data= CE.comparison, aes(x = initial_biomass, y = CE, colour = Nutrients)) +
  geom_point(position = "jitter") + 
  xlab("Initial biomass (g)") + 
  ylab("Neighbour effect") +
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p(height) = 0.407")), hjust = 1.1, vjust = 1.1, size = 5) +
  theme_classic(base_size = 18)
CE.vs.initial.biomass

#### put it together ####
CE.initialsize.multiplot <- cowplot::plot_grid(
  CE.vs.initial.biomass + theme(legend.position = "none"),
  CE.vs.initial.height + theme(legend.position = "none"),
  labels = 'auto')
CE.initialsize.multiplot

#extract legend
legend <- get_legend(
  # create some space to the left of the legend
  CE.vs.initial.biomass + theme(legend.box.margin = margin(0, 0, 0, 12))
)

#add legend to plot
all_plot <- plot_grid(CE.initialsize.multiplot, legend, rel_widths = c(3, .6))
all_plot

ggsave(filename = "figures/Competitive effect and initial plant size.png", 
       all_plot,
       width = 15,
       height = 6,
       units = "in"
)
