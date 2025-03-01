## Meso plants initial/final size relationships
## Author: Emily H
## Created: Novemeber 14, 2024
## Last edited: February 17, 2025

#install.packages("tidyverse")
#install.packages("vegan")
#install.packages("car")
#install.packages("cowplot")
#install.packages("lme4")
#install.packages("lmerTest")
#devtools::install_github("cttobin/ggthemr")

library(vegan)
library(tidyverse)
library(car)
library(cowplot)
library(ggthemr)
library(lme4)
library(lmerTest)

ggthemr("fresh")
palette_colors <- ggthemr::swatch()

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

# read in initial size measures and biomass estimates for meso plants
initial.bio <- read_rds("output/meso initial biomasses.rds")
str(initial.bio)

# read in final size measures 
final.bio <- read.csv("data/2022 Mesocosm AB and BG biomasses-2.csv") %>%
  dplyr::select(Pot, Block, Nutrients, Position.in.pot, Plant.identity, Neighbour.identity, Dead.at.harvest.,
                Aboveground.biomass..g., Belowground.biomass..g., Total.biomass..g.) %>%
  dplyr::rename(AB.biomass = Aboveground.biomass..g.,
                BG.biomass = Belowground.biomass..g.,
                total.biomass = Total.biomass..g.,
                pot.position = Position.in.pot) %>%
  subset(Plant.identity != "") %>% #removes cells where the died or was not found at harvest
  mutate(Neighbour.identity = ifelse(Neighbour.identity %in% "", "alone", Neighbour.identity)) %>%
  # clarify which plants where grown alone
  mutate(Dead.alive = if_else(Dead.at.harvest. == "Y",1,0)) %>%
  dplyr::select(!Dead.at.harvest.)
str(final.bio)

#### CV for species at start and finish ####
## join initial biomass data with meso df 
dat <- final.bio %>%
  left_join(initial.bio, by = c("Pot", "pot.position")) %>%
  dplyr::select(!c(Species, Block.y)) %>% #drop duplicated columns
  rename(Block = Block.x) %>%
  filter(Plant.identity != "Hel hoo") %>% #remove hel hoo data because I don't have initial biomasses for those plants
  filter(Neighbour.identity != "Hel hoo") %>%
  mutate(initial_biomass = if_else(initial_biomass < 0, NA_real_, initial_biomass)) %>% #replace negative initial.bio values with NA
  mutate(wneighbour = if_else(Neighbour.identity == "alone", "alone", "with neighbour"))
str(dat)

## check on an outlier
brome <- dat %>%
  filter(Plant.identity == "Bro ine") %>%
  filter(Nutrients == "Unfertilized") %>%
  filter(Neighbour.identity == "alone") %>%
  summarise(mean.initial.bio = mean(initial_biomass, na.rm = TRUE),
            sd.initial.bio = sd(initial_biomass, na.rm = TRUE),
            cv.initial.bio = sd.initial.bio/mean.initial.bio,
            mean.initial.height = mean(Height, na.rm = TRUE),
            sd.initial.height = sd(Height, na.rm = TRUE),
            cv.initial.height = mean.initial.height/sd.initial.height) 

## calculate cv for size measures
cvs <- dat %>% 
  group_by(Plant.identity, wneighbour, Nutrients) %>%
  summarise(mean.initial.bio = mean(initial_biomass, na.rm = TRUE),
            sd.initial.bio = sd(initial_biomass, na.rm = TRUE),
            cv.initial.bio = sd.initial.bio/mean.initial.bio,
            mean.initial.height = mean(Height, na.rm = TRUE),
            sd.initial.height = sd(Height, na.rm = TRUE),
            cv.initial.height = mean.initial.height/sd.initial.height,
            mean.final.bio = mean(total.biomass, na.rm = TRUE),
            sd.final.bio = sd(total.biomass, na.rm = TRUE),
            cv.final.bio = sd.final.bio/mean.final.bio)
str(cvs)

initial.cvs <- cvs %>%
  dplyr::select(Plant.identity, wneighbour, Nutrients, mean.initial.bio, sd.initial.bio, cv.initial.bio) %>%
  rename(mean = mean.initial.bio,
         sd = sd.initial.bio,
         cv = cv.initial.bio) %>%
  mutate(Time = "initial biomass", .before = mean)
str(initial.cvs)

final.cvs <- cvs %>%
  dplyr::select(Plant.identity, wneighbour, Nutrients, mean.final.bio, sd.final.bio, cv.final.bio) %>%
  rename(mean = mean.final.bio,
         sd = sd.final.bio,
         cv = cv.final.bio) %>%
  mutate(Time = "final biomass", .before = mean)
str(final.cvs)

cvs.forplot <- rbind(initial.cvs, final.cvs) %>%
  filter(!(Plant.identity == "Ast fal")) %>%
  filter(!(Plant.identity == "Ast lae")) %>% 
  filter(!(Plant.identity == "Fra vir")) #remove species grown only in the unfertilized treatment
str(cvs.forplot)


## plot
cv.boxplot <- ggplot(cvs.forplot, aes(x=Plant.identity, y=cv, fill = Nutrients)) + 
  geom_boxplot() +
  facet_wrap(~Time) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_classic(base_size = 18)
cv.boxplot

ggsave(filename = "final and initial biomass cvs.png", 
       cv.boxplot,
       path = "figures/",
       width = 16,
       height = 8,
       units = "in"
)

#### what % of final size was initial size? ####
size.comparison <- dat %>%
  mutate(percent.final.size = initial_biomass/total.biomass*100) %>% 
  filter(!(Plant.identity == "Geu tri" & #remove a crazy outlier
             wneighbour == "with neighbour" & 
             Nutrients == "Unfertilized" & 
             Pot == 97)) %>%
  group_by(Plant.identity, wneighbour, Nutrients) %>%
  summarise(mean.percent.final.size = mean(percent.final.size),
            se.percent.final.size = sd(percent.final.size)/length(percent.final.size)) %>% 
  unite(Treatment, c("Nutrients", "wneighbour")) %>% 
  filter(!(Plant.identity == "Ast fal")) %>%
  filter(!(Plant.identity == "Ast lae")) %>% 
  filter(!(Plant.identity == "Fra vir")) #remove species grown only in the unfertilized treatment
str(size.comparison)  

## plot
size.comparison.plots <- ggplot(size.comparison, aes(x = Plant.identity, y = mean.percent.final.size, fill = Treatment)) +
  geom_col(position = "dodge") +
  geom_errorbar(position = "dodge", aes(ymin = mean.percent.final.size-se.percent.final.size, ymax = mean.percent.final.size+se.percent.final.size)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  ylab("Mean % of final size comprised of initial size") +
  theme_classic(base_size = 20)
size.comparison.plots ## come back to this

ggsave(filename = "initial percent of final size.png", 
       size.comparison.plots,
       path = "figures/",
       width = 16,
       height = 8,
       units = "in"
)
