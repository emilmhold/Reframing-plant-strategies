## Clean mesocosm biomass data
## Author: Emily H
## Created: July 15. 2024
## Last edited: December 17, 2024

#install.packages("tidyverse")
#install.packages("vegan")

library(vegan)
library(tidyverse)

#setwd("C:/Users/emilyh/Documents/GitHub/Mesocosm_networks")
setwd("~/Documents/GitHub/Mesocosm networks/Mesocosm_networks")

options(scipen = 999)

#### import data ####
## import mesocosm data
meso1 <- read.csv("data/2022 Mesocosm AB and BG biomasses-2.csv") %>%
  dplyr::select(Pot, Block, Nutrients, Position.in.pot, Plant.identity, Neighbour.identity, Dead.at.harvest.,
                Aboveground.biomass..g., Belowground.biomass..g., Total.biomass..g.) %>%
  dplyr::rename(AB.biomass = Aboveground.biomass..g.,
                BG.biomass = Belowground.biomass..g.,
                total.biomass = Total.biomass..g.,
                pot.position = Position.in.pot) %>%
  subset(Plant.identity != "") %>% #removes cells where no plant grew (these pots were transferred to the alone treatment)
  mutate(Neighbour.identity = ifelse(Neighbour.identity %in% "", "alone", Neighbour.identity)) %>%
  # clarify which plants where grown alone
  mutate(Dead.alive = if_else(Dead.at.harvest. == "Y",1,0)) %>%
  dplyr::select(!Dead.at.harvest.) %>% 
  filter(Plant.identity != "Hel hoo" & Neighbour.identity != "Hel hoo") %>% #remove hel hoo data because I don't have initial biomasses for those plants
  filter(Plant.identity != "Ast lae" & Neighbour.identity != "Ast lae") %>% #remove species with data only in the unfertilized treatment
  filter(Plant.identity != "Ast fal" & Neighbour.identity != "Ast fal") %>% #remove species with data only in the unfertilized treatment
  filter(Plant.identity != "Fra vir" & Neighbour.identity != "Fra vir")  #remove species with data only in the unfertilized treatment
str(meso1) #1141 plants

## import initial biomass data
initial.biomass <- read_rds("output/meso initial biomasses.rds") %>%
  rename(Plant.identity = Species) %>%
  dplyr::select(Pot, Plant.identity, pot.position, initial_biomass)
str(initial.biomass)

## join initial biomass data with meso df 
meso <- meso1 %>%
  left_join(initial.biomass, by = c("Pot", "pot.position")) %>%
  filter(Plant.identity.x != "Hel hoo") %>% #remove hel hoo data because I don't have initial biomasses for those plants
  dplyr::select(!Plant.identity.y) %>%
  rename(Plant.identity = Plant.identity.x) %>%
  mutate(biomass.change = total.biomass - initial_biomass) %>% # calculate change in biomass
  filter(biomass.change > 0) #exclude negative values
str(meso) #955 plants

##export df
write_rds(meso1, "output/meso biomass changes.rds")

## create df of alone plants
alone.df <- meso1 %>%
  filter(Neighbour.identity == "alone") %>% #select alone plants - 154 alone plants
  group_by(Plant.identity, Nutrients) %>% 
  summarise(alone.replicates = length(total.biomass),
            #alone.mean.ab.biomass = mean (AB.biomass, na.rm = TRUE),
            #alone.se.ab.biomass = sd(AB.biomass)/sqrt(length(AB.biomass)),
            #alone.mean.bg.biomass = mean (BG.biomass, na.rm = TRUE),
            #alone.se.bg.biomass = sd(BG.biomass)/sqrt(length(BG.biomass)),
            alone.mean.total.biomass = mean (total.biomass, na.rm = TRUE),
            alone.se.total.biomass = sd(total.biomass)/sqrt(length(total.biomass))) #calculate spp mean biomass values
str(alone.df)

## create df of neighbour plants
neighbour.df <- meso1 %>%
  filter(Neighbour.identity != "alone") %>% #select plants grown with neighbours
  group_by(Plant.identity, Neighbour.identity, Nutrients) %>%
  summarise(replicates = length(total.biomass),
            #mean.ab.biomass = mean (AB.biomass),
            #se.ab.biomass = sd(AB.biomass)/sqrt(length(AB.biomass)),
            #mean.bg.biomass = mean (BG.biomass),
            #se.bg.biomass = sd(BG.biomass)/sqrt(length(BG.biomass)),
            mean.total.biomass = mean(total.biomass),
            se.total.biomass = sd(total.biomass)/sqrt(length(total.biomass))) #calculate spp mean biomass values
str(neighbour.df)

#### export dfs ####
write_rds(meso, "output/meso_new.rds")
write_rds(alone.df, "output/alone_df_new.rds")
write_rds(neighbour.df, "output/neighbour_df_new.rds")
