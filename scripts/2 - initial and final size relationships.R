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

#setwd("C:/Users/emilyh/Documents/GitHub/Mesocosm_networks")
setwd("~/Documents/GitHub/Mesocosm networks/Mesocosm_networks")

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

#### correlations among initial size and CE ####
##import data
CE.comparison <- read_rds("output/competitive effect df.rds")
str(CE.comparison)

#check distribution of variables
hist(CE.comparison$CE)
CE.comparison$logCE <- log(CE.comparison$CE)
hist(CE.comparison$logCE)

hist(CE.comparison$initial.height)

hist(CE.comparison$initial_biomass)
CE.comparison$loginitial_biomass <- log(CE.comparison$initial_biomass)
hist(CE.comparison$loginitial_biomass)

#correlation test
cor.test(CE.comparison$initial.height, CE.comparison$CE, method = "pearson") #p-value = 0.33
cor.test(CE.comparison$initial.height, CE.comparison$initial_biomass, method = "pearson") #p-value = 0.48
cor.test(CE.comparison$initial_biomass, CE.comparison$CE, method = "pearson") #p-value = 0.57

#model
initial.height.vs.CE <- lmer(CE ~ initial.height*Nutrients + (1|Plant.identity) + (1|Block), data = CE.comparison)
summary(initial.height.vs.CE)
Anova(initial.height.vs.CE)
qqnorm(resid(initial.height.vs.CE))
qqline(resid(initial.height.vs.CE))
shapiro.test(resid(initial.height.vs.CE))

initial.biomass.vs.CE <- lmer(CE ~ initial_biomass*Nutrients + (1|Plant.identity) + (1|Block), data = CE.comparison)
summary(initial.biomass.vs.CE)
Anova(initial.biomass.vs.CE)
qqnorm(resid(initial.biomass.vs.CE))
qqline(resid(initial.biomass.vs.CE))
shapiro.test(resid(initial.biomass.vs.CE))

#plot
CE.initial.height.plot <- ggplot(data = CE.comparison, aes(x = initial.height, y = CE, colour = Nutrients)) +
  geom_point() + 
  geom_smooth(method = lm, se = FALSE, fullrange = FALSE) +
  xlab("Initial height (cm)") +
  ylab("Neighbour effect value") +
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p(height) = 0.870")), 
           hjust = 1, vjust = 1, size = 5) +
  theme_classic(base_size = 18)
CE.initial.height.plot

CE.initial.biomass.plot <- ggplot(data = CE.comparison, aes(x = initial_biomass, y = CE, colour = Nutrients)) +
  geom_point() + 
  geom_smooth(method = lm, se = FALSE, fullrange = FALSE) +
  xlab("Initial biomass estimate (g)") +
  ylab(" ") +
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p(biomass) = 0.064")), 
           hjust = 1, vjust = 1, size = 5) +
  theme_classic(base_size = 18)
CE.initial.biomass.plot

#### correlations among initial size and CR ####
##rename species column in initial.bio
initial.bio <- initial.bio %>% rename(Plant.identity = Species)

##import data
CR.comparison <- read_rds("output/NNR df.rds") %>%
  left_join(initial.bio, by = c("Pot", "Plant.identity", "pot.position")) %>%
  dplyr::select(!Block.y) %>%
  dplyr::rename(initial.height = Height,
                Block = Block.x) %>%
  filter(initial_biomass > 0) #remove impossible numbers
str(CR.comparison)

#check distribution of variables
hist(CR.comparison$NNR)
hist(CR.comparison$initial.height)
hist(CR.comparison$initial_biomass)
CR.comparison$loginitial_biomass <- log(CR.comparison$initial_biomass)
hist(CR.comparison$loginitial_biomass)

#correlation test
cor.test(CR.comparison$initial.height, CR.comparison$NNR, method = "pearson") #p-value = 0.01
cor.test(CR.comparison$initial.height, CR.comparison$initial_biomass, method = "pearson") #p-value = 0.944
cor.test(CR.comparison$initial_biomass, CR.comparison$NNR, method = "pearson") #p-value < 0.001

#model
initial.height.vs.CR <- lmer(NNR ~ initial.height*Nutrients + (1|Plant.identity) + (1|Block), data = CR.comparison)
summary(initial.height.vs.CR)
Anova(initial.height.vs.CR)
qqnorm(resid(initial.height.vs.CR))
qqline(resid(initial.height.vs.CR))
shapiro.test(resid(initial.height.vs.CR))

initial.biomass.vs.CR <- lmer(NNR ~ initial_biomass*Nutrients + (1|Plant.identity) + (1|Block), data = CR.comparison)
summary(initial.biomass.vs.CR)
Anova(initial.biomass.vs.CR)
qqnorm(resid(initial.biomass.vs.CR))
qqline(resid(initial.biomass.vs.CR))
shapiro.test(resid(initial.biomass.vs.CR))

#plot
CR.initial.height.plot <- ggplot(data = CR.comparison, aes(x = initial.height, y = NNR, colour = Nutrients)) +
  geom_point() + 
  geom_smooth(method = lm, se = FALSE, fullrange = FALSE) +
  xlab("Initial height (cm)") +
  ylab("Neighbour response value") +
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p(height) < 0.001")), 
           hjust = 1, vjust = 1, size = 5) +
  theme_classic(base_size = 18)
CR.initial.height.plot

CR.initial.biomass.plot <- ggplot(data = CR.comparison, aes(x = initial_biomass, y = NNR, colour = Nutrients)) +
  geom_point() + 
  geom_smooth(method = lm, se = FALSE, fullrange = FALSE) +
  xlab("Initial biomass estimate (g)") +
  ylab(" ") +
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p(interaction) < 0.001")), 
           hjust = 1, vjust = 1, size = 5) +
  theme_classic(base_size = 18)
CR.initial.biomass.plot

#### put plots together ####
initial.size.plots <- cowplot::plot_grid(
  CE.initial.height.plot + theme(legend.position = "none"), CE.initial.biomass.plot + theme(legend.position = "none"),
  CR.initial.height.plot + theme(legend.position = "none"), CR.initial.biomass.plot + theme(legend.position = "none"),
  labels = 'auto')
initial.size.plots

#extract legend
legend <- get_legend(
  # create some space to the left of the legend
  CE.initial.height.plot + theme(legend.box.margin = margin(0, 0, 0, 12))
)

#add legend to plot
initial.size.all_plot <- plot_grid(initial.size.plots, 
                                   legend, 
                                   rel_widths = c(3, .5),
                                   nrow = 1)
initial.size.all_plot

#export
ggsave(filename = "initial size and comp plot.png", 
       initial.size.all_plot,
       path = "figures/",
       width = 10,
       height = 8,
       units = "in"
)
