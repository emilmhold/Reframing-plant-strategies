## Meso plants initial/final size relationships
## Author: Emily H
## Created: March 1, 2025
## Last edited: March 1, 2025

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

#### correlations among initial size and CE ####
##import data
CE.comparison <- read_rds("output/CE df.rds")
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
