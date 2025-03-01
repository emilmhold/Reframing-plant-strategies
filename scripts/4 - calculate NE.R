## Calculate neighbour effect
## Author: Emily H
## Created: November 15, 2024
## Last edited: March 1, 2025

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

