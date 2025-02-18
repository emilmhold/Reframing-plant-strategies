## Relationships between neighbour effect and response
## Author: Emily H
## Created: December 16, 2024
## Last edited: February 17, 2025

#install.packages("tidyverse")
#install.packages("vegan")
#install.packages("reshape2")
#install.packages("lme4")
#install.packages("lmerTest")
#install.packages("car")
#install.packages("cowplot")
#install.packages("broom")
#devtools::install_github("cttobin/ggthemr")


library(vegan)
library(tidyverse)
library(reshape2)
library(lme4)
library(lmerTest)
library(car)
library(ggthemr)
library(cowplot)
library(broom)


#setwd("C:/Users/emilyh/Documents/GitHub/Mesocosm_networks")
setwd("~/Documents/GitHub/Mesocosm networks/Mesocosm_networks")

##prep colour palettes
ggthemr("fresh")
palette_colors <- ggthemr::swatch()  #extract the colors from the palette
custom_colors <- palette_colors[c(3, 4)] #select the 3rd and 4th colors for a later plot
# Define a new palette with more colors
new_palette <- define_palette(
  swatch = c("#111111", "#65ADC2", "#233B43", "#E84646", "#C29365", "#362C21", "#316675", "#168E7F", "#109B37", "#F2B134", "#8E5C9E", "#476A6F", "#9D4E3C", "#5E4C3A"),
  gradient = c(lower = "#132B43", upper = "#56B1F7"),
  background = "white"
)

# Set the new palette
ggthemr(new_palette)

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
  full_join(NE.df, by = c("Plant.identity", "Neighbour.identity", "Nutrients", "Pot", "pot.position"), relationship = "many-to-many")
str(competitive.df) 

## summarize for species
competitive.summary <- competitive.df %>%
  group_by(Plant.identity, Neighbour.identity, Nutrients) %>%
  summarise(mean.NR = mean(NNR, na.rm = TRUE),
            se.NR = sd(NNR, na.rm = TRUE)/sqrt(length(NNR)),
            mean.NE = mean(CE, na.rm = TRUE),
            se.NE = sd(CE, na.rm = TRUE)/sqrt(length(CE)))
str(competitive.summary)

#### interaction counts plots ####
## find prevalence of positive vs. negative interactions for NR
NR.interaction.counts <- competitive.df %>%
  group_by(Nutrients) %>%
  summarise(Facilitation = sum(NNR>0, na.rm = TRUE),
            Competition = sum(NNR<0, na.rm = TRUE),
            total.interactions = length(NNR)) %>%
  mutate(percent.pos = Facilitation/total.interactions*100,
         percent.neg = Competition/total.interactions*100) %>%
  pivot_longer(cols = c(Facilitation, Competition), names_to = "Interaction.type", values_to = "Interaction.count") %>%
  mutate(metric = "Neighbour response")
print(NR.interaction.counts)

## find prevalence of positive vs. negative interactions for NE
NE.interation.counts <- competitive.df %>%
  group_by(Nutrients) %>%
  summarise(Facilitation = sum(CE>0, na.rm = TRUE),
            Competition = sum(CE<0, na.rm = TRUE),
            total.interactions = length(CE)) %>%
  mutate(percent.pos = Facilitation/total.interactions*100,
         percent.neg = Competition/total.interactions*100) %>%
  pivot_longer(cols = c(Facilitation, Competition), names_to = "Interaction.type", values_to = "Interaction.count") %>%
  mutate(metric = "Neighbour effect")
print(NE.interation.counts)

## join dfs
interaction.counts <- rbind(NE.interation.counts, NR.interaction.counts)

#find min and max NR values
min(competitive.df$NNR, na.rm = TRUE)
max(competitive.df$NNR, na.rm = TRUE)
min(competitive.df$CE, na.rm = TRUE)
max(competitive.df$CE, na.rm = TRUE)

##plot
interaction.counts.plot <-ggplot(dat = interaction.counts, aes(x=Interaction.type, y = Interaction.count, fill = Nutrients)) + 
  geom_col(position = "stack") + 
  xlab("Interaction type") +
  ylab("Count of interactions") +
  facet_grid(~metric) +
  theme_classic(base_size = 18)
interaction.counts.plot

## export
ggsave(filename = "Interaction counts.png", 
       interaction.counts.plot,
       path = "figures/",
       width = 8,
       height = 6,
       units = "in"
)

#### interactions by species plots ####
## count positive vs. negative interactions for NR by spp
NR.spp.interaction.counts <- competitive.df %>%
  group_by(Plant.identity, Nutrients) %>%
  summarise(Facilitation = sum(NNR>0, na.rm = TRUE),
            Competition = sum(NNR<0, na.rm = TRUE)) %>%
  pivot_longer(cols = c(Facilitation, Competition), names_to = "Interaction.type", values_to = "Interaction.count") %>%
  mutate(metric = "Neighbour response")
str(NR.spp.interaction.counts)

## find prevalence of positive vs. negative interactions for NE
NE.spp.interation.counts <- competitive.df %>%
  group_by(Plant.identity, Nutrients) %>%
  summarise(Facilitation = sum(CE>0, na.rm = TRUE),
            Competition = sum(CE<0, na.rm = TRUE)) %>%
  pivot_longer(cols = c(Facilitation, Competition), names_to = "Interaction.type", values_to = "Interaction.count") %>%
  mutate(metric = "Neighbour effect")
str(NE.spp.interation.counts)

## join dfs
spp.interaction.counts <- rbind(NE.spp.interation.counts, NR.spp.interaction.counts)

##plot
spp.interaction.counts.plot <-ggplot(dat = spp.interaction.counts, aes(x=Plant.identity, y = Interaction.count, fill = Interaction.type)) + 
  geom_col(position = "dodge") + 
  xlab("Species identity") +
  ylab("Count of interactions") +
  labs(fill = "Interaction type") +
  facet_grid(Nutrients ~ metric, switch = "y") +
  theme(strip.placement = "outside") +
  scale_fill_manual(values = custom_colors) +
  scale_x_discrete(limits = c("Gai ari", "Tra dub", "Hed alp", "Bou gra", "Bro ine", "Ely tra",
                              "Koe mac", "Poa pra", "Sti com", "Sti vir", "Ane pat", "Geu tri", "Pot arg"),
                   guide = guide_axis(angle = 90)) +
  theme_classic(base_size = 18)
spp.interaction.counts.plot

## export
ggsave(filename = "Interactions by spp.png", 
       spp.interaction.counts.plot,
       path = "figures/",
       width = 8,
       height = 6,
       units = "in"
)

#### join spp interaction counts and overall interaction counts plots ####
interactions.counts.plot.new <- cowplot::plot_grid(interaction.counts.plot, spp.interaction.counts.plot,
                                                   nrow = 2,
                                                   labels = "auto")
interactions.counts.plot.new

## export
ggsave(filename = "combined interaction counts plots.png", 
       interactions.counts.plot.new,
       path = "figures/",
       width = 8,
       height = 12,
       units = "in"
)

#### NE vs NR models ####
## check distributions
hist(competitive.df$NNR)
hist(competitive.df$CE)

#LMMs
NE.vs.NR <- lmer(mean.NR ~ mean.NE*Nutrients + (1|Plant.identity) , data = competitive.summary)
summary(NE.vs.NR)
print(Anova(NE.vs.NR), digits = 4)
qqnorm(resid(NE.vs.NR))
qqline(resid(NE.vs.NR))
shapiro.test(resid(NE.vs.NR))

NR.vs.NE <- lmer(mean.NE ~ mean.NR*Nutrients + (1|Plant.identity), data = competitive.summary)
summary(NR.vs.NE)
print(Anova(NR.vs.NE), digits = 6)
qqnorm(resid(NR.vs.NE))
qqline(resid(NR.vs.NE))
shapiro.test(resid(NR.vs.NE))

#### NE vs NR plot ####
NE.vs.NR.by.spp.plot <- ggplot(data = competitive.summary, aes(x = mean.NE, y = mean.NR)) + 
  geom_point(aes(colour = Plant.identity, shape = Nutrients), size = 2) + 
  geom_smooth(method = lm, se = FALSE, fullrange = FALSE, aes(linetype = Nutrients)) +
  xlab("Neighbour effect") +
  ylab("Neighbour response") +
  labs(colour = "Species identity") +
  annotate("text", x = 3.5, y = -4.5, 
           label = "bold('p(NE x Nutrients) < 0.001')", parse = TRUE, hjust = 1, vjust = 0, size = 5) +
  theme_classic(base_size = 18)
NE.vs.NR.by.spp.plot

##export
ggsave(filename = "NE vs NR by spp plot.png", 
       NE.vs.NR.by.spp.plot,
       path = "figures/",
       width = 8,
       height = 6,
       units = "in"
)

#### mean NE vs. facilitation counts ####
## calculate spp mean NE values
mean.NE <- competitive.df %>%
  group_by(Plant.identity, Nutrients) %>%
  summarise(spp.mean.NE = mean(CE, na.rm = TRUE))
str(mean.NE)

##join with interaction counts
mean.NE.plus.counts <- mean.NE %>%
  full_join(spp.interaction.counts, by = c("Plant.identity", "Nutrients")) %>%
  filter(Interaction.type == "Facilitation") ##include only facilitation counts

## split into NE/NR dfs for analyses
mean.NE.plus.NE.counts <- mean.NE.plus.counts %>% filter(metric == 'Neighbour effect')
mean.NE.plus.NR.counts <- mean.NE.plus.counts %>% filter(metric == 'Neighbour response')

## models
mean.NE.NR.mod <- lmer(Interaction.count ~ spp.mean.NE*Nutrients + (1|Plant.identity), data = mean.NE.plus.NR.counts)
qqnorm(resid(mean.NE.NR.mod))
qqline(resid(mean.NE.NR.mod))
shapiro.test(resid(mean.NE.NR.mod))
summary(mean.NE.NR.mod)
print(Anova(mean.NE.NR.mod), digits = 6)

mean.NE.NE.mod <- lmer(Interaction.count ~ spp.mean.NE*Nutrients + (1|Plant.identity), data = mean.NE.plus.NE.counts)
qqnorm(resid(mean.NE.NE.mod))
qqline(resid(mean.NE.NE.mod))
shapiro.test(resid(mean.NE.NE.mod))
summary(mean.NE.NE.mod)
print(Anova(mean.NE.NE.mod), digits = 6)

## plot
mean.NE.plots <- ggplot(mean.NE.plus.counts, aes(x = spp.mean.NE, y = Interaction.count)) +
  geom_point(aes(colour = Plant.identity, shape = Nutrients)) +
  geom_smooth(data = subset(mean.NE.plus.counts, metric == "Neighbour effect"), 
              method = "lm", se = FALSE, fullrange = FALSE, colour = "black",
              aes(linetype = Nutrients)) +
  facet_wrap(~metric) +
  ylim(0,24) +
  ylab("Count of facilitative interactions") +
  xlab("Species mean neighbour effect") +
  labs(colour = "Species identity") +
  theme_classic(base_size = 20)
mean.NE.plots

##export
ggsave(filename = "mean NE vs facilitation.png", 
       mean.NE.plots,
       path = "figures/",
       width = 8,
       height = 6,
       units = "in"
)

#### make table S5 ####
##reformat spp interaction counts
spp.interaction.counts.new <- spp.interaction.counts %>%
  pivot_wider(names_from = Interaction.type, values_from = Interaction.count) #expand columns again

##calculate means
summary_stats <- competitive.df %>%
  pivot_longer(cols = NNR:CE, names_to = "index.name", values_to = "index.value") %>%
  mutate(index.name = replace(index.name, index.name == 'CE', 'Neighbour effect')) %>%
  mutate(index.name = replace(index.name, index.name == 'NNR', 'Neighbour response')) %>%
  group_by(Plant.identity, Nutrients, index.name) %>%
  summarize(mean_value = mean(index.value, na.rm = TRUE),
            se_value = sd(index.value, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

##join with summary stats and format for export
tables5 <- summary_stats %>%
  rename(metric = index.name) %>%
  full_join(spp.interaction.counts.new, by = c("Plant.identity", "Nutrients", "metric")) %>%
  mutate(mean_se = paste0(round(mean_value, 3), " ± ", round(se_value, 3))) %>% #join mean and se columns and round to three decimal places
  dplyr::select(Plant.identity, Nutrients, metric, mean_se, Competition, Facilitation) %>%
  rename('Species code' = Plant.identity,
         'Interaction type' = metric,
         'Mean interaction value ± SE' = mean_se,
         'Count of competitive interactions' = Competition,
         'Count of facilitative interactions' = Facilitation) #rename columns

## export
write_csv(tables5, "output/tableS5.csv")
