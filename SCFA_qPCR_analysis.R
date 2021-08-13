#Set up environment----

setwd("/Users/elizabethmallott/Dropbox/Projects/Gut_microbiome/Caatinga_marmosets")

library(tidyverse)
library(nlme)
library(multcomp)
library(car)

#Import data----

butyrate = read_csv("SCFA_caatinga_marmosets_results.csv")

#Linear models

butyrate$Group = as.factor(butyrate$Group)
butyrate$Season = as.factor(butyrate$Season)
butyrate$Sex = as.factor(butyrate$Sex)
butyrate$Preservative = as.factor(butyrate$Preservative)

bcoa = lme(fixed = BCoA ~ Season + Sex + Preservative, data = butyrate, 
           random = ~1|Group)
summary(bcoa)
Anova(bcoa)

bcoa_full = lm(BCoA ~ Season + Group + Sex + Preservative, 
            data = butyrate)
summary(bcoa_full)
Anova(bcoa_full)
summary(glht(bcoa_full,linfct=mcp(Group="Tukey")))

