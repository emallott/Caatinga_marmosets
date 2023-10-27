#Set up environment----
setwd("/Users/mallott/Dropbox/Projects/Gut_microbiome/Caatinga_marmosets/diet_data")

library(tidyverse)
library(ggpubr)
library(vegan)
library(car)
library(multcomp)

#Import invert data ----
inverts_raw = read_tsv("invert_results_clean_marmoset.txt")
invert = inverts_raw %>% 
  group_by(TaxID) %>% 
  summarize(across(everything(), sum)) %>% 
  gather(SampleID, count, 2:68) %>% 
  spread(TaxID, count)
metadata = read_csv("caatinga_metadata_diet.csv")

invert_meta = invert %>% left_join(metadata, by = "SampleID") 
invert_filtered = invert_meta %>% 
  mutate(total = rowSums(invert_meta[,2:82])) %>% 
  filter(total > 0) %>% 
  select_if(~ !is.numeric(.) || sum(.) != 0) 
write_csv(invert_filtered, "invert_marmoset_filtered.csv")
invert_filtered = read_csv("invert_marmoset_filtered.csv")
invert_filtered_family = read_csv("invert_marmoset_filtered_family.csv") #after re-annotating and collapsing at family level

inverts_raw_div = read_tsv("invert_results_clean_marmoset_div.txt")
invert_meta_div = inverts_raw_div %>% 
  gather(SampleID, count, 2:68) %>% 
  spread(ID, count) %>% 
  left_join(metadata, by = "SampleID") 
invert_filtered_div = invert_meta_div %>% 
  mutate(total = rowSums(invert_meta_div[,2:1369])) %>% 
  filter(total > 0) %>% 
  select_if(~ !is.numeric(.) || sum(.) != 0) 


#Import plant data ----

plants_raw = read_tsv("plant_results_clean_marmoset.txt")
plant = plants_raw %>% 
  group_by(TaxID) %>% 
  summarize(across(everything(), sum)) %>% 
  gather(SampleID, count, 2:68) %>% 
  spread(TaxID, count)
metadata = read_csv("caatinga_metadata_diet.csv")

plant_meta = plant %>% left_join(metadata, by = "SampleID") 
plant_filtered = plant_meta %>% 
  mutate(total = rowSums(plant_meta[,2:245])) %>% 
  filter(total > 0) %>% 
  select_if(~ !is.numeric(.) || sum(.) != 0) 
write_csv(plant_filtered, "plant_marmoset_filtered.csv")
plant_filtered = read_csv("plant_marmoset_filtered.csv")
plant_filtered_family = read_csv("plant_marmoset_filtered_family.csv") #after re-annotating and collapsing at family level

plants_raw_div = read_tsv("plant_results_clean_marmoset_div.txt")
plant_meta_div = plants_raw_div %>% 
  gather(SampleID, count, 2:68) %>% 
  spread(ID, count) %>% 
  left_join(metadata, by = "SampleID") 
plant_filtered_div = plant_meta_div %>% 
  mutate(total = rowSums(plant_meta_div[,2:5322])) %>% 
  filter(total > 0) %>% 
  select_if(~ !is.numeric(.) || sum(.) != 0)

#Import vert data ----

verts_raw = read_tsv("newvert_results_clean_marmoset.txt")
vert = verts_raw %>% 
  group_by(TaxID) %>% 
  summarize(across(everything(), sum)) %>% 
  gather(SampleID, count, 2:67) %>% 
  spread(TaxID, count)
metadata = read_csv("caatinga_metadata_diet.csv")

vert_meta = vert %>% left_join(metadata, by = "SampleID") 
vert_filtered = vert_meta %>% 
  mutate(total = rowSums(vert_meta[,2:54])) %>% 
  filter(total > 0) %>% 
  select_if(~ !is.numeric(.) || sum(.) != 0) 
write_csv(vert_filtered, "newvert_marmoset_filtered.csv")
vert_filtered = read_csv("newvert_marmoset_filtered.csv")
vert_filtered_family = read_csv("newvert_marmoset_filtered_family_noprimate.csv") #after re-annotating, collapsing at family level, and removing reads assigned to primates or higher levels (eg, unassigned)

verts_raw_div = read_tsv("newvert_results_clean_marmoset_div.txt")
verts_raw_div_noprim = verts_raw_div %>% 
  filter(TaxID != "Alouatta" & TaxID != "Amniota" & TaxID != "Amniota" &
           TaxID != "Aotus" & TaxID != "Aotus azarai" & TaxID != "Ateles" &
           TaxID != "Ateles belzebuth" & TaxID != "Atelinae" & TaxID != "Boreoeutheria" &
           TaxID != "Bracyteles arachnoides" & TaxID != "Callicebinae" & TaxID != "Callithrix" &
           TaxID != "Callithrix jacchus" & TaxID != "Callithrix pygmaea" & TaxID != "Catarrhini" &
           TaxID != "Cebus" & TaxID != "Cercopithecidae" & TaxID != "Chiropotes" &
           TaxID != "Euarchontoglires" & TaxID != "Euteleostomi" & TaxID != "Eutheria" &
           TaxID != "Hominidae" & TaxID != "Hominoidea" & TaxID != "Lagothrix lagotricha" &
           TaxID != "Leontopithecus rosalia" & TaxID != "Metatheria" & TaxID != "Pan troglodytes" &
           TaxID != "Pithecia pithecia" & TaxID != "Pitheciinae" & TaxID != "Platyrrhini" &
           TaxID != "Plecturocebus" & TaxID != "Plecturocebus cupreus" & TaxID != "Primates" &
           TaxID != "Saguinus midas" & TaxID != "Saimiri" & TaxID != "Simiiformes" & 
           TaxID != "Theria" & TaxID != "<NA>")

vert_meta_div = verts_raw_div_noprim %>% 
  dplyr::select(!TaxID) %>% 
  gather(SampleID, count, 2:67) %>% 
  spread(ID, count) %>% 
  left_join(metadata, by = "SampleID") 
vert_filtered_div = vert_meta_div %>% 
  mutate(total = rowSums(vert_meta_div[,2:1657])) %>% 
  filter(total > 0) %>% 
  select_if(~ !is.numeric(.) || sum(.) != 0)
         

#Calculate frequencies per season and per group ----

invert_family_season = invert_filtered_family %>% 
  group_by(Season) %>% 
  summarize(across(Anostostomatidae:`Unclassified Lepidoptera`,
                   ~sum(. > 0)/n(), .names = "{col}.freq"))
write_csv(invert_family_season, "invert_season_family_freq.csv")

invert_family_group = invert_filtered_family %>% 
  group_by(Group) %>% 
  summarize(across(Anostostomatidae:`Unclassified Lepidoptera`,
                   ~sum(. > 0)/n(), .names = "{col}.freq"))
write_csv(invert_family_group, "invert_season_group_freq.csv")

invert_family = invert_filtered_family %>% 
  summarize(across(Anostostomatidae:`Unclassified Lepidoptera`,
                   ~sum(. > 0)/n(), .names = "{col}.freq"))
write_csv(invert_family, "invert_family_freq_total.csv")

plant_family_season = plant_filtered_family %>% 
  group_by(Season) %>% 
  summarize(across(Anacardiaceae:Vochysiaceae,
                   ~sum(. > 0)/n(), .names = "{col}.freq"))
write_csv(plant_family_season, "plant_season_family_freq.csv")

plant_family_group = plant_filtered_family %>% 
  group_by(Group) %>% 
  summarize(across(Anacardiaceae:Vochysiaceae,
                   ~sum(. > 0)/n(), .names = "{col}.freq"))
write_csv(plant_family_group, "plant_group_family_freq.csv")

vert_family_season = vert_filtered_family %>% 
  group_by(Season) %>% 
  summarize(across(Columbidae:`Unassigned Aves`,
                   ~sum(. > 0)/n(), .names = "{col}.freq"))
write_csv(vert_family_season, "vert_season_family_freq.csv")

vert_family_group = vert_filtered_family %>% 
  group_by(Group) %>% 
  summarize(across(Columbidae:`Unassigned Aves`,
                   ~sum(. > 0)/n(), .names = "{col}.freq"))
write_csv(vert_family_group, "vert_group_family_freq.csv")

#Calculate and plot species richness -----

invert_family_richness = apply(invert_filtered_family[,2:26] > 0, 1, sum)
invert_family_richness = as.data.frame(invert_family_richness)
invert_family_richness = invert_family_richness %>% 
  cbind(invert_filtered_family[,27:45])

invert_family_richness$Group = as.factor(invert_family_richness$Group)
invert_family_richness$Season = as.factor(invert_family_richness$Season)

invert_rich = lm(invert_family_richness ~ Season + Group + Preservative, data = invert_family_richness)
Anova(invert_rich)
summary(glht(invert_rich,linfct=mcp(Season="Tukey")))

richness_season_summary = invert_family_richness %>% 
  group_by(Season) %>% 
  summarize(average = mean(invert_family_richness))

inv_rich_season = ggboxplot(invert_family_richness, x = "Season",
                     y = "invert_family_richness", color = "Season",
                     add = "jitter", palette = c("gray", "black"),
                     ylab = "Family richness",
                     title = " ",
                     add.params = list(fill = "white"))
inv_rich_season = ggpar(inv_rich_season, legend = "right", font.y = 16,
                        font.legend = 16, font.ytickslab = 14,
                        font.main = 17, ylim = c(0,24)) + rremove("xlab") + 
  rremove("x.text") + rremove("legend.title") + 
  stat_compare_means(label = "p.signif", 
                     comparisons = list(c("Wet", "Dry")),
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.05, 0.5, 1),
                                        symbols = c("*", "*", "*", "*", "*")))

inv_rich_group = ggboxplot(invert_family_richness, x = "Group",
                            y = "invert_family_richness", color = "Group",
                            add = "jitter", palette = "Set1",
                            ylab = "Family richness",
                           title = "Arthopods",
                            add.params = list(fill = "white"))
inv_rich_group = ggpar(inv_rich_group, legend = "right", font.y = 16,
                       font.legend = 16, font.ytickslab = 14,
                       font.main = 17, ylim = c(0,24)) + rremove("xlab") + 
  rremove("x.text") + rremove("legend.title")


tiff(file="inv_rich.tif", res=300, width=10, height=3, units="in")
ggarrange(inv_rich_season, inv_rich_group, ncol = 2, nrow=1, 
          common.legend = F, legend = "right", align = "h", 
          labels = c("A", "B"), widths = c(2,4))

dev.off()

plant_family_richness = apply(plant_filtered_family[,2:46] > 0, 1, sum)
plant_family_richness = as.data.frame(plant_family_richness)
plant_family_richness = plant_family_richness %>% 
  cbind(plant_filtered_family[,47:65])

plant_family_richness$Group = as.factor(plant_family_richness$Group)
plant_family_richness$Season = as.factor(plant_family_richness$Season)

plant_rich = lm(plant_family_richness ~ Season + Group + Preservative, data = plant_family_richness)
Anova(plant_rich)
summary(glht(plant_rich,linfct=mcp(Group="Tukey")))


pla_rich_season = ggboxplot(plant_family_richness, x = "Season",
                            y = "plant_family_richness", color = "Season",
                            add = "jitter", palette = c("gray", "black"),
                            ylab = "Family richness",
                            title = " ",
                            add.params = list(fill = "white"))
pla_rich_season = ggpar(pla_rich_season, legend = "right", font.y = 16,
                        font.legend = 16, font.ytickslab = 14,
                        font.main = 17, ylim = c(0,24)) + rremove("xlab") + 
  rremove("x.text") + rremove("legend.title")

pla_rich_group = ggboxplot(plant_family_richness, x = "Group",
                           y = "plant_family_richness", color = "Group",
                           add = "jitter", palette = "Set1",
                           ylab = "Family richness",
                           title = "Plants",
                           add.params = list(fill = "white"))
pla_rich_group = ggpar(pla_rich_group, legend = "right", font.y = 16,
                       font.legend = 16, font.ytickslab = 14,
                       font.main = 17, ylim = c(0,24)) + rremove("xlab") + 
  rremove("x.text") + rremove("legend.title") + 
  stat_compare_means(comparisons = c("Princess","House"), label = "p.signif")


tiff(file="pla_rich.tif", res=300, width=10, height=3, units="in")
ggarrange(pla_rich_season, pla_rich_group, ncol = 2, nrow=1, 
          common.legend = F, legend = "right", align = "h", 
          labels = c("A", "B"), widths = c(2,4))

dev.off()

vert_family_richness = apply(vert_filtered_family[,2:5] > 0, 1, sum)
vert_family_richness = as.data.frame(vert_family_richness)
vert_family_richness = vert_family_richness %>% 
  cbind(vert_filtered_family[,6:24])

vert_family_richness$Group = as.factor(vert_family_richness$Group)
vert_family_richness$Season = as.factor(vert_family_richness$Season)

vert_rich = lm(vert_family_richness ~ Season + Group + Preservative, data = vert_family_richness)
Anova(vert_rich)
summary(glht(vert_rich,linfct=mcp(Group="Tukey")))

vert_rich_season = ggboxplot(vert_family_richness, x = "Season",
                            y = "vert_family_richness", color = "Season",
                            add = "jitter", palette = c("gray", "black"),
                            ylab = "Family richness",
                            title = " ",
                            add.params = list(fill = "white"))
vert_rich_season = ggpar(vert_rich_season, legend = "right", font.y = 16,
                         font.legend = 16, font.ytickslab = 14,
                         font.main = 17, ylim = c(0,24)) + rremove("xlab") + 
  rremove("x.text") + rremove("legend.title")

vert_rich_group = ggboxplot(vert_family_richness, x = "Group",
                           y = "vert_family_richness", color = "Group",
                           add = "jitter", palette = "Set1",
                           ylab = "Family richness",
                           title = "Vertebrates",
                           add.params = list(fill = "white"))
vert_rich_group = ggpar(vert_rich_group, legend = "right", font.y = 16,
                        font.legend = 16, font.ytickslab = 14,
                        font.main = 17, ylim = c(0,24)) + rremove("xlab") + 
  rremove("x.text") + rremove("legend.title")


tiff(file="vert_rich.tif", res=300, width=10, height=3, units="in")
ggarrange(vert_rich_season, vert_rich_group, ncol = 2, nrow=1, 
          common.legend = F, legend = "right", align = "h", 
          labels = c("A", "B"), widths = c(2,4))

dev.off()

#Calculate shannon diversity ----

plant_shannon = as.data.frame(diversity(plant_filtered_div[,2:1876], index = "shannon"))
invert_shannon = as.data.frame(diversity(invert_filtered_div[,2:465], index = "shannon"))
vert_shannon = as.data.frame(diversity(vert_filtered_div[,2:787], index = "shannon"))

plant_shannon = plant_shannon %>% 
  rename("shannon_plant" = "diversity(plant_filtered_div[, 2:1876], index = \"shannon\")")
invert_shannon = invert_shannon %>% 
  rename("shannon_invert" = "diversity(invert_filtered_div[, 2:465], index = \"shannon\")")
vert_shannon = vert_shannon %>% 
  rename("shannon_vert" = "diversity(vert_filtered_div[, 2:787], index = \"shannon\")")

plant_div = cbind(plant_filtered_div[,1877:1895], plant_shannon)
invert_div = cbind(invert_filtered_div[,466:484], invert_shannon)
vert_div = cbind(vert_filtered_div[,788:806], vert_shannon)

alpha_plant = lm(shannon_plant ~ Season + Group + Preservative, data = plant_div)
Anova(alpha_plant)

invert_div$Group = as.factor(invert_div$Group)
invert_div$Season = as.factor(invert_div$Season)
invert_div$Preservative = as.factor(invert_div$Preservative)

alpha_invert = lm(shannon_invert ~ Season + Group + Preservative, data = invert_div)
Anova(alpha_invert)
summary(glht(alpha_invert,linfct=mcp(Season="Tukey")))
summary(glht(alpha_invert,linfct=mcp(Preservative="Tukey")))

shannon_season_summary = invert_div %>% 
  group_by(Season) %>% 
  summarize(average = mean(shannon_invert))

shannon_preservative_summary = invert_div %>% 
  group_by(Preservative) %>% 
  summarize(average = mean(shannon_invert))

alpha_vert = lm(shannon_vert ~ Season + Group + Preservative, data = vert_div)
Anova(alpha_vert)

plant_shan_season = ggboxplot(plant_div, x = "Season",
                             y = "shannon_plant", color = "Season",
                             add = "jitter", palette = c("gray", "black"),
                             ylab = "Shannon diversity",
                             add.params = list(fill = "white"))
plant_shan_season = ggpar(plant_shan_season, legend = "right", font.y = 16,
                          font.legend = 16, font.ytickslab = 14,
                          font.main = 17, ylim = c(0,4)) + rremove("xlab") + 
  rremove("x.text") + rremove("legend.title")

invert_shan_season = ggboxplot(invert_div, x = "Season",
                              y = "shannon_invert", color = "Season",
                              add = "jitter", palette = c("gray", "black"),
                              ylab = "Shannon diversity",
                              add.params = list(fill = "white"))
invert_shan_season = ggpar(invert_shan_season, legend = "right", font.y = 16,
                           font.legend = 16, font.ytickslab = 14,
                           font.main = 17, ylim = c(0,4)) + rremove("xlab") + 
  rremove("x.text") + rremove("legend.title") + 
  stat_compare_means(label = "p.signif", 
                     comparisons = list(c("Wet", "Dry")),
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.05, 0.5, 1),
                                        symbols = c("*", "*", "*", "*", "*")))


vert_shan_season = ggboxplot(vert_div, x = "Season",
                               y = "shannon_vert", color = "Season",
                               add = "jitter", palette = c("gray", "black"),
                               ylab = "Shannon diversity",
                               add.params = list(fill = "white"))
vert_shan_season = ggpar(vert_shan_season, legend = "right", font.y = 16,
                         font.legend = 16, font.ytickslab = 14,
                         font.main = 17, ylim = c(0,4)) + rremove("xlab") + 
  rremove("x.text") + rremove("legend.title")

tiff(file="shannon_div_season.tif", res=300, width=10, height=3, units="in")
ggarrange(plant_shan_season, invert_shan_season, vert_shan_season, ncol = 3, nrow=1, 
          common.legend = T, legend = "right", align = "h", 
          labels = c("A", "B", "C"))

dev.off()

plant_shan_group = ggboxplot(plant_div, x = "Group",
                              y = "shannon_plant", color = "Group",
                              add = "jitter", palette = "Set1",
                              ylab = "Shannon diversity",
                              add.params = list(fill = "white"))
plant_shan_group = ggpar(plant_shan_group, legend = "right", font.y = 16,
                         font.legend = 16, font.ytickslab = 14,
                         font.main = 17, ylim = c(0,4)) + rremove("xlab") + 
  rremove("x.text") + rremove("legend.title")

invert_shan_group = ggboxplot(invert_div, x = "Group",
                               y = "shannon_invert", color = "Group",
                               add = "jitter", palette = "Set1",
                               ylab = "Shannon diversity",
                               add.params = list(fill = "white"))
invert_shan_group = ggpar(invert_shan_group, legend = "right", font.y = 16,
                          font.legend = 16, font.ytickslab = 14,
                          font.main = 17, ylim = c(0,4)) + rremove("xlab") + 
  rremove("x.text") + rremove("legend.title") 

vert_shan_group = ggboxplot(vert_div, x = "Group",
                             y = "shannon_vert", color = "Group",
                             add = "jitter", palette = "Set1",
                             ylab = "Shannon diversity",
                             add.params = list(fill = "white"))
vert_shan_group = ggpar(vert_shan_group, legend = "right", font.y = 16,
                        font.legend = 16, font.ytickslab = 14,
                        font.main = 17, ylim = c(0,4)) + rremove("xlab") + 
  rremove("x.text") + rremove("legend.title")

tiff(file="shannon_div_group.tif", res=300, width=10, height=9, units="in")
ggarrange(plant_shan_group, invert_shan_group, vert_shan_group, ncol = 1, nrow=3, 
          common.legend = T, legend = "right", align = "h", 
          labels = c("A", "B", "C"))

dev.off()

#Diversity graphs ----

tiff(file = "diversity_combo.tif", res = 300, width = 20, height = 12,
     units = "in")
plot1 = ggarrange(inv_rich_season, invert_shan_season, vert_rich_season, 
                  vert_shan_season,pla_rich_season, plant_shan_season,  
                  ncol = 2, nrow = 3, common.legend = T, legend = "right",
                  align = "hv", labels = c("B", " ", "D", " ",
                                           "F", " "))
plot2 = ggarrange(inv_rich_group, invert_shan_group, vert_rich_group, 
                  vert_shan_group, pla_rich_group, plant_shan_group, 
                  ncol = 2, nrow = 3, common.legend = T, legend = "right",
                  align = "hv", labels = c("A", " ", "C", " ", 
                                           "E", " "))
div_plot = ggarrange(plot2, plot1, ncol = 2, nrow = 1, widths = c(3,2),
          align = "hv")
div_plot
dev.off()

setEPS()
postscript(file="diversity_combo.eps", width=20, height=12, paper = "special")
div_plot
dev.off()

#Proportion of reads to each dietary category ----

prop = read_csv("diet_proportion.csv")
prop$Group = as.factor(prop$Group)

prop_plant = lm(Plant_relab ~ Season + Group + Preservative, data = prop)
Anova(prop_plant)
summary(glht(prop_plant,linfct=mcp(Group="Tukey")))

plant_prop_season = ggboxplot(prop, x = "Season",
                              y = "Plant_relab", color = "Season",
                              add = "jitter", palette = c("gray", "black"),
                              ylab = "Relative abundance",
                              title = "Plants",
                              add.params = list(fill = "white"))
plant_prop_season = ggpar(plant_prop_season, legend = "right") + rremove("xlab") + 
  rremove("x.text") + rremove("legend.title")

prop_invert = lm(Invert_relab ~ Season + Group + Preservative, data = prop)
Anova(prop_invert)
summary(glht(prop_invert,linfct=mcp(Group="Tukey")))

invert_prop_season = ggboxplot(prop, x = "Season",
                              y = "Invert_relab", color = "Season",
                              add = "jitter", palette = c("gray", "black"),
                              ylab = "Relative abundance",
                              title = "Invertebrates",
                              add.params = list(fill = "white"))
invert_prop_season = ggpar(invert_prop_season, legend = "right") + rremove("xlab") + 
  rremove("x.text") + rremove("legend.title")

prop_vert = lm(Vert_relab ~ Season + Group + Preservative, data = prop)
Anova(prop_vert)
summary(glht(prop_vert,linfct=mcp(Group="Tukey")))

vert_prop_season = ggboxplot(prop, x = "Season",
                               y = "Vert_relab", color = "Season",
                               add = "jitter", palette = c("gray", "black"),
                               ylab = "Relative abundance",
                             title = "Vertebrates",
                               add.params = list(fill = "white"))
vert_prop_season = ggpar(vert_prop_season, legend = "right") + rremove("xlab") + 
  rremove("x.text") + rremove("legend.title")

plant_prop_group = ggboxplot(prop, x = "Group",
                              y = "Plant_relab", color = "Group",
                              add = "jitter", palette = "Set1",
                              ylab = "Relative abundance",
                              title = "Plants",
                              add.params = list(fill = "white"))
plant_prop_group = ggpar(plant_prop_group, legend = "right") + rremove("xlab") + 
  rremove("x.text") + rremove("legend.title")

invert_prop_group = ggboxplot(prop, x = "Group",
                               y = "Invert_relab", color = "Group",
                               add = "jitter", palette = "Set1",
                               ylab = "Relative abundance",
                              title = "Invertebrates",
                               add.params = list(fill = "white"))
invert_prop_group = ggpar(invert_prop_group, legend = "right") + rremove("xlab") + 
  rremove("x.text") + rremove("legend.title")

vert_prop_group = ggboxplot(prop, x = "Group",
                             y = "Vert_relab", color = "Group",
                             add = "jitter", palette = "Set1",
                             ylab = "Relative abundance",
                            title = "Vertebrates",
                             add.params = list(fill = "white"))
vert_prop_group = ggpar(vert_prop_group, legend = "right") + rremove("xlab") + 
  rremove("x.text") + rremove("legend.title")


tiff(file="diet_prop_season.tif", res=300, width=10, height=3, units="in")
ggarrange(plant_prop_season, invert_prop_season, vert_prop_season, ncol = 3, nrow=1, 
          common.legend = T, legend = "right", align = "h", 
          labels = c("A", "B", "C"))

dev.off()

tiff(file="diet_prop_group.tif", res=300, width=10, height=9, units="in")
ggarrange(plant_prop_group, invert_prop_group, vert_prop_group, ncol = 1, nrow=3, 
          common.legend = T, legend = "right", align = "v", 
          labels = c("A", "B", "C"))

dev.off()

#Diet Proportion Plots----

##Period----

season_all = read_csv("sankey_plots_season.csv")
season_all$Family <- factor(season_all$Family, 
                            levels = unique(season_all$Family))

prop_season = ggplot(season_all, aes(x=Period, y=Family)) +
  geom_point(aes(size = Proportion, color = `Diet component`)) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black", angle = 90),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_text(face = "bold",size = 14, color = "black"),
    legend.text = element_text(face = "plain",size = 12, color = "black")) +
  labs(color="Diet component") +
  scale_color_manual(values=c("#2171b5","#807dba","#67000d"))+
  scale_size_continuous(range = c(-0.4,3))

##Group----

group_all = read_csv("sankey_plots_group.csv")
group_all$Family <- factor(group_all$Family, 
                            levels = unique(group_all$Family))

prop_group = ggplot(group_all, aes(x=Group, y=Family)) +
  geom_point(aes(size = Proportion, color = `Diet component`)) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 10, color = "black", angle = 90),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(face = "bold",size = 14, color = "black"),
        legend.text = element_text(face = "plain",size = 12, color = "black")) +
  labs(color="Diet component") +
  scale_color_manual(values=c("#2171b5","#807dba","#67000d"))+
  scale_size_continuous(range = c(-0.4,3))

##Combo----
prop_plot = ggarrange(prop_season, prop_group,  
                  ncol = 2, common.legend = T, legend = "right",
                  align = "hv", widths = c(1.5, 2))
prop_plot

setEPS()
postscript(file="proportion_combo.eps", width=8, height=10, paper = "special")
prop_plot
dev.off()

