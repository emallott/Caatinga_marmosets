#Set up environment----

setwd("/Users/mallott/Dropbox/Projects/Gut_microbiome/Caatinga_marmosets/16S_results")

#Import data----
unweighted = as.dist(read.table("unweighted-distance-matrix.tsv", header = T))
weighted = as.dist(read.table("weighted-distance-matrix.tsv", header = T))
metadata = read.csv("caatinga_metadata_r_8956.csv", header=T)

unweighted_nounk = as.dist(read.table("unweighted-distance-matrix-nounk.tsv", header = T))
weighted_nounk = as.dist(read.table("weighted-distance-matrix-nounk.tsv", header = T))
metadata_nounk = read.csv("caatinga_metadata_r_8956_nounk.csv", header=T)

unweighted_monthcomp = as.dist(read.table("unweighted-distance-matrix-monthcomp.tsv", header = T))
weighted_monthcomp = as.dist(read.table("weighted-distance-matrix-monthcomp.tsv", header = T))
metadata_monthcomp = read.csv("caatinga_metadata_r_8956-monthcomp.csv", header=T)

unweighted_dryonly = as.dist(read.table("unweighted-distance-matrix-dryonly.tsv", header = T))
weighted_dryonly = as.dist(read.table("weighted-distance-matrix-dryonly.tsv", header = T))
metadata_dryonly = read.csv("caatinga_metadata_r_8956-dryonly.csv", header=T)

#Permanovas----
library(vegan)

set.seed(1018)

adonis2(unweighted~Season+Group+Age+Sex+Preservative, data=metadata, 
        by = "margin", permutations = 5000)
adonis2(unweighted~Season+Human_food+Age+Sex+Preservative, data=metadata, 
        by = "margin", permutations = 5000)
adonis2(unweighted~Season+Domestic_animal+Age+Sex+Preservative, data=metadata, 
        by = "margin", permutations = 5000)
adonis2(unweighted~Year+Season+Group+Age+Sex+Preservative, data=metadata, 
        by = "margin", permutations = 5000)
anova(betadisper(unweighted, group = metadata$Year))
anova(betadisper(unweighted, group = metadata$Season))
anova(betadisper(unweighted, group = metadata$Group))
anova(betadisper(unweighted, group = metadata$Human_food))
anova(betadisper(unweighted, group = metadata$Domestic_animal))
anova(betadisper(unweighted, group = metadata$Age))
anova(betadisper(unweighted, group = metadata$Sex))
anova(betadisper(unweighted, group = metadata$Preservative))

adonis2(weighted~Season+Group+Age+Sex+Preservative, data=metadata, 
        by = "margin", permutations = 5000)
adonis2(weighted~Season+Human_food+Age+Sex+Preservative, data=metadata, 
        by = "margin", permutations = 5000)
adonis2(weighted~Season+Domestic_animal+Age+Sex+Preservative, data=metadata, 
        by = "margin", permutations = 5000)
adonis2(weighted~Year+Season+Group+Age+Sex+Preservative, data=metadata, 
        by = "margin", permutations = 5000)
anova(betadisper(weighted, group = metadata$Year))
anova(betadisper(weighted, group = metadata$Season))
anova(betadisper(weighted, group = metadata$Group))
anova(betadisper(weighted, group = metadata$Human_food))
anova(betadisper(weighted, group = metadata$Domestic_animal))
anova(betadisper(weighted, group = metadata$Age))
anova(betadisper(weighted, group = metadata$Sex))
anova(betadisper(weighted, group = metadata$Preservative))

adonis2(unweighted_nounk~Season+Group+Age+Sex+Preservative, data=metadata_nounk, 
        by = "margin", permutations = 5000)
adonis2(unweighted_nounk~Season+Human_food+Age+Sex+Preservative, data=metadata_nounk, 
        by = "margin", permutations = 5000)
adonis2(unweighted_nounk~Season+Domestic_animal+Age+Sex+Preservative, data=metadata_nounk, 
        by = "margin", permutations = 5000)
adonis2(unweighted_nounk~Year+Season+Group+Age+Sex+Preservative, data=metadata_nounk, 
        by = "margin", permutations = 5000)
anova(betadisper(unweighted_nounk, group = metadata_nounk$Season))
anova(betadisper(unweighted_nounk, group = metadata_nounk$Group))
anova(betadisper(unweighted_nounk, group = metadata_nounk$Human_food))
anova(betadisper(unweighted_nounk, group = metadata_nounk$Domestic_animal))
anova(betadisper(unweighted_nounk, group = metadata_nounk$Age))
anova(betadisper(unweighted_nounk, group = metadata_nounk$Sex))
anova(betadisper(unweighted_nounk, group = metadata_nounk$Preservative))
anova(betadisper(unweighted_nounk, group = metadata_nounk$Year))

adonis2(weighted_nounk~Season+Group+Age+Sex+Preservative, data=metadata_nounk, 
        by = "margin", permutations = 5000)
adonis2(weighted_nounk~Season+Human_food+Age+Sex+Preservative, data=metadata_nounk, 
        by = "margin", permutations = 5000)
adonis2(weighted_nounk~Season+Domestic_animal+Age+Sex+Preservative, data=metadata_nounk, 
        by = "margin", permutations = 5000)
adonis2(weighted_nounk~Year+Season+Group+Age+Sex+Preservative, data=metadata_nounk, 
        by = "margin", permutations = 5000)
anova(betadisper(weighted_nounk, group = metadata_nounk$Season))
anova(betadisper(weighted_nounk, group = metadata_nounk$Group))
anova(betadisper(weighted_nounk, group = metadata_nounk$Human_food))
anova(betadisper(weighted_nounk, group = metadata_nounk$Domestic_animal))
anova(betadisper(weighted_nounk, group = metadata_nounk$Age))
anova(betadisper(weighted_nounk, group = metadata_nounk$Sex))
anova(betadisper(weighted_nounk, group = metadata_nounk$Preservative))
anova(betadisper(weighted_nounk, group = metadata_nounk$Year))

adonis2(unweighted_monthcomp ~ Year + Season + Group + Preservative, 
        data=metadata_monthcomp, 
        by = "margin", permutations = 4999)
pairwise.adonis(unweighted_monthcomp, metadata_monthcomp$SeasonYear)
adonis2(unweighted_monthcomp ~ Year + Season + Group, 
        data=metadata_monthcomp, 
        by = "margin", permutations = 4999)
pairwise.adonis(weighted_monthcomp, metadata_monthcomp$SeasonYear)


adonis2(unweighted_dryonly ~ Year + Group + Preservative, 
        data=metadata_dryonly, 
        by = "margin", permutations = 4999)
pairwise.adonis(unweighted_monthcomp, metadata_monthcomp$SeasonYear)
adonis2(unweighted_dryonly ~ Year + Group + Preservative, 
        data=metadata_dryonly, 
        by = "margin", permutations = 4999)
pairwise.adonis(weighted_monthcomp, metadata_monthcomp$SeasonYear)

#Alpha diversity----
faith = read.table("faithpd.tsv", header=T)
otus = read.table("observed_features.tsv", header = T)
shannon = read.table("shannon.tsv", header = T)

library(tidyverse)
library(nlme)
library(multcomp)
library(car)

alpha = inner_join(metadata, faith, by = "SampleID") %>% inner_join(otus, by = "SampleID") %>% inner_join(shannon, by = "SampleID")
alpha$Group = as.factor(alpha$Group)
alpha$Human_food = as.factor(alpha$Human_food)
alpha$Domestic_animal = as.factor(alpha$Domestic_animal)
alpha$Season = as.factor(alpha$Season)
alpha$Age = as.factor(alpha$Age)
alpha$Sex = as.factor(alpha$Sex)
alpha$AgeSex = as.factor(alpha$AgeSex)
alpha$Preservative = as.factor(alpha$Preservative)
alpha$Year = as.factor(alpha$Year)

alpha_nounk = alpha %>% 
  filter(Sex != "Unknown" & Age != "Unknown")

f = lme(fixed=faith_pd~Season + Age + Sex + Preservative, data=alpha, random = ~1|Group)
summary(f)
Anova(f)
summary(glht(f,linfct=mcp(Age="Tukey")))

f_full = lm(faith_pd ~ Season+Group+Age+Sex+Preservative, 
            data = alpha)
summary(f_full)
Anova(f_full)
summary(glht(f_full,linfct=mcp(Season="Tukey")))
summary(glht(f_full,linfct=mcp(Preservative="Tukey")))

f_full_nounk = lm(faith_pd ~ Season+Group+Age+Sex+Preservative, 
            data = alpha_nounk)
summary(f_full_nounk)
Anova(f_full_nounk)
summary(glht(f_full_nounk,linfct=mcp(Season="Tukey")))
summary(glht(f_full_nounk,linfct=mcp(Preservative="Tukey")))

f_full_nounk_year = lm(faith_pd ~ Year+Season+Group+Age+Sex+Preservative, 
                  data = alpha_nounk)
summary(f_full_nounk_year)
Anova(f_full_nounk_year)
summary(glht(f_full_nounk_year,linfct=mcp(Season="Tukey")))

faith_preservative_summary = alpha_nounk %>% 
  group_by(Preservative) %>% 
  summarize(average = mean(faith_pd))
faith_season_summary = alpha_nounk %>% 
  group_by(Season) %>% 
  summarize(average = mean(faith_pd))

f_full_nounk_agesex = lm(faith_pd ~ Season+Group+AgeSex+Preservative, 
                  data = alpha_nounk)
summary(f_full_nounk_agesex)
Anova(f_full_nounk_agesex)
summary(glht(f_full_nounk_agesex,linfct=mcp(Season="Tukey")))
summary(glht(f_full_nounk_agesex,linfct=mcp(Preservative="Tukey")))

f_full_human = lm(faith_pd ~ Season+Human_food+Age+Sex+Preservative, 
            data = alpha)
summary(f_full_human)
Anova(f_full_human)

f_full_animal = lm(faith_pd ~ Season+Domestic_animal+Age+Sex+Preservative, 
                  data = alpha)
summary(f_full_animal)
Anova(f_full_animal)
summary(glht(f_full_animal,linfct=mcp(Season="Tukey")))
summary(glht(f_full_animal,linfct=mcp(Preservative="Tukey")))

o = lme(fixed=observed_features~Season + Age + Sex + Preservative, 
        data=alpha, random = ~1|Group)
summary(o)
Anova(o)

o_full = lm(observed_features ~ Season+Group+Age+Sex+Preservative, 
            data = alpha)
summary(o_full)
Anova(o_full)
summary(glht(o_full,linfct=mcp(Group="Tukey")))
summary(glht(o_full,linfct=mcp(Sex="Tukey")))

o_full_nounk = lm(observed_features ~ Season+Group+Age+Sex+Preservative, 
            data = alpha_nounk)
summary(o_full_nounk)
Anova(o_full_nounk)
summary(glht(o_full_nounk,linfct=mcp(Group="Tukey")))

o_full_nounk_year = lm(observed_features ~ Year+Season+Group+Age+Sex+Preservative, 
                  data = alpha_nounk)
summary(o_full_nounk_year)
Anova(o_full_nounk_year)
summary(glht(o_full_nounk_year,linfct=mcp(Group="Tukey")))

o_full_nounk_agesex = lm(observed_features ~ Season+Group+AgeSex+Preservative, 
                  data = alpha_nounk)
summary(o_full_nounk_agesex)
Anova(o_full_nounk_agesex)
summary(glht(o_full_nounk_agesex,linfct=mcp(Group="Tukey")))
summary(glht(o_full_nounk_agesex,linfct=mcp(Sex="Tukey")))

o_human = lm(observed_features ~ Season+Human_food+Age+Sex+Preservative, 
            data = alpha)
summary(o_human)
Anova(o_human)

o_full_animal = lm(observed_features ~ Season+Domestic_animal+Age+Sex+Preservative, 
            data = alpha)
summary(o_full_animal)
Anova(o_full_animal)

s = lme(fixed=shannon_entropy~Season + Age + Sex + Preservative, 
        data=alpha, random = ~1|Group)
summary(s)
Anova(s)

s_full = lm(shannon_entropy ~ Season+Group+Age+Sex+Preservative, 
            data = alpha)
summary(s_full)
Anova(s_full)
summary(glht(s_full,linfct=mcp(Group="Tukey")))
summary(glht(s_full,linfct=mcp(Sex="Tukey")))

s_full_nounk = lm(shannon_entropy ~ Season+Group+Age+Sex+Preservative, 
            data = alpha_nounk)
summary(s_full_nounk)
Anova(s_full_nounk)
summary(glht(s_full_nounk,linfct=mcp(Group="Tukey")))

s_full_nounk_year = lm(shannon_entropy ~ Year+Season+Group+Age+Sex+Preservative, 
                  data = alpha_nounk)
summary(s_full_nounk_year)
Anova(s_full_nounk_year)
summary(glht(s_full_nounk_year,linfct=mcp(Group="Tukey")))

s_full_human = lm(shannon_entropy ~ Season+Human_food+Age+Sex+Preservative, 
            data = alpha)
summary(s_full_human)
Anova(s_full_human)
summary(glht(s_full_human,linfct=mcp(Preservative="Tukey")))

s_full_animal = lm(shannon_entropy ~ Season+Domestic_animal+Age+Sex+Preservative, 
            data = alpha)
summary(s_full_animal)
Anova(s_full_animal)
summary(glht(s_full_animal,linfct=mcp(Group="Tukey")))
summary(glht(s_full_animal,linfct=mcp(Sex="Tukey")))

library(ggpubr)

plot1 = ggboxplot(alpha_nounk, x = "Group", 
                  y = "faith_pd", color = "Group", 
                  palette = "Set1", add = "jitter", 
                  add.params = list(fill = "white"), 
                  ylab = "Faith's Phylogenetic Diversity") 
plot1 = ggpar(plot1, legend = "right", font.y = 16,
              font.legend = 16, font.ytickslab = 14) + 
  rremove("xlab") + 
  rremove("x.text") + rremove("x.ticks") + 
  rremove("legend.title")

plot2 = ggboxplot(alpha_nounk, x = "Group", 
                  y = "observed_features", color = "Group", 
                  palette = "Set1", add = "jitter", 
                  add.params = list(fill = "white"), 
                  ylab = "Observed ASVs") 
plot2 = ggpar(plot2, legend = "right", font.y = 16,
              font.legend = 16, font.ytickslab = 14) + 
  rremove("xlab") +  
  rremove("x.ticks")  + rremove("x.text") + 
  rremove("legend.title")

plot3 = ggboxplot(alpha_nounk, x = "Group", 
                  y = "shannon_entropy", color = "Group", 
                  palette = "Set1", add = "jitter", 
                  add.params = list(fill = "white"), 
                  ylab = "Shannon Diversity") 
plot3 = ggpar(plot3, legend = "right", font.y = 16,
              font.legend = 16, font.ytickslab = 14) + 
  rremove("xlab") +  
  rremove("x.ticks")  + rremove("x.text") + 
  rremove("legend.title") + 
  stat_compare_means(label = "p.signif", 
                     comparisons = list(c("Princess", "House")),
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.05, 0.5, 1),
                                        symbols = c("*", "*", "*", "*", "ns")))

plot4 = ggboxplot(alpha_nounk, x = "Season", 
                  y = "faith_pd", color = "Season", 
                  palette = c("black","darkgray"), add = "jitter", 
                  add.params = list(fill = "white"), 
                  ylab = "Faith's Phylogenetic Diversity") 
plot4 = ggpar(plot4, legend = "right", font.y = 16,
              font.legend = 16, font.ytickslab = 14) + 
  rremove("xlab") +  
  rremove("x.text") + rremove("x.ticks") + 
  rremove("legend.title") + 
  stat_compare_means(label = "p.signif", 
                     comparisons = list(c("Dry", "Wet")),
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.05, 0.8, 1),
                                        symbols = c("*", "*", "*", "*", "ns")))

plot5 = ggboxplot(alpha_nounk, x = "Season", 
                  y = "observed_features", color = "Season", 
                  palette = c("black","darkgray"), add = "jitter", 
                  add.params = list(fill = "white"), 
                  ylab = "Observed ASVs") 
plot5 = ggpar(plot5, legend = "right", font.y = 16,
              font.legend = 16, font.ytickslab = 14) + 
  rremove("xlab") + 
  rremove("x.ticks")  + rremove("x.text") + 
  rremove("legend.title")

plot6 = ggboxplot(alpha_nounk, x = "Season", 
                  y = "shannon_entropy", color = "Season", 
                  palette = c("black","darkgray"), add = "jitter", 
                  add.params = list(fill = "white"), 
                  ylab = "Shannon Diversity") 
plot6 = ggpar(plot6, legend = "right", font.y = 16,
              font.legend = 16, font.ytickslab = 14) + 
  rremove("xlab") + 
  rremove("x.ticks")  + rremove("x.text") + 
  rremove("legend.title")

plot7 = ggboxplot(alpha_nounk, x = "Preservative", 
                  y = "faith_pd", color = "Preservative", 
                  palette = "Dark2", add = "jitter", 
                  add.params = list(fill = "white"), 
                  ylab = "Faith's Phylogenetic Diversity") 
plot7 = ggpar(plot7, legend = "right", font.y = 16,
              font.legend = 16, font.ytickslab = 14) + 
  rremove("xlab") +  
  rremove("x.text") + rremove("x.ticks") + 
  rremove("legend.title") + 
  stat_compare_means(label = "p.signif", 
                     comparisons = list(c("Ethanol", "RNAlater")),
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.05, 0.5, 1),
                                        symbols = c("*", "*", "*", "*", "ns")))

plot8 = ggboxplot(alpha_nounk, x = "Preservative", 
                  y = "observed_features", color = "Preservative", 
                  palette = "Dark2", add = "jitter", 
                  add.params = list(fill = "white"), 
                  ylab = "Observed ASVs") 
plot8 = ggpar(plot8, legend = "right", font.y = 16,
              font.legend = 16, font.ytickslab = 14) + 
  rremove("xlab") + 
  rremove("x.ticks")  + rremove("x.text") + 
  rremove("legend.title")

plot9 = ggboxplot(alpha_nounk, x = "Preservative", 
                  y = "shannon_entropy", color = "Preservative", 
                  palette = "Dark2", add = "jitter", 
                  add.params = list(fill = "white"), 
                  ylab = "Shannon Diversity") 
plot9 = ggpar(plot9, legend = "right", font.y = 16,
              font.legend = 16, font.ytickslab = 14) + 
  rremove("xlab") + 
  rremove("x.ticks")  + rremove("x.text") + 
  rremove("legend.title")

tiff(file="alpha_combined.tif", res=300, width=15, height=9, units="in")
alpha_plot = ggarrange(ggarrange(plot1, plot2, plot3, nrow = 1, ncol = 3, 
          common.legend = T, align = "h", legend = "right"),
          ggarrange(ggarrange(plot4, plot5, plot6, 
                    nrow = 1, ncol = 3, 
                    common.legend = T, align = "h", 
                    legend = "right"),
                    ggarrange(plot7, plot8, plot9, 
                              nrow = 1, ncol = 3, 
                              common.legend = T, align = "h", 
                              legend = "right"),
                    nrow = 1, ncol = 2,
                    align = "h", labels = c("B", "C")),
          nrow = 2, ncol = 1, align = "hv",
          labels = c("A", ""))
dev.off()

setEPS()
postscript(file="alpha_combined.eps", width=15, height=9, paper = "special")
alpha_plot
dev.off()

#Phyla differences----
library(tidyverse)
library(phyloseq)
library(ANCOMBC)

phyla = read_tsv("phyla-table-full.tsv")
metadata_full = read_csv("caatinga_metadata.csv")
taxonomy = phyla %>% dplyr::select(Phyla)

phyla_matrix = phyla %>% column_to_rownames("Phyla") %>% as.matrix()
phyla_phylo = otu_table(phyla_matrix, taxa_are_rows = T)
meta_phylo = metadata_full %>% column_to_rownames("SampleID") %>% 
  sample_data()

phyla_p = phyloseq(phyla_phylo, meta_phylo)

###Season ----

season = ancombc(phyla_p, "Season + Group + Age + Sex + Preservative",
                 p_adj_method = "fdr", lib_cut = 10000,
                 group = "Season")
res_s_df = data.frame(
  Species = row.names(season$res$beta),
  beta = unlist(season$res$beta),
  se = unlist(season$res$se),
  W = unlist(season$res$W),
  p_val = unlist(season$res$p_val),
  q_val = unlist(season$res$q_val),
  diff_abn = unlist(season$res$diff_abn))

fdr_ancom_s <- res_s_df %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Season")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_s, "Diff_abund_season_phyla.csv")

###Group global ----

group = ancombc(phyla_p, "Season + Group + Age + Sex + Preservative",
                p_adj_method = "fdr", lib_cut = 10000,
                group = "Group", global = T)

res_g = group$res_global %>% 
  filter(diff_abn == "TRUE") %>% 
  rownames_to_column(var = "Phyla") %>% 
  inner_join(taxonomy)

write_csv(res_g, "Diff_abund_group_phyla.csv")

###Group pairs ----

####Alg vs. Coq ----
meta_algcoq = metadata_full %>% 
  filter(Group == "Algaroba" | Group == "Coqueiro") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_algcoq = phyloseq(phyla_phylo, meta_algcoq)
group_algcoq = ancombc(phyla_p_algcoq, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_algcoq = data.frame(
  Species = row.names(group_algcoq$res$beta),
  beta = unlist(group_algcoq$res$beta),
  se = unlist(group_algcoq$res$se),
  W = unlist(group_algcoq$res$W),
  p_val = unlist(group_algcoq$res$p_val),
  q_val = unlist(group_algcoq$res$q_val),
  diff_abn = unlist(group_algcoq$res$diff_abn))

fdr_ancom_g_algcoq <- res_g_df_algcoq %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_algcoq, "Diff_abund_group_phyla_algcoq.csv")

####Alg vs. Cow ----
meta_algcow = metadata_full %>% 
  filter(Group == "Algaroba" | Group == "Cow") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_algcow = phyloseq(phyla_phylo, meta_algcow)
group_algcow = ancombc(phyla_p_algcow, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_algcow = data.frame(
  Species = row.names(group_algcow$res$beta),
  beta = unlist(group_algcow$res$beta),
  se = unlist(group_algcow$res$se),
  W = unlist(group_algcow$res$W),
  p_val = unlist(group_algcow$res$p_val),
  q_val = unlist(group_algcow$res$q_val),
  diff_abn = unlist(group_algcow$res$diff_abn))

fdr_ancom_g_algcow <- res_g_df_algcow %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_algcow, "Diff_abund_group_phyla_algcow.csv")

####Alg vs. F ----
meta_algf = metadata_full %>% 
  filter(Group == "Algaroba" | Group == "F group") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_algf = phyloseq(phyla_phylo, meta_algf)
group_algf = ancombc(phyla_p_algf, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_algf = data.frame(
  Species = row.names(group_algf$res$beta),
  beta = unlist(group_algf$res$beta),
  se = unlist(group_algf$res$se),
  W = unlist(group_algf$res$W),
  p_val = unlist(group_algf$res$p_val),
  q_val = unlist(group_algf$res$q_val),
  diff_abn = unlist(group_algf$res$diff_abn))

fdr_ancom_g_algf <- res_g_df_algf %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_algf, "Diff_abund_group_phyla_algf.csv")

####Alg vs. House ----
meta_alghou = metadata_full %>% 
  filter(Group == "Algaroba" | Group == "House") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_alghou = phyloseq(phyla_phylo, meta_alghou)
group_alghou = ancombc(phyla_p_alghou, "Group",
                     p_adj_method = "fdr", lib_cut = 10000,
                     group = "Group")

res_g_df_alghou = data.frame(
  Species = row.names(group_alghou$res$beta),
  beta = unlist(group_alghou$res$beta),
  se = unlist(group_alghou$res$se),
  W = unlist(group_alghou$res$W),
  p_val = unlist(group_alghou$res$p_val),
  q_val = unlist(group_alghou$res$q_val),
  diff_abn = unlist(group_alghou$res$diff_abn))

fdr_ancom_g_alghou <- res_g_df_alghou %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_alghou, "Diff_abund_group_phyla_alghou.csv")

####Alg vs. Key ----
meta_algkey = metadata_full %>% 
  filter(Group == "Algaroba" | Group == "Key") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_algkey = phyloseq(phyla_phylo, meta_algkey)
group_algkey = ancombc(phyla_p_algkey, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_algkey = data.frame(
  Species = row.names(group_algkey$res$beta),
  beta = unlist(group_algkey$res$beta),
  se = unlist(group_algkey$res$se),
  W = unlist(group_algkey$res$W),
  p_val = unlist(group_algkey$res$p_val),
  q_val = unlist(group_algkey$res$q_val),
  diff_abn = unlist(group_algkey$res$diff_abn))

fdr_ancom_g_algkey <- res_g_df_algkey %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_algkey, "Diff_abund_group_phyla_algkey.csv")

####Alg vs. Princess ----
meta_algprin = metadata_full %>% 
  filter(Group == "Algaroba" | Group == "Princess") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_algprin = phyloseq(phyla_phylo, meta_algprin)
group_algprin = ancombc(phyla_p_algprin, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_algprin = data.frame(
  Species = row.names(group_algprin$res$beta),
  beta = unlist(group_algprin$res$beta),
  se = unlist(group_algprin$res$se),
  W = unlist(group_algprin$res$W),
  p_val = unlist(group_algprin$res$p_val),
  q_val = unlist(group_algprin$res$q_val),
  diff_abn = unlist(group_algprin$res$diff_abn))

fdr_ancom_g_algprin <- res_g_df_algprin %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_algprin, "Diff_abund_group_phyla_algprin.csv")

####Alg vs. Road ----
meta_algroa = metadata_full %>% 
  filter(Group == "Algaroba" | Group == "Road") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_algroa = phyloseq(phyla_phylo, meta_algroa)
group_algroa = ancombc(phyla_p_algroa, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_algroa = data.frame(
  Species = row.names(group_algroa$res$beta),
  beta = unlist(group_algroa$res$beta),
  se = unlist(group_algroa$res$se),
  W = unlist(group_algroa$res$W),
  p_val = unlist(group_algroa$res$p_val),
  q_val = unlist(group_algroa$res$q_val),
  diff_abn = unlist(group_algroa$res$diff_abn))

fdr_ancom_g_algroa <- res_g_df_algroa %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_algroa, "Diff_abund_group_phyla_algroa.csv")

####Coqueiro vs. Cow ----
meta_coqcow = metadata_full %>% 
  filter(Group == "Coqueiro" | Group == "Cow") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_coqcow = phyloseq(phyla_phylo, meta_coqcow)
group_coqcow = ancombc(phyla_p_coqcow, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_coqcow = data.frame(
  Species = row.names(group_coqcow$res$beta),
  beta = unlist(group_coqcow$res$beta),
  se = unlist(group_coqcow$res$se),
  W = unlist(group_coqcow$res$W),
  p_val = unlist(group_coqcow$res$p_val),
  q_val = unlist(group_coqcow$res$q_val),
  diff_abn = unlist(group_coqcow$res$diff_abn))

fdr_ancom_g_coqcow <- res_g_df_coqcow %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_coqcow, "Diff_abund_group_phyla_coqcow.csv")

####Coqueiro vs. F ----
meta_coqf = metadata_full %>% 
  filter(Group == "Coqueiro" | Group == "F group") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_coqf = phyloseq(phyla_phylo, meta_coqf)
group_coqf = ancombc(phyla_p_coqf, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_coqf = data.frame(
  Species = row.names(group_coqf$res$beta),
  beta = unlist(group_coqf$res$beta),
  se = unlist(group_coqf$res$se),
  W = unlist(group_coqf$res$W),
  p_val = unlist(group_coqf$res$p_val),
  q_val = unlist(group_coqf$res$q_val),
  diff_abn = unlist(group_coqf$res$diff_abn))

fdr_ancom_g_coqf <- res_g_df_coqf %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_coqf, "Diff_abund_group_phyla_coqf.csv")

####Coqueiro vs. House ----
meta_coqhou = metadata_full %>% 
  filter(Group == "Coqueiro" | Group == "House") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_coqhou = phyloseq(phyla_phylo, meta_coqhou)
group_coqhou = ancombc(phyla_p_coqhou, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_coqhou = data.frame(
  Species = row.names(group_coqhou$res$beta),
  beta = unlist(group_coqhou$res$beta),
  se = unlist(group_coqhou$res$se),
  W = unlist(group_coqhou$res$W),
  p_val = unlist(group_coqhou$res$p_val),
  q_val = unlist(group_coqhou$res$q_val),
  diff_abn = unlist(group_coqhou$res$diff_abn))

fdr_ancom_g_coqhou <- res_g_df_coqhou %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_coqhou, "Diff_abund_group_phyla_coqhou.csv")

####Coqueiro vs. Key ----
meta_coqkey = metadata_full %>% 
  filter(Group == "Coqueiro" | Group == "Key") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_coqkey = phyloseq(phyla_phylo, meta_coqkey)
group_coqkey = ancombc(phyla_p_coqkey, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_coqkey = data.frame(
  Species = row.names(group_coqkey$res$beta),
  beta = unlist(group_coqkey$res$beta),
  se = unlist(group_coqkey$res$se),
  W = unlist(group_coqkey$res$W),
  p_val = unlist(group_coqkey$res$p_val),
  q_val = unlist(group_coqkey$res$q_val),
  diff_abn = unlist(group_coqkey$res$diff_abn))

fdr_ancom_g_coqkey <- res_g_df_coqkey %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_coqkey, "Diff_abund_group_phyla_coqkey.csv")

####Coqueiro vs. Princess ----
meta_coqprin = metadata_full %>% 
  filter(Group == "Coqueiro" | Group == "Princess") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_coqprin = phyloseq(phyla_phylo, meta_coqprin)
group_coqprin = ancombc(phyla_p_coqprin, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_coqprin = data.frame(
  Species = row.names(group_coqprin$res$beta),
  beta = unlist(group_coqprin$res$beta),
  se = unlist(group_coqprin$res$se),
  W = unlist(group_coqprin$res$W),
  p_val = unlist(group_coqprin$res$p_val),
  q_val = unlist(group_coqprin$res$q_val),
  diff_abn = unlist(group_coqprin$res$diff_abn))

fdr_ancom_g_coqprin <- res_g_df_coqprin %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_coqprin, "Diff_abund_group_phyla_coqprin.csv")

####Coqueiro vs. Road ----
meta_coqroa = metadata_full %>% 
  filter(Group == "Coqueiro" | Group == "Road") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_coqroa = phyloseq(phyla_phylo, meta_coqroa)
group_coqroa = ancombc(phyla_p_coqroa, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_coqroa = data.frame(
  Species = row.names(group_coqroa$res$beta),
  beta = unlist(group_coqroa$res$beta),
  se = unlist(group_coqroa$res$se),
  W = unlist(group_coqroa$res$W),
  p_val = unlist(group_coqroa$res$p_val),
  q_val = unlist(group_coqroa$res$q_val),
  diff_abn = unlist(group_coqroa$res$diff_abn))

fdr_ancom_g_coqroa <- res_g_df_coqroa %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_coqroa, "Diff_abund_group_phyla_coqroa.csv")

####Cow vs. F ----
meta_cowf = metadata_full %>% 
  filter(Group == "Cow" | Group == "F group") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_cowf = phyloseq(phyla_phylo, meta_cowf)
group_cowf = ancombc(phyla_p_cowf, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_cowf = data.frame(
  Species = row.names(group_cowf$res$beta),
  beta = unlist(group_cowf$res$beta),
  se = unlist(group_cowf$res$se),
  W = unlist(group_cowf$res$W),
  p_val = unlist(group_cowf$res$p_val),
  q_val = unlist(group_cowf$res$q_val),
  diff_abn = unlist(group_cowf$res$diff_abn))

fdr_ancom_g_cowf <- res_g_df_cowf %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_cowf, "Diff_abund_group_phyla_cowf.csv")

####Cow vs. House ----
meta_cowhou = metadata_full %>% 
  filter(Group == "Cow" | Group == "House") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_cowhou = phyloseq(phyla_phylo, meta_cowhou)
group_cowhou = ancombc(phyla_p_cowhou, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_cowhou = data.frame(
  Species = row.names(group_cowhou$res$beta),
  beta = unlist(group_cowhou$res$beta),
  se = unlist(group_cowhou$res$se),
  W = unlist(group_cowhou$res$W),
  p_val = unlist(group_cowhou$res$p_val),
  q_val = unlist(group_cowhou$res$q_val),
  diff_abn = unlist(group_cowhou$res$diff_abn))

fdr_ancom_g_cowhou <- res_g_df_cowhou %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_cowhou, "Diff_abund_group_phyla_cowhou.csv")

####Cow vs. Key ----
meta_cowkey = metadata_full %>% 
  filter(Group == "Cow" | Group == "Key") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_cowkey = phyloseq(phyla_phylo, meta_cowkey)
group_cowkey = ancombc(phyla_p_cowkey, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_cowkey = data.frame(
  Species = row.names(group_cowkey$res$beta),
  beta = unlist(group_cowkey$res$beta),
  se = unlist(group_cowkey$res$se),
  W = unlist(group_cowkey$res$W),
  p_val = unlist(group_cowkey$res$p_val),
  q_val = unlist(group_cowkey$res$q_val),
  diff_abn = unlist(group_cowkey$res$diff_abn))

fdr_ancom_g_cowkey <- res_g_df_cowkey %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_cowkey, "Diff_abund_group_phyla_cowkey.csv")

####Cow vs. Princess ----
meta_cowprin = metadata_full %>% 
  filter(Group == "Cow" | Group == "Princess") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_cowprin = phyloseq(phyla_phylo, meta_cowprin)
group_cowprin = ancombc(phyla_p_cowprin, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_cowprin = data.frame(
  Species = row.names(group_cowprin$res$beta),
  beta = unlist(group_cowprin$res$beta),
  se = unlist(group_cowprin$res$se),
  W = unlist(group_cowprin$res$W),
  p_val = unlist(group_cowprin$res$p_val),
  q_val = unlist(group_cowprin$res$q_val),
  diff_abn = unlist(group_cowprin$res$diff_abn))

fdr_ancom_g_cowprin <- res_g_df_cowprin %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_cowprin, "Diff_abund_group_phyla_cowprin.csv")

####Cow vs. Road ----
meta_cowroa = metadata_full %>% 
  filter(Group == "Cow" | Group == "Road") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_cowroa = phyloseq(phyla_phylo, meta_cowroa)
group_cowroa = ancombc(phyla_p_cowroa, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_cowroa = data.frame(
  Species = row.names(group_cowroa$res$beta),
  beta = unlist(group_cowroa$res$beta),
  se = unlist(group_cowroa$res$se),
  W = unlist(group_cowroa$res$W),
  p_val = unlist(group_cowroa$res$p_val),
  q_val = unlist(group_cowroa$res$q_val),
  diff_abn = unlist(group_cowroa$res$diff_abn))

fdr_ancom_g_cowroa <- res_g_df_cowroa %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_cowroa, "Diff_abund_group_phyla_cowroa.csv")

####F vs. House ----
meta_fhou = metadata_full %>% 
  filter(Group == "F group" | Group == "House") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_fhou = phyloseq(phyla_phylo, meta_fhou)
group_fhou = ancombc(phyla_p_fhou, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_fhou = data.frame(
  Species = row.names(group_fhou$res$beta),
  beta = unlist(group_fhou$res$beta),
  se = unlist(group_fhou$res$se),
  W = unlist(group_fhou$res$W),
  p_val = unlist(group_fhou$res$p_val),
  q_val = unlist(group_fhou$res$q_val),
  diff_abn = unlist(group_fhou$res$diff_abn))

fdr_ancom_g_fhou <- res_g_df_fhou %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_fhou, "Diff_abund_group_phyla_fhou.csv")

####F vs. Key ----
meta_fkey = metadata_full %>% 
  filter(Group == "F group" | Group == "Key") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_fkey = phyloseq(phyla_phylo, meta_fkey)
group_fkey = ancombc(phyla_p_fkey, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_fkey = data.frame(
  Species = row.names(group_fkey$res$beta),
  beta = unlist(group_fkey$res$beta),
  se = unlist(group_fkey$res$se),
  W = unlist(group_fkey$res$W),
  p_val = unlist(group_fkey$res$p_val),
  q_val = unlist(group_fkey$res$q_val),
  diff_abn = unlist(group_fkey$res$diff_abn))

fdr_ancom_g_fkey <- res_g_df_fkey %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_fkey, "Diff_abund_group_phyla_fkey.csv")

####F vs. Princess ----
meta_fprin = metadata_full %>% 
  filter(Group == "F group" | Group == "Princess") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_fprin = phyloseq(phyla_phylo, meta_fprin)
group_fprin = ancombc(phyla_p_fprin, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_fprin = data.frame(
  Species = row.names(group_fprin$res$beta),
  beta = unlist(group_fprin$res$beta),
  se = unlist(group_fprin$res$se),
  W = unlist(group_fprin$res$W),
  p_val = unlist(group_fprin$res$p_val),
  q_val = unlist(group_fprin$res$q_val),
  diff_abn = unlist(group_fprin$res$diff_abn))

fdr_ancom_g_fprin <- res_g_df_fprin %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_fprin, "Diff_abund_group_phyla_fprin.csv")

####F vs. Road ----
meta_froa = metadata_full %>% 
  filter(Group == "F group" | Group == "Road") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_froa = phyloseq(phyla_phylo, meta_froa)
group_froa = ancombc(phyla_p_froa, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_froa = data.frame(
  Species = row.names(group_froa$res$beta),
  beta = unlist(group_froa$res$beta),
  se = unlist(group_froa$res$se),
  W = unlist(group_froa$res$W),
  p_val = unlist(group_froa$res$p_val),
  q_val = unlist(group_froa$res$q_val),
  diff_abn = unlist(group_froa$res$diff_abn))

fdr_ancom_g_froa <- res_g_df_froa %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_froa, "Diff_abund_group_phyla_froa.csv")

####House vs. Key ----
meta_houkey = metadata_full %>% 
  filter(Group == "House" | Group == "Key") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_houkey = phyloseq(phyla_phylo, meta_houkey)
group_houkey = ancombc(phyla_p_houkey, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_houkey = data.frame(
  Species = row.names(group_houkey$res$beta),
  beta = unlist(group_houkey$res$beta),
  se = unlist(group_houkey$res$se),
  W = unlist(group_houkey$res$W),
  p_val = unlist(group_houkey$res$p_val),
  q_val = unlist(group_houkey$res$q_val),
  diff_abn = unlist(group_houkey$res$diff_abn))

fdr_ancom_g_houkey <- res_g_df_houkey %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_houkey, "Diff_abund_group_phyla_houkey.csv")

####House vs. Princess ----
meta_houprin = metadata_full %>% 
  filter(Group == "House" | Group == "Princess") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_houprin = phyloseq(phyla_phylo, meta_houprin)
group_houprin = ancombc(phyla_p_houprin, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_houprin = data.frame(
  Species = row.names(group_houprin$res$beta),
  beta = unlist(group_houprin$res$beta),
  se = unlist(group_houprin$res$se),
  W = unlist(group_houprin$res$W),
  p_val = unlist(group_houprin$res$p_val),
  q_val = unlist(group_houprin$res$q_val),
  diff_abn = unlist(group_houprin$res$diff_abn))

fdr_ancom_g_houprin <- res_g_df_houprin %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_houprin, "Diff_abund_group_phyla_houprin.csv")

####House vs. Road ----
meta_houroa = metadata_full %>% 
  filter(Group == "House" | Group == "Road") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_houroa = phyloseq(phyla_phylo, meta_houroa)
group_houroa = ancombc(phyla_p_houroa, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_houroa = data.frame(
  Species = row.names(group_houroa$res$beta),
  beta = unlist(group_houroa$res$beta),
  se = unlist(group_houroa$res$se),
  W = unlist(group_houroa$res$W),
  p_val = unlist(group_houroa$res$p_val),
  q_val = unlist(group_houroa$res$q_val),
  diff_abn = unlist(group_houroa$res$diff_abn))

fdr_ancom_g_houroa <- res_g_df_houroa %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_houroa, "Diff_abund_group_phyla_houroa.csv")

####Key vs. Princess ----
meta_keyprin = metadata_full %>% 
  filter(Group == "Key" | Group == "Princess") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_keyprin = phyloseq(phyla_phylo, meta_keyprin)
group_keyprin = ancombc(phyla_p_keyprin, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_keyprin = data.frame(
  Species = row.names(group_keyprin$res$beta),
  beta = unlist(group_keyprin$res$beta),
  se = unlist(group_keyprin$res$se),
  W = unlist(group_keyprin$res$W),
  p_val = unlist(group_keyprin$res$p_val),
  q_val = unlist(group_keyprin$res$q_val),
  diff_abn = unlist(group_keyprin$res$diff_abn))

fdr_ancom_g_keyprin <- res_g_df_keyprin %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_keyprin, "Diff_abund_group_phyla_keyprin.csv")

####Key vs. Road ----
meta_keyroa = metadata_full %>% 
  filter(Group == "Key" | Group == "Road") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_keyroa = phyloseq(phyla_phylo, meta_keyroa)
group_keyroa = ancombc(phyla_p_keyroa, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_keyroa = data.frame(
  Species = row.names(group_keyroa$res$beta),
  beta = unlist(group_keyroa$res$beta),
  se = unlist(group_keyroa$res$se),
  W = unlist(group_keyroa$res$W),
  p_val = unlist(group_keyroa$res$p_val),
  q_val = unlist(group_keyroa$res$q_val),
  diff_abn = unlist(group_keyroa$res$diff_abn))

fdr_ancom_g_keyroa <- res_g_df_keyroa %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_keyroa, "Diff_abund_group_phyla_keyroa.csv")

####Princess vs. Road ----
meta_prinroa = metadata_full %>% 
  filter(Group == "Princess" | Group == "Road") %>% 
  column_to_rownames("SampleID") %>% 
  sample_data()
phyla_p_prinroa = phyloseq(phyla_phylo, meta_prinroa)
group_prinroa = ancombc(phyla_p_prinroa, "Group",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "Group")

res_g_df_prinroa = data.frame(
  Species = row.names(group_prinroa$res$beta),
  beta = unlist(group_prinroa$res$beta),
  se = unlist(group_prinroa$res$se),
  W = unlist(group_prinroa$res$W),
  p_val = unlist(group_prinroa$res$p_val),
  q_val = unlist(group_prinroa$res$q_val),
  diff_abn = unlist(group_prinroa$res$diff_abn))

fdr_ancom_g_prinroa <- res_g_df_prinroa %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Phyla"))

write_csv(fdr_ancom_g_prinroa, "Diff_abund_group_phyla_prinroa.csv")



#Phyla graphs----
library(ggpubr)

phyla_t = read_tsv("phyla-table-full-t.txt")
phyla_meta = metadata_full %>% inner_join(phyla_t)

bact = ggboxplot(phyla_meta, x = "Group", 
                  y = "Bacteroidetes", color = "Group", 
                  palette = "Set1", add = "jitter", 
                  add.params = list(fill = "white"), 
                  ylab = "Relative Abundance", 
                 title = "Bacteroidetes") 
bact = ggpar(bact, legend = "right", font.y = 16,
             font.legend = 16, font.ytickslab = 14) + rremove("xlab") + 
  rremove("x.text") + rremove("x.ticks") + 
  rremove("legend.title") + 
  stat_compare_means(label = "p.signif", 
                     comparisons = list(c("Cow", "Road")),
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.05, 0.5, 1),
                                        symbols = c("*", "*", "*", "*", "*")))

cyan = ggboxplot(phyla_meta, x = "Group", 
                 y = "Cyanobacteria", color = "Group", 
                 palette = "Set1", add = "jitter", 
                 add.params = list(fill = "white"), 
                 ylab = "Relative Abundance", 
                 title = "Cyanobacteria") 
cyan = ggpar(cyan, legend = "right", font.y = 16,
             font.legend = 16, font.ytickslab = 14) + rremove("xlab") + 
  rremove("x.text") + rremove("x.ticks") + 
  rremove("legend.title") + 
  stat_compare_means(label = "p.signif", 
                     comparisons = list(c("Algaroba", "House"),
                                        c("Coqueiro", "House"),
                                        c("Cow", "House"),
                                        c("Cow", "Road"),
                                        c("F group", "House"),
                                        c("House", "Key"),
                                        c("House", "Princess")),
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.05, 0.5, 1),
                                        symbols = c("*", "*", "*", "*", "*")))

firm = ggboxplot(phyla_meta, x = "Group", 
                 y = "Firmicutes", color = "Group", 
                 palette = "Set1", add = "jitter", 
                 add.params = list(fill = "white"), 
                 ylab = "Relative Abundance", 
                 title = "Firmicutes") 
firm = ggpar(firm, legend = "right", font.y = 16,
             font.legend = 16, font.ytickslab = 14) + rremove("xlab") + 
  rremove("x.text") + rremove("x.ticks") + 
  rremove("legend.title") + 
  stat_compare_means(label = "p.signif", 
                     comparisons = list(c("Algaroba", "Cow"),
                                        c("Algaroba", "Key"),
                                        c("Algaroba", "Princess"),
                                        c("Coqueiro", "House"),
                                        c("Cow", "House"),
                                        c("F group", "House"),
                                        c("House", "Key"),
                                        c("House", "Princess")),
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.05, 0.5, 1),
                                        symbols = c("*", "*", "*", "*", "*")))

fuso = ggboxplot(phyla_meta, x = "Group", 
                 y = "Fusobacteria", color = "Group", 
                 palette = "Set1", add = "jitter", 
                 add.params = list(fill = "white"), 
                 ylab = "Relative Abundance", 
                 title = "Fusobacteria") 
fuso = ggpar(fuso, legend = "right", font.y = 16,
             font.legend = 16, font.ytickslab = 14) + rremove("xlab") + 
  rremove("x.text") + rremove("x.ticks") + 
  rremove("legend.title") + 
  stat_compare_means(label = "p.signif", 
                     comparisons = list(c("Algaroba", "Key"),
                                        c("Algaroba", "Princess"),
                                        c("Coqueiro", "Cow"),
                                        c("Cow", "House"),
                                        c("Cow", "Key"),
                                        c("Cow", "Princess"),
                                        c("Cow", "Road"),
                                        c("House", "Key")),
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.05, 0.5, 1),
                                        symbols = c("*", "*", "*", "*", "*")))

prot = ggboxplot(phyla_meta, x = "Group", 
                 y = "Proteobacteria", color = "Group", 
                 palette = "Set1", add = "jitter", 
                 add.params = list(fill = "white"), 
                 ylab = "Relative Abundance", 
                 title = "Proteobacteria") 
prot = ggpar(prot, legend = "right", font.y = 16,
             font.legend = 16, font.ytickslab = 14) + rremove("xlab") + 
  rremove("x.text") + rremove("x.ticks") + 
  rremove("legend.title") + 
  stat_compare_means(label = "p.signif", 
                     comparisons = list(c("Algaroba", "Cow"), 
                                        c("Algaroba", "F group"),
                                        c("Algaroba", "Key"),
                                        c("Algaroba", "Princess"),
                                        c("Coqueiro", "Cow"),
                                        c("Cow", "House"),
                                        c("Cow", "Key"),
                                        c("Cow", "Princess"),
                                        c("Cow", "Road"),
                                        c("F group", "Key")),
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.05, 0.5, 1),
                                        symbols = c("*", "*", "*", "*", "*")))

tm7 = ggboxplot(phyla_meta, x = "Group", 
                 y = "TM7", color = "Group", 
                 palette = "Set1", add = "jitter", 
                 add.params = list(fill = "white"), 
                 ylab = "Relative Abundance", 
                 title = "TM7") 
tm7 = ggpar(tm7, legend = "right", font.y = 16,
            font.legend = 16, font.ytickslab = 14) + rremove("xlab") + 
  rremove("x.text") + rremove("x.ticks") + 
  rremove("legend.title")

tiff(file="phyla_combined.tif", res=300, width=15, height=10, units="in")
phyla_graph = ggarrange(bact, cyan, firm, fuso, prot, tm7, nrow = 2, ncol = 3, 
          common.legend = T, align = "h", legend = "right")
phyla_graph
dev.off()

setEPS()
postscript(file="phyla_combined.eps", width=15, height=10, paper = "special")
phyla_graph
dev.off()

#Family differences ----
library(tidyverse)
library(phyloseq)
library(ANCOMBC)

family = read_tsv("family-table-full.tsv")
metadata_full = read_csv("caatinga_metadata.csv")
taxonomy = family %>% dplyr::select(Family)

family_matrix = family %>% column_to_rownames("Family") %>% as.matrix()
family_phylo = otu_table(family_matrix, taxa_are_rows = T)
meta_phylo = metadata_full %>% column_to_rownames("SampleID") %>% 
  sample_data()

family_p = phyloseq(family_phylo, meta_phylo)

season = ancombc(family_p, "Season + Group + Age + Sex + Preservative",
                 p_adj_method = "fdr", lib_cut = 10000,
                 group = "Season")
res_s_df = data.frame(
  Species = row.names(season$res$beta),
  beta = unlist(season$res$beta),
  se = unlist(season$res$se),
  W = unlist(season$res$W),
  p_val = unlist(season$res$p_val),
  q_val = unlist(season$res$q_val),
  diff_abn = unlist(season$res$diff_abn))

fdr_ancom_s <- res_s_df %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Season")) %>% 
  inner_join(taxonomy, by = c("Species" = "Family"))

write_csv(fdr_ancom_s, "Diff_abund_season_family.csv")

group = ancombc(family_p, "Season + Group + Age + Sex + Preservative",
                p_adj_method = "fdr", lib_cut = 10000,
                group = "Group", global = T)

res_g = group$res_global %>% 
  filter(diff_abn == "TRUE") %>% 
  rownames_to_column(var = "Family") %>% 
  inner_join(taxonomy)

res_g_df = data.frame(
  Species = row.names(season$res$beta),
  beta = unlist(season$res$beta),
  se = unlist(season$res$se),
  W = unlist(season$res$W),
  p_val = unlist(season$res$p_val),
  q_val = unlist(season$res$q_val),
  diff_abn = unlist(season$res$diff_abn))

fdr_ancom_g <- res_g_df %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Family"))

write_csv(res_g, "Diff_abund_group_family.csv")

#Genera differences ----
library(tidyverse)
library(phyloseq)
library(ANCOMBC)

genus = read_tsv("genus-table-full.tsv")
metadata_full = read_csv("caatinga_metadata.csv")
taxonomy = genus %>% dplyr::select(Genus)

genus_matrix = genus %>% column_to_rownames("Genus") %>% as.matrix()
genus_phylo = otu_table(genus_matrix, taxa_are_rows = T)
meta_phylo = metadata_full %>% column_to_rownames("SampleID") %>% 
  sample_data()

genus_p = phyloseq(genus_phylo, meta_phylo)

season = ancombc(genus_p, "Season + Group + Age + Sex + Preservative",
                 p_adj_method = "fdr", lib_cut = 10000,
                 group = "Season")
res_s_df = data.frame(
  Species = row.names(season$res$beta),
  beta = unlist(season$res$beta),
  se = unlist(season$res$se),
  W = unlist(season$res$W),
  p_val = unlist(season$res$p_val),
  q_val = unlist(season$res$q_val),
  diff_abn = unlist(season$res$diff_abn))

fdr_ancom_s <- res_s_df %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Season")) %>% 
  inner_join(taxonomy, by = c("Species" = "Genus"))

write_csv(fdr_ancom_s, "Diff_abund_season_genus.csv")

group = ancombc(genus_p, "Season + Group + Age + Sex + Preservative",
                p_adj_method = "fdr", lib_cut = 10000,
                group = "Group", global = T)

res_g = group$res_global %>% 
  filter(diff_abn == "TRUE") %>% 
  rownames_to_column(var = "Genus") %>% 
  inner_join(taxonomy)

res_g_df = data.frame(
  Species = row.names(season$res$beta),
  beta = unlist(season$res$beta),
  se = unlist(season$res$se),
  W = unlist(season$res$W),
  p_val = unlist(season$res$p_val),
  q_val = unlist(season$res$q_val),
  diff_abn = unlist(season$res$diff_abn))

fdr_ancom_g <- res_g_df %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "Genus"))

write_csv(res_g, "Diff_abund_group_genus.csv")

#ASVS differences ---- 
library(tidyverse)
library(phyloseq)
library(ANCOMBC)

asvs = read_tsv("feature-table-full.tsv")
metadata_full = read_csv("caatinga_metadata.csv")
taxonomy = asvs %>% dplyr::select(ASV:taxonomy)

asvs_matrix = asvs %>% dplyr::select(-taxonomy) %>% 
  column_to_rownames("ASV") %>% as.matrix()
asvs_phylo = otu_table(asvs_matrix, taxa_are_rows = T)
meta_phylo = metadata_full %>% column_to_rownames("SampleID") %>% 
  sample_data()

asvs_p = phyloseq(asvs_phylo, meta_phylo)

season = ancombc(asvs_p, "Season + Group + Age + Sex + Preservative",
              p_adj_method = "fdr", lib_cut = 10000,
              group = "Season")
res_s_df = data.frame(
  Species = row.names(season$res$beta),
  beta = unlist(season$res$beta),
  se = unlist(season$res$se),
  W = unlist(season$res$W),
  p_val = unlist(season$res$p_val),
  q_val = unlist(season$res$q_val),
  diff_abn = unlist(season$res$diff_abn))

fdr_ancom_s <- res_s_df %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Season")) %>% 
  inner_join(taxonomy, by = c("Species" = "ASV"))

write_csv(fdr_ancom_s, "Diff_abund_season_asvs.csv")

group = ancombc(asvs_p, "Season + Group + Age + Sex + Preservative",
                 p_adj_method = "fdr", lib_cut = 10000,
                 group = "Group", global = T)

res_g = group$res_global %>% 
  filter(diff_abn == "TRUE") %>% 
  rownames_to_column(var = "ASV") %>% 
  inner_join(taxonomy)

res_g_df = data.frame(
  Species = row.names(season$res$beta),
  beta = unlist(season$res$beta),
  se = unlist(season$res$se),
  W = unlist(season$res$W),
  p_val = unlist(season$res$p_val),
  q_val = unlist(season$res$q_val),
  diff_abn = unlist(season$res$diff_abn))

fdr_ancom_g <- res_g_df %>%
  dplyr::filter(q_val < 0.05) %>% 
  rownames_to_column(var = "Comparison") %>% 
  dplyr::filter(str_detect(Comparison, "^Group")) %>% 
  inner_join(taxonomy, by = c("Species" = "ASV"))

write_csv(res_g, "Diff_abund_group_asvs.csv")

#NMDS plots----
library(ggplot2)
library(ggtext)
library(vegan)
library(cowplot)

mds_otus_weighted<-metaMDS(weighted_nounk, k=2, trymax=499)
mds_otus_weighted_points<-mds_otus_weighted$points
mds_otus_weighted_points2<-merge(x=mds_otus_weighted_points, y = metadata_nounk, 
                                 by.x = "row.names", by.y = "SampleID")


mds_otus_unweighted<-metaMDS(unweighted_nounk, k=2, trymax=499)
mds_otus_unweighted_points<-mds_otus_unweighted$points
mds_otus_unweighted_points2<-merge(x=mds_otus_unweighted_points, y = metadata_nounk, 
                                   by.x = "row.names", by.y = "SampleID")


w_taxa <- ggplot(mds_otus_weighted_points2, 
                 aes(x = MDS1, y = MDS2, color = Group, 
                     shape = Period)) +
  geom_point(size=3) + scale_color_brewer(palette = 'Set1') +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
        legend.key=element_blank()) + 
  theme(axis.title.x=element_text(size=rel(2)), 
        axis.title.y=element_text(size=rel(2)),
        plot.title = element_text(size=rel(3)),
        legend.title = element_text(size=rel(2)),
        legend.text = element_text(size = rel(1.8))) + 
  ggtitle("Weighted UniFrac") +
  stat_ellipse(aes(x = MDS1, y = MDS2, group = Period, 
                   linetype = Period), 
               type = "t", level = 0.9) + 
  scale_linetype_manual(values = c(1,2)) +
  annotate(geom = "richtext", fill = NA, label.color = NA,
           label = "Period: p = 0.167, R<sup>2</sup> = 2.2%<br>
           <b>Group: p < 0.001, R<sup>2</sup> = 27.9%</b><br>
           Age: p = 0.147, R<sup>2</sup> = 6.0%<br>
           Sex: p = 0.164, R<sup>2</sup> = 2.1%<br>
           Preservative: p = 0.129, R<sup>2</sup> = 2.5%", 
           x = -Inf, y = Inf, size = 5,
           hjust = 0, vjust = 1)
w_taxa

uw_taxa <- ggplot(mds_otus_unweighted_points2, 
                  aes(x = MDS1, y = MDS2, color = Group, 
                      shape = Period)) +
  geom_point(size=3) + scale_color_brewer(palette = 'Set1') +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
        legend.key=element_blank()) + 
  theme(axis.title.x=element_text(size=rel(2)), 
        axis.title.y=element_text(size=rel(2)),
        plot.title = element_text(size=rel(3)),
        legend.title = element_text(size=rel(2)),
        legend.text = element_text(size = rel(1.8))) + 
  ggtitle("Unweighted UniFrac") +
  stat_ellipse(aes(x = MDS1, y = MDS2, group = Period, linetype = Period), 
               type = "t", level = 0.9) + 
  scale_linetype_manual(values = c(1,2)) +
  annotate(geom = "richtext", fill = NA, label.color = NA,
           label = "<b>Period: p = 0.011, R<sup>2</sup> = 4.0%<br>
           Group: p < 0.001, R<sup>2</sup> = 20.3%</b><br>
           Age: p = 0.211, R<sup>2</sup> = 5.2%<br>
           Sex: p = 0.162, R<sup>2</sup> = 1.9%<br>
           <b>Preservative: p < 0.001, R<sup>2</sup> = 6.6%</b>", 
           x = -Inf, 
           y = Inf, 
           size = 5, hjust = 0, vjust = 1)
uw_taxa

tiff(file = "nmds_plot_combined_taxa_forpub.tif", 
     res = 300, width = 18, height = 8, units="in")
legend1 = get_legend(w_taxa + theme(legend.box.margin = margin(0, 0, 0, 12)))
col1 = plot_grid(w_taxa + theme(legend.position = "none"),
                 nrow = 1, ncol = 1, labels = c('A'), 
                 label_size = 20)
col2 = plot_grid(uw_taxa + theme(legend.position = "none"),
                 nrow = 1, ncol = 1, labels = c('B'),
                 label_size = 20)
nmds = plot_grid(col1, col2, legend1, 
          nrow = 1, ncol = 3, rel_widths = c(1.5, 1.5, 0.75), 
          align = "hv", axis = "t")
nmds
dev.off()

setEPS()
postscript(file="nmds_plot_combined_taxa_forpub.eps", width=18, height=8, paper = "special")
nmds
dev.off()

#Taxa plots----
library(ggpubr)
#ggplot is great, but the syntax is confusing. ggpubr wraps ggplot in easier-to-understand syntax

tiff(file="firm_group.tif", res=150, width=8, height=4, units="in")
plot1 = ggviolin(phyla, x = "Group", 
                 y = "Relative_Firmicutes", fill = "Group", 
                 palette = "Set1", add = "boxplot", 
                 add.params = list(fill = "white"), ylab = "Relative Firmicutes") 
ggpar(plot1, legend = "right") + rremove("xlab") + rremove("x.text") + rremove("legend.title")
dev.off()

library(limma)
sv = read.csv('Season_Venn.csv', header=TRUE)
a= vennCounts(sv)
vennDiagram(a,include='both')
