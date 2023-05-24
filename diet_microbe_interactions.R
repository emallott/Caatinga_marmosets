#Set up environment----
setwd("/Users/mallott/Dropbox/Projects/Gut_microbiome/Caatinga_marmosets/")

library(vegan)
library(tidyverse)

#16S-diet associations----

##Mantel tests----

metadata = read.csv("16S_results/caatinga_metadata_r_10000.csv", header=T)

bray_16S = as.dist(read.table("16S_results/bray-distance-matrix.tsv", header = T))

inverts = read.csv("diet_data/invert_marmoset_filtered_family.csv", header=T)
invert_meta = inner_join(metadata, inverts, by = c("SampleID" = "ID"))
invert_relab = decostand(invert_meta[,17:41], method = "total", MARGIN = 1)
bray_invert = vegdist(invert_relab, method = "bray")

plants = read.csv("diet_data/plant_marmoset_filtered_family.csv", header=T)
plant_meta = inner_join(metadata, plants, by = c("SampleID" = "ID"))
plant_relab = decostand(plant_meta[,17:61], method = "total", MARGIN = 1)
bray_plant = vegdist(plant_relab, method = "bray")

###16S-inverts----
mantel(bray_invert, bray_16S, method = "spearman", permutations = 499)

###16S-plants----
mantel(bray_plant, bray_16S, method = "spearman", permutations = 499)

##Diversity correlations----

faith_16S = read.table("16S_results/faithpd.tsv", header=T)
otus_16S = read.table("16S_results/observed-otus.tsv", header = T)
shannon_16S = read.table("16S_results/shannon.tsv", header = T)
alpha = inner_join(metadata, faith_16S, by = "SampleID") %>% 
  inner_join(otus_16S, by = "SampleID") %>% 
  inner_join(shannon_16S, by = "SampleID")

metadata_diet = read_csv("diet_data/caatinga_metadata_diet.csv")

invert_filtered_family = read_csv("diet_data/invert_marmoset_filtered_family.csv")
invert_family_richness = apply(invert_filtered_family[,2:26] > 0, 1, sum)
invert_family_richness = as.data.frame(invert_family_richness)
invert_family_richness = invert_family_richness %>% 
  cbind(invert_filtered_family[,27:45]) %>%
  dplyr::select(1:2)

plant_filtered_family = read_csv("diet_data/plant_marmoset_filtered_family.csv")
plant_family_richness = apply(plant_filtered_family[,2:46] > 0, 1, sum)
plant_family_richness = as.data.frame(plant_family_richness)
plant_family_richness = plant_family_richness %>% 
  cbind(plant_filtered_family[,47:65]) %>%
  dplyr::select(1:3)

inverts_raw_div = read_tsv("diet_data/invert_results_clean_marmoset_div.txt")
invert_meta_div = inverts_raw_div %>% 
  gather(SampleID, count, 2:68) %>% 
  spread(ID, count) %>% 
  left_join(metadata_diet, by = "SampleID") 
invert_filtered_div = invert_meta_div %>% 
  mutate(total = rowSums(invert_meta_div[,2:1369])) %>% 
  filter(total > 0) %>% 
  select_if(~ !is.numeric(.) || sum(.) != 0)

plants_raw_div = read_tsv("diet_data/plant_results_clean_marmoset_div.txt")
plant_meta_div = plants_raw_div %>% 
  gather(SampleID, count, 2:68) %>% 
  spread(ID, count) %>% 
  left_join(metadata_diet, by = "SampleID") 
plant_filtered_div = plant_meta_div %>% 
  mutate(total = rowSums(plant_meta_div[,2:5322])) %>% 
  filter(total > 0) %>% 
  select_if(~ !is.numeric(.) || sum(.) != 0)

plant_shannon = as.data.frame(diversity(plant_filtered_div[,2:1876], index = "shannon"))
invert_shannon = as.data.frame(diversity(invert_filtered_div[,2:465], index = "shannon"))

plant_shannon = plant_shannon %>% 
  rename("shannon_plant" = "diversity(plant_filtered_div[, 2:1876], index = \"shannon\")")
invert_shannon = invert_shannon %>% 
  rename("shannon_invert" = "diversity(invert_filtered_div[, 2:465], index = \"shannon\")")

plant_div = cbind(plant_filtered_div[,1877:1895], plant_shannon) %>% 
  dplyr::select(ID, shannon_plant)
invert_div = cbind(invert_filtered_div[,466:484], invert_shannon) %>% 
  dplyr::select(ID, shannon_invert)

alpha_comp_richness = alpha %>% 
  inner_join(invert_family_richness, by = c("SampleID" = "ID"))  %>% 
  inner_join(plant_family_richness, by = c("SampleID" = "ID"))

alpha_comp_shannon = alpha %>% 
  inner_join(plant_div[], by = c("SampleID" = "ID")) %>% 
  inner_join(invert_div, by = c("SampleID" = "ID")) 

###16S-inverts----
cor.test(alpha_comp_richness$observed_otus, 
         alpha_comp_richness$invert_family_richness, 
         method = "spearman")

cor.test(alpha_comp_shannon$shannon_entropy, 
         alpha_comp_shannon$shannon_invert, 
         method = "spearman")

###16S-plants----
cor.test(alpha_comp_richness$observed_otus, 
         alpha_comp_richness$plant_family_richness, 
         method = "spearman")

cor.test(alpha_comp_shannon$shannon_entropy, 
         alpha_comp_shannon$shannon_plant, 
         method = "spearman")

###Correlation scatterplots----

plot.invert.s = ggscatter(alpha_comp_shannon, x = "shannon_invert",
                          y = "shannon_entropy",
                          add = "reg.line", conf.int = T,
                          xlab = "Shannon Diversity - Invertebrates",
                          ylab = "Shannon Diversity - Microbes")

plot.plant.s = ggscatter(alpha_comp_shannon, x = "shannon_plant",
                         y = "shannon_entropy", add = "reg.line", 
                         conf.int = T,
                         xlab = "Shannon Diversity - Plants",
                         ylab = "Shannon Diversity - Microbes")

tiff(file="diet_microb_shannon_corr.tif", res=300, width=8, height=4, units="in")
ggarrange(plot.plant.s, plot.invert.s,  
          ncol = 2, nrow=1, 
          align = "h")
dev.off()

#Shotgun-diet associations----

##Mantel tests----

metadata_gene = read.csv("shotgun/caatinga_shotgun_gf_metadata.csv", header = T)
metadata_pa = read.csv("shotgun/caatinga_shotgun_metadata.csv", header = T)

bray_gene = as.dist(read.table("shotgun/bray_genefamilies_unstrat_mantel.txt", header = T))
bray_path = as.dist(read.table("shotgun/bray_pathabund_unstrat_mantel.txt", header = T))

inverts = read.csv("diet_data/invert_marmoset_filtered_family.csv", header=T)
invert_meta = inner_join(metadata_gene, inverts, by = c("Monkey" = "ID"))
invert_relab = decostand(invert_meta[,18:42], method = "total", MARGIN = 1)
bray_invert = vegdist(invert_relab, method = "bray")

plants = read.csv("diet_data/plant_marmoset_filtered_family.csv", header=T)
plant_meta = inner_join(metadata_gene, plants, by = c("Monkey" = "ID"))
plant_relab = decostand(plant_meta[,18:62], method = "total", MARGIN = 1)
bray_plant = vegdist(plant_relab, method = "bray")

###Function-inverts----
mantel(bray_invert, bray_path, method = "spearman", permutations = 499)
mantel(bray_invert, bray_gene, method = "spearman", permutations = 499)

###Function-plants----
mantel(bray_plant, bray_path, method = "spearman", permutations = 499)
mantel(bray_plant, bray_gene, method = "spearman", permutations = 499)

##Diversity correlations----

shannon_gene = read.table("shotgun/shannon_genefam_unstrat.tsv", header = T)
evenness_gene = read.table("shotgun/evenness_genefam_unstrat.tsv", header = T)
obs_otus_gene = read.table("shotgun/observed_otus_genefam_unstrat.tsv", header = T)

shannon_pa = read.table("shotgun/shannon_pathabund_unstrat.tsv", header = T)
evenness_pa = read.table("shotgun/evenness_pathabund_unstrat.tsv", header = T)
obs_otus_pa = read.table("shotgun/observed_otus_pathabund_unstrat.tsv", header = T)

alpha_gene = inner_join(metadata_gene, shannon_gene, by = "SampleID") %>% 
  inner_join(evenness_gene, by = "SampleID") %>% 
  inner_join(obs_otus_gene, by = "SampleID")
alpha_gene$Group = as.factor(alpha_gene$Group)

alpha_pa = inner_join(metadata_pa, shannon_pa, by = "SampleID") %>% 
  inner_join(evenness_pa, by = "SampleID") %>% 
  inner_join(obs_otus_pa, by = "SampleID")

invert_filtered_family = read_csv("diet_data/invert_marmoset_filtered_family.csv")
invert_family_richness = apply(invert_filtered_family[,2:26] > 0, 1, sum)
invert_family_richness = as.data.frame(invert_family_richness)
invert_family_richness = invert_family_richness %>% 
  cbind(invert_filtered_family[,27:45]) %>%
  dplyr::select(1:2)

plant_filtered_family = read_csv("diet_data/plant_marmoset_filtered_family.csv")
plant_family_richness = apply(plant_filtered_family[,2:46] > 0, 1, sum)
plant_family_richness = as.data.frame(plant_family_richness)
plant_family_richness = plant_family_richness %>% 
  cbind(plant_filtered_family[,47:65]) %>%
  dplyr::select(1:3)

inverts_raw_div = read_tsv("diet_data/invert_results_clean_marmoset_div.txt")
invert_meta_div = inverts_raw_div %>% 
  gather(SampleID, count, 2:68) %>% 
  spread(ID, count) %>% 
  mutate(Monkey = str_extract(SampleID, pattern = "[^_]*")) %>% 
  left_join(metadata_gene, by = c("Monkey" = "Monkey1"))
invert_filtered_div = invert_meta_div %>% 
  mutate(total = rowSums(invert_meta_div[,2:1369])) %>% 
  filter(total > 0) %>% 
  select_if(~ !is.numeric(.) || sum(.) != 0)

plants_raw_div = read_tsv("diet_data/plant_results_clean_marmoset_div.txt")
plant_meta_div = plants_raw_div %>% 
  gather(SampleID, count, 2:68) %>% 
  spread(ID, count)  %>% 
  mutate(Monkey = str_extract(SampleID, pattern = "[^_]*")) %>% 
  left_join(metadata_gene, by = c("Monkey" = "Monkey1"))
plant_filtered_div = plant_meta_div %>% 
  mutate(total = rowSums(plant_meta_div[,2:5322])) %>% 
  filter(total > 0) %>% 
  select_if(~ !is.numeric(.) || sum(.) != 0)

plant_shannon = as.data.frame(diversity(plant_filtered_div[,2:1876], index = "shannon"))
invert_shannon = as.data.frame(diversity(invert_filtered_div[,2:465], index = "shannon"))

plant_shannon = plant_shannon %>% 
  rename("shannon_plant" = "diversity(plant_filtered_div[, 2:1876], index = \"shannon\")")
invert_shannon = invert_shannon %>% 
  rename("shannon_invert" = "diversity(invert_filtered_div[, 2:465], index = \"shannon\")")

plant_div = cbind(plant_filtered_div[,1877:1894], plant_shannon) %>% 
  dplyr::select(Monkey, shannon_plant)
invert_div = cbind(invert_filtered_div[,466:483], invert_shannon) %>% 
  dplyr::select(Monkey, shannon_invert)

alpha_gene_comp_richness = alpha_gene %>% 
  inner_join(invert_family_richness, by = c("Monkey" = "ID")) %>% 
  inner_join(plant_family_richness, by = c("Monkey" = "ID"))

alpha_path_comp_richness = alpha_pa %>% 
  inner_join(invert_family_richness, by = c("Monkey" = "ID")) %>% 
  inner_join(plant_family_richness, by = c("Monkey" = "ID"))

###Function-inverts----
cor.test(alpha_gene_comp_richness$observed_features, 
         alpha_gene_comp_richness$invert_family_richness, 
         method = "spearman")

cor.test(alpha_path_comp_richness$observed_features, 
         alpha_path_comp_richness$invert_family_richness, 
         method = "spearman")

###Function-plants----
cor.test(alpha_gene_comp_richness$observed_features, 
         alpha_gene_comp_richness$plant_family_richness, 
         method = "spearman")

cor.test(alpha_path_comp_richness$observed_features, 
         alpha_path_comp_richness$plant_family_richness, 
         method = "spearman")


