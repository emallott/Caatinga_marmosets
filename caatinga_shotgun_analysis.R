#Set up environment ----

setwd("/Users/mallott/Dropbox/Projects/Gut_microbiome/Caatinga_marmosets/shotgun")

#Permanovas ----

library(vegan)
library(pairwiseAdonis)

bray_gene = as.dist(read.table("bray_genefamilies_unstrat.tsv", header = T))
jaccard_gene = as.dist(read.table("jaccard_genefamilies_unstrat.tsv", header = T))
bray_pa = as.dist(read.table("bray_pathabund_unstrat.tsv", header = T))
jaccard_pa = as.dist(read.table("jaccard_pathabund_unstrat.tsv", header = T))
metadata_gene = read.csv("caatinga_shotgun_gf_metadata.csv", header = T)
metadata_pa = read.csv("caatinga_shotgun_metadata.csv", header = T)

adonis2(bray_gene ~ Period + Group + Sex + Preservative, 
        data=metadata_gene, 
        by = "margin", permutations = 4999)
adonis2(bray_gene ~ Season + Human_food + Sex + Preservative, 
        data=metadata_gene, 
        by = "margin", permutations = 4999)
adonis2(bray_gene ~ Season + Domestic_animal + Sex + Preservative, 
        data=metadata_gene, 
        by = "margin", permutations = 4999)
adonis2(bray_gene ~ Year + Season + Group + Sex + Preservative, 
        data=metadata_gene, 
        by = "margin", permutations = 4999)

anova(betadisper(bray_gene, group = metadata_gene$Season))
anova(betadisper(bray_gene, group = metadata_gene$Group))
anova(betadisper(bray_gene, group = metadata_gene$Human_food))
anova(betadisper(bray_gene, group = metadata_gene$Domestic_animal))
anova(betadisper(bray_gene, group = metadata_gene$Sex))
anova(betadisper(bray_gene, group = metadata_gene$Preservative))
anova(betadisper(bray_gene, group = metadata_gene$Year))

pairwise.adonis(bray_gene, factors = metadata_gene$Group, 
                perm = 4999, 
                p.adjust.m='holm')

adonis2(jaccard_gene ~ Season + Group + Sex +  Preservative, 
        data=metadata_gene, 
        by = "margin", permutations = 4999)
adonis2(jaccard_gene ~ Season + Human_food + Sex +  Preservative, 
        data=metadata_gene, 
        by = "margin", permutations = 4999)
adonis2(jaccard_gene ~ Season + Domestic_animal + Sex +  Preservative, 
        data=metadata_gene, 
        by = "margin", permutations = 4999)
adonis2(jaccard_gene ~ Year + Season + Group + Sex +  Preservative, 
        data=metadata_gene, 
        by = "margin", permutations = 4999)

anova(betadisper(jaccard_gene, group = metadata_gene$Season))
anova(betadisper(jaccard_gene, group = metadata_gene$Group))
anova(betadisper(jaccard_gene, group = metadata_gene$Human_food))
anova(betadisper(jaccard_gene, group = metadata_gene$Domestic_animal))
anova(betadisper(jaccard_gene, group = metadata_gene$Sex))
anova(betadisper(jaccard_gene, group = metadata_gene$Preservative))
anova(betadisper(jaccard_gene, group = metadata_gene$Year))

pairwise.adonis(jaccard_gene, factors = metadata_gene$Group, 
                perm = 4999, 
                p.adjust.m='holm')

adonis2(bray_pa ~ Season + Group + Sex +  Preservative, 
        data=metadata_pa, 
        by = "margin", permutations = 4999)
adonis2(bray_pa ~ Season + Human_food + Sex +  Preservative, 
        data=metadata_pa, 
        by = "margin", permutations = 4999)
adonis2(bray_pa ~ Season + Domestic_animal + Sex +  Preservative, 
        data=metadata_pa, 
        by = "margin", permutations = 4999)
adonis2(bray_pa ~ Year + Season + Group + Sex +  Preservative, 
        data=metadata_pa, 
        by = "margin", permutations = 4999)

anova(betadisper(bray_pa, group = metadata_pa$Season))
anova(betadisper(bray_pa, group = metadata_pa$Group))
anova(betadisper(bray_pa, group = metadata_pa$Human_food))
anova(betadisper(bray_pa, group = metadata_pa$Domestic_animal))
anova(betadisper(bray_pa, group = metadata_pa$Sex))
anova(betadisper(bray_pa, group = metadata_pa$Preservative))
anova(betadisper(bray_pa, group = metadata_pa$Year))

pairwise.adonis(bray_pa, factors = metadata_pa$Group, 
                perm = 4999, 
                p.adjust.m='holm')

adonis2(jaccard_pa ~ Season + Group + Sex +  Preservative, 
        data=metadata_pa, 
        by = "margin", permutations = 4999)
adonis2(jaccard_pa ~ Season + Human_food + Sex +  Preservative, 
        data=metadata_pa, 
        by = "margin", permutations = 4999)
adonis2(jaccard_pa ~ Season + Domestic_animal + Sex +  Preservative, 
        data=metadata_pa, 
        by = "margin", permutations = 4999)
adonis2(jaccard_pa ~ Year + Season + Group + Sex +  Preservative, 
        data=metadata_pa, 
        by = "margin", permutations = 4999)

anova(betadisper(jaccard_pa, group = metadata_pa$Season))
anova(betadisper(jaccard_pa, group = metadata_pa$Group))
anova(betadisper(jaccard_pa, group = metadata_pa$Human_food))
anova(betadisper(jaccard_pa, group = metadata_pa$Domestic_animal))
anova(betadisper(jaccard_pa, group = metadata_pa$Sex))
anova(betadisper(jaccard_pa, group = metadata_pa$Preservative))
anova(betadisper(jaccard_pa, group = metadata_pa$Year))

pairwise.adonis(jaccard_pa, factors = metadata_pa$Group, 
                perm = 4999, 
                p.adjust.m='holm')

#NMDS plots----
library(ggplot2)
library(ggtext)

set.seed(1018)
mds_otus_bray_gene<-metaMDS(bray_gene, k=2, trymax=499)
mds_otus_bray_gene_points<-mds_otus_bray_gene$points
mds_otus_bray_gene_points2<-merge(x=mds_otus_bray_gene_points, y = metadata_gene, 
                                  by.x = "row.names", by.y = "SampleID")
tiff(file="nmds_plot_bray_genefam_unstrat.tif", res=300, width=8, height=6, units="in")
braygf <- ggplot(mds_otus_bray_gene_points2, aes(x = MDS1, y = MDS2, 
                                                 color = Group, shape = Period)) +
  geom_point(size=3) + scale_color_manual(values = c("#7D3560", "#148F77", "#098BD9", 
                                                     "#97CE2F", "#616161", "#FCB076", "#E784C1")) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
        legend.key=element_blank()) + 
  theme(axis.title.x=element_text(size=rel(1)), 
        axis.title.y=element_text(size=rel(1)),
        plot.title = element_text(size=rel(2)),
        legend.title = element_text(size=rel(1.5)),
        legend.text = element_text(size = rel(1))) + 
  ggtitle("Bray-Curtis:\nGene families") +
  stat_ellipse(aes(x = MDS1, y = MDS2, group = Period, 
                   linetype = Period), 
               type = "t", level = 0.9) + 
  scale_linetype_manual(values = c(1,2)) +
  annotate(geom = "richtext", fill = NA, label.color = NA,
           label = "Period: p = 0.534, 
           R<sup>2</sup> = 2.1%<br>Group: p = 0.434, 
           R<sup>2</sup> = 21.7%<br>Sex: p = 0.339, 
           R<sup>2</sup> = 8.2%<br>Preservative: p = 0.283, 
           R<sup>2</sup> = 4.4%", 
           x = -Inf, y = -Inf, size = 5,
           hjust = 0, vjust = 0)
braygf
dev.off()

mds_otus_jaccard_gene<-metaMDS(jaccard_gene, k=2, trymax=499)
mds_otus_jaccard_gene_points<-mds_otus_jaccard_gene$points
mds_otus_jaccard_gene_points2<-merge(x=mds_otus_jaccard_gene_points, y = metadata_gene, 
                                     by.x = "row.names", by.y = "SampleID")
tiff(file="nmds_plot_jaccard_genefam_unstrat.tif", res=300, width=8, height=6, units="in")
jaccgf <- ggplot(mds_otus_jaccard_gene_points2, aes(x = MDS1, y = MDS2, 
                                                    color = Group, shape = Period)) +
  geom_point(size=3) + scale_color_manual(values = c("#7D3560", "#148F77", "#098BD9", 
                                                     "#97CE2F", "#616161", "#FCB076", "#E784C1")) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
        legend.key=element_blank()) + 
  theme(axis.title.x=element_text(size=rel(1)), 
        axis.title.y=element_text(size=rel(1)),
        plot.title = element_text(size=rel(2)),
        legend.title = element_text(size=rel(1.8)),
        legend.text = element_text(size = rel(1.4))) + 
  ggtitle("Gene families") +
  stat_ellipse(aes(x = MDS1, y = MDS2, group = Period, linetype = Period), 
               type = "t", level = 0.9) + 
  scale_linetype_manual(values = c(1,2)) +
  annotate(geom = "richtext", fill = NA, label.color = NA,
           label = "Period: p = 0.073, 
           R<sup>2</sup> = 7.4%<br><b>Group: p = 0.013, 
           R<sup>2</sup> = 41.0%<br></b>Sex: p = 0.849,
           R<sup>2</sup> = 3.9%<br>Preservative: 
           p = 0.311, R<sup>2</sup> = 3.6%", 
           x = -Inf, y = Inf, size = 5,
           hjust = 0, vjust = 1)
jaccgf
dev.off()

mds_otus_bray_pa<-metaMDS(bray_pa, k=2, trymax=499)
mds_otus_bray_pa_points<-mds_otus_bray_pa$points
mds_otus_bray_pa_points2<-merge(x=mds_otus_bray_pa_points, y = metadata_pa, 
                                by.x = "row.names", by.y = "SampleID")
tiff(file="nmds_plot_bray_pathabund_unstrat.tif", res=300, width=8, height=6, units="in")
braypa <- ggplot(mds_otus_bray_pa_points2, aes(x = MDS1, y = MDS2, 
                                               color = Group, shape = Period)) +
  geom_point(size=3) + scale_color_manual(values = c("#7D3560", "#148F77", "#098BD9", 
                                                     "#97CE2F", "#616161", "#FCB076", "#E784C1")) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
        legend.key=element_blank()) + 
  theme(axis.title.x=element_text(size=rel(2)), 
        axis.title.y=element_text(size=rel(2)),
        plot.title = element_text(size=rel(3)),
        legend.title = element_text(size=rel(2)),
        legend.text = element_text(size = rel(1.8))) + 
  ggtitle("Bray-Curtis:\nPathway abundance") +
  stat_ellipse(aes(x = MDS1, y = MDS2, group = Period, linetype = Period), 
               type = "t", level = 0.9) + 
  scale_linetype_manual(values = c(1,2)) +
  annotate(geom = "richtext", fill = NA, label.color = NA,
           label = "Period: p = 0.442, 
           R<sup>2</sup> = 2.2%<br>Group: p = 0.427, 
           R<sup>2</sup> = 21.0%<br>Sex: p = 0.321, 
           R<sup>2</sup> = 8.1%<br>Preservative: 
           p = 0.222, R<sup>2</sup> = 5.2%", 
           x = -Inf, y = -Inf, size = 5,
           hjust = 0, vjust = 0)
braypa
dev.off()

mds_otus_jaccard_pa<-metaMDS(jaccard_pa, k=2, trymax=499)
mds_otus_jaccard_pa_points<-mds_otus_jaccard_pa$points
mds_otus_jaccard_pa_points2<-merge(x=mds_otus_jaccard_pa_points, y = metadata_pa, 
                                   by.x = "row.names", by.y = "SampleID")
tiff(file="nmds_plot_jaccard_pathabund_unstrat.tif", res=300, width=8, height=6, units="in")
jaccpa <- ggplot(mds_otus_jaccard_pa_points2, aes(x = MDS1, y = MDS2, 
                                                  color = Group, shape = Period)) +
  geom_point(size=3) + scale_color_manual(values = c("#7D3560", "#148F77", "#098BD9",
                                                     "#97CE2F", "#616161", "#FCB076", "#E784C1")) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
        legend.key=element_blank()) + 
  theme(axis.title.x=element_text(size=rel(1)), 
        axis.title.y=element_text(size=rel(1)),
        plot.title = element_text(size=rel(2)),
        legend.title = element_text(size=rel(1.5)),
        legend.text = element_text(size = rel(1))) + 
  ggtitle("Pathway Abundance") +
  stat_ellipse(aes(x = MDS1, y = MDS2, group = Period, linetype = Period), 
               type = "t", level = 0.9) + 
  scale_linetype_manual(values = c(1,2)) +
  annotate(geom = "richtext", fill = NA, label.color = NA,
           label = "Period: p = 0.098, 
           R<sup>2</sup> = 7.1%<br><b>Group: p = 0.030, 
           R<sup>2</sup> = 40.6%<br></b>Sex: p = 0.880,
           R<sup>2</sup> = 3.1%<br>Preservative: 
           p = 0.248, R<sup>2</sup> = 4.2%", 
           x = -Inf, y = Inf, size = 5,
           hjust = 0, vjust = 1)
jaccpa
dev.off()

library(cowplot)

tiff(file = "nmds_plot_combined_shotgun.tif", res = 300, width = 18, height = 14, units="in")
legend1 = get_legend(braygf + theme(legend.box.margin = margin(0, 0, 0, 12)))
col1 = plot_grid(braygf + theme(legend.position = "none"),
                 braypa + theme(legend.position = "none"),
                 nrow = 2, ncol = 1, labels = c('A', 'C'), 
                 label_size = 20)
col2 = plot_grid(jaccgf + theme(legend.position = "none"),
                 jaccpa + theme(legend.position = "none"),
                 nrow = 2, ncol = 1, labels = c('B', 'D'),
                 label_size = 20)
plot_grid(col1, col2, legend1, 
          nrow = 1, ncol = 3, rel_widths = c(1.5, 1.5, 0.75), align = "hv", 
          axis = "t")
dev.off()

tiff(file = "nmds_plot_bray_combined_shotgun.tif", 
     res = 300, width = 18, height = 8, units="in")
legend1 = get_legend(braygf + theme(legend.box.margin = margin(0, 0, 0, 12)))
col1 = plot_grid(braygf + theme(legend.position = "none"),
                 nrow = 1, ncol = 1, labels = c('A'), 
                 label_size = 20)
col2 = plot_grid(braypa + theme(legend.position = "none"),
                 nrow = 1, ncol = 1, labels = c('B'),
                 label_size = 20)
plot_grid(col1, col2, legend1, 
          nrow = 1, ncol = 3, rel_widths = c(1.5, 1.5, 0.75), 
          align = "hv", axis = "t")
dev.off()

#Alpha diversity models ----

library(tidyverse)
library(nlme)
library(multcomp)
library(car)

shannon_gene = read.table("shannon_genefam_unstrat.tsv", header = T)
evenness_gene = read.table("evenness_genefam_unstrat.tsv", header = T)
obs_otus_gene = read.table("observed_otus_genefam_unstrat.tsv", header = T)

shannon_pa = read.table("shannon_pathabund_unstrat.tsv", header = T)
evenness_pa = read.table("evenness_pathabund_unstrat.tsv", header = T)
obs_otus_pa = read.table("observed_otus_pathabund_unstrat.tsv", header = T)

alpha_gene = inner_join(metadata_gene, shannon_gene, by = "SampleID") %>% 
  inner_join(evenness_gene, by = "SampleID") %>% 
  inner_join(obs_otus_gene, by = "SampleID")
alpha_gene$Group = as.factor(alpha_gene$Group)

alpha_pa = inner_join(metadata_pa, shannon_pa, by = "SampleID") %>% 
  inner_join(evenness_pa, by = "SampleID") %>% 
  inner_join(obs_otus_pa, by = "SampleID")

sg <- lm(shannon_entropy ~ Season + Group + Sex + Preservative, data = alpha_gene)
summary(sg)
Anova(sg)
summary(glht(sg,linfct=mcp(Season="Tukey")))

sgh <- lm(shannon_entropy ~ Season + Human_food + Sex + Preservative, data = alpha_gene)
summary(sgh)
Anova(sgh)

sga <- lm(shannon_entropy ~ Season + Domestic_animal + Sex + Preservative, data = alpha_gene)
summary(sga)
Anova(sga)

sgy <- lm(shannon_entropy ~ Year + Season + Group + Sex + Preservative, data = alpha_gene)
summary(sgy)
Anova(sgy)
summary(glht(sgysg,linfct=mcp(Season="Tukey")))

eg <- lm(pielou_evenness ~ Season + Group + Sex + Preservative, data = alpha_gene)
summary(eg)
Anova(eg)
summary(glht(eg,linfct=mcp(Season="Tukey")))

egh <- lm(pielou_evenness ~ Season + Human_food + Sex + Preservative, data = alpha_gene)
summary(egh)
Anova(egh)

ega <- lm(pielou_evenness ~ Season + Domestic_animal + Sex + Preservative, data = alpha_gene)
summary(ega)
Anova(ega)

egy <- lm(pielou_evenness ~ Year + Season + Group + Sex + Preservative, data = alpha_gene)
summary(egy)
Anova(egy)

og <- lm(observed_features ~ Season + Group + Sex + Preservative, data = alpha_gene)
summary(og)
Anova(og)
summary(glht(og,linfct=mcp(Group="Tukey")))

ogh <- lm(observed_features ~ Season + Human_food + Sex + Preservative, data = alpha_gene)
summary(ogh)
Anova(ogh)

oga <- lm(observed_features ~ Season + Domestic_animal + Sex + Preservative, data = alpha_gene)
summary(oga)
Anova(oga)

ogy <- lm(observed_features ~ Year + Season + Group + Sex + Preservative, data = alpha_gene)
summary(ogy)
Anova(ogy)

sp <- lm(shannon_entropy ~ Season + Group + Sex + Preservative, data = alpha_pa)
summary(sp)
Anova(sp)
summary(glht(sp,linfct=mcp(Season="Tukey")))

sph <- lm(shannon_entropy ~ Season + Human_food + Sex + Preservative, data = alpha_pa)
summary(sph)
Anova(sph)

spa <- lm(shannon_entropy ~ Season + Domestic_animal + Sex + Preservative, data = alpha_pa)
summary(spa)
Anova(spa)

spy <- lm(shannon_entropy ~ Year + Season + Group + Sex + Preservative, data = alpha_pa)
summary(spy)
Anova(spy)

ep <- lm(pielou_evenness ~ Season + Group + Sex + Preservative, data = alpha_pa)
summary(ep)
Anova(ep)
summary(glht(ep,linfct=mcp(Season="Tukey")))

eph <- lm(pielou_evenness ~ Season + Human_food + Sex + Preservative, data = alpha_pa)
summary(eph)
Anova(eph)

epa <- lm(pielou_evenness ~ Season + Domestic_animal + Sex + Preservative, data = alpha_pa)
summary(epa)
Anova(epa)

epy <- lm(pielou_evenness ~ Year + Season + Group + Sex + Preservative, data = alpha_pa)
summary(epy)
Anova(epy)

op <- lm(observed_features ~ Season + Group + Sex + Preservative, data = alpha_pa)
summary(op)
Anova(op)
summary(glht(op,linfct=mcp(Season="Tukey")))

oph <- lm(observed_features ~ Season + Human_food + Sex + Preservative, data = alpha_pa)
summary(oph)
Anova(oph)

opa <- lm(observed_features ~ Season + Domestic_animal + Sex + Preservative, data = alpha_pa)
summary(opa)
Anova(opa)
summary(glht(opa,linfct=mcp(Season="Tukey")))

opy <- lm(observed_features ~ Year + Season + Group + Sex + Preservative, data = alpha_pa)
summary(opy)
Anova(opy)

#Diversity plots ----
library(ggpubr)

plot1 = ggboxplot(alpha_gene, x = "Group", 
                  y = "observed_features", color = "Group", 
                  palette = c("#7D3560", "#148F77", "#098BD9",
                              "#97CE2F", "#616161", "#FCB076", "#E784C1"), 
                  add = "jitter", 
                  add.params = list(fill = "white"), 
                  ylab = "Observed Gene Families") 
plot1 = ggpar(plot1, legend = "right", font.y = 16,
              font.legend = 16, font.ytickslab = 14) + rremove("xlab") + 
  rremove("x.text") + rremove("x.ticks") + 
  rremove("legend.title") + rremove("legend") + 
  stat_compare_means(label = "p.signif", 
                     comparisons = list(c("Princess", "Key")),
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.05, 0.5, 1),
                                        symbols = c("*", "*", "*", "*", "ns")))


legend1 = get_legend(jaccgf + theme(legend.box.margin = margin(0, 0, 0, 12)))
col1 = plot_grid(jaccgf + theme(legend.position = "none"),
                 nrow = 1, ncol = 1)
col2 = plot_grid(jaccpa + theme(legend.position = "none"),
                 nrow = 1, ncol = 1)

tiff(file="div_combined.tif", res=300, width=20, height=6, units="in")
shot_div = plot_grid(plot1, col1, col2, legend1, 
          nrow = 1, ncol = 4, rel_widths = c(1.5, 1.5, 1.5, 0.75), align = "hv", 
          axis = "t", labels = c("A", "B", "", ""), label_size = 28)
shot_div
dev.off()

setEPS()
postscript(file="div_combined.eps", width=20, height=6, paper = "special")
shot_div
dev.off()


#Pathway abundance GLMMs----

library(glmmTMB)
library(multcomp)
library(car)
library(tidyverse)
library(fdrtool)

glht_glmmTMB <- function (model, ..., component="cond") {
  glht(model, ...,
       coef. = function(x) fixef(x)[[component]],
       vcov. = function(x) vcov(x)[[component]],
       df = NULL)
}
modelparm.glmmTMB <- function (model, coef. = function(x) fixef(x)[[component]],
                               vcov. = function(x) vcov(x)[[component]],
                               df = NULL, component="cond", ...) {
  multcomp:::modelparm.default(model, coef. = coef., vcov. = vcov.,
                               df = df, ...)
}

path_abund = read.csv("caatinga_pathabundance_relab_unstratified.csv", 
                      header = T)
metadata_pa = read.csv("caatinga_shotgun_metadata.csv", header = T)

pa = metadata_pa %>% inner_join(path_abund, by = "SampleID")

season_estimatematrix = mat.or.vec(459,2)
season_pvaluematrix = mat.or.vec(459,2)
group_estimatematrix = mat.or.vec(459,2)
group_pvaluematrix = mat.or.vec(459,2)
sex_estimatematrix = mat.or.vec(459,2)
sex_pvaluematrix = mat.or.vec(459,2)
preserv_estimatematrix = mat.or.vec(459,2)
preserv_pvaluematrix = mat.or.vec(459,2)

for(i in 17:479) {
  variable = pa[,i]
  b = try(glmmTMB(variable ~ Season + Group + Sex + Preservative, 
                  data = pa, family = nbinom2))
  anova = Anova(b)
  season_estimatematrix[i-16,2] = anova[1,1]
  season_pvaluematrix[i-16,2] = anova[1,3]
  group_estimatematrix[i-16,2] = anova[2,1]
  group_pvaluematrix[i-16,2] = anova[2,3]
  sex_estimatematrix[i-16,2] = anova[3,1]
  sex_pvaluematrix[i-16,2] = anova[3,3]
  preserv_estimatematrix[i-16,2] = anova[4,1]
  preserv_pvaluematrix[i-16,2] = anova[4,3]
  season_estimatematrix[i-16,1] = names(pa)[i]
  season_pvaluematrix[i-16,1] = names(pa)[i]
  group_estimatematrix[i-16,1] = names(pa)[i]
  group_pvaluematrix[i-16,1] = names(pa)[i]
  sex_estimatematrix[i-16,1] = names(pa)[i]
  sex_pvaluematrix[i-16,1] = names(pa)[i]
  preserv_estimatematrix[i-16,1] = names(pa)[i]
  preserv_pvaluematrix[i-16,1] = names(pa)[i]
}

pa_season = bind_cols(as.data.frame(season_estimatematrix[,1:2]), 
                     as.data.frame(season_pvaluematrix[,2]))
write.csv(pa_season, "pa_season.csv")
pa_season_nona = filter(pa_season, pa_season[,3] != "NaN") 
pa_season_nona[,3] = as.numeric(as.character(pa_season_nona[,3]))
pa_season_corrected = bind_cols(pa_season_nona, 
                               as.data.frame(fdrtool(pa_season_nona[,3], 
                                                     statistic = "pvalue", plot = F)))
write.csv(pa_season_corrected, "pa_season_corrected.csv")

pa_group = bind_cols(as.data.frame(group_estimatematrix[,1:2]), 
                      as.data.frame(group_pvaluematrix[,2]))
write.csv(pa_group, "pa_group.csv")
pa_group_nona = filter(pa_group, pa_group[,3] != "NaN") 
pa_group_nona[,3] = as.numeric(as.character(pa_group_nona[,3]))
pa_group_corrected = bind_cols(pa_group_nona, 
                                as.data.frame(fdrtool(pa_group_nona[,3], 
                                                      statistic = "pvalue", plot = F)))
write.csv(pa_group_corrected, "pa_group_corrected.csv")

pa_sex = bind_cols(as.data.frame(sex_estimatematrix[,1:2]), 
                     as.data.frame(sex_pvaluematrix[,2]))
write.csv(pa_sex, "pa_sex.csv")
pa_sex_nona = filter(pa_sex, pa_sex[,3] != "NaN") 
pa_sex_nona[,3] = as.numeric(as.character(pa_sex_nona[,3]))
pa_sex_corrected = bind_cols(pa_sex_nona, 
                               as.data.frame(fdrtool(pa_sex_nona[,3], 
                                                     statistic = "pvalue", plot = F)))
write.csv(pa_sex_corrected, "pa_sex_corrected.csv")


pa_preserv = bind_cols(as.data.frame(preserv_estimatematrix[,1:2]), 
                      as.data.frame(preserv_pvaluematrix[,2]))
write.csv(pa_preserv, "pa_preserv.csv")
pa_preserv_nona = filter(pa_preserv, pa_preserv[,3] != "NaN") 
pa_preserv_nona[,3] = as.numeric(as.character(pa_preserv_nona[,3]))
pa_preserv_corrected = bind_cols(pa_preserv_nona, 
                                as.data.frame(fdrtool(pa_preserv_nona[,3], 
                                                      statistic = "pvalue", plot = F)))
write.csv(pa_preserv_corrected, "pa_preserv_corrected.csv")

a = glmmTMB(X1CMET2.PWY..N10.formyl.tetrahydrofolate.biosynthesis ~ Season + Group + Preservative, 
            data = pa, family = nbinom2)
summary(a)
Anova(a)
summary(glht(a,linfct=mcp(Season="Tukey")))

#Gene family GLMMs----

library(glmmTMB)
library(multcomp)
library(car)
library(tidyverse)
library(fdrtool)

glht_glmmTMB <- function (model, ..., component="cond") {
  glht(model, ...,
       coef. = function(x) fixef(x)[[component]],
       vcov. = function(x) vcov(x)[[component]],
       df = NULL)
}
modelparm.glmmTMB <- function (model, coef. = function(x) fixef(x)[[component]],
                               vcov. = function(x) vcov(x)[[component]],
                               df = NULL, component="cond", ...) {
  multcomp:::modelparm.default(model, coef. = coef., vcov. = vcov.,
                               df = df, ...)
}

gene_fam = read.csv("caatinga_genefamilies_ko_cpm_unstratified.csv", 
                      header = T)
metadata_gf = read.csv("caatinga_shotgun_gf_metadata.csv", header = T)

gf = metadata_gf %>% inner_join(gene_fam, by = "SampleID")

season_estimatematrix = mat.or.vec(5522,2)
season_pvaluematrix = mat.or.vec(5522,2)
group_estimatematrix = mat.or.vec(5522,2)
group_pvaluematrix = mat.or.vec(5522,2)
sex_estimatematrix = mat.or.vec(5522,2)
sex_pvaluematrix = mat.or.vec(5522,2)
preserv_estimatematrix = mat.or.vec(5522,2)
preserv_pvaluematrix = mat.or.vec(5522,2)

for(i in 17:5538) {
  variable = gf[,i]
  b = try(glmmTMB(variable ~ Season + Group + Sex + Preservative, 
                  data = gf, family = nbinom2))
  anova = Anova(b)
  season_estimatematrix[i-16,2] = anova[1,1]
  season_pvaluematrix[i-16,2] = anova[1,3]
  group_estimatematrix[i-16,2] = anova[2,1]
  group_pvaluematrix[i-16,2] = anova[2,3]
  sex_estimatematrix[i-16,2] = anova[3,1]
  sex_pvaluematrix[i-16,2] = anova[3,3]
  preserv_estimatematrix[i-16,2] = anova[4,1]
  preserv_pvaluematrix[i-16,2] = anova[4,3]
  season_estimatematrix[i-16,1] = names(gf)[i]
  season_pvaluematrix[i-16,1] = names(gf)[i]
  group_estimatematrix[i-16,1] = names(gf)[i]
  group_pvaluematrix[i-16,1] = names(gf)[i]
  sex_estimatematrix[i-16,1] = names(gf)[i]
  sex_pvaluematrix[i-16,1] = names(gf)[i]
  preserv_estimatematrix[i-16,1] = names(gf)[i]
  preserv_pvaluematrix[i-16,1] = names(gf)[i]
}

gf_season = bind_cols(as.data.frame(season_estimatematrix[,1:2]), 
                      as.data.frame(season_pvaluematrix[,2]))
write.csv(gf_season, "gf_season.csv")
gf_season_nona = filter(gf_season, gf_season[,3] != "NaN") 
gf_season_nona[,3] = as.numeric(as.character(gf_season_nona[,3]))
gf_season_corrected = bind_cols(gf_season_nona, 
                                as.data.frame(fdrtool(gf_season_nona[,3], 
                                                      statistic = "pvalue", plot = F)))
write.csv(gf_season_corrected, "gf_season_corrected.csv")

gf_group = bind_cols(as.data.frame(group_estimatematrix[,1:2]), 
                     as.data.frame(group_pvaluematrix[,2]))
write.csv(gf_group, "gf_group.csv")
gf_group_nona = filter(gf_group, gf_group[,3] != "NaN") 
gf_group_nona[,3] = as.numeric(as.character(gf_group_nona[,3]))
gf_group_corrected = bind_cols(gf_group_nona, 
                               as.data.frame(fdrtool(gf_group_nona[,3], 
                                                     statistic = "pvalue", plot = F)))
write.csv(gf_group_corrected, "gf_group_corrected.csv")

gf_sex = bind_cols(as.data.frame(sex_estimatematrix[,1:2]), 
                   as.data.frame(sex_pvaluematrix[,2]))
write.csv(gf_sex, "gf_sex.csv")
gf_sex_nona = filter(gf_sex, gf_sex[,3] != "NaN") 
gf_sex_nona[,3] = as.numeric(as.character(gf_sex_nona[,3]))
gf_sex_corrected = bind_cols(gf_sex_nona, 
                             as.data.frame(fdrtool(gf_sex_nona[,3], 
                                                   statistic = "pvalue", plot = F)))
write.csv(gf_sex_corrected, "gf_sex_corrected.csv")


gf_preserv = bind_cols(as.data.frame(preserv_estimatematrix[,1:2]), 
                       as.data.frame(preserv_pvaluematrix[,2]))
write.csv(gf_preserv, "gf_preserv.csv")
gf_preserv_nona = filter(gf_preserv, gf_preserv[,3] != "NaN") 
gf_preserv_nona[,3] = as.numeric(as.character(gf_preserv_nona[,3]))
gf_preserv_corrected = bind_cols(gf_preserv_nona, 
                                 as.data.frame(fdrtool(gf_preserv_nona[,3], 
                                                       statistic = "pvalue", plot = F)))
write.csv(gf_preserv_corrected, "gf_preserv_corrected.csv")

#Taxonomy permanovas ----

setwd("/Users/elizabethmallott/Dropbox/Projects/Gut_microbiome/Caatinga_marmosets/shotgun")

library(vegan)

species = read.csv("merge_abundance_table_species_t.csv", header = T)
metadata = read.csv("caatinga_shotgun_taxa_metadata.csv", header = T)

bray_species = vegdist(species[2:24], method = "bray")
jaccard_species = vegdist(species[2:24], method = "jaccard")

adonis2(bray_species ~ Season + Group + Preservative, data=metadata, 
        by = "margin", permutations = 5000)

adonis2(jaccard_species ~ Season + Group + Preservative, data=metadata, 
        by = "margin", permutations = 5000)
