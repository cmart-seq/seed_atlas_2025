suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.At.tair.db))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(RColorBrewer))
library(harmony)
set.seed(123)

setwd("../04_manual_annotation/review_embryo_subclustering")

######## reading in data ######## 
#loading all of the timepoint annotated seurats
dap3 <- readRDS("../../04_manual_annotation/outputs/DAP3_wcze_subs/DAP3_wcze_subs_annotated.rds") 
DefaultAssay(dap3) <- "RNA"
dap3 <- JoinLayers(dap3)
dap5 <- readRDS("../../04_manual_annotation/outputs/DAP5_wcze_subs/DAP5_wcze_subs_annotated.rds")
DefaultAssay(dap5) <- "RNA"
dap5 <- JoinLayers(dap5)
dap7 <- readRDS("../../04_manual_annotation/outputs/DAP7/DAP7_annotated.rds")
DefaultAssay(dap7) <- "RNA"
dap7 <- JoinLayers(dap7)

#loading the integrated, module-scored atlas object
atlas <- readRDS("../../06_atlas_merging/outputs/ATLAS_integrated_annotated_rpca.rds")  
DefaultAssay(atlas) <- "RNA"
atlas <- JoinLayers(atlas)

#loading the embryo-only data, reading in the pseudotimed object: 
embryo <- readRDS("../../05_across_timepoints/01_merging/EMB/EMB_tp_integrated_cc_reg_harmony_30PCs.rds")
DefaultAssay(embryo) <- "RNA"
embryo <- JoinLayers(embryo)

#Finally, reading in the DE for the atlas L2 annotations, not timed. This will be presented next to the proportionality 
atlas_l2_markers <- read_csv("../../07_differential_expression/outputs/ATLAS/ATLAS_level_2_annotation_markers.csv")

######## setting color keys ######## 
#adding the colors for the level 1 and 2 clusters

color_key1 = c("Embryo" = "#57A15D", 
               "Endosperm" = "#D44A90", 
               "Seed coat" = "#96ADD0", 
               "Funiculus" = "#f1c40f", 
               "Ovule" = "#a6acaf")

level_1_order <- c("Embryo",
                   "Endosperm",
                   "Seed coat", 
                   "Funiculus", 
                   "Ovule")

color_key2 = c("EMB" = "#57A15D", 
               "PEN" = "#ceaab5", 
               "MCE" = "#97576b", 
               "CZE" = "#ce0d48",
               "FUN" = "#f1c40f",
               "OVL" = "#a6acaf", 
               "ii1" = "#aed6f1", 
               "ii1'" = "#5dade2",
               "ii1'/ii2"= "#5dade2",
               "ii2" = "#21618c",
               "oi1" = "#bb8fce", 
               "oi2" = "#6c3483", 
               "CZSC" = "#088F8F",
               "CPT" = "#044747")

level_2_order = c("EMB",
                  "PEN",
                  "MCE",
                  "CZE",
                  "FUN",
                  "OVL",
                  "ii1",
                  "ii1'",
                  "ii1'/ii2",
                  "ii2",
                  "oi1",
                  "oi2",
                  "CZSC",
                  "CPT")


#subclustering individual embryos
#plotting nodine markers, making module scores
nodine <- read_csv("inputs/nodine_markers.csv")
#filetering the nodine lab markers for specific groups
#"s" means means that the gene strongly labels this subtype
upper_protoderm <- nodine[nodine$`upper protoderm?` == "s", ]
shoot_meristem_initials <- nodine[nodine$`shoot meristem initials?`== "s", ]
vascular_initials <- nodine[nodine$`vascular initials?`== "s", ]
quiescent_center_initials <- nodine[nodine$`quiescent center initials?`== "s", ]
suspensor <- nodine[nodine$`suspensor`== "s", ]
lower_protoderm <- nodine[nodine$`lower protoderm?`== "s", ]
upper_inner_periphery <- nodine[nodine$`upper inner periphery?`== "s", ]
ground_tissue_initials <- nodine[nodine$`ground tissue initials?`== "s", ]
columella_initials <- nodine[nodine$`columella initials?`== "s", ]

#performing enrichment analysis
atlas <- AddModuleScore(
  object = atlas,
  features = list(aux_iaa = upper_protoderm$`TAIR ID`),
  name = 'upper_protoderm'
)

atlas <- AddModuleScore(
  object = atlas,
  features = list(aux_iaa = shoot_meristem_initials$`TAIR ID`),
  name = 'shoot_meristem_initials'
)

atlas <- AddModuleScore(
  object = atlas,
  features = list(aux_iaa = vascular_initials$`TAIR ID`),
  name = 'vascular_initials'
)

atlas <- AddModuleScore(
  object = atlas,
  features = list(aux_iaa = quiescent_center_initials$`TAIR ID`),
  name = 'quiescent_center_initials'
)

atlas <- AddModuleScore(
  object = atlas,
  features = list(aux_iaa = suspensor$`TAIR ID`),
  name = 'suspensor'
)

atlas <- AddModuleScore(
  object = atlas,
  features = list(aux_iaa = lower_protoderm$`TAIR ID`),
  name = 'lower_protoderm'
)

atlas <- AddModuleScore(
  object = atlas,
  features = list(aux_iaa = upper_inner_periphery$`TAIR ID`),
  name = 'upper_inner_periphery'
)

atlas <- AddModuleScore(
  object = atlas,
  features = list(aux_iaa = ground_tissue_initials$`TAIR ID`),
  name = 'ground_tissue_initials'
)

atlas <- AddModuleScore(
  object = atlas,
  features = list(aux_iaa = columella_initials$`TAIR ID`),
  name = 'columella_initials'
)

#with the analysis performed at the atlas level, subsetting to embryo datasets
embryo <- subset(atlas, subset = level_1_annotation== "Embryo" )
dap3_embryo <- subset(embryo, subset = timepoint== "3DAP" )
dap5_embryo <- subset(embryo, subset = timepoint== "5DAP")
dap7_embryo <- subset(embryo, subset = timepoint== "7DAP")

#clustering to isolate upper/lower protoderm
dap3_embryo <- FindVariableFeatures(dap3_embryo, selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = VariableFeatures(dap3_embryo), vars.to.regress = c("S.Score", "G2M.Score")) %>%
  RunPCA(npcs=20) %>%
  RunUMAP(dims = 1:20) %>%
  FindNeighbors(dims = 1:20, verbose = F) %>%
  FindClusters(resolution = 1.2, cluster.name = "ccreg_clusters") 
png(paste0("outputs/dap3_embryo_subclustered.png"), width = 5, height = 5, units = "in",res = 240, bg = "transparent")
DimPlot(dap3_embryo)
dev.off()

#clustering to isolate upper/lower protoderm
dap5_embryo  <- FindVariableFeatures(dap5_embryo, selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = VariableFeatures(dap5_embryo), vars.to.regress = c("S.Score", "G2M.Score")) %>%
  RunPCA(npcs=20) %>%
  RunHarmony(reduction = "pca", group.by.vars = "bio_rep", 
             dims.use = 1:20, reduction.save = 'harmony', plot_convergence = F) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 1.5, cluster.name = "ccreg_clusters") %>% 
  RunUMAP(reduction = "harmony", dims = 1:20, reduction.name = "umap.harmony")

png(paste0("outputs/dap5_embryo_subclustered.png"), width = 5, height = 5, units = "in",res = 240, bg = "transparent")
DimPlot(dap5_embryo, reduction = "umap.harmony")
dev.off()

#now plotting enrichment of these modscores in the new clusters
all_modscores <- c("upper_protoderm1", "shoot_meristem_initials1", "vascular_initials1", "quiescent_center_initials1", "suspensor1", 
                   "lower_protoderm1", "upper_inner_periphery1", "ground_tissue_initials1", "columella_initials1")
all_embryos <- list("DAP3" = dap3_embryo, "DAP5" = dap5_embryo)

#ID-ing the suspensor: those in the 90th percentile of modscore enrichment
sus_cells <-WhichCells(dap3_embryo, expression = suspensor1 > quantile(dap3_embryo$suspensor1, 0.9))
png(paste0("outputs/dap3_embryo_subclustered_for_sus.png"), width = 5, height = 5, units = "in",res = 240, bg = "transparent")
DimPlot(dap3_embryo, cells.highlight = sus_cells, reduction = "umap")
dev.off()

#transferring the annotations to the embryo object
dap3_embryo$ccreg_clusters_old <- as.character(dap3_embryo$ccreg_clusters)
dap3_embryo$ccreg_clusters <- as.character(dap3_embryo$ccreg_clusters)
dap3_embryo$ccreg_clusters[colnames(dap3_embryo) %in% sus_cells] <- "5"

#renaming the embryo subclusters, DAP3
embryo_annotations <- read_csv("inputs/Embryo_subclusters.csv")[,1:6]
colnames(embryo_annotations) <- c("ccreg_clusters", "level_1_annotation",
                                  "level_2_annotation", "level_3_annotation_abbr",
                                  "level_3_annotation_full", "timepoint")
embryo_anno_dap3 <- filter(embryo_annotations, timepoint == "DAP3")
embryo_anno_dap5 <- filter(embryo_annotations, timepoint == "DAP5")

#making all new plotting IDs for L3 in DAP3
#making a "plotting ID" for the level 3 annotations
plotting_ID3_table <- embryo_anno_dap3 %>% dplyr::select(c(level_1_annotation, level_2_annotation, level_3_annotation_abbr, level_3_annotation_full)) %>% 
  group_by(level_1_annotation) %>% 
  arrange(level_2_annotation) %>% 
  arrange(level_1_annotation) %>% 
  ungroup() %>%
  unique()

plotting_ID3_table$level_3_plotting <- seq(1, length(unique(plotting_ID3_table$level_3_annotation_full)), 1)
plotting_ID3_table <- dplyr::select(plotting_ID3_table, c(level_3_annotation_full, level_3_plotting))

#merging and saving with the old metadata
embryo_anno_dap3 <- left_join(embryo_anno_dap3, plotting_ID3_table, by ="level_3_annotation_full")

#now, updating the rest of the 3 DAP annotations
#full level 3
level_3_annotation_vector_full <- embryo_anno_dap3$level_3_annotation_full
names(level_3_annotation_vector_full) <- embryo_anno_dap3$ccreg_clusters
#full level 3 updating IDs, adding full names
dap3_embryo <- SetIdent(dap3_embryo, value = "ccreg_clusters")
dap3_embryo <- RenameIdents(object = dap3_embryo, level_3_annotation_vector_full)
dap3_embryo$level_3_annotation_full <- Idents(dap3_embryo)

#abbr level 3
level_3_annotation_vector_abbr <- embryo_anno_dap3$level_3_annotation_abbr
names(level_3_annotation_vector_abbr) <- embryo_anno_dap3$level_3_annotation_full
#abbr level 3 updating IDs, adding full names
dap3_embryo <- RenameIdents(object = dap3_embryo, level_3_annotation_vector_abbr)
dap3_embryo$level_3_annotation_abbr <- Idents(dap3_embryo)

#plotting level 3
level_3_annotation_vector_plotting <- embryo_anno_dap3$level_3_plotting
names(level_3_annotation_vector_plotting) <- embryo_anno_dap3$level_3_annotation_full
#abbr level 3 updating IDs, adding full names
dap3_embryo <- SetIdent(dap3_embryo, value = "level_3_annotation_full")
dap3_embryo <- RenameIdents(object = dap3_embryo, level_3_annotation_vector_plotting)
dap3_embryo$level_3_plotting <- Idents(dap3_embryo)

#adding the plotting number to the level 2 annotations
dap3_embryo[[]] <-mutate(dap3_embryo[[]], level_3_annotation_full = paste0(level_3_plotting, ". ",level_3_annotation_full) )
dap3_embryo[[]] <-mutate(dap3_embryo[[]], level_3_annotation_abbr = paste0(level_3_plotting, ". ",level_3_annotation_abbr) )

#factoring for the right order, level 3
level_3_factor_full <- dap3_embryo[[]] %>% dplyr::select(c(level_3_annotation_full, level_3_plotting))
level_3_factor_full$level_3_plotting <- as.numeric(as.character(level_3_factor_full$level_3_plotting))
level_3_factor_full <-   level_3_factor_full %>%
  arrange(level_3_plotting) %>% 
  pull(level_3_annotation_full) %>% unique()

dap3_embryo$level_3_annotation_full <- factor(dap3_embryo$level_3_annotation_full, levels = level_3_factor_full)

level_3_factor_abbr <- dap3_embryo[[]] %>% dplyr::select(c(level_3_annotation_abbr, level_3_plotting))%>% arrange(level_3_plotting) %>% pull(level_3_annotation_abbr) %>% unique()
dap3_embryo$level_3_annotation_abbr <- factor(dap3_embryo$level_3_annotation_abbr, levels = level_3_factor_abbr)

#plotting and saving
#new 
pdf("outputs/dap3_embryo_subclusters_annotated.pdf")
DimPlot(dap3_embryo, group.by = "level_3_annotation_full")
DimPlot(dap3_embryo, group.by = "level_3_annotation_abbr")
DimPlot(dap3_embryo, group.by = "level_3_plotting")
dev.off()

#renaming the embryo subclusters, DAP5

#making a "plotting ID" for the level 3 annotations
plotting_ID3_table <- embryo_anno_dap5 %>% dplyr::select(c(level_1_annotation, level_2_annotation, level_3_annotation_abbr, level_3_annotation_full)) %>% 
  group_by(level_1_annotation) %>% 
  arrange(level_2_annotation) %>% 
  arrange(level_1_annotation) %>% 
  ungroup() %>%
  unique()

plotting_ID3_table$level_3_plotting <- seq(1, length(unique(plotting_ID3_table$level_3_annotation_full)), 1)
plotting_ID3_table <- dplyr::select(plotting_ID3_table, c(level_3_annotation_full, level_3_plotting))

#merging and saving with the old metadata
embryo_anno_dap5 <- left_join(embryo_anno_dap5, plotting_ID3_table, by ="level_3_annotation_full")

#new, marker based IDs
#full level 3
level_3_annotation_vector_full <- embryo_anno_dap5$level_3_annotation_full
names(level_3_annotation_vector_full) <- embryo_anno_dap5$ccreg_clusters
#full level 3 updating IDs, adding full names
dap5_embryo <- SetIdent(dap5_embryo, value = "ccreg_clusters")
dap5_embryo <- RenameIdents(object = dap5_embryo, level_3_annotation_vector_full)
dap5_embryo$level_3_annotation_full <- Idents(dap5_embryo)

#abbr level 3
level_3_annotation_vector_abbr <- embryo_anno_dap5$level_3_annotation_abbr
names(level_3_annotation_vector_abbr) <- embryo_anno_dap5$level_3_annotation_full
#abbr level 3 updating IDs, adding full names
dap5_embryo <- RenameIdents(object = dap5_embryo, level_3_annotation_vector_abbr)
dap5_embryo$level_3_annotation_abbr <- Idents(dap5_embryo)

#plotting level 3
level_3_annotation_vector_plotting <- embryo_anno_dap5$level_3_plotting
names(level_3_annotation_vector_plotting) <- embryo_anno_dap5$level_3_annotation_full
#abbr level 3 updating IDs, adding full names
dap5_embryo <- SetIdent(dap5_embryo, value = "level_3_annotation_full")
dap5_embryo <- RenameIdents(object = dap5_embryo, level_3_annotation_vector_plotting)
dap5_embryo$level_3_plotting <- Idents(dap5_embryo)

#adding the plotting number to the level 2 annotations
dap5_embryo[[]] <-mutate(dap5_embryo[[]], level_3_annotation_full = paste0(level_3_plotting, ". ",level_3_annotation_full) )
dap5_embryo[[]] <-mutate(dap5_embryo[[]], level_3_annotation_abbr = paste0(level_3_plotting, ". ",level_3_annotation_abbr) )

#factoring for the right order, level 3
level_3_factor_full <- dap5_embryo[[]] %>% dplyr::select(c(level_3_annotation_full, level_3_plotting))
level_3_factor_full$level_3_plotting <- as.numeric(as.character(level_3_factor_full$level_3_plotting))
level_3_factor_full <-   level_3_factor_full %>%
  arrange(level_3_plotting) %>% 
  pull(level_3_annotation_full) %>% unique()

dap5_embryo$level_3_annotation_full <- factor(dap5_embryo$level_3_annotation_full, levels = level_3_factor_full)

level_3_factor_abbr <- dap5_embryo[[]] %>% dplyr::select(c(level_3_annotation_abbr, level_3_plotting))%>% arrange(level_3_plotting) %>% pull(level_3_annotation_abbr) %>% unique()
dap5_embryo$level_3_annotation_abbr <- factor(dap5_embryo$level_3_annotation_abbr, levels = level_3_factor_abbr)

#plotting and saving
#new 
pdf("outputs/dap5_embryo_subclusters_annotated.pdf")
DimPlot(dap5_embryo, group.by = "level_3_annotation_full", reduction = "umap.harmony")
DimPlot(dap5_embryo, group.by = "level_3_annotation_abbr", reduction = "umap.harmony")
DimPlot(dap5_embryo, group.by = "level_3_plotting", reduction = "umap.harmony")
dev.off()

#saving objects
saveRDS(dap3_embryo, "outputs/dap3_embryo_subclustered.rds")
saveRDS(dap5_embryo, "outputs/dap5_embryo_subclustered.rds")








