#This script takes all of the timepoint-annotated datasets and merges them into a single atlas dataset
#example usage: 
#Rscript atlas_merging_rpca.R --wd . --dataset ATLAS --n_VarGenes 4000

library(Seurat)
library(stringr)
library(dplyr)
library(argparse)
library(scCustomize)
library(STACAS)
library(ggplot2)
library(grid)

#reproducible randomness
set.seed(123)

#getting the args
parser <- ArgumentParser()

parser$add_argument("--wd", help = "working directory") 
parser$add_argument("--dataset", help = "sample name, of the form TIMEPOINT (e.g. DAP5)") 
parser$add_argument("--n_VarGenes", help = "number of variable genes to detect", type = "double")

# Parse the arguments
args <- parser$parse_args()

#setting the directory
working = args$wd
setwd(working)

#getting the dataset name
dataname = args$dataset

#making the output directory if it doesnt exist
if (!dir.exists("outputs")) {
  dir.create("outputs", recursive = TRUE)
}


if (!dir.exists(file.path("outputs/figures"))) {
  dir.create(file.path("outputs/figures"), recursive = TRUE)
}

outpath = file.path(args$wd, "outputs")

#loading all of the seurats
dap3 <- readRDS("../04_manual_annotation/outputs/DAP3_wcze_subs/DAP3_wcze_subs_annotated.rds")
DefaultAssay(dap3) <- "RNA"
dap3 <- JoinLayers(dap3)
dap5 <- readRDS("../04_manual_annotation/outputs/DAP5_wcze_subs/DAP5_wcze_subs_annotated.rds")
DefaultAssay(dap5) <- "RNA"
dap5 <- JoinLayers(dap5)
dap7 <- readRDS("../04_manual_annotation/outputs/DAP7/DAP7_annotated.rds")
DefaultAssay(dap7) <- "RNA"
dap7 <- JoinLayers(dap7)

#labeling the ii1'/ii2 layer in the DAP7 dataset ii1', since that layer undergoes less PCD than 
#ii2, so DAP7 ii1'/ii2 might be more similar to the ii1' layers of other time points
#https://bmcplantbiol.biomedcentral.com/articles/10.1186/1471-2229-8-35/figures/2

dap7[[]]$level_2_annotation <- gsub("ii1'/ii2", "ii1'", dap7[[]]$level_2_annotation)

# Merging
seu3_5 <- merge(dap3, y = dap5)
atlas <- merge(seu3_5, y = dap7)

DefaultAssay(atlas) <- "RNA"  #next step is DE/clustering
atlas <- JoinLayers(atlas)

#adding the colors for the level 1 and 2 clusters

color_key_level1 = c("Embryo" = "#57A15D", 
               "Endosperm" = "#D44A90", 
               "Seed coat" = "#96ADD0", 
               "Funiculus" = "#f1c40f", 
               "Ovule" = "#a6acaf")

level_1_order <- c("Embryo",
             "Endosperm",
             "Seed coat", 
             "Funiculus", 
             "Ovule")

atlas$level_1_annotation <- factor(atlas$level_1_annotation, levels = level_1_order)

color_key_level2 = c("EMB" = "#57A15D", 
               "PEN" = "#ceaab5", 
               "MCE" = "#97576b", 
               "CZE" = "#ce0d48",
               "FUN" = "#f1c40f",
               "OVL" = "#a6acaf", 
               "ii1" = "#aed6f1", 
               "ii1'" = "#5dade2", 
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
            "ii2",
            "oi1",
            "oi2",
            "CZSC",
            "CPT")

atlas$level_2_annotation <- factor(atlas$level_2_annotation, levels = level_2_order)

#re-analyzing 
atlas <- NormalizeData(atlas) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = args$n_VarGenes) %>%
  ScaleData(features = rownames(atlas)) %>%
  RunPCA(npcs=80) %>%
  RunUMAP(dims = 1:80)


#adding the timepoints to the current annotations
atlas@meta.data$level_2_annotation_timed <- sapply(1:nrow(atlas@meta.data), function(x){return(paste(atlas@meta.data$timepoint[x], atlas@meta.data$level_2_annotation[x]))})
#updating the DAP7 inner integuments to be more accurate, "DAP7 ii1'" ->"DAP7 ii1'/ii2"
atlas@meta.data$level_2_annotation_timed <- gsub("DAP7 ii1'", "DAP7 ii1'/ii2", atlas@meta.data$level_2_annotation_timed)
#now, adding a level 2 timed plotting number for easier reading
#making a "plotting ID" for the level 3 annotations
l2_timed_table <- atlas[[]] %>% dplyr::select(c(timepoint, level_2_annotation, level_2_annotation_timed)) %>% 
  group_by(timepoint) %>% 
  arrange(level_2_annotation) %>% 
  arrange(timepoint) %>% 
  ungroup() %>%
  unique()

l2_timed_table$level_2_timed_plotting <- seq(1, nrow(l2_timed_table), 1)
l2_timed_table <- dplyr::select(l2_timed_table, c(level_2_annotation_timed, level_2_timed_plotting))
l2_timed_table$level_2_annotation_timed_final <- sapply(1:nrow(l2_timed_table), function(x){return(paste0(l2_timed_table$level_2_timed_plotting[x],". ", l2_timed_table$level_2_annotation_timed[x]))})

#renaming
level_2_timed_vector <- l2_timed_table$level_2_annotation_timed_final
names(level_2_timed_vector) <- l2_timed_table$level_2_annotation_timed
atlas <- SetIdent(atlas, value = "level_2_annotation_timed")
atlas <- RenameIdents(object = atlas, level_2_timed_vector)
atlas$level_2_annotation_timed_final <- Idents(atlas)

#plotting number only
level_2_timed_num_vector <- l2_timed_table$level_2_timed_plotting
names(level_2_timed_num_vector) <- l2_timed_table$level_2_annotation_timed
atlas <- SetIdent(atlas, value = "level_2_annotation_timed")
atlas <- RenameIdents(object = atlas, level_2_timed_num_vector)
atlas$level_2_timed_plotting <- Idents(atlas)

#now just level 1:
atlas@meta.data$level_1_annotation_timed <- sapply(1:nrow(atlas@meta.data), function(x){return(paste(atlas@meta.data$timepoint[x], atlas@meta.data$level_1_annotation[x]))})

#a couple extra steps for the level 3 timing
level_3_cleaned <- sapply(1:nrow(atlas@meta.data), function(x){return(sub("^\\d+\\.\\s*", "", atlas@meta.data$level_3_annotation[x]))})
atlas@meta.data$level_3_annotation_timed <- sapply(1:nrow(atlas@meta.data), function(x){return(paste(atlas@meta.data$timepoint[x], level_3_cleaned[x]))})

#saving the pre integration object
saveRDS(atlas, "outputs/ATLAS_merged_annotated_orig.rds")
#atlas <- readRDS("outputs/ATLAS_merged_annotated.rds")

dp <- 540
pdf("outputs/figures/ATLAS_merged_figures_orig.pdf")
#level 1
atlas <- SetIdent(atlas, value = "level_1_annotation")
DimPlot(object = atlas, cols = color_key_level1, label = TRUE, raster = TRUE, raster.dpi = c(dp, dp)) +
  theme(legend.position = "none")
key_plot1 <- DimPlot(object = atlas, cols = color_key_level1) 
legend1 <- cowplot::get_legend(key_plot1)
grid.newpage()
grid.draw(legend1)

#level 2
atlas <- SetIdent(atlas, value = "level_2_annotation")
DimPlot(object = atlas, cols = color_key_level2, label = FALSE, raster = TRUE, raster.dpi = c(dp, dp)) +
  theme(legend.position = "none")
key_plot2 <- DimPlot(object = atlas, cols = color_key_level2) 
legend2 <- cowplot::get_legend(key_plot2)
grid.newpage()
grid.draw(legend2)

# now, timepoints
DimPlot(object = atlas, group.by = "timepoint", raster = TRUE, raster.dpi = c(dp, dp)) +
  theme(legend.position = "none")
key_plot3 <- DimPlot(object = atlas, group.by = "timepoint") 
legend3 <- cowplot::get_legend(key_plot3)
grid.newpage()
grid.draw(legend3)

# now, timepoints
DimPlot(object = atlas, group.by = "bio_rep", raster = TRUE, raster.dpi = c(dp, dp)) +
  theme(legend.position = "none")
key_plot4 <- DimPlot(object = atlas, group.by = "bio_rep")
legend4 <- cowplot::get_legend(key_plot4)
grid.newpage()
grid.draw(legend4)

FeaturePlot(atlas, "S.Score", raster = TRUE, raster.dpi = c(dp, dp))
FeaturePlot(atlas, "G2M.Score", raster = TRUE, raster.dpi = c(dp, dp))

#level one timepoints
DimPlot(object = atlas, group.by = "level_1_annotation_timed", label = TRUE, raster = TRUE, raster.dpi = c(dp, dp))+
  theme(legend.position = "none")
key_plot5 <- DimPlot(object = atlas, group.by = "level_1_annotation_timed", label = TRUE)
legend5 <- cowplot::get_legend(key_plot5)
grid.newpage()
grid.draw(legend5)


#level two timepoints
DimPlot(object = atlas, group.by = "level_2_timed_plotting", label = TRUE, raster = TRUE, raster.dpi = c(dp, dp)) +
  theme(legend.position = "none")
key_plot6 <- DimPlot(object = atlas, group.by = "level_2_annotation_timed_final", label = FALSE)
legend6 <- cowplot::get_legend(key_plot6)
grid.newpage()
grid.draw(legend6)

dev.off()

#now integration
#splitting for integration

atlas <- readRDS("outputs/ATLAS_merged_annotated_orig.rds")
DefaultAssay(atlas)<- "RNA"
atlas.split <- atlas %>% SplitObject(split.by = "timepoint") 

#from https://ucdavis-bioinformatics-training.github.io/2023-December-Single-Cell-RNA-Seq-Analysis/data_analysis/08-integration
atlas.split <- lapply(atlas.split, function(sce){
  sce = NormalizeData(sce)
  sce = FindVariableFeatures(sce, selection.method = "vst", nfeatures = args$n_VarGenes)
})

features <- SelectIntegrationFeatures(object.list = atlas.split)

atlas.split <- lapply(atlas.split,function(sce){
  sce = ScaleData(sce, features = rownames(sce))
  RunPCA(sce, features = features, npcs = 80)
})

anchors <- FindIntegrationAnchors(object.list = atlas.split, anchor.features = features, reduction = "rpca")

atlas <- IntegrateData(anchorset = anchors)

DefaultAssay(atlas) <- "integrated"
#keep in mind - just use the integrated slot for clustering and viz

atlas <- ScaleData(atlas, assay="integrated", features = rownames(atlas)) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = args$n_VarGenes) %>%
  RunPCA(assay="integrated", npcs = 80) %>%
  RunUMAP(dims = 1:80, reduction.name = "umap.rpca.dr") %>%
  FindNeighbors(dims = 1:80, verbose = F) %>%
  FindClusters(resolution = 0.8, cluster.name = "clusters_integrated") 

saveRDS(atlas, "outputs/ATLAS_integrated_annotated_rpca.rds")

dp <- 540
pdf("outputs/figures/ATLAS_integrated_figures_rpca.pdf")
#level 1
atlas <- SetIdent(atlas, value = "level_1_annotation")
DimPlot(object = atlas, cols = color_key_level1, label = TRUE, raster = TRUE, raster.dpi = c(dp, dp)) +
  theme(legend.position = "none")
key_plot1 <- DimPlot(object = atlas, cols = color_key_level1) 
legend1 <- cowplot::get_legend(key_plot1)
grid.newpage()
grid.draw(legend1)

#level 2
atlas <- SetIdent(atlas, value = "level_2_annotation")
DimPlot(object = atlas, cols = color_key_level2, label = FALSE, raster = TRUE, raster.dpi = c(dp, dp)) +
  theme(legend.position = "none")
key_plot2 <- DimPlot(object = atlas, cols = color_key_level2) 
legend2 <- cowplot::get_legend(key_plot2)
grid.newpage()
grid.draw(legend2)

# now, timepoints
DimPlot(object = atlas, group.by = "timepoint", raster = TRUE, raster.dpi = c(dp, dp)) +
  theme(legend.position = "none")
key_plot3 <- DimPlot(object = atlas, group.by = "timepoint") 
legend3 <- cowplot::get_legend(key_plot3)
grid.newpage()
grid.draw(legend3)

# now, timepoints
DimPlot(object = atlas, group.by = "bio_rep", raster = TRUE, raster.dpi = c(dp, dp)) +
  theme(legend.position = "none")
key_plot4 <- DimPlot(object = atlas, group.by = "bio_rep")
legend4 <- cowplot::get_legend(key_plot4)
grid.newpage()
grid.draw(legend4)

FeaturePlot(atlas, "S.Score", raster = TRUE, raster.dpi = c(dp, dp))
FeaturePlot(atlas, "G2M.Score", raster = TRUE, raster.dpi = c(dp, dp))

#level one timepoints
DimPlot(object = atlas, group.by = "level_1_annotation_timed", label = TRUE, raster = TRUE, raster.dpi = c(dp, dp))+
  theme(legend.position = "none")
key_plot5 <- DimPlot(object = atlas, group.by = "level_1_annotation_timed", label = TRUE)
legend5 <- cowplot::get_legend(key_plot5)
grid.newpage()
grid.draw(legend5)

#removing/rebuilding old metadata
atlas <- readRDS("outputs/ATLAS_integrated_annotated_rpca.rds")
head(atlas[[]])
atlas@meta.data <- atlas@meta.data[, !(colnames(atlas@meta.data) %in% c("level_3_annotation_timed", 
                                                                       "optimal_clusters_res_1.1", 
                                                                       "optimal_clusters_res_1.5_wbasal_cyst_and_nod", 
                                                                       "original_clusters_res_1.5", "optimal_clusters_res_1.5_wbasal_cyst", 
                                                                       "RNA_snn_res.3", "optimal_clusters_res_1.5", "round2_clusters_integrated", 
                                                                       "round1_clusters_integrated", "RNA_snn_res.0.8", "seurat_clusters", 
                                                                       "round1_clusters_merged", "round2_clusters_merged", "optimal_clusters_res_1.4", 
                                                                       "optimal_clusters_res_1.2_wbasal_cyst", "optimal_clusters_res_1.4_wbasal_cyst", 
                                                                       "original_clusters_res_1.4", "clusters_integrated"))]
#adding a "timed" level 3
atlas$level_3_annotation_full_timed <- paste0(atlas$timepoint, " ", atlas$level_3_annotation_full)
#reordering
atlas@meta.data <- atlas@meta.data[, c("orig.ident", "bio_rep", "timepoint", "nCount_RNA", "nFeature_RNA", "percent.chlor", "percent.mito", 
                                       "percent.cc", "S.Score", "G2M.Score", "Phase", 
                                       "doublet_score", "doublet_class", "outlier_RNA", "outlier_GENES", 
                                       "level_1_annotation", "level_1_annotation_timed", 
                                       "level_2_annotation", "level_2_annotation_timed", "level_2_timed_plotting", "level_2_annotation_timed_final", 
                                       "level_3_annotation_full" , "level_3_annotation_full_timed", "level_3_annotation_abbr",  "level_3_plotting", 
                                       "Nodule_modscore", "Cyst_modscore", "NodLike_modscore")]
saveRDS(atlas, "outputs/ATLAS_integrated_annotated_rpca.rds")

