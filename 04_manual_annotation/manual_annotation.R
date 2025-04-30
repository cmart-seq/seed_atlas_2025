#This script assigns annotations to clusters based on a metadata sheet. The metadata pairs the optimal clustering resolution
#with a level 1, 2, and 3, annotation based on published markers or in situ validation

#in Martin et al 2025, this step was performed twice for 3 and 5 DAP in order to annotate the CZE clusters that were identified in a subclustering analysis
#example usage: 
#Rscript manual_annotation.R  --wd . --dataset DAP7  --harmony FALSE --seupath ../03_clustering/outputs/DAP7/DAP7_clustered.rds --mpath inputs/DAP7_level_3_manual_metadata.csv

library(tidyr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(scCustomize)
library(argparse)
library(readr)
library(stringr)
library(org.At.tair.db)

set.seed(54321)

#getting the args
parser <- ArgumentParser()

parser$add_argument("--wd", help = "working directory")
parser$add_argument("--dataset", help = "sample name, of the form TIMEPOINT (one of DAP3, DAP5, DAP7, or ATLAS)")
parser$add_argument("--harmony", help = "if data was harmony-integrated") # 
parser$add_argument("--seupath", help = "path to seurat object") # 
parser$add_argument("--mpath", help = "path to metadata") # 
  
#parse the arguments
args <- parser$parse_args()

dataname <- args$dataset

#orient to the working directory
setwd(args$wd)

#making the output directory if it doesnt exist
#outputs
if (!dir.exists("outputs")) {
  dir.create("outputs", recursive = TRUE)
}

#one for the sample
if (!dir.exists(file.path("outputs",dataname))) {
  dir.create(file.path("outputs", dataname), recursive = TRUE)
}

#one for figures
if (!dir.exists(file.path("outputs", dataname, "figures"))) {
  dir.create(file.path("outputs", dataname, "figures"), recursive = TRUE)
}

outpath = file.path(args$wd, "outputs", dataname)

#checking if integrated
red <- ifelse(args$harmony, "umap.harmony", "umap")

#loading the seurat
seu <- readRDS(args$seupath)
DefaultAssay(seu) <- "RNA"

#updating the "timepoint" metadata so it fits with convention
conv_timepoint <- function(origID){
  orig_timepoint <- str_split(origID, "_")[[1]][1]
  new_timepoint <- paste0(substr(orig_timepoint,4,4), "DAP")
}

seu$timepoint <- sapply(seu$orig.ident, conv_timepoint)

#getting the resolution
res <- grep("optimal_clusters_res_",colnames(seu[[]]), value = TRUE)

seu <- SetIdent(seu, value = res)

#making a dimplot of optimal clusters
pdf(paste0(outpath, "/figures/", dataname, "_optimal_clusters.pdf"))
DimPlot(seu, label = TRUE, reduction = red, alpha = 0.3) +ggtitle(res)
dev.off()

#making color keys for the unified levels, level 1 and 2

color_key1 = c("Embryo" = "#57A15D", 
               "Endosperm" = "#D44A90", 
               "Seed coat" = "#96ADD0", 
               "Funiculus" = "#f1c40f", 
               "Ovule" = "#a6acaf")

levels1 <- c("Embryo",
             "Endosperm",
             "Seed coat", 
             "Funiculus", 
             "Ovule")

color_key2 = c("EMB" = "#57A15D", 
               "PEN" = "#F17DB1", 
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
              
levels2 = c("EMB",
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

#loading the cited genes

cited_genes <- read_csv("marker_citations.csv")

#also making sure that all of the genes in the signature are present in the object
cited_genes <- filter(cited_genes, TAIR %in% rownames(seu))

#setting up outputs
if (!dir.exists(paste0(outpath, "/figures/genes/"))) {
  dir.create(paste0(outpath, "/figures/genes/"), recursive = TRUE)
}

#plotting all of the genes
for (i in na.omit(unique(cited_genes$`general expression pattern`))){
  pdf(paste0(outpath, "/figures/genes/", sub(" ", "_", i), "_markers.pdf"))
  genes <- filter(cited_genes, `general expression pattern` == i)
 for(q in 1:nrow(genes)){
    plot(FeaturePlot(seu, genes$TAIR[q], reduction = red)+
           ggtitle(paste0(genes$TAIR[q], ", ", genes$symbol[q]))) 
  }
  dev.off()
}

#inspecting the plots, manually labeling
#I will label level 3 first (most granular)
#the below IDs were identified using the plots generated in the automated annotation and the known
#gene plotting above

if(dataname == "DAP3" | dataname == "DAP3_wcze_subs" ){
  
  #reading in the manual metadata
  metadata <- read_csv(args$mpath) %>% filter_all(any_vars( . != ""))
  
  #making a "plotting ID" for the level 3 annotations
  plotting_ID3_table <- metadata %>% dplyr::select(c(level_1_annotation, level_2_annotation, level_3_annotation_abbr, level_3_annotation_full)) %>% 
    group_by(level_1_annotation) %>% 
    arrange(level_2_annotation) %>% 
    arrange(level_1_annotation) %>% 
    ungroup() %>%
    unique()
  
  plotting_ID3_table$level_3_plotting <- seq(1, length(unique(plotting_ID3_table$level_3_annotation_full)), 1)
  plotting_ID3_table <- dplyr::select(plotting_ID3_table, c(level_3_annotation_full, level_3_plotting))
  
  #merging and saving with the old metadata
  metadata <- left_join(metadata, plotting_ID3_table, by ="level_3_annotation_full")
  
  #new, marker based IDs
  #full level 3
  level_3_annotation_vector_full <- metadata$level_3_annotation_full
  names(level_3_annotation_vector_full) <- metadata[[res]]
  #full level 3 updating IDs, adding full names
  seu <- SetIdent(seu, value = res)
  seu <- RenameIdents(object = seu, level_3_annotation_vector_full)
  seu$level_3_annotation_full <- Idents(seu)
  
  #abbr level 3
  level_3_annotation_vector_abbr <- metadata$level_3_annotation_abbr
  names(level_3_annotation_vector_abbr) <- metadata[[res]]
  #abbr level 3 updating IDs, adding full names
  seu <- SetIdent(seu, value = res)
  seu <- RenameIdents(object = seu, level_3_annotation_vector_abbr)
  seu$level_3_annotation_abbr <- Idents(seu)
  
  #plotting level 3
  level_3_annotation_vector_plotting <- metadata$level_3_plotting
  names(level_3_annotation_vector_plotting) <- metadata[[res]]
  #abbr level 3 updating IDs, adding full names
  seu <- SetIdent(seu, value = res)
  seu <- RenameIdents(object = seu, level_3_annotation_vector_plotting)
  seu$level_3_plotting <- Idents(seu)
  
  #adding the plotting number to the level 2 annotations
  seu[[]] <-mutate(seu[[]], level_3_annotation_full = paste0(level_3_plotting, ". ",level_3_annotation_full) )
  seu[[]] <-mutate(seu[[]], level_3_annotation_abbr = paste0(level_3_plotting, ". ",level_3_annotation_abbr) )
  
  #factoring for the right order, level 3
  level_3_factor_full <- seu[[]] %>% dplyr::select(c(level_3_annotation_full, level_3_plotting))
  level_3_factor_full$level_3_plotting <- as.numeric(as.character(level_3_factor_full$level_3_plotting))
  level_3_factor_full <-   level_3_factor_full %>%
    arrange(level_3_plotting) %>% 
    pull(level_3_annotation_full) %>% unique()
  
  seu$level_3_annotation_full <- factor(seu$level_3_annotation_full, levels = level_3_factor_full)
  
  level_3_factor_abbr <- seu[[]] %>% dplyr::select(c(level_3_annotation_abbr, level_3_plotting))%>% arrange(level_3_plotting) %>% pull(level_3_annotation_abbr) %>% unique()
  seu$level_3_annotation_abbr <- factor(seu$level_3_annotation_abbr, levels = level_3_factor_abbr)

  #factoring for the right order, level 2
  
  #level 2 
  level_2_annotation_vector <- metadata$level_2_annotation
  names(level_2_annotation_vector) <- metadata[[res]]
  #abbr level 3 updating IDs, adding full names
  seu <- SetIdent(seu, value = res)
  seu <- RenameIdents(object = seu, level_2_annotation_vector)
  seu$level_2_annotation <- Idents(seu)
  
  seu$level_2_annotation <- factor(seu$level_2_annotation, levels = levels2)
  
  #level 1 
  level_1_annotation_vector <- metadata$level_1_annotation
  names(level_1_annotation_vector) <- metadata[[res]]
  #abbr level 3 updating IDs, adding full names
  seu <- SetIdent(seu, value = res)
  seu <- RenameIdents(object = seu, level_1_annotation_vector)
  seu$level_1_annotation <- Idents(seu)
  
  #plotting
  seu <- SetIdent(seu, value = "level_3_annotation_full")
  
  #getting the wordy legend
  p1 <- DimPlot(seu, label = T, repel = T, reduction = red, alpha = 0.3) 
  key <- get_legend(p1)
  
  #plotting without the legend
  seu <- SetIdent(seu, value = "level_3_plotting")
  p2 <- DimPlot(seu, label = T, repel = T, reduction = red, alpha = 0.3) + theme(legend.position = "none")
  p3 <- DimPlot(seu, repel = T, raster = TRUE, reduction = red, alpha = 0.3) + theme(legend.position = "none")
  
  pdf(paste0(outpath, "/figures/", dataname, "_level_3_cluster_key.pdf"), width = 12, height = 8)
  plot(plot_grid(key, ncol = 1))
  dev.off()
  
  pdf(paste0(outpath, "/figures/", dataname, "_level_3_clusters.pdf"), width = 8, height = 8)
  plot(p2)
  dev.off()
  
  png(paste0(outpath, "/figures/", dataname, "_level_3_clusters.png"), width = 8, height = 8, units = "in", res = 300)
  plot(p2)
  dev.off()
  
  #now level 2
  seu <- SetIdent(seu, value = "level_2_annotation")
  
  pdf(paste0(outpath, "/figures/", dataname, "_level_2_clusters.pdf"), width = 8, height = 8)
  plot(DimPlot(seu, label = T, repel = T, cols = color_key2, reduction = red, alpha = 0.3))
  dev.off()
  
  png(paste0(outpath, "/figures/", dataname, "_level_2_clusters.png"), width = 8, height = 8, units = "in", res = 300)
  plot(DimPlot(seu, label = TRUE, repel = TRUE, cols = color_key2, reduction = red, alpha = 0.3))
  dev.off()
  
  #level 1
  seu <- SetIdent(seu, value = "level_1_annotation")
  
  pdf(paste0(outpath, "/figures/", dataname, "_level_1_clusters.pdf"), width = 8, height = 8)
  plot(DimPlot(seu, label = T, repel = T, cols = color_key1, reduction = red, alpha = 0.3))
  dev.off()
  
  png(paste0(outpath, "/figures/", dataname, "_level_1_clusters.png"), width = 8, height = 8, units = "in", res = 300)
  plot(DimPlot(seu, label = TRUE, repel = TRUE, cols = color_key1, reduction = red, alpha = 0.3))
  dev.off()
  
  saveRDS(seu, paste0(outpath, "/", dataname, "_annotated.rds"))
  
}else if(dataname == "DAP5" | dataname == "DAP5_wcze_subs" ){
  #reading in the manual metadata
  metadata <- read_csv(args$mpath) %>% filter_all(any_vars( . != ""))
  
  #making a "plotting ID" for the level 3 annotations
  plotting_ID3_table <- metadata %>% dplyr::select(c(level_1_annotation, level_2_annotation, level_3_annotation_abbr, level_3_annotation_full)) %>% 
    group_by(level_1_annotation) %>% 
    arrange(level_2_annotation) %>% 
    arrange(level_1_annotation) %>% 
    ungroup() %>%
    unique()
  
  plotting_ID3_table$level_3_plotting <- seq(1, length(unique(plotting_ID3_table$level_3_annotation_full)), 1)
  plotting_ID3_table <- dplyr::select(plotting_ID3_table, c(level_3_annotation_full, level_3_plotting))
  
  #merging and saving with the old metadata
  metadata <- left_join(metadata, plotting_ID3_table, by ="level_3_annotation_full")
  
  #new, marker based IDs
  #full level 3
  level_3_annotation_vector_full <- metadata$level_3_annotation_full
  names(level_3_annotation_vector_full) <- metadata[[res]]
  #full level 3 updating IDs, adding full names
  seu <- SetIdent(seu, value = res)
  seu <- RenameIdents(object = seu, level_3_annotation_vector_full)
  seu$level_3_annotation_full <- Idents(seu)
  
  #abbr level 3
  level_3_annotation_vector_abbr <- metadata$level_3_annotation_abbr
  names(level_3_annotation_vector_abbr) <- metadata[[res]]
  #abbr level 3 updating IDs, adding full names
  seu <- SetIdent(seu, value = res)
  seu <- RenameIdents(object = seu, level_3_annotation_vector_abbr)
  seu$level_3_annotation_abbr <- Idents(seu)
  
  #plotting level 3
  level_3_annotation_vector_plotting <- metadata$level_3_plotting
  names(level_3_annotation_vector_plotting) <- metadata[[res]]
  #abbr level 3 updating IDs, adding full names
  seu <- SetIdent(seu, value = res)
  seu <- RenameIdents(object = seu, level_3_annotation_vector_plotting)
  seu$level_3_plotting <- Idents(seu)
  
  #adding the plotting number to the level 2 annotations
  seu[[]] <-mutate(seu[[]], level_3_annotation_full = paste0(level_3_plotting, ". ",level_3_annotation_full) )
  seu[[]] <-mutate(seu[[]], level_3_annotation_abbr = paste0(level_3_plotting, ". ",level_3_annotation_abbr) )
  
  #factoring for the right order, level 3
  level_3_factor_full <- seu[[]] %>% dplyr::select(c(level_3_annotation_full, level_3_plotting))
  level_3_factor_full$level_3_plotting <- as.numeric(as.character(level_3_factor_full$level_3_plotting))
  level_3_factor_full <-   level_3_factor_full %>%
    arrange(level_3_plotting) %>% 
    pull(level_3_annotation_full) %>% unique()
  
  seu$level_3_annotation_full <- factor(seu$level_3_annotation_full, levels = level_3_factor_full)
  
  level_3_factor_abbr <- seu[[]] %>% dplyr::select(c(level_3_annotation_abbr, level_3_plotting))%>% arrange(level_3_plotting) %>% pull(level_3_annotation_abbr) %>% unique()
  seu$level_3_annotation_abbr <- factor(seu$level_3_annotation_abbr, levels = level_3_factor_abbr)
  
  #factoring for the right order, level 2
  
  
  #level 2 
  level_2_annotation_vector <- metadata$level_2_annotation
  names(level_2_annotation_vector) <- metadata[[res]]
  #abbr level 3 updating IDs, adding full names
  seu <- SetIdent(seu, value = res)
  seu <- RenameIdents(object = seu, level_2_annotation_vector)
  seu$level_2_annotation <- Idents(seu)
  
  seu$level_2_annotation <- factor(seu$level_2_annotation, levels = levels2)
  
  #level 1 
  level_1_annotation_vector <- metadata$level_1_annotation
  names(level_1_annotation_vector) <- metadata[[res]]
  #abbr level 3 updating IDs, adding full names
  seu <- SetIdent(seu, value = res)
  seu <- RenameIdents(object = seu, level_1_annotation_vector)
  seu$level_1_annotation <- Idents(seu)
  
  #plotting
  seu <- SetIdent(seu, value = "level_3_annotation_full")
  
  #getting the wordy legend
  p1 <- DimPlot(seu, label = T, repel = T, reduction = red, alpha = 0.3) 
  key <- get_legend(p1)
  
  #plotting without the legend
  seu <- SetIdent(seu, value = "level_3_plotting")
  p2 <- DimPlot(seu, label = T, repel = T, reduction = red, alpha = 0.3) + theme(legend.position = "none")
  p3 <- DimPlot(seu, repel = T, raster = TRUE, reduction = red, alpha = 0.3) + theme(legend.position = "none")
  
  pdf(paste0(outpath, "/figures/", dataname, "_level_3_cluster_key.pdf"), width = 12, height = 8)
  plot(plot_grid(key, ncol = 1))
  dev.off()
  
  pdf(paste0(outpath, "/figures/", dataname, "_level_3_clusters.pdf"), width = 8, height = 8)
  plot(p2)
  dev.off()
  
  png(paste0(outpath, "/figures/", dataname, "_level_3_clusters.png"), width = 8, height = 8, units = "in", res = 300)
  plot(p2)
  dev.off()
  
  #now level 2
  seu <- SetIdent(seu, value = "level_2_annotation")
  
  pdf(paste0(outpath, "/figures/", dataname, "_level_2_clusters.pdf"), width = 8, height = 8)
  plot(DimPlot(seu, label = T, repel = T, cols = color_key2, reduction = red, alpha = 0.3))
  dev.off()
  
  png(paste0(outpath, "/figures/", dataname, "_level_2_clusters.png"), width = 8, height = 8, units = "in", res = 300)
  plot(DimPlot(seu, label = TRUE, repel = TRUE, cols = color_key2, reduction = red, alpha = 0.3))
  dev.off()
  
  #level 1
  seu <- SetIdent(seu, value = "level_1_annotation")
  
  pdf(paste0(outpath, "/figures/", dataname, "_level_1_clusters.pdf"), width = 8, height = 8)
  plot(DimPlot(seu, label = T, repel = T, cols = color_key1, reduction = red, alpha = 0.3))
  dev.off()
  
  png(paste0(outpath, "/figures/", dataname, "_level_1_clusters.png"), width = 8, height = 8, units = "in", res = 300)
  plot(DimPlot(seu, label = TRUE, repel = TRUE, cols = color_key1, reduction = red, alpha = 0.3))
  dev.off()
  
  saveRDS(seu, paste0(outpath, "/", dataname, "_annotated.rds"))
  
}else if(dataname == "DAP7"){
  #need to update for the inner integument layers
  
  color_key2 = c("EMB" = "#57A15D", 
                 "PEN" = "#ceaab5", 
                 "MCE" = "#97576b", 
                 "CZE" = "#ce0d48",
                 "FUN" = "#f1c40f",
                 "OVL" = "#a6acaf", 
                 "ii1" = "#aed6f1", 
                 "ii1'/ii2" = "#21618c",
                 "oi1" = "#bb8fce", 
                 "oi2" = "#6c3483", 
                 "CZSC" = "#088F8F")
  
  levels2 = c("EMB",
              "PEN",
              "MCE",
              "CZE",
              "FUN",
              "OVL",
              "ii1",
              "ii1'/ii2",
              "oi1",
              "oi2",
              "CZSC")
  
  #reading in the manual metadata
  metadata <- read_csv(args$mpath) %>% filter_all(any_vars( . != ""))
  
  #making a "plotting ID" for the level 3 annotations
  plotting_ID3_table <- metadata %>% dplyr::select(c(level_1_annotation, level_2_annotation, level_3_annotation_abbr, level_3_annotation_full)) %>% 
    group_by(level_1_annotation) %>% 
    arrange(level_2_annotation) %>% 
    arrange(level_1_annotation) %>% 
    ungroup() %>%
    unique()
  
  plotting_ID3_table$level_3_plotting <- seq(1, length(unique(plotting_ID3_table$level_3_annotation_full)), 1)
  plotting_ID3_table <- dplyr::select(plotting_ID3_table, c(level_3_annotation_full, level_3_plotting))
  
  #merging and saving with the old metadata
  metadata <- left_join(metadata, plotting_ID3_table, by ="level_3_annotation_full")
  
  #new, marker based IDs
  #full level 3
  level_3_annotation_vector_full <- metadata$level_3_annotation_full
  names(level_3_annotation_vector_full) <- metadata[[res]]
  #full level 3 updating IDs, adding full names
  seu <- SetIdent(seu, value = res)
  seu <- RenameIdents(object = seu, level_3_annotation_vector_full)
  seu$level_3_annotation_full <- Idents(seu)
  
  #abbr level 3
  level_3_annotation_vector_abbr <- metadata$level_3_annotation_abbr
  names(level_3_annotation_vector_abbr) <- metadata[[res]]
  #abbr level 3 updating IDs, adding full names
  seu <- SetIdent(seu, value = res)
  seu <- RenameIdents(object = seu, level_3_annotation_vector_abbr)
  seu$level_3_annotation_abbr <- Idents(seu)
  
  #plotting level 3
  level_3_annotation_vector_plotting <- metadata$level_3_plotting
  names(level_3_annotation_vector_plotting) <- metadata[[res]]
  #abbr level 3 updating IDs, adding full names
  seu <- SetIdent(seu, value = res)
  seu <- RenameIdents(object = seu, level_3_annotation_vector_plotting)
  seu$level_3_plotting <- Idents(seu)
  
  #adding the plotting number to the level 2 annotations
  seu[[]] <-mutate(seu[[]], level_3_annotation_full = paste0(level_3_plotting, ". ",level_3_annotation_full) )
  seu[[]] <-mutate(seu[[]], level_3_annotation_abbr = paste0(level_3_plotting, ". ",level_3_annotation_abbr) )
  
  #factoring for the right order, level 3
  level_3_factor_full <- seu[[]] %>% dplyr::select(c(level_3_annotation_full, level_3_plotting))
  level_3_factor_full$level_3_plotting <- as.numeric(as.character(level_3_factor_full$level_3_plotting))
  level_3_factor_full <-   level_3_factor_full %>%
    arrange(level_3_plotting) %>% 
    pull(level_3_annotation_full) %>% unique()
  
  seu$level_3_annotation_full <- factor(seu$level_3_annotation_full, levels = level_3_factor_full)
  
  level_3_factor_abbr <- seu[[]] %>% dplyr::select(c(level_3_annotation_abbr, level_3_plotting))%>% arrange(level_3_plotting) %>% pull(level_3_annotation_abbr) %>% unique()
  seu$level_3_annotation_abbr <- factor(seu$level_3_annotation_abbr, levels = level_3_factor_abbr)
  
  #factoring for the right order, level 2
  
  
  #level 2 
  level_2_annotation_vector <- metadata$level_2_annotation
  names(level_2_annotation_vector) <- metadata[[res]]
  #abbr level 3 updating IDs, adding full names
  seu <- SetIdent(seu, value = res)
  seu <- RenameIdents(object = seu, level_2_annotation_vector)
  seu$level_2_annotation <- Idents(seu)
  
  seu$level_2_annotation <- factor(seu$level_2_annotation, levels = levels2)
  
  #level 1 
  level_1_annotation_vector <- metadata$level_1_annotation
  names(level_1_annotation_vector) <- metadata[[res]]
  #abbr level 3 updating IDs, adding full names
  seu <- SetIdent(seu, value = res)
  seu <- RenameIdents(object = seu, level_1_annotation_vector)
  seu$level_1_annotation <- Idents(seu)
  
  #plotting
  seu <- SetIdent(seu, value = "level_3_annotation_full")
  
  #getting the wordy legend
  p1 <- DimPlot(seu, label = T, repel = T, reduction = red, alpha = 0.3) 
  key <- get_legend(p1)
  
  #plotting without the legend
  seu <- SetIdent(seu, value = "level_3_plotting")
  p2 <- DimPlot(seu, label = T, repel = T, reduction = red, alpha = 0.3) + theme(legend.position = "none")
  p3 <- DimPlot(seu, repel = T, raster = TRUE, reduction = red, alpha = 0.3) + theme(legend.position = "none")
  
  pdf(paste0(outpath, "/figures/", dataname, "_level_3_cluster_key.pdf"), width = 12, height = 8)
  plot(plot_grid(key, ncol = 1))
  dev.off()
  
  pdf(paste0(outpath, "/figures/", dataname, "_level_3_clusters.pdf"), width = 8, height = 8)
  plot(p2)
  dev.off()
  
  png(paste0(outpath, "/figures/", dataname, "_level_3_clusters.png"), width = 8, height = 8, units = "in", res = 300)
  plot(p2)
  dev.off()
  
  #now level 2
  seu <- SetIdent(seu, value = "level_2_annotation")
  
  pdf(paste0(outpath, "/figures/", dataname, "_level_2_clusters.pdf"), width = 8, height = 8)
  plot(DimPlot(seu, label = T, repel = T, cols = color_key2, reduction = red, alpha = 0.3))
  dev.off()
  
  png(paste0(outpath, "/figures/", dataname, "_level_2_clusters.png"), width = 8, height = 8, units = "in", res = 300)
  plot(DimPlot(seu, label = TRUE, repel = TRUE, cols = color_key2, reduction = red, alpha = 0.3))
  dev.off()
  
  #level 1
  seu <- SetIdent(seu, value = "level_1_annotation")
  
  pdf(paste0(outpath, "/figures/", dataname, "_level_1_clusters.pdf"), width = 8, height = 8)
  plot(DimPlot(seu, label = T, repel = T, cols = color_key1, reduction = red, alpha = 0.3))
  dev.off()
  
  png(paste0(outpath, "/figures/", dataname, "_level_1_clusters.png"), width = 8, height = 8, units = "in", res = 300)
  plot(DimPlot(seu, label = TRUE, repel = TRUE, cols = color_key1, reduction = red, alpha = 0.3))
  dev.off()
  
  
  saveRDS(seu, paste0(outpath, "/", dataname, "_annotated.rds"))

}
