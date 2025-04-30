#this script performs harmony integration for the level 2 annotations that are present at all timepoints to prepare for pseudotime analysis

library(tidyverse)
library(Seurat)
library(scran)
library(bluster)
library(argparse)
library(cluster)
library(dplyr)
library(cowplot)
library(harmony)
library(org.At.tair.db)

set.seed(54321)

setwd("05_across_timepoints/01_merging")

#making the output directory if it doesnt exist
#reading in data
dap3 <- readRDS("../../04_manual_annotation/outputs/DAP3_wcze_subs/DAP3_wcze_subs_annotated.rds") %>% SetIdent(value = "level_2_annotation")
dap5 <- readRDS("../../04_manual_annotation/outputs/DAP5_wcze_subs/DAP5_wcze_subs_annotated.rds") %>% SetIdent(value = "level_2_annotation")
DefaultAssay(dap5) <- "RNA"
dap5[['integrated']] <- NULL
dap7 <- readRDS("../../04_manual_annotation/outputs/DAP7/DAP7_annotated.rds") %>% SetIdent(value = "level_2_annotation")

#updating the ii1'/ii2 in DAP7 so that it can be processed with the other ii1' clusters
dap7 <- RenameIdents(dap7, "ii1'/ii2" = "ii1'")
dap7$level_2_annotation <- Idents(dap7)

shared_idents <- c("CZE", "PEN","EMB", "MCE", "ii1'","CPT", "ii2","oi1","oi2")

for(i in shared_idents){
  #in all three: EMB, PEN,  MCE,  CZE, ii1, ii1', oi1, oi2, CZSC
  if(i %in% unique(dap3$level_2_annotation) && i %in% unique(dap5$level_2_annotation) && i %in% unique(dap7$level_2_annotation)){
    subs3 <- subset(dap3, subset = level_2_annotation == i)
    subs5 <- subset(dap5, subset = level_2_annotation == i)
    subs7 <- subset(dap7, subset = level_2_annotation == i)
    #merging
    seu3_5 <- merge(subs3, y = subs5)
    seu <- merge(seu3_5, y = subs7)
    
  #Just in DAP3 and DAP5: FUN, ii2, CPT
  } else if(i %in% unique(dap3$level_2_annotation) && i %in% unique(dap5$level_2_annotation)){
    subs3 <- subset(dap3, subset = level_2_annotation == i)
    subs5 <- subset(dap5, subset = level_2_annotation == i)
    #merging
    seu <- merge(subs3, y = subs5)
  }

  #now performing the integration by timepoint
  DefaultAssay(seu)<- "RNA"
  seu <- JoinLayers(seu)
  
  #getting ready
  pseudo_prep <- function(seurat_object, res){
    seurat_object <- NormalizeData(seurat_object) %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 2500) %>%
      #using the cell cycle scoring from the orignal analysis, regressing out to hopefully get a simpler pseudotime trajectory
      ScaleData(vars.to.regress = c("S.Score", "G2M.Score"), assay = "RNA") %>%
      RunPCA(npcs = 30) %>%
      FindNeighbors(dims = 1:30, verbose = F) %>%
      FindClusters(resolution = res) %>% 
      RunUMAP(dims = 1:30, verbose = F)
    return(seurat_object)
  }
  
  seu <- pseudo_prep(seu, 0.5)

  #running Harmony across timepoints
  seu <- RunHarmony(object = seu, reduction = "pca", group.by.vars = "timepoint", 
                               dims.use = 1:30, reduction.save = 'harmony', plot_convergence = F) %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
    FindClusters(resolution = 0.5, cluster.name = "harmony_clusters") %>% 
    RunUMAP(reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
  
  seu$level_3_annotation_abbr_timed <- paste0(seu$timepoint,"_", seu$level_3_annotation_abbr)

  if (!dir.exists(i)) {
    dir.create(i, recursive = TRUE)
  }
  
  #saving the object
  saveRDS(seu, paste0(i,"/",i,"_tp_integrated_cc_reg_harmony.rds"))

  pdf(paste0(i,"/",i,"_tp_integrated_cc_reg_harmony.pdf"))
  #timepoint
  plot(DimPlot(seu, reduction = "umap.harmony", group.by = c("timepoint"), label.size = 2))
  plot(DimPlot(seu, reduction = "umap.harmony", group.by = c("harmony_clusters"), label.size = 2))
  #plotting with previous annotations
  seu <- SetIdent(seu, value = "level_3_annotation_abbr_timed")
  plot(DimPlot(seu, label = T, repel = T, reduction = "umap.harmony") + theme(legend.position = "none"))
  dev.off()
}


#also doing this for the level 1 annotations
#performing analysis on everything but the ovule and funiculus
shared_idents <- c("Seed coat", "Embryo", "Endosperm")
for(i in shared_idents){
  #in all three: 
  #i = "Endosperm"
  subs3 <- subset(dap3, subset = level_1_annotation == i)
  subs5 <- subset(dap5, subset = level_1_annotation == i)
  subs7 <- subset(dap7, subset = level_1_annotation == i)
  #merging
  seu3_5 <- merge(subs3, y = subs5)
  seu <- merge(seu3_5, y = subs7)
  
  #now performing the integration by timepoint
  DefaultAssay(seu)<- "RNA"
  seu <- JoinLayers(seu)
  
  #getting ready
  pseudo_prep <- function(seurat_object, res){
    seurat_object <- NormalizeData(seurat_object) %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 2500) %>%
      #using the cell cycle scoring from the orignal analysis, regressing out to hopefully get a simpler pseudotime trajectory
      ScaleData(vars.to.regress = c("S.Score", "G2M.Score"), assay = "RNA") %>%
      RunPCA(npcs = 30) %>%
      FindNeighbors(dims = 1:30, verbose = F) %>%
      FindClusters(resolution = res) %>% 
      RunUMAP(dims = 1:30, verbose = F)
    return(seurat_object)
  }
  
  seu <- pseudo_prep(seu, 0.5)
  
  #running Harmony across timepoints
  seu <- RunHarmony(object = seu, reduction = "pca", group.by.vars = "timepoint", 
                    dims.use = 1:30, reduction.save = 'harmony', plot_convergence = F) %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
    FindClusters(resolution = 0.5, cluster.name = "harmony_clusters") %>% 
    #for MCE, had to reduce to 0.1
    RunUMAP(reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
  
  seu$level_3_annotation_abbr_timed <- paste0(seu$timepoint,"_", seu$level_3_annotation_abbr)
  
  if (!dir.exists(i)) {
    dir.create(i, recursive = TRUE)
  }
  
  #saving the object
  saveRDS(seu, paste0(i,"/",i,"_tp_integrated_cc_reg_harmony.rds"))
  
  pdf(paste0(i,"/",i,"_tp_integrated_cc_reg_harmony.pdf"))
  #timepoint
  plot(DimPlot(seu, reduction = "umap.harmony", group.by = c("timepoint"), label.size = 2))
  plot(DimPlot(seu, reduction = "umap.harmony", group.by = c("harmony_clusters"), label.size = 2))
  #plotting with previous annotations
  seu <- SetIdent(seu, value = "level_3_annotation_abbr_timed")
  plot(DimPlot(seu, label = T, repel = T, reduction = "umap.harmony") + theme(legend.position = "none"))
  dev.off()
}





