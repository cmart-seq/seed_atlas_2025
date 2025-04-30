#selecting tissues that can undergo level 2 pseudotime analysis, inspected the outputs of the 01 merging step to determine the root cluster

library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
library(dplyr)
library(readr)
library(pheatmap)
setwd("05_across_timepoints/02_pseudotime/")

set.seed(5432)

tissues <- c("CZE", "EMB", "MCE")
initial_cluster <- c(6, 5, 1)

for(i in 1:length(tissues)){
  dataname <- tissues[i]
  root_cluster <- initial_cluster[[i]]
  
  #reading in  data
  seu <- readRDS(paste0("../01_merging/",dataname, "/", dataname, "_tp_integrated_cc_reg_harmony.rds"))
  DefaultAssay(seu) <- "RNA"
  seu <- JoinLayers(seu)
  seu <- SetIdent(seu, value = "harmony_clusters")
  
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
  
  #converting for monocle3
  cds <- as.cell_data_set(seu)
  
  #Assigning partitions, only one
  recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
  names(recreate.partitions) <- cds@colData@rownames
  recreate.partitions <- as.factor(recreate.partitions)
  cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
  
  #Assigning cluster information
  list.cluster <- seu@active.ident
  cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
  
  #Assigning UMAP coordinates
  cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- seu@reductions$umap.harmony@cell.embeddings
  
  #plotting
  cluster.before.traj <-plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F, 
                                   group_label_size = 5) + theme(legend.position = "right")
  
  pdf(paste0("outputs/", dataname, "/figures/", dataname, "_initial_clusters.pdf"))
  plot(cluster.before.traj)
  plot(DimPlot(seu, group.by = "timepoint", reduction = "umap.harmony"))
  dev.off()
  
  #only one partition
  cds <- learn_graph(cds, use_partition = F)
  
  pdf(paste0("outputs/", dataname, "/figures/", dataname, "_after_learn_graph.pdf"))
  plot(plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
             label_branch_points = T, label_roots = T, label_leaves = F,
             group_label_size = 5))
  plot(DimPlot(seu, group.by = "timepoint", reduction = "umap.harmony"))
  dev.off()
  
  cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, clusters(cds) %in% root_cluster]))
  pdf(paste0("outputs/", dataname, "/figures/", dataname, "_ordered_cells.pdf"))
  plot(plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
             label_branch_points = T, label_roots = F, label_leaves = F))
  dev.off()
  
  head(pseudotime(cds), 10)
  
  cds$monocle3_pseudotime <- pseudotime(cds)
  data.pseudo <- as.data.frame(colData(cds))
  
  pdf(paste0("outputs/", dataname, "/figures/", dataname, "_clusters_by_pseudo.pdf"))
  plot(ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime), fill = seurat_clusters)) + geom_boxplot())
  dev.off()
  
  #Find genes that change as a function of pseudotime
  deg <- graph_test(cds, neighbor_graph = "principal_graph")
  deg_pass <- deg %>% arrange(q_value) %>% filter(status == "OK")
  deg_pass$gene <- rownames(deg_pass)
  write.csv(deg_pass, paste0("outputs/", dataname, "/", dataname, "_genes_in_pseudo.csv"))
  
  good_genes <- rownames(deg_pass)[1:4]
  
  #plotting a few
  pdf(paste0("outputs/", dataname, "/figures/", dataname, "_clusters_by_pseudo.pdf"))
  plot(FeaturePlot(seu, features = good_genes, reduction = "umap.harmony"))
  dev.off()
  
  #Add pseudotime values into the seuratobject
  seu$level_2_pseudotime <- pseudotime(cds)
  pdf(paste0("outputs/", dataname, "/figures/", dataname, "_pseudo_on_seurat.pdf"))
  plot(FeaturePlot(seu, features = "level_2_pseudotime", reduction = "umap.harmony"))
  dev.off()
  
  #plot some genes in pseudotime
  pdf(paste0("outputs/", dataname, "/figures/", dataname, "_genes_in_pseudo.pdf"))
  my_genes <-good_genes
  cds_subset <- cds[my_genes,]
  plot(plot_genes_in_pseudotime(cds_subset, color_cells_by = "monocle3_pseudotime", label_by_short_name = FALSE))
  dev.off()
  
  #saving seu
  saveRDS(seu, paste0("outputs/", dataname, "/", dataname, "_pseudo.rds"))
  #saving cds
  saveRDS(cds, paste0("outputs/", dataname, "/", dataname, "_pseudo_cds.rds"))
  
  #saving just cells and pseudo
  write_csv(data.frame(cells = rownames(seu[[]]),
                       level_2_pseudotime = seu[[]]$level_2_pseudotime), paste0("outputs/", dataname, "/", dataname, "_cells_in_pseudo.csv"))
  
}

