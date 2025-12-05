#performing pseudotime on 5 DAP PEN and CZE to define the apical-basal CZE trajectory
library(monocle3)
library(Seurat)
library(dplyr)
library(SeuratWrappers)
library(ggplot2)
library(readr)
library(harmony)
setwd("02_pseudotime/czen_spectrum_pseudo")

set.seed(123)

#loading the data
dap5 <- readRDS("../../../04_manual_annotation/outputs/DAP5_wcze_subs/DAP5_wcze_subs_annotated.rds")

#getting ready, not integrating
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

#preparing the 5DAP object for this analysis
cze_pen_dap5 <- subset(dap5, subset = level_2_annotation %in% c("CZE", "PEN"))
cze_pen_dap5 <- pseudo_prep(cze_pen_dap5, 0.5)

pdf("cze_pen_dap5_preharmony.pdf")
plot(DimPlot(cze_pen_dap5, label = TRUE, reduction = "umap.harmony", group.by = c("bio_rep"), label.size = 2)) +ggtitle("CZE and PEN 5DAP")
plot(DimPlot(cze_pen_dap5, label = TRUE, reduction = "umap.harmony", group.by = c("level_3_annotation_abbr"), label.size = 2))+ggtitle("CZE and PEN 5DAP")
dev.off()

#running Harmony across bio_reps, 5DAP cze and pen
cze_pen_dap5 <- RunHarmony(object = cze_pen_dap5, reduction = "pca", group.by.vars = "bio_rep", 
                         dims.use = 1:30, reduction.save = 'harmony', plot_convergence = F) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.5, cluster.name = "harmony_clusters") %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

pdf("cze_pen_dap5_postharmony.pdf")
plot(DimPlot(cze_pen_dap5, label = TRUE, reduction = "umap.harmony", group.by = c("bio_rep"), label.size = 2)) +ggtitle("CZE and PEN 5DAP")
plot(DimPlot(cze_pen_dap5, label = TRUE, reduction = "umap.harmony", group.by = c("harmony_clusters"), label.size = 2))+ggtitle("CZE and PEN 5DAP")
plot(DimPlot(cze_pen_dap5, label = TRUE, reduction = "umap.harmony", group.by = c("level_3_annotation_abbr"), label.size = 2))+ggtitle("CZE and PEN 5DAP")
dev.off()

#performing pseudotime analysis
chalazal_pseudotime <- function(seu, root, name){
  #setting to integrated clusters
  seu <- SetIdent(seu, value = "harmony_clusters")
  
  #converting for monocle3
  cds <- as.cell_data_set(seu)
  
  #Assigning partitions, only one
  recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
  names(recreate.partitions) <- cds@colData@rownames
  recreate.partitions <- as.factor(recreate.partitions)
  cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
  
  #Assign cluster information
  list.cluster <- seu@active.ident
  cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
  
  #Assign UMAP coordinates
  cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- seu@reductions$umap.harmony@cell.embeddings
  
  #plotting
  cluster.before.traj <-plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F, 
                                   group_label_size = 5) + theme(legend.position = "right")
  
  pdf(paste0(name, "_pre_trajectory.pdf"))
  plot(cluster.before.traj)
  plot(DimPlot(seu, group.by = "timepoint", reduction = "umap.harmony"))
  dev.off()
  
  #only one partition
  cds <- learn_graph(cds, use_partition = T)
  
  pdf(paste0(name, "_after_learn_graph.pdf"))
  plot(plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
                  label_branch_points = T, label_roots = T, label_leaves = F,
                  group_label_size = 5))
  plot(DimPlot(seu, group.by = "timepoint", reduction = "umap.harmony"))
  dev.off()
  
  cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, clusters(cds) == root]))
  pdf(paste0(name, "_ordered_cells.pdf"))
  plot(plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
                  label_branch_points = T, label_roots = F, label_leaves = F))
  dev.off()
  
  cds$monocle3_pseudotime <- pseudotime(cds)
  data.pseudo <- as.data.frame(colData(cds))
  
  pdf(paste0(name, "_clusters_by_pseudo.pdf"))
  plot(ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime), fill = seurat_clusters)) + geom_boxplot())
  dev.off()
  
  #Find genes that change as a function of pseudotime
  deg <- graph_test(cds, neighbor_graph = "principal_graph")
  deg_pass <- deg %>% arrange(q_value) %>% filter(status == "OK")
  write.csv(deg_pass, paste0(name, "_genes_in_pseudo.csv"))
  
  good_genes <- rownames(deg_pass)[1:4]
  
  #plotting a few
  pdf(paste0(name, "_genes_by_pseudo.pdf"))
  plot(FeaturePlot(seu, features = good_genes, reduction = "umap.harmony"))
  dev.off()
  
  #Add pseudotime values into the seuratobject
  seu$level_2_pseudotime <- pseudotime(cds)
  pdf(paste0(name, "_pseudo_on_seurat.pdf"))
  plot(FeaturePlot(seu, features = "level_2_pseudotime", reduction = "umap.harmony" ))
  dev.off()
  
  #plot some genes in pseudotime
  pdf(paste0(name, "_good_genes_pseudo.pdf"))
  my_genes <-good_genes
  cds_subset <- cds[my_genes,]
  plot(plot_genes_in_pseudotime(cds_subset, color_cells_by = "monocle3_pseudotime", label_by_short_name = FALSE))
  dev.off()
  
  #saving seu
  saveRDS(seu, paste0(name, "_pseudo_seu.rds"))
  #saving cds
  saveRDS(cds, paste0(name, "_pseudo_cds.rds"))
  
  #saving just cells and pseudo
  write_csv(data.frame(cells = rownames(seu[[]]),
                       level_2_pseudotime = seu[[]]$level_2_pseudotime), paste0(name, "_cells_in_pseudo.csv"))
  
}

chalazal_pseudotime(cze_pen_dap5, 11, "cze_pen_dap5")  






