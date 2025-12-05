#performing pseudotime on 5 DAP PEN and CZE to define the apical-basal CZE trajectory
library(monocle3)
library(Seurat)
library(dplyr)
library(SeuratWrappers)
library(ggplot2)
library(readr)
library(harmony)

set.seed(123)

#loading the data
dap5 <- readRDS("../../../04_manual_annotation/review_embryo_subclustering/outputs/dap5_annotation_reviews.rds")
dap3 <- readRDS("../../../04_manual_annotation/review_embryo_subclustering/outputs/dap3_annotation_reviews.rds")

#getting ready, not integrating
pseudo_prep <- function(seurat_object, res, npcs){
  seurat_object <- NormalizeData(seurat_object) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2500) %>%
    #using the cell cycle scoring from the orignal analysis, regressing out to hopefully get a simpler pseudotime trajectory
    ScaleData(vars.to.regress = c("S.Score", "G2M.Score"), assay = "RNA") %>%
    RunPCA(npcs = npcs) %>%
    FindNeighbors(dims = 1:npcs, verbose = F) %>%
    FindClusters(resolution = res) %>% 
    RunUMAP(dims = 1:npcs, verbose = F)
  return(seurat_object)
}

#preparing the 5DAP object for this analysis
cze_pen_dap5 <- subset(dap5, subset = level_2_annotation %in% c("CZE", "PEN"))
cze_pen_dap3 <- subset(dap3, subset = level_2_annotation %in% c("CZE", "PEN"))

cze_pen_dap <- merge(cze_pen_dap3, y = cze_pen_dap5)
cze_pen_dap <- JoinLayers(cze_pen_dap)
cze_pen_dap <- pseudo_prep(cze_pen_dap, 0.5, 30)

pdf("elbow_plot35.pdf")
ElbowPlot(cze_pen_dap, reduction = "pca", ndims = 30)
dev.off()

#rerunning with fewer PCs
cze_pen_dap <- pseudo_prep(cze_pen_dap, 0.5, 20)

pdf("cze_pen_dap35_preharmony.pdf")
plot(DimPlot(cze_pen_dap, label = TRUE, reduction = "umap", group.by = c("bio_rep"), label.size = 2)) +ggtitle("CZE and PEN 5DAP")
plot(DimPlot(cze_pen_dap, label = TRUE, reduction = "umap", group.by = c("level_3_annotation_abbr"), label.size = 2))+ggtitle("CZE and PEN 5DAP")
dev.off()

#running Harmony across bio_reps, 5DAP cze and pen
cze_pen_dap <- RunHarmony(object = cze_pen_dap, reduction = "pca", group.by.vars = c("timepoint"), 
                         dims.use = 1:10, reduction.save = 'harmony', plot_convergence = F) %>%
  FindNeighbors(reduction = "harmony", dims = 1:10) %>% 
  FindClusters(resolution = 0.5, cluster.name = "harmony_clusters") %>% 
  RunUMAP(reduction = "harmony", dims = 1:10, reduction.name = "umap.harmony")

pdf("cze_pen_dap35_postharmony.pdf")
plot(DimPlot(cze_pen_dap, label = TRUE, reduction = "umap.harmony", group.by = c("timepoint"), label.size = 2)) +ggtitle("CZE and PEN 5DAP")
plot(DimPlot(cze_pen_dap, label = TRUE, reduction = "umap.harmony", group.by = c("harmony_clusters"), label.size = 2))+ggtitle("CZE and PEN 5DAP")
plot(DimPlot(cze_pen_dap, label = TRUE, reduction = "umap.harmony", group.by = c("level_3_annotation_abbr"), label.size = 2))+ggtitle("CZE and PEN 5DAP")
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

#anchoring in just the 3 DAP PEN
#anchoring to PEN
chalazal_pseudotime(cze_pen_dap, 0, "cze_pen_3dappen")  

#anchoring to basal cyst
chalazal_pseudotime(cze_pen_dap, 9, "cze_pen_dap35")  

#reversing
cze_pen_dap$rev_clusters <- sapply(cze_pen_dap$harmony_clusters, function(x){ifelse(x %in% c(11, 9, 6,0), x, 1)})
cze_pen_dap <- SetIdent(cze_pen_dap, value = "rev_clusters")
chalazal_pseudotime(cze_pen_dap, 1, "cze_pen_daprev35")  

pdf("RALFL3_WYO_35.pdf")
FeaturePlot(cze_pen_dap, "AT1G23147", reduction = "umap.harmony") + scale_color_viridis_c(option = "rocket", direction = -1)
FeaturePlot(cze_pen_dap, "AT3G49307", reduction = "umap.harmony") + scale_color_viridis_c(option = "rocket", direction = -1)
dev.off()





