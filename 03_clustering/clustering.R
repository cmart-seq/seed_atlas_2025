#this script performs a parameter sweep of the Seurat FindClusters function, then evaluates the clustering using 
#bluster functions that calculate silohuette width and neighborhood purity
#After a clustering resolution is found, the Picard 2021 CZE markers are used in a GSEA to aid cluster annotation, 
#and an initial differential expression analysis is performed

#example usage: 
#Rscript clustering.R --wd . --integrated FALSE --timepoint DAP3 --dataname DAP3 --low_res 1 --high_res 2 --inc 0.1 --cze_marks inputs


library(tidyverse)
library(Seurat)
library(scran)
library(bluster)
library(argparse)
library(clustree)
library(cluster)
library(org.At.tair.db)

set.seed(123)

#getting the args
parser <- ArgumentParser()

parser$add_argument("--wd", help = "working directory")
parser$add_argument("--integrated", help = "is the data integrated? (TRUE/FALSE)")
parser$add_argument("--timepoint", help = "sample name, of the form TIMEPOINT (e.g. DAP5)")
parser$add_argument("--dataname", help = "prefix for outputs")
parser$add_argument("--low_res", help = "lowest resolution", type = "double")
parser$add_argument("--high_res", help = "max resolution", type = "double")
parser$add_argument("--inc", help = "increments for resolution scan", type = "double")
parser$add_argument("--cze_marks", help = "CZE markers path")

# Parse the arguments
args <- parser$parse_args()

#storing args 
wd <- args$wd
integrated <- args$integrated
timepoint <- args$timepoint
dataname <- args$dataname
low_res <- args$low_res
high_res <- args$high_res
inc <- args$inc
cze_marks <- args$cze_marks

#orienting to the right directory
setwd(wd)

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

outpath = file.path(wd, "outputs", dataname)

#loading the seurat
if(integrated){
  seu <- readRDS(paste0("../02_merge_libraries_filtering_and_batch_effects/outputs/", timepoint, "/", timepoint, "_integrated_round2.rds")) 
} else{
  seu <- readRDS(paste0("../02_merge_libraries_filtering_and_batch_effects/outputs/", timepoint, "/", timepoint, "_merged_round2.rds"))
}


#loading the CZE markers for the enrichment analysis
colettes_nodule_marks_E12 <- read.csv(file.path(cze_marks, "colettes_nodule_marks_E12.csv"))
colettes_cyst_marks_E13 <- read.csv(file.path(cze_marks, "colettes_cyst_marks_E13.csv"))
colettes_nodule_like_marks_E14 <- read.csv(file.path(cze_marks, "colettes_nodule_like_marks_E14.csv"))

#setting up dataframes that will hold metrics from the clustering parameter sweep
resolutions <- seq(low_res, high_res, by = inc)  # Define a sequence of resolutions to test
metrics <- data.frame(Resolution = numeric(), NumClusters = numeric(), AvgSilhouetteWidth = numeric(), 
                      NeighborPurity = numeric(), ClusterRMSD = numeric(), pct_bel_zero = numeric(), area_bel_zero = numeric())

purities<- data.frame(purity = numeric(),   
                      maximum = numeric(),
                      resolution = numeric())

#resolution parameter sweep
seu_scan <- seu
for (res in resolutions) {

  seu_scan <- FindClusters(seu_scan, resolution = res)
  clusters <- Idents(seu_scan)
  
  #getting the PCA embeddings used for clustering for scoring
  data <- Embeddings(seu_scan, reduction = "pca")
  
  #silhouette scores
  sil_scores <- approxSilhouette(data, clusters)
  
  #average silhouette width
  avg_sil_width <- mean(sil_scores[, 3])  #third column contains silhouette widths
  
  #neighbor purity, this is what we will use to determine the optimal clustering
  purity <- neighborPurity(data, clusters)
  purity$resolution <- res
  purities <- rbind(purities, data.frame(purity))
  avg_purity <- mean(purity$purity)
  
  #cluster RMSD
  rmsd <- clusterRMSD(data, clusters)
  avg_rmsd <- mean(rmsd)
  
  #plotting the silhouettes 
  seu_scan$sil <- sil_scores[, 3]
  
  #convert silhouette object to data frame
  sil_df <- as.data.frame(seu_scan[[]]) %>% dplyr::select(paste0("RNA_snn_res.",res), sil)
  
  colnames(sil_df) <- c("cluster", "sil")
  
  #getting the "area" under the zero line
  area_below_zero <- sum(abs(sil_df$sil[sil_df$sil < 0]))
  
  #getting ncells under the zero line
  pct_cells_below_zero <- nrow(sil_df[sil_df$sil < 0,])/nrow(sil_df)
  
  # Store the results
  metrics <- rbind(metrics, data.frame(Resolution = res, NumClusters = length(unique(clusters)), 
                                       AvgSilhouetteWidth = avg_sil_width, NeighborPurity = avg_purity, 
                                       ClusterRMSD = avg_rmsd, pct_bel_zero = pct_cells_below_zero, 
                                       area_bel_zero = area_below_zero))
  
  # Ordering samples within clusters by silhouette width
  sil_df <- sil_df[order(sil_df$cluster, -sil_df$sil), ]
  sil_df$order <- seq_len(nrow(sil_df))  # Add sequential order for x-axis
  pdf(paste0(outpath,"/figures/",dataname,"_sils_res_",res,".pdf"), width = 12, height = 6)
  plot(ggplot(sil_df, aes(x = order, y = sil, color = cluster)) +
    geom_bar(stat = "identity", position = "dodge")+
    geom_hline(yintercept = mean(sil_df$sil), color = "red", linetype = "dashed") +
    scale_x_continuous(breaks = NULL) +  # Hide x-axis labels (optional)
    #scale_color_manual(values = c("red", "green", "blue")) +  # Customize cluster colors
    labs(
      title = paste0("Silhouette Plot, res = ",res, ", pct below 0 = ", pct_cells_below_zero, ", area below zero = ", area_below_zero),
      x = "Samples (Grouped by Cluster)",
      y = "Silhouette Width",
      color = "Cluster"
    ) +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank()))  # Remove vertical gridlines
  
  dev.off()
}

#plotting the relationships between the clusters with increasing resolution
pdf(paste0(outpath,"/figures/",dataname,"_clustree.pdf"))
  plot(clustree(seu_scan, prefix = "RNA_snn_res."))
dev.off()

#find the when the neighborhood purity drops off
#calculating slopes between consecutive points
slopes <- diff(metrics$NeighborPurity) / diff(metrics$Resolution)

#finding the index of the maximum slope
max_slope_index <- which.max(abs(slopes))

#extracting the point to the left of the line segment with the greatest slope
point_to_left <- metrics[max_slope_index, ]

#this is our "best resolution" that we will use for the most refined level of annotation (level 3)
best_resolution <- point_to_left$Resolution

#this is our "best resolution" that we will use for the most refined level of annotation (level 3)
write.table(metrics, file = paste0(outpath, "/", dataname, "_all_metrics.txt"), sep = "\t", row.names = FALSE)
metrics <- read.delim(paste0(outpath, "/", dataname, "_all_metrics.txt"))

#plotting 
#the number of clusters vs. average silhouette width
pdf(paste0(outpath,"/figures/",dataname,"_nclusters_v_avgsil.pdf"))
ggplot(metrics, aes(x = NumClusters, y = AvgSilhouetteWidth)) +
  geom_point() +
  geom_line() +
  labs(title = "Number of Clusters vs. Average Silhouette Width",
       x = "Number of Clusters",
       y = "Average Silhouette Width") +
  theme_minimal()
dev.off()

#the resolution vs. average silhouette width
pdf(paste0(outpath,"/figures/",dataname,"_res_v_avgsil.pdf"))
ggplot(metrics, aes(x = Resolution, y = AvgSilhouetteWidth)) +
  geom_point() +
  geom_line() +
  labs(title = "Resolution vs. Average Silhouette Width",
       x = "Resolution",
       y = "Average Silhouette Width") +
  geom_vline(xintercept = best_resolution,              # Vertical line at x = 5
             color = "red",               # Line color
             linetype = "dashed",         # Line type
             size = 1.5)+                  # Line thickness
  theme_minimal()
dev.off()

#the resolution vs. neighbor purity, with a vline at the best resolution
pdf(paste0(outpath,"/figures/",dataname,"_res_v_purity.pdf"))
ggplot(metrics, aes(x = Resolution, y = NeighborPurity)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = best_resolution,              # Vertical line at x = 5
             color = "red",               # Line color
             linetype = "dashed",         # Line type
             size = 1.5)+                  # Line thickness
  labs(title = "Resolution vs. Neighbor Purity",
       x = "Resolution",
       y = "Neighbor Purity") +
  theme_minimal()
dev.off()

#clusters vs. neighbor purity
pdf(paste0(outpath,"/figures/",dataname,"_clusters_v_purity.pdf"))
ggplot(metrics, aes(x = NumClusters, y = NeighborPurity)) +
  geom_point() +
  geom_line() +
  labs(title = "NumClusters vs. Neighbor Purity",
       x = "NumClusters",
       y = "Neighbor Purity") +
  theme_minimal()
dev.off()

#the resolution vs. cluster RMSD
pdf(paste0(outpath,"/figures/",dataname,"_res_v_rmsd.pdf"))
ggplot(metrics, aes(x = Resolution, y = ClusterRMSD)) +
  geom_point() +
  geom_line() +
  labs(title = "Resolution vs. Cluster RMSD",
       x = "Resolution",
       y = "Cluster RMSD") +
  geom_vline(xintercept = best_resolution,              # Vertical line at x = 5
             color = "red",               # Line color
             linetype = "dashed",         # Line type
             size = 1.5)+ 
  theme_minimal()
dev.off()


#now, re-doing the analysis with the best cluster resolution
if(integrated){
  seu <-FindNeighbors(seu, dims = 1:ncol(seu@reductions$pca), verbose = F, reduction = "harmony") %>%
    FindClusters(resolution = best_resolution, cluster.name = paste0("optimal_clusters_res_",best_resolution)) 
  
  seu <-SetIdent(seu, value = paste0("optimal_clusters_res_",best_resolution))
  
  #plotting the final clustering
  pdf(paste0(outpath,"/figures/",dataname,"_optimal_clustering_integrated.pdf"))
  plot(DimPlot(seu, label = TRUE, reduction = "umap.harmony"))
  dev.off()
  
} else{
  seu <-FindNeighbors(seu, dims = 1:ncol(seu@reductions$pca), verbose = F) %>%
    FindClusters(resolution = best_resolution, cluster.name = paste0("optimal_clusters_res_",best_resolution))

  seu <-SetIdent(seu, value = paste0("optimal_clusters_res_",best_resolution))
  
  #plotting the final clustering
  pdf(paste0(outpath,"/figures/",dataname,"_optimal_clustering_merged.pdf"))
  plot(DimPlot(seu, label = TRUE))
  dev.off()
  
}

#performing DE for annotation

DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)
seu <-SetIdent(seu, value = paste0("optimal_clusters_res_",best_resolution))
optimal_cluster_markers <- FindAllMarkers(seu)
optimal_cluster_markers <- mutate(optimal_cluster_markers, diff = pct.1 - pct.2) 
optimal_cluster_markers$timepoint <-dataname
#adding the gene symbol
all_genes <- unique(optimal_cluster_markers$gene)
gene_mapping <- AnnotationDbi::select(org.At.tair.db,
                                      keys = all_genes,
                                      columns = c("SYMBOL"),
                                      keytype = "TAIR")
gene_mapping <- nest(gene_mapping, symbol = SYMBOL)
optimal_cluster_markers <- left_join(optimal_cluster_markers, gene_mapping, by = c("gene" = "TAIR"))
optimal_cluster_markers$symbol <- sapply(optimal_cluster_markers$symbol, function(x){paste0(unlist(x), collapse = ",")})

#saving
write.csv(optimal_cluster_markers, paste0(outpath,"/", dataname,"_optimal_cluster_markers.csv"))


#%cc gene stats
cc_grouping<- seu[[]] %>%
  dplyr::select(c(paste0("optimal_clusters_res_",best_resolution), "Phase")) %>%
  group_by(eval(parse(text = paste0("optimal_clusters_res_",best_resolution))), Phase) %>%
  summarise(CellCount = n(), .groups = "drop")

colnames(cc_grouping) <- c("bestres", "Phase", "CellCount")
cc_grouping %>%
  group_by(bestres) %>%
  mutate(Percentage = (CellCount / sum(CellCount)) * 100)%>%
  slice_max(Percentage, n = 1) %>% write.csv(paste0(outpath,"/", dataname,"_cc_phase.csv"))

#top 10 best markers
markers <- read.csv(paste0(outpath,"/", dataname,"_optimal_cluster_markers.csv"))

best_genes <- markers %>%
  mutate(gene_symbol = paste0(gene, " [", symbol, "]")) %>% 
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  dplyr::slice_head(n=10) %>% 
  ungroup() %>%
  dplyr::select(cluster, gene_symbol) %>% 
  nest(genes = gene_symbol)

#saving
best_genes$genes <- sapply(best_genes$genes, function(x){ paste0(unlist(x), collapse = ", ")})
write.csv(best_genes, paste0(outpath,"/", dataname,"_optimal_cluster_markers_top10.csv"))


#performing the CZE subtype enrichment analysis
#performing a module score analysis
#colettes_nodule_marks_E12, colettes_cyst_marks_E13, colettes_nodule_like_marks_E14

#Adding module score to the Seurat object
colettes_nodule_marks_E12 <- intersect(colettes_nodule_marks_E12$TAIR, rownames(seu))
colettes_cyst_marks_E13 <- intersect(colettes_cyst_marks_E13$TAIR, rownames(seu))
colettes_nodule_like_marks_E14 <- intersect(colettes_nodule_like_marks_E14$TAIR, rownames(seu))

seu <- AddModuleScore(object = seu,
                      features = list(colettes_nodule_marks_E12,
                                      colettes_cyst_marks_E13,
                                      colettes_nodule_like_marks_E14),
                                   name = "chalazal")

#The module score will be stored in the metadata with the prefix "ModuleScore1"
#renaming for clarity
colnames(seu[[]])[grep("chalazal1", colnames(seu[[]]))] <- "Nodule_modscore"
colnames(seu[[]])[grep("chalazal2", colnames(seu[[]]))] <- "Cyst_modscore"
colnames(seu[[]])[grep("chalazal3", colnames(seu[[]]))] <- "NodLike_modscore"

#getting rid of the old metadata
seu@meta.data <- seu@meta.data[, !grepl("chalazal", colnames(seu@meta.data))]

res <- grep("optimal_clusters_res_", colnames(seu[[]]), value = TRUE)
seu <-SetIdent(seu, value = res)

pdf(paste0(outpath,"/figures/",dataname,"_CZE_nodule_enrichment.pdf"), width = 10, height = 5)
VlnPlot(seu, "Nodule_modscore") +ggtitle("nodule")
dev.off()

pdf(paste0(outpath,"/figures/",dataname,"_CZE_cyst_enrichment.pdf"), width = 10, height = 5)
VlnPlot(seu, "Cyst_modscore")+ggtitle("cyst")
dev.off()

pdf(paste0(outpath,"/figures/",dataname,"_CZE_nodlike_enrichment.pdf"), width = 10, height = 5)
VlnPlot(seu, "NodLike_modscore")+ggtitle("nodule-like")
dev.off()

saveRDS(seu, paste0(outpath,"/", dataname,"_clustered.rds"))


