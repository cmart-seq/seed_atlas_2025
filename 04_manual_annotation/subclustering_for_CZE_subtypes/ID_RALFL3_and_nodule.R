#below is the code used to identify CZE subtypes in the 3 and 5 DAP datasets

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

setwd("04_manual_annotation/subclustering_for_CZE_subtypes")

#reading in data
dap3 <- readRDS("../../04_manual_annotation/outputs/DAP3/DAP3_annotated.rds")
DefaultAssay(dap3) <- "RNA"
dap3 <- JoinLayers(dap3)
dap5 <- readRDS("../../04_manual_annotation/outputs/DAP5/DAP5_annotated.rds")
DefaultAssay(dap5) <- "RNA"
dap5 <- JoinLayers(dap5)

#plotting de novo clusters

png(paste0("3DAP_de_novo_clusters_umap.png"), width = 4, height = 4, res = 480, units = "in")
DimPlot(dap3, group.by = "seurat_clusters", label = TRUE)+theme(legend.position = "none") +ggtitle("de novo clusters")
dev.off()

png(paste0("5DAP_de_novo_clusters_umap.png"), width = 4, height = 4, res = 480, units = "in")
DimPlot(dap5, group.by = "seurat_clusters", reduction = "umap.harmony", label = TRUE)+theme(legend.position = "none")+ggtitle("de novo clusters")
dev.off()

#plotting RALFL3 and WYO 
png(paste0("RALFL3_3DAP_de_novo_clusters_umap.png"), width = 4, height = 4, res = 480, units = "in")
FeaturePlot(dap3, "AT1G23147", alpha = 0.5)+ggtitle("RALFL3")
dev.off()

png(paste0("WYO_3DAP_de_novo_clusters_umap.png"), width = 4, height = 4, res = 480, units = "in")
FeaturePlot(dap3, "AT3G49307", alpha = 0.5)+ggtitle("WYO")
dev.off()

png(paste0("RALFL3_5DAP_de_novo_clusters_umap.png"), width = 4, height = 4, res = 480, units = "in")
FeaturePlot(dap5, "AT1G23147", reduction = "umap.harmony", alpha = 0.5)+ggtitle("RALFL3")
dev.off()

png(paste0("WYO_5DAP_de_novo_clusters_umap.png"), width = 4, height = 4, res = 480, units = "in")
FeaturePlot(dap5, "AT3G49307", alpha = 0.5, reduction = "umap.harmony" )+ggtitle("WYO")
dev.off()

png(paste0("3DAP_cyst_module_score.png"), width = 4, height = 4, res = 480, units = "in")
FeaturePlot(dap3, "Cyst_modscore", cols = c("lightgrey", "#ce0d48"), alpha = 0.5)+ggtitle("Cyst module score")
dev.off()

png(paste0("5DAP_cyst_module_score.png"), width = 4, height = 4, res = 480, units = "in")
FeaturePlot(dap5, "Cyst_modscore", cols = c("lightgrey", "#ce0d48"),reduction = "umap.harmony", alpha = 0.5)+ggtitle("Cyst module score")
dev.off()

#subclustering out the 3DAP RALFL3+ nuclei in 3 and 5 dap. Subsetting both datasets to just the endosperm and whichever cluster has RALFL3
dap3_sub <- subset(dap3, subset = level_2_annotation %in% c("CZE", "MCE", "PEN", "CPT", "ii1"))
dap5_sub <- subset(dap5, subset = level_2_annotation %in% c("CZE", "MCE", "PEN", "oi2"))

#re-analyzing
dap3_sub <- NormalizeData(dap3_sub) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(dap3_sub)) %>%
  RunPCA(npcs=30) %>%
  RunUMAP(dims = 1:30) %>%
  FindNeighbors(dims = 1:30, verbose = F) %>%
  FindClusters(resolution = 2)

png(paste0("DAP3_subclustered.png"), width = 5, height = 5, res = 480, units = "in")
DimPlot(dap3_sub, label = TRUE) +theme(legend.position = "none")
dev.off()

#performing a ralf+ labeling
png("RALFL3_expressing.png", width = 5, height = 5, res = 480, units = "in")
DimPlot(dap3_sub,cells.highlight = WhichCells(dap3_sub, expression = AT1G23147 > 2))
dev.off()

#performing a ralf+ labeling
png("RALFL3_expressing_dap3.png", width = 5, height = 5, res = 480, units = "in")
FeaturePlot(dap3_sub,"AT1G23147", order = TRUE)
dev.off()

dap5_sub <- NormalizeData(dap5_sub) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
  ScaleData(features = rownames(dap5_sub)) %>%
  RunPCA(npcs=30)%>%
  RunHarmony(reduction = "pca", group.by.vars = "bio_rep", 
                  dims.use = 1:30, reduction.save = 'harmony', plot_convergence = F) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 1.5) %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

png("DAP5_subclustered.png", width = 5, height = 5, res = 480, units = "in")
DimPlot(dap5_sub, label = TRUE, reduction = "umap.harmony")+theme(legend.position = "none")
dev.off()

#performing a ralf+ labeling
png("RALFL3_expressing_dap5.png", width = 5, height = 5, res = 480, units = "in")
FeaturePlot(dap5_sub,"AT1G23147", order = TRUE, reduction = "umap.harmony")
dev.off()

saveRDS(dap3_sub, "dap3_RALFL3_subclustered.rds")
dap3_sub <- readRDS("dap3_RALFL3_subclustered.rds")
saveRDS(dap5_sub, "dap5_RALFL3_subclustered.rds")
dap5_sub <- readRDS("dap5_RALFL3_subclustered.rds")

#cze subtype markers
cz_genes <- c("AT3G49307", "AT1G23147", "AT3G03260", "AT2G44240", "AT5G10440","AT1G02580", "AT1G27040")
names(cz_genes) <- c("AT3G49307", "RALFL3", "HDG8","AT2G44240","CYCD4;2","MEA", "NPF4.5")

dap3_sub$RNA_snn_res.2 <- factor(dap3_sub$RNA_snn_res.2, levels = seq(0, 24, 1))

#plotting the markers validated by HCR and in situ
png("RALFL3_markers_3DAP_subclustered.png", width = 9.5, height = 7, res = 480, units = "in")
p<- DotPlot(dap3_sub, cz_genes, scale = FALSE)+ 
  scale_colour_gradientn(colours=brewer.pal(11,"YlOrRd"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        panel.grid.major = element_line(color = "grey"), panel.grid.minor = element_line(color = "grey"))

# Extract plot data
dot_data <- p$data

# Overlay points with black outlines
p + geom_point(data = dot_data, aes(size = pct.exp), shape = 21, fill = NA, color = "black", stroke = 0.5)

dev.off()

dap5_sub$RNA_snn_res.1.5 <- factor(dap5_sub$RNA_snn_res.1.5, levels = seq(0, 24, 1))

png("RALFL3_markers_5DAP_subclustered.png", width = 9.5, height = 7, res = 480, units = "in")
p<- DotPlot(dap5_sub, cz_genes, scale = FALSE)+ 
  scale_colour_gradientn(colours=brewer.pal(11,"YlOrRd"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        panel.grid.major = element_line(color = "grey"), panel.grid.minor = element_line(color = "grey"))

# Extract plot data
dot_data <- p$data

# Overlay points with black outlines
p + geom_point(data = dot_data, aes(size = pct.exp), shape = 21, fill = NA, color = "black", stroke = 0.5)

dev.off()

#isolating the clusters of interest
#3DAP: 
dap3_ralfl3 <- WhichCells(dap3_sub, ident = 22)
dap3_nod <- WhichCells(dap3_sub, ident = 18)
dap3_nodlike <- WhichCells(dap3_sub, ident = 17) 

#5DAP: 
dap5_ralfl3 <- WhichCells(dap5_sub, ident = 23)
dap5_nodlike <- WhichCells(dap5_sub, ident = 18) 
dap5_nod <- WhichCells(dap5_sub, ident = 20)
dap5_cyst <- WhichCells(dap5_sub, ident = 19)

#renaming the cell types in the DAP3 object. I will send it through the re-annotation after this
dap3 <- SetIdent(dap3, value = "optimal_clusters_res_1.3")
dap3 <- SetIdent(dap3, cells = dap3_ralfl3, value = 31)
dap3 <- SetIdent(dap3, cells = dap3_nodlike, value = 32)
dap3 <- SetIdent(dap3, cells = dap3_nod, value = 33)
dap3$optimal_clusters_res_1.3_czesubtypes <- Idents(dap3)

#renaming the cell types in the DAP5 object.
dap5 <- SetIdent(dap5, value = "optimal_clusters_res_1.2")
dap5 <- SetIdent(dap5, cells = dap5_ralfl3, value = 29)
dap5 <- SetIdent(dap5, cells = dap5_nodlike, value = 30)
dap5 <- SetIdent(dap5, cells = dap5_nod, value = 31)
dap5 <- SetIdent(dap5, cells = dap5_cyst, value = 32)
dap5$optimal_clusters_res_1.2_czesubtypes <- Idents(dap5)

#updating the seurats
colnames(dap3[[]])[which(colnames(dap3[[]]) == "optimal_clusters_res_1.3")] = "original_clusters_res_1.3"
dap3$optimal_clusters_res_1.3 <- NULL
colnames(dap5[[]])[which(colnames(dap5[[]]) == "optimal_clusters_res_1.2")] = "original_clusters_res_1.2"
dap5$optimal_clusters_res_1.2 <- NULL

saveRDS(dap3, "DAP3_clustered_wcze_subs.rds")
saveRDS(dap5, "DAP5_clustered_wcze_subs.rds")

#plotting to make sure the marker genes are correct

#checking
pdf("DAP3_new_cze_types_highlighted.pdf", width = 7, height = 6)
DimPlot(dap3, group.by = "optimal_clusters_res_1.3_czesubtypes", label = T)
DimPlot(dap3, cells.highlight = dap3_ralfl3)
DimPlot(dap3, cells.highlight = dap3_nodlike)
DimPlot(dap3, cells.highlight = dap3_nod)
dev.off()

pdf("DAP5_new_cze_types_highlighted.pdf", width = 5, height = 6)
DimPlot(dap5, group.by = "optimal_clusters_res_1.2_czesubtypes")
DimPlot(dap5, cells.highlight = dap5_ralfl3)
DimPlot(dap5, cells.highlight = dap5_nodlike)
DimPlot(dap5, cells.highlight = dap5_nod)
DimPlot(dap5, cells.highlight = dap5_cyst)
dev.off()


