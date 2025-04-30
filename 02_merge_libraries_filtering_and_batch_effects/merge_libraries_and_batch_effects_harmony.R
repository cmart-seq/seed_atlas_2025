#this script merges libraries by timepoint and is run in two rounds in order to filter low-quality clusters and tune dimensionality reduction parameters
#round 1: merging libraries, adding metadata and setting variable genes, optionally integrating, plotting QCs
#round 2: filtering low quality clusters based on round 1, setting variable genes, optionally integrating, plotting QCs

#In Martin et al 2025, the only dataset that required batch correction is DAP5, all other timepoints did not show clustering based on replicate

#example usage:
#Rscript merge_libraries_and_batch_effects_harmony.R --wd . --round 1 --merge_only TRUE --dataset DAP3 --libraries DAP3_1a DAP3_1b DAP3_2 --n_VarGenes 3000
#Rscript merge_libraries_and_batch_effects_harmony.R --wd . --round 2 --merge_only TRUE --dataset DAP3 --n_VarGenes 3000 --ndims 60

library(Seurat)
library(stringr)
library(dplyr)
library(argparse)
library(scCustomize)
library(ggplot2)
library(harmony)

#reproducible randomness
set.seed(123)

#getting the args
parser <- ArgumentParser()
parser$add_argument("--wd", help = "working directory")
parser$add_argument("--round", help = "before (1) or after (2) inspection of initial QCs", type = "double")
parser$add_argument("--merge_only", help = "Perform merging, not integration (TRUE/FALSE)")
parser$add_argument("--int_variable", help = "meta data column that will be used to split the object for integration")
parser$add_argument("--dataset", help = "sample name, of the form TIMEPOINT (e.g. DAP5)")
parser$add_argument("--libraries", nargs = "+", help = "list of libraries for normalization and integration")
parser$add_argument("--n_VarGenes", help = "number of variable genes to detect", type = "double")
parser$add_argument("--ndims", help = "number of PCs to generate based on jackstraw analysis", type = "double", nargs = "+")

#parse the arguments
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

if (!dir.exists(file.path("outputs", dataname))) {
  dir.create(file.path("outputs", dataname), recursive = TRUE)
}

if (!dir.exists(file.path("outputs", dataname, "figures"))) {
  dir.create(file.path("outputs", dataname, "figures"), recursive = TRUE)
}

outpath = file.path(args$wd, "outputs", dataname)

##### 1) round 1: merging libraries, adding metadata setting variable genes, optionally integrating, plotting QCs ##### 

if(args$round == 1 ){
  #merging libraries
  libraries = args$libraries
  
  seu_list = list()
  
  for (i in 1:length(libraries)){
    seu_list[i] <- readRDS(paste0("../01_per_library_QCs/outputs/", libraries[i], "/", libraries[i] ,"_soupx_doublets.rds"))
  }
  #merging, setting to RNA at this step.  
  seu = merge(seu_list[[1]], y = seu_list[-1])
  DefaultAssay(seu) <- "RNA"
  seu <- JoinLayers(seu)
  
  #adding some metadata
  #biological replicate
  seu$bio_rep <- sapply(seu$orig.ident, function(x){
    q <- str_split(x, pattern = "_")[[1]][2]
    q <- gsub("[^0-9.-]", "", q)
    t <- str_split(x, pattern = "_")[[1]][1]
    return(paste0(t, "_", q))
  })
  
  #timepoint
  seu$timepoint <- sapply(seu$orig.ident, function(x){
    q <- str_split(x, pattern = "_")[[1]][1]
    return(q)
  })
  
  #Beginning with the RNA slot
    seu <- NormalizeData(seu) %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = args$n_VarGenes) %>%
      ScaleData(features = rownames(seu)) %>%
      RunPCA(npcs=30) %>%
      RunUMAP(dims = 1:30) %>%
      FindNeighbors(dims = 1:30, verbose = F) %>%
      FindClusters(resolution = 0.8, cluster.name = "round1_clusters_merged")  
  
  saveRDS(seu, paste0(outpath,"/", dataname,"_merged_only_round1.rds"))
  seu<-readRDS(paste0(outpath,"/", dataname,"_merged_only_round1.rds")) 
  
  #making the variable genes plot to make sure the n_VarGenes is reasonable, this is for the RNA slot.
  pdf(file = paste0(outpath,"/figures/",dataname,"_var_genes_merged_only_round1.pdf"))
  plot(VariableFeaturePlot_scCustom(seurat_object = seu, num_features = 30, repel = TRUE))
  dev.off()
  
  pdf(paste0(outpath,"/figures/",dataname,"_merged_only_clusters_round1.pdf"))
  plot(DimPlot(seu, label = T, repel = T)) 
  plot(DimPlot(seu, label = T, repel = T, group.by = "round1_clusters_merged")) 
  plot(DimPlot(seu, label = T, repel = T, group.by = "bio_rep")) 
  dev.off()
  
  #now, plotting QCs on the filtered clusters to identify obviously low quality cells 
  pdf(paste0(outpath,"/figures/",dataname,"_cluster_QCs_merged_only_round1.pdf"), width = 12, height = 6)
  
  plot(VlnPlot(object = seu, features = 'nCount_RNA', alpha = 0.1))
  plot(VlnPlot(object = seu, features = 'nFeature_RNA', alpha = 0.1))
  plot(VlnPlot(object = seu, features = 'percent.chlor', alpha = 0.1))
  plot(VlnPlot(object = seu, features = 'percent.mito', alpha = 0.1))
  dev.off()
  
  if(args$merge_only == FALSE){
    
    #splitting for integration
    DefaultAssay(seu)<- "RNA"
    
    #running Harmony across timepoints
    seu <- RunHarmony(object = seu, reduction = "pca", group.by.vars = args$int_variable, 
                        dims.use = 1:30, reduction.save = 'harmony', plot_convergence = F) %>%
      FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
      FindClusters(resolution = 0.8, cluster.name = "round1_clusters_integrated") %>% 
      RunUMAP(reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
    
    pdf(paste0(outpath,"/figures/",dataname,"_integrated_clusters_initial.pdf"))
    plot(DimPlot(seu, label = T, repel = T, reduction = "umap.harmony") )
    plot(DimPlot(seu, label = T, repel = T, group.by = "round1_clusters_integrated", reduction = "umap.harmony") )
    plot(DimPlot(seu, label = T, repel = T, group.by = args$int_variable, reduction = "umap.harmony") )
    dev.off()
    
    #saving 
    seu <- SetIdent(seu, value = "round1_clusters_integrated")
    saveRDS(seu, paste0(outpath,"/", dataname,"_integrated_round1.rds"))
    
    #now, plotting QCs on the filtered clusters to identify obviously low quality cells
    pdf(paste0(outpath,"/figures/",dataname,"_cluster_QCs_integrated_round1.pdf"), width = 12, height = 6)
    
    plot(VlnPlot(object = seu, features = 'nCount_RNA', alpha = 0.1))
    plot(VlnPlot(object = seu, features = 'nFeature_RNA', alpha = 0.1))
    plot(VlnPlot(object = seu, features = 'percent.chlor', alpha = 0.1))
    plot(VlnPlot(object = seu, features = 'percent.mito', alpha = 0.1))
    dev.off()
    
    
  }
  
  #performing jackstraw to select appropriate PCs
  ndi <- dim(seu[["pca"]])[2]
  
  #testing
  seu <- JackStraw(seu, num.replicate = 70, verbose = TRUE, dims = ndi)
  
  seu <- ScoreJackStraw(seu, dims = 1:ndi)
  
  #plotting
  n_20s <- floor(ndi/20)
  remainders<- ndi-(n_20s*20)
  
  pdf(paste0(paste0(outpath,"/figures/",dataname,"_nPCs_elbow_jackstraw.pdf")))
  plot(ElbowPlot(seu, ndims = ndi))
  for(i in 1:n_20s){
    b <- i*20-20+1
    e <- i*20
    plot(JackStrawPlot(seu, dims = seq(b, e, 1)))
  }
  b <- ndi - remainders + 1
  e <- ndi
  plot(JackStrawPlot(seu, dims = seq(b, e, 1)))
  dev.off()
  
  ##### 2) round 2: filtering low quality clusters based on round 1, setting variable genes, optionally integrating, plotting QCs##### 
  
  } else if(args$round == 2 ){
  #filtering low quality cells/clusters and performing integration
    library(stringr)
    shortname <- str_split(dataname, "_")[[1]][1]
    shortpath <- file.path(args$wd, "outputs", shortname)
  #loading the seurat
  if(args$merge_only == FALSE){
    seu <- readRDS(paste0(shortpath,"/", shortname,"_integrated_round1.rds")) 
  } else{
    seu <- readRDS(paste0(shortpath,"/", shortname,"_merged_only_round1.rds"))
  }
    
  if(sum(grepl("round1_clusters_integrated", colnames(seu[[]]))) == 1){
    #if the dataset was integrated, set the integrated clusters as the ident
    seu <- SetIdent(seu, value = "round1_clusters_integrated")
  } else {
    seu <- SetIdent(seu, value = "round1_clusters_merged")
  }
  
  
  #beginning by re-analyzing the RNA assay
  DefaultAssay(seu) <- "RNA"
  seu <- JoinLayers(seu)
  
  #analyzing the merged object, scaling and regressing out cc genes if need be
    seu <- NormalizeData(seu) %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = args$n_VarGenes) %>%
      ScaleData(features = rownames(seu)) %>%
      RunPCA(npcs=args$ndims) %>%
      RunUMAP(dims = 1:args$ndims) %>%
      FindNeighbors(dims = 1:args$ndims, verbose = F) %>%
      FindClusters(resolution = 0.8, cluster.name = "round2_clusters_merged")
    
  saveRDS(seu, paste0(outpath,"/", dataname,"_merged_round2.rds"))
  seu<-readRDS(paste0(outpath,"/", dataname,"_merged_round2.rds")) 
  
  pdf(paste0(outpath,"/figures/",dataname,"_merged_clusters_round2.pdf"))
  plot(DimPlot(seu, label = T, repel = T)) 
  if(args$merge_only == FALSE){
    plot(DimPlot(seu, label = T, repel = T, group.by = args$int_variable))
  } 
  plot(DimPlot(seu, label = T, repel = T, group.by = "round2_clusters_merged")) 
  plot(DimPlot(seu, label = T, repel = T, group.by = "bio_rep")) 
  dev.off()
  
  if(args$merge_only == FALSE){
    #splitting for integration
    DefaultAssay(seu)<- "RNA"
    seu <- JoinLayers(seu) 
    #running Harmony across timepoints
    seu <- RunHarmony(object = seu, reduction = "pca", group.by.vars = args$int_variable, 
                      dims.use = 1:args$ndims, reduction.save = 'harmony', plot_convergence = F) %>%
      FindNeighbors(reduction = "harmony", dims = 1:args$ndims) %>% 
      FindClusters(resolution = 0.8, cluster.name = "round2_clusters_integrated") %>% 
      RunUMAP(reduction = "harmony", dims = 1:args$ndims, reduction.name = "umap.harmony")
    
    pdf(paste0(outpath,"/figures/",dataname,"_integrated_clusters_round2.pdf"))
    plot(DimPlot(seu, label = T, repel = T, reduction = "umap.harmony") )
    plot(DimPlot(seu, label = T, repel = T, group.by = "round2_clusters_integrated", reduction = "umap.harmony") )
    plot(DimPlot(seu, label = T, repel = T, group.by = args$int_variable, reduction = "umap.harmony") )
    dev.off()
    
    #saving 
    seu <- SetIdent(seu, value = "round2_clusters_integrated")
    saveRDS(seu, paste0(outpath,"/", dataname,"_integrated_round2.rds"))
    
  }

  #now, plotting QCs on the filtered clusters
  pdf(paste0(outpath,"/figures/",dataname,"_cluster_QCs_round2.pdf"), width = 12, height = 6)
  
  plot(VlnPlot(object = seu, features = 'nCount_RNA', alpha = 0.1))
  plot(VlnPlot(object = seu, features = 'nFeature_RNA', alpha = 0.1))
  plot(VlnPlot(object = seu, features = 'percent.chlor', alpha = 0.1))
  plot(VlnPlot(object = seu, features = 'percent.mito', alpha = 0.1))
  
  dev.off()
  
}
