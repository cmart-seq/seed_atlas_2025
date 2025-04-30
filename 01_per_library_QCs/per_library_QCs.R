#this code performs the initial QCs and filtering steps of individual libraries after alignment with CellRanger

#example usage: 
#Rscript per_library_QCs.R --wd . --sample DAP3_2 --tenXlibname DAP3_colcol_3 --ccGenes inputs/Menges2005_Mphase_Sphase.csv

library(Seurat)
library(dplyr)
library(SoupX)
library(ggplot2)
library(scCustomize)
library(argparse)
library(optparse)

#reproducible randomness
set.seed(123)

#getting the args
parser <- ArgumentParser()

#add arguments
parser$add_argument("--wd", help = "path to directory containing the 'libraries' and 'output' directories")
parser$add_argument("--sample", help = "sample name, of the form TIMEPOINT_REP (e.g. DAP5_3)")
parser$add_argument("--tenXlibname", help = "name of the root directory for 10X cellranger outs")
parser$add_argument("--ccGenes", help = "list of cell cycle genes")

#parse the command line arguments
args <- parser$parse_args()

#set working directory, with "library" and "output" directories
setwd(file.path(args$wd))

#getting the dataset name
dataname = args$sample

##### helper functions: #####

vanilla_seurat <- function(seurat_object){
  seurat_object <- NormalizeData(seurat_object) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30, verbose = F) %>%
    FindClusters(resolution = 0.8) %>% #using defaults for preprocessing, this will be updated later
    RunUMAP(dims = 1:30, verbose = F)
  return(seurat_object)
}

is_outlier = function(seu_obj, metric, nmads){
  M = seu_obj@meta.data[,metric]
  mad = median(abs(M-mean(M)))
  results = sapply(M, function(x){ifelse( x + nmads*mad < median(M) || x - nmads*mad > median(M), TRUE, FALSE)})
  return(results)
}

##### 1) reading in the data, pointing to the working directory, setting up outputs #####

#making the output directory if it doesn't exist
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
filt.matrix <- Read10X_h5(paste0("libraries/",args$tenXlibname ,"/outs/filtered_feature_bc_matrix.h5"),use.names = T)
seu  <- CreateSeuratObject(counts = filt.matrix)
raw.matrix <- Read10X_h5(paste0("libraries/",args$tenXlibname,"/outs/raw_feature_bc_matrix.h5"),use.names = T)


##### 2) modeling ambient RNA with SoupX ##### 
#create a soup channel
soup.channel  <- SoupChannel(raw.matrix, filt.matrix)

#initializing the object
seu <- vanilla_seurat(seu)

#transferring data for soupX
meta    <- seu@meta.data
umap    <- seu@reductions$umap@cell.embeddings
soup.channel  <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
soup.channel  <- setDR(soup.channel, umap)
#running soupx
soup.channel  <- autoEstCont(soup.channel, doPlot = FALSE)
adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)
DropletUtils:::write10xCounts(paste0(outpath,"/soupx_corrected"), adj.matrix)

#making a new seurat object, removing the uninformative genes, those not detected in at least 10 cells. Removing very low gene cells
mtx = Read10X(data.dir = paste0(outpath,"/soupx_corrected"), gene.column = 1)
seu = CreateSeuratObject(counts = mtx, project = dataname, min.cells = 10, min.features = 250)

##### 3) finding gene and RNA outliers by getting the MADs for each covariate in each cell ##### 
#mark cells as outliers if they differ by 5 MADs
#median absolute deviations: 
#   median(abs(Xi-Xavg))

#adding outlier status to the seurat object
seu@meta.data[,"outlier_RNA"] <- is_outlier(seu, "nCount_RNA", 5)
seu@meta.data[,"outlier_GENES"] <- is_outlier(seu, "nFeature_RNA", 5)

#plotting outliers
pdf(file = paste0(outpath,"/figures/",dataname,"_outliers.pdf"))
FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "outlier_RNA") + ggtitle("RNA outliers")
FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "outlier_GENES") + ggtitle("Gene outliers")
invisible(dev.off())

#filtering the dataset
remove = colnames(seu)[unique(c(which(seu@meta.data$outlier_RNA),which(seu@meta.data$outlier_GENES)))]
orig = ncol(seu)
seu = subset(seu, cells = remove, invert = TRUE)

#rerunning vanilla
seu <- vanilla_seurat(seu)
DefaultAssay(seu) <- "RNA"

#adding the chlor and mito genes as percentage gene sets
seu[["percent.chlor"]] <- PercentageFeatureSet(seu, pattern = "^ATCG")
seu[["percent.mito"]] <- PercentageFeatureSet(seu, pattern = "^ATMG")
#getting cc genes, assiging phase
cc_genes = read.csv(args$ccGenes)
cc_genes$TAIR <-  toupper(cc_genes$TAIR)
s.genes <- intersect(rownames(seu), cc_genes[grep("Sphase", cc_genes$phase), ]$TAIR)
g2m.genes <- intersect(rownames(seu), cc_genes[grep("Mphase", cc_genes$phase), ]$TAIR)

#performing a cell cycle scoring 
seu = CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
seu[["percent.cc"]] <- PercentageFeatureSet(seu, features = union(s.genes, g2m.genes))

#making the variable genes plot to make sure the n_VarGenes is reasonable, this is for the RNA slot.
pdf(file = paste0(outpath,"/figures/",dataname,"_var_genes.pdf"))
VariableFeaturePlot_scCustom(seurat_object = seu, num_features = 30, repel = TRUE)
dev.off()

##### 4) plotting general QCs ##### 
#violins
pdf(file = paste0(outpath,"/figures/",dataname,"_contamination.pdf"), width = 10, height = 6)
plot(VlnPlot(seu, features = c("nFeature_RNA"), 
             ncol = 3, alpha = 0.1, group.by = "seurat_clusters"))
plot(VlnPlot(seu, features = c("nCount_RNA"), 
             ncol = 3, alpha = 0.1, group.by = "seurat_clusters"))
plot(VlnPlot(seu, features = c("percent.mito"), 
             ncol = 3, alpha = 0.1, group.by = "seurat_clusters"))
plot(VlnPlot(seu, features = c("percent.chlor"), 
             ncol = 3, alpha = 0.1, group.by = "seurat_clusters"))
invisible(dev.off())

#scatterplot
pdf(file = paste0(outpath,"/figures/",dataname,"_nRNA_nGene.pdf"))
plot(FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
invisible(dev.off())

pdf(file = paste0(outpath,"/figures/",dataname,"_nRNA_mito.pdf"))
plot(FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mito"))
invisible(dev.off())

#saving
saveRDS(seu, paste0(outpath,"/",dataname,"soupx.rds"))





