library(Seurat)
library(scDblFinder)
library(dplyr)
library(SoupX)
library(ggplot2)
library(scCustomize)
library(argparse)
library(optparse)
library(stringr)

#reproducible randomness
set.seed(123)


#cd /lab/solexa_gehring/carly/snRNA_seed_dev/clean_analysis/01_per_library_QCs
#sbatch -p 20 --mem=24gb --mail-type=ALL --wrap 'Rscript score_doublets.R --wd . --seu_path outputs/DAP3_1a/DAP3_1asoupx.rds --sample DAP3_1a --rm_clusters'


#getting the args
parser <- ArgumentParser()

parser$add_argument("--wd", help = "working directory") #"/lab/solexa_gehring/carly/snRNA_seed_dev/clean_analysis/02_merge_libraries_filtering_and_batch_effects"
parser$add_argument("--seu_path", help = "seurat") #
parser$add_argument("--sample", help = "prefix") #
parser$add_argument("--rm_clusters", help = "clusters to drop", nargs = "+", type = "double") #

# Parse the arguments
args <- parser$parse_args()

#setting the directory
working = args$wd
setwd(working)
dataname = args$sample

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

#reading in the data 
seu <- readRDS(args$seu_path)

#removing visibly low quality clusters
if(!is.null(args$rm_clusters)){
  seu <- SetIdent(seu, value = "seurat_clusters")
  seu <- subset(seu, subset = seurat_clusters %in% args$rm_clusters, invert = TRUE) 
}

#paste0(outpath,"/",dataname,"_soupx_doublets.rds")
##### Doublet detection on filtered dataset with scDblFinder ##### 
#preliminary vanilla seurat analysis for doublet detection
capture <- capture.output(sce <- scDblFinder(as.SingleCellExperiment(seu), clusters=FALSE), type = "message")
threshhold <- as.numeric(substring(capture[8], 17))

#saving scores
seu$doublet_score = sce$scDblFinder.score
seu$doublet_class = sce$scDblFinder.class

#removing doublets
seu <- subset(seu, subset = doublet_score < threshhold)

#saving object
saveRDS(seu, paste0(outpath,"/",dataname,"_soupx_doublets.rds"))
seu<-readRDS(paste0(outpath,"/",dataname,"_soupx_doublets.rds"))





