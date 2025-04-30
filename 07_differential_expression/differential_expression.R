#this script performs differential expression analysis using the Seurat FindAllMarkers function for all levels within a timepoint dataset
#example usage: 
#Rscript differential_expression.R --wd .  --sample DAP3

library(Seurat)
library(dplyr)
library(readr)
library(argparse)
library(org.At.tair.db)
library(tidyr)

set.seed(123)

#getting the args
parser <- ArgumentParser()

#add arguments
parser$add_argument("--wd", help = "path to directory containing the 'libraries' and 'output' directories") # /lab/solexa_gehring/carly/snRNA_seed_dev/clean_analysis/07_differential_expression/
parser$add_argument("--sample", help = "sample name, of the form TIMEPOINT_REP (e.g. DAP5_3)") #

#Parse the command line arguments
args <- parser$parse_args()

#Set working directory, with "library" and "output" directories
setwd(file.path(args$wd))

#getting the dataset name
dataname = args$sample

#example usage: 
#sbatch -p 20 --mem=24gb --mail-type=ALL --job-name DAP3 --wrap "Rscript differential_expression.R --wd .  --sample DAP3"
#sbatch -p 20 --mem=24gb --mail-type=ALL --job-name DAP5 --wrap "Rscript differential_expression.R --wd .  --sample DAP5"
#sbatch -p 20 --mem=24gb --mail-type=ALL --job-name DAP7 --wrap "Rscript differential_expression.R --wd .  --sample DAP7"
#sbatch -p 20 --mem=24gb --mail-type=ALL --job-name ATLAS --wrap "Rscript differential_expression.R --wd .  --sample ATLAS"

##### main code: #####

#1) reading in the data, pointing to the working directory, setting up outputs

setwd(args$wd)

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

#reading in the data 
if(dataname == "ATLAS"){
  seu <- readRDS("../06_atlas_merging/outputs/ATLAS_integrated_annotated_rpca.rds") 
  DefaultAssay(seu) <- "RNA"
  seu <- JoinLayers(seu)
} else if(dataname == "DAP3"){
  seu <- readRDS("../04_manual_annotation/outputs/DAP3_wcze_subs/DAP3_wcze_subs_annotated.rds")
  DefaultAssay(seu) <- "RNA"
  seu <- JoinLayers(seu)
} else if (dataname == "DAP5") {
  seu <- readRDS("../04_manual_annotation/outputs/DAP5_wcze_subs/DAP5_wcze_subs_annotated.rds") 
  DefaultAssay(seu) <- "RNA"
  seu <- JoinLayers(seu)
}else{
  seu <- readRDS(paste0("../04_manual_annotation/outputs/", dataname, "/", dataname, "_annotated.rds")) 
  DefaultAssay(seu) <- "RNA"
  seu <- JoinLayers(seu)
}

#cleaning up the names
#gsub("/", "_", string)

#performing DE

if(dataname == "ATLAS"){
  for(i in c("level_1_annotation","level_2_annotation", "level_1_annotation_timed", "level_2_annotation_timed_final")){
    #cleaning up names
    seu$level_2_annotation <- sapply(seu$level_2_annotation, function(x){gsub("/", "_", x)})
    seu$level_2_annotation_timed_final <- sapply(seu$level_2_annotation_timed_final, function(x){gsub("/", "_", x)})
    
    #finding markers
    seu <- SetIdent(seu, value = i)
    all_marks <- FindAllMarkers(seu)
    all_marks$diff <- all_marks$pct.1 - all_marks$pct.2
    
    #adding the symbol
    all_genes <- unique(all_marks$gene)
    gene_mapping <- AnnotationDbi::select(org.At.tair.db,
                                          keys = all_genes,
                                          columns = c("SYMBOL"),
                                          keytype = "TAIR")
    gene_mapping <- nest(gene_mapping, symbol = SYMBOL)
    all_marks <- left_join(all_marks, gene_mapping, by = c("gene" = "TAIR"))
    all_marks$symbol <- sapply(all_marks$symbol, function(x){paste0(unlist(x), collapse = ",")})
    
    write_csv(all_marks, paste0("outputs/", dataname,"/", dataname,"_",i,"_markers.csv"))
    #getting the top markers for each cluster, plotting
    for(x in unique(all_marks$cluster)){
      clust_sub <- filter(all_marks, cluster == x) %>% 
        filter(avg_log2FC > 0) %>% 
        arrange(desc(diff))
      pdf(paste0("outputs/", dataname, "/figures/", dataname,"_cluster_",i, "_", x, "_top_marks.pdf"))
      for(y in 1:10){
        plot(FeaturePlot(seu, clust_sub$gene[y]))
      }
      dev.off()
    }
  }
  
} else{
  for(i in c("level_1_annotation","level_2_annotation","level_3_annotation_full")){
    #cleaning up names
    seu$level_2_annotation <- sapply(seu$level_2_annotation, function(x){gsub("/", "_", x)})
    seu$level_3_annotation_full <- sapply(seu$level_3_annotation_full, function(x){gsub("/", "_", x)})
    
    #finding markers
    seu <- SetIdent(seu, value = i)
    all_marks <- FindAllMarkers(seu)
    all_marks$diff <- all_marks$pct.1 - all_marks$pct.2
    
    #adding the symbol
    all_genes <- unique(all_marks$gene)
    gene_mapping <- AnnotationDbi::select(org.At.tair.db,
                                          keys = all_genes,
                                          columns = c("SYMBOL"),
                                          keytype = "TAIR")
    gene_mapping <- nest(gene_mapping, symbol = SYMBOL)
    all_marks <- left_join(all_marks, gene_mapping, by = c("gene" = "TAIR"))
    all_marks$symbol <- sapply(all_marks$symbol, function(x){paste0(unlist(x), collapse = ",")})
    
    
    write_csv(all_marks, paste0("outputs/", dataname,"/", dataname,"_",i,"_markers.csv"))
    #getting the top markers for each cluster, plotting
    for(x in unique(all_marks$cluster)){
      clust_sub <- filter(all_marks, cluster == x) %>% 
        filter(avg_log2FC > 0) %>% 
        arrange(desc(diff))
      pdf(paste0("outputs/", dataname, "/figures/", dataname,"_cluster_",i, "_", x, "_top_marks.pdf"))
      for(y in 1:10){
        plot(FeaturePlot(seu, clust_sub$gene[y]))
      }
      dev.off()
    }
  }
}
