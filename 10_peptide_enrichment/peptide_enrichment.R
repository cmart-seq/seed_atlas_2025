#performing GSEA for peptide families in all timepoint objects and the atlas dataset

library(Seurat)
library(dplyr)
library(readr)
library(clusterProfiler)
library(org.At.tair.db)
library(stringr)
library(tidyr)
library(ggplot2)
library(seqinr)
library(cowplot)

set.seed(123)

setwd("/lab/solexa_gehring/carly/snRNA_seed_dev/clean_analysis/10_peptide_enrichment")
#reading in data
dap3 <- readRDS("../09_signalling_transport_gene_enrichment/outputs/DAP3_merged_annotated_sigmods.rds") 
DefaultAssay(dap3) <- "RNA"
dap3 <- JoinLayers(dap3)

dap5 <- readRDS("../09_signalling_transport_gene_enrichment/outputs/DAP5_merged_annotated_sigmods.rds")
DefaultAssay(dap5) <- "RNA"
dap5 <- JoinLayers(dap5)

dap7 <- readRDS("../09_signalling_transport_gene_enrichment/outputs/DAP7_merged_annotated_sigmods.rds")
DefaultAssay(dap7) <- "RNA"
dap7 <- JoinLayers(dap7)

atlas <- readRDS(paste0("../09_signalling_transport_gene_enrichment/outputs/ATLAS_merged_annotated_sigmods.rds"))  
DefaultAssay(atlas) <- "RNA"
atlas <- JoinLayers(atlas)

#color schemes
color_key1 = c("Embryo" = "#57A15D", 
               "Endosperm" = "#D44A90", 
               "Seed coat" = "#96ADD0", 
               "Funiculus" = "#f1c40f", 
               "Ovule" = "#a6acaf")

levels1 <- c("Embryo",
             "Endosperm",
             "Seed coat", 
             "Funiculus", 
             "Ovule")

color_key2 = c("EMB" = "#57A15D", 
               "PEN" = "#ceaab5", 
               "MCE" = "#97576b", 
               "CZE" = "#ce0d48",
               "FUN" = "#f1c40f",
               "OVL" = "#a6acaf", 
               "ii1" = "#aed6f1", 
               "ii1'" = "#5dade2", 
               "ii2" = "#21618c",
               "oi1" = "#bb8fce", 
               "oi2" = "#6c3483", 
               "CZSC" = "#088F8F", 
               "CPT" = "#044747")

levels2 = c("EMB",
            "PEN",
            "MCE",
            "CZE",
            "FUN",
            "OVL",
            "ii1",
            "ii1'",
            "ii2",
            "oi1",
            "oi2",
            "CZSC", 
            "CPT")

color_key2_d7 = c("EMB" = "#57A15D", 
                  "PEN" = "#ceaab5", 
                  "MCE" = "#97576b", 
                  "CZE" = "#ce0d48",
                  "FUN" = "#f1c40f",
                  "OVL" = "#a6acaf", 
                  "ii1" = "#aed6f1", 
                  "ii1'/ii2" = "#5dade2", 
                  "oi1" = "#bb8fce", 
                  "oi2" = "#6c3483", 
                  "CZSC" = "#088F8F", 
                  "CPT" = "#044747")

levels2_d7 = c("EMB",
               "PEN",
               "MCE",
               "CZE",
               "FUN",
               "OVL",
               "ii1",
               "ii1'/ii2",
               "ii2",
               "oi1",
               "oi2",
               "CZSC", 
               "CPT")

#making a level_3 color scheme for all three levels
#level 3, DAP3
obj_list <- c(dap3, dap5, dap7)
names(obj_list) <-c("DAP3", "DAP5", "DAP7")

for(i in 1:3){
  lvl3_labels <- sort(unique(obj_list[[i]][[]]$level_3_annotation_full))
  
  if(names(obj_list)[i] == "DAP7"){
    color_key2_tbl <- tibble(level_2_annotation = names(color_key2_d7), color = color_key2_d7)
  }else{
    color_key2_tbl <- tibble(level_2_annotation = names(color_key2), color = color_key2)
  }
  
  
  color_key_l3 <- tibble(level_2_annotation = obj_list[[i]][[]]$level_2_annotation, 
                         level_3_annotation_full = obj_list[[i]][[]]$level_3_annotation_full,
                         level_3_plotting = obj_list[[i]][[]]$level_3_plotting) %>% unique() %>%
    left_join(color_key2_tbl, by = "level_2_annotation")
  
  
  color_key_l3_anno <- color_key_l3$color
  names(color_key_l3_anno) <- color_key_l3$level_3_annotation_full
  
  write.csv(data.frame(tissue = names(color_key_l3_anno),
                       color = color_key_l3_anno), paste0("inputs/",names(obj_list[i]), "_level_3_anno_color_scheme.csv"))
  
  color_key_l3_plotting <- color_key_l3$color
  names(color_key_l3_plotting) <- color_key_l3$level_3_plotting
  
  write.csv(data.frame(tissue = names(color_key_l3_plotting),
                       color = color_key_l3_plotting), paste0("inputs/", names(obj_list[i]), "_level_3_plotting_color_scheme.csv"))
}

#peptide data
#getting the peptides from the signalP analysis
pep_fasta <- read.fasta("inputs/signalP/signalP_results/processed_entries.fasta")
head(names(pep_fasta))
#making them locus IDs
pep_genes <- sapply(names(pep_fasta), function(x){str_split(x, "[.]")[[1]][1]}) %>% unname() %>% unique()

#getting the ghorbani peptides
peps_erv346 <- read_csv("inputs/erv346_Supplementary_Data/all_peptide_annotations.csv")
#filtering to arabidopsis
peps_erv346 <- mutate(peps_erv346, species = sapply(peps_erv346$`Gene ID`, function(x){species = strsplit(x, "_")[[1]][1]})) %>%
  filter(species == "tair10")

#filtering, only keeping the primary transcript
peps_erv346 <- mutate(peps_erv346, iso = sapply(peps_erv346$`Gene ID`, function(x){species = strsplit(x, "[_.]")[[1]][3]}))%>%
  filter(iso == 1)

peps_erv346 <- mutate(peps_erv346, tair = sapply(peps_erv346$`Gene ID`, function(x){species = strsplit(x, "[_.]")[[1]][2]}))

#saving
write_csv(peps_erv346, "inputs/erv346_Supplementary_Data/all_peptide_annotations_simplified.csv")

#filtering ghorbani based on signalp
peps_erv346_sp6 <- filter(peps_erv346, tair %in% pep_genes)
#filtering ghorbani/signalP based on atlas expression
peps_erv346_sp6 <- filter(peps_erv346_sp6, tair %in% rownames(atlas))
#updating the SSP annotation
peps_erv346_sp6$SSP <- replace_na(peps_erv346_sp6$SSP, "no_family")
dim(peps_erv346_sp6)

#annotating ssp families that were not included in ghorbani
def_fasta <- read.fasta("inputs/DEFL_TAIR_proteins.fa")
defs <- unname(sapply(names(def_fasta), function(x){str_split(x, "[.]")[[1]][1]}))

ltp_fasta <- read.fasta("inputs/LTP_TAIR_proteins.fa")
ltp <- union(unname(sapply(names(ltp_fasta), function(x){str_split(x, "[.]")[[1]][1]})), peps_erv346_sp6[peps_erv346_sp6$SSP == "LTP",]$tair)

lcr_fasta <- read.fasta("inputs/LCR_TAIR_proteins.fa")
lcr <- unname(sapply(names(lcr_fasta), function(x){str_split(x, "[.]")[[1]][1]}))

pmei_fasta <- read.fasta("inputs/PMEI_TAIR_proteins.fa")
pmei <- unname(sapply(names(pmei_fasta), function(x){str_split(x, "[.]")[[1]][1]}))

#LTPs, LCRs, DEFLs
peps_erv346_sp6$SSP <- sapply(1:length(peps_erv346_sp6$SSP), 
                                     function(x){
                                       if(peps_erv346_sp6$tair[x] %in% defs) {
                                         return("DEFL")
                                       } else if(peps_erv346_sp6$tair[x] %in% ltp) {
                                         return("LTP")
                                       } else if(peps_erv346_sp6$tair[x] %in% lcr) {
                                         return("LCR")
                                       } else if(peps_erv346_sp6$tair[x] %in% pmei) {
                                         return("PMEI")
                                       } else {
                                         return(peps_erv346_sp6$SSP[x])
                                       }
                                     })

#manually curated SSP
manual_fams <- read_csv("inputs/manual_curation.csv") %>% dplyr::select(Locus, manual_classification)
colnames(manual_fams) <- c("tair", "SSP")

#adding the updated SSP annotations 
peps_erv346_sp6_man_fams <- filter(peps_erv346_sp6, tair %in% manual_fams$tair) %>% 
  dplyr::select(!SSP) %>%
  left_join(manual_fams, by = "tair")

#updating the no_fam_genes in the ghorbani table
peps_erv346_sp6 <- filter(peps_erv346_sp6, !tair %in% manual_fams$tair) %>%
  rbind(peps_erv346_sp6_man_fams)

#generating gene lists
thionin <- peps_erv346_sp6[peps_erv346_sp6$SSP == "THIONIN",]$tair

pmei <- peps_erv346_sp6[peps_erv346_sp6$SSP == "PMEI",]$tair

cx8 <- peps_erv346_sp6[peps_erv346_sp6$SSP == "CX8",]$tair

diri <- peps_erv346_sp6[peps_erv346_sp6$SSP == "DIR",]$tair

si <- peps_erv346_sp6[peps_erv346_sp6$SSP == "SPH",]$tair

#saving the peptide annotation table
write.csv(peps_erv346_sp6, "outputs/annotated_peptide_table_ghorbani_SP6_curation.csv")

#ATH1 genes, from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL198#:~:text=The%20GeneChip%C2%AE%20Arabidopsis%20ATH1,represented%20on%20the%20ATH1%20array.
ath1_raw <- read_delim("inputs/GPL198-17390.txt", delim = "\t")
cleaned_genes <- unique(na.omit(unname(unlist(sapply(ath1_raw$AGI, function(x){str_trim(str_split(x, "\\///")[[1]])})))))


##### performing a module score analysis of SSPs ####

peptide_enrichment_analysis <- function(seurat_object, timepoint){
  
  #making the output directory if it doesnt exist
  if (!dir.exists("outputs")) {
    dir.create("outputs", recursive = TRUE)
  }
  
  if (!dir.exists(file.path("outputs", timepoint))) {
    dir.create(file.path("outputs", timepoint), recursive = TRUE)
  }
  
  if (!dir.exists(file.path("outputs", timepoint, "figures"))) {
    dir.create(file.path("outputs", timepoint, "figures"), recursive = TRUE)
  }
  
  #all peps
  peps_erv346_sp6 <- filter(peps_erv346_sp6, tair %in% pep_genes) %>% filter(tair %in% rownames(seurat_object))
  num_peps <- length(unique(pep_genes))
  
  #assigned peps
  asn_SSPs <- filter(peps_erv346_sp6, SSP != "no_family")%>% 
    pull(tair) %>% unique()
  
  #RALFLs
  ralfs <- filter(peps_erv346_sp6, SSP == "RALFL") %>% pull(tair)
  num_ralfls <- length(ralfs)
  
  #GASA
  gasas <- filter(peps_erv346_sp6, SSP == "GASA") %>% pull(tair)
  num_gasas <- length(gasas)
  
  #SCRL
  scrls <- filter(peps_erv346_sp6, SSP == "SCRL") %>% pull(tair)
  num_scrls <- length(scrls)
  
  #CLE
  cles <- filter(peps_erv346_sp6, SSP == "CLE") %>% pull(tair)
  num_cles <- length(cles)
  
  #TPD
  tpds <- filter(peps_erv346_sp6, SSP == "TPD") %>% pull(tair)
  num_tpds <- length(tpds)
  
  #EPF
  epfs <- filter(peps_erv346_sp6, SSP == "EPF") %>% pull(tair)
  num_tpds <- length(epfs)
  
  #PSY
  psys <- filter(peps_erv346_sp6, SSP == "PSY") %>% pull(tair)
  num_psys <- length(psys)
  
  #PSK
  psks <- filter(peps_erv346_sp6, SSP == "PSK") %>% pull(tair)
  num_psks <- length(psks)
  
  #CEP
  ceps <- filter(peps_erv346_sp6, SSP == "CEP") %>% pull(tair)
  num_ceps <- length(ceps)
  
  #EPF
  epfs <- filter(peps_erv346_sp6, SSP == "EPF") %>% pull(tair)
  num_epfs <- length(epfs)
  
  #GLV
  glvs <- filter(peps_erv346_sp6, SSP == "GLV") %>% pull(tair)
  num_glvs <- length(glvs)
  
  #IDA/IDL
  idas <- filter(peps_erv346_sp6, SSP == "IDA/IDL") %>% pull(tair)
  num_idas <- length(idas)
  
  #PIP/PIPL
  pippipls <- filter(peps_erv346_sp6, SSP == "PIP/PIPL") %>% pull(tair)
  num_pippipls <- length(pippipls)
  
  #LTP
  ltp <- intersect(ltp, rownames(seurat_object))
  num_ltp <- length(ltp)
  
  #LCR
  lcr <- intersect(lcr, rownames(seurat_object))
  num_lcr <- length(lcr)
  
  #DEFL
  defs <- intersect(defs, rownames(seurat_object))
  num_defs <- length(defs)
  
  #THIONIN
  thionins <- intersect(thionin, rownames(seurat_object))
  num_thionins <- length(thionins)
  
  #PMEI
  pmeis <- intersect(pmei, rownames(seurat_object))
  num_pmeis <- length(pmeis)
  
  #CX8
  cx8s <- intersect(cx8, rownames(seurat_object))
  num_cx8s <- length(cx8s)
  
  #DIR
  diris <- intersect(diri, rownames(seurat_object))
  num_diris <- length(diris)
  
  #SI
  sis <- intersect(si, rownames(seurat_object))
  num_sis <- length(sis)

  #listing for looping
  pep_sets <- list(pippipls, idas, glvs, epfs, ceps, psks, psys, ralfs, gasas, cles, scrls, tpds, 
                   ltp, lcr, defs, thionins, pmeis, cx8s, diris, sis)
  names(pep_sets) <- c("PIPPIPL_modsc", "IDA_modsc", "GLV_modsc", "EPF_modsc", "CEP_modsc", "PSK_modsc", "PSY_modsc", 
                       "RALFL_modsc", "GASA_modsc", "CLE_modsc", "SCRL_modsc", "TPD_modsc", 
                       "LTP_modsc", "LCR_modsc", "DEFL_modsc", 
                       "Thionin_modsc", "PMEI_modsc", "CX8_modsc", "DIR_modsc", "SI_modsc")
  
  #reporting the number of peps in each category
  write.csv(data.frame(sapply(pep_sets, length)), paste0("outputs/", timepoint,"/",timepoint,"_pepfam_mod_numbers.csv"))
  
  for(i in 1:length(pep_sets)){
    #performing the module scoring SCRL
    if(length(pep_sets[[i]])>1){
      features <- list(c(pep_sets[[i]]))
      seurat_object <- AddModuleScore(
        object = seurat_object,
        features = features,
        name = names(pep_sets[i])
      )
    }
  }
  
  #performing the module scoring, all peptides
  pep_features <- list(c(unique(pep_genes)))
  seurat_object <- AddModuleScore(
    object = seurat_object,
    features = pep_features,
    name = 'All_SSP_modsc'
  )
  
  #performing the module scoring, all assigned peptides
  asn_pep_features <- list(c(unique(asn_SSPs)))
  seurat_object <- AddModuleScore(
    object = seurat_object,
    features = asn_pep_features,
    name = 'Asn_SSP_modsc'
  )
  
  
  #saving the new object
  saveRDS(seurat_object, paste0("outputs/", timepoint, "/",timepoint, "_peps.rds"))
  
  return(seurat_object)
  
  
}

dap3 <- peptide_enrichment_analysis(dap3, "3DAP")
dap5 <- peptide_enrichment_analysis(dap5, "5DAP")
dap7 <- peptide_enrichment_analysis(dap7, "7DAP")
atlas <- peptide_enrichment_analysis(atlas, "ATLAS")


#and a plotting function for the new mod scores
plotting_pep_mods <- function(seurat_object, timepoint, level, order, color_key){
  seurat_object <- SetIdent(seurat_object, value = level) 
  seurat_object@active.ident <- factor(x = seurat_object@active.ident, levels = order)
  
  if (!dir.exists(file.path("outputs", timepoint, "figures/pngs"))) {
    dir.create(file.path("outputs", timepoint, "figures/pngs"), recursive = TRUE)
  }

  features_list <- intersect(colnames(seurat_object[[]]), c("All_SSP_modsc1", "Asn_SSP_modsc1","IDA_modsc1", "GLV_modsc1", 
                                                            "EPF_modsc1", "CEP_modsc1", "PSK_modsc1", "PSY_modsc1", 
                                                            "RALFL_modsc1", "GASA_modsc1", "CLE_modsc1", "SCRL_modsc1", 
                                                            "TPD_modsc1", "LTP_modsc1", "LCR_modsc1", "DEFL_modsc1", 
                                                            "Thionin_modsc1", "PMEI_modsc1", "CX8_modsc1", "DIR_modsc1", "SI_modsc1"))
  
  for(y in features_list){
    png(paste0("outputs/", timepoint, "/figures/pngs/",y, "_", timepoint, "_", level ,"_pep_mod_scores_vlnplots.png"),
        width     = 3.25,
        height    = 3.25,
        units     = "in",
        res       = 1200,
        pointsize = 4)
    plot(VlnPlot(object = seurat_object, features = y, alpha = 0.02) +
           geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
           scale_fill_manual(values = color_key)+
           ggtitle(paste0(y))+
           theme(legend.position = "none")+
           theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                 axis.text.y = element_text(size = 10), 
                 legend.position = "none"))
    dev.off()
    
    png(paste0("outputs/", timepoint, "/figures/pngs/",y, "_", timepoint, "_", level ,"_pep_mod_scores_boxplots.png"),
        width     = 3.25,
        height    = 3.25,
        units     = "in",
        res       = 1200,
        pointsize = 4)
    plot(ggplot(seurat_object[[]], aes(x = eval(parse(text=level)), y = eval(parse(text=y)), fill = eval(parse(text=level)))) +
           geom_boxplot(outlier.shape = NA) +
           scale_fill_manual(values = color_key) +  # Use the color key for the fill colors
           #coord_cartesian(ylim = c(0, 2.5)) +        # Set limits without truncating data
           theme_classic()+
           geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
           theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                 axis.text.y = element_text(size = 10), 
                 legend.position = "none"))
    dev.off()
  }
  
  pdf(paste0("outputs/", timepoint, "/figures/",timepoint, "_", level ,"_pep_mod_scores_vlnplots.pdf"), width = 4, height = 4)
  features_list <- intersect(colnames(seurat_object[[]]), c("All_SSP_modsc1", "Asn_SSP_modsc1","IDA_modsc1", "GLV_modsc1", 
                                                            "EPF_modsc1", "CEP_modsc1", "PSK_modsc1", "PSY_modsc1", 
                                                            "RALFL_modsc1", "GASA_modsc1", "CLE_modsc1", "SCRL_modsc1", 
                                                            "TPD_modsc1", "LTP_modsc1", "LCR_modsc1", "DEFL_modsc1", 
                                                            "Thionin_modsc1", "PMEI_modsc1", "CX8_modsc1", "DIR_modsc1", "SI_modsc1"))
  
  for(y in features_list){
    plot(VlnPlot(object = seurat_object, features = y, alpha = 0.02) +
           geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
           scale_fill_manual(values = color_key)+
           ggtitle(paste0(y))+
           theme(legend.position = "none"))
    
    plot(ggplot(seurat_object[[]], aes(x = eval(parse(text=level)), y = eval(parse(text=y)), fill = eval(parse(text=level)))) +
           geom_boxplot(outlier.shape = NA) +
           scale_fill_manual(values = color_key) +  # Use the color key for the fill colors
           #coord_cartesian(ylim = c(0, 2.5)) +        # Set limits without truncating data
           theme_classic()+
           geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
           theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
                 axis.text.y = element_text(size = 18),
                 legend.position = "none"))
  }
  dev.off()
}

#3DAP
plotting_pep_mods(dap3, "3DAP", "level_2_annotation", levels2, color_key2)
plotting_pep_mods(dap3, "3DAP", "level_1_annotation", levels1, color_key1)

#5DAP
plotting_pep_mods(dap5, "5DAP", "level_2_annotation", levels2, color_key2)
plotting_pep_mods(dap5, "5DAP", "level_1_annotation", levels1, color_key1)

#7DAP
plotting_pep_mods(dap7, "7DAP", "level_2_annotation", levels2_d7, color_key2_d7)
plotting_pep_mods(dap7, "7DAP", "level_1_annotation", levels1, color_key1)

#ATLAS
plotting_pep_mods(atlas, "ATLAS", "level_2_annotation", levels2, color_key2)
plotting_pep_mods(atlas, "ATLAS", "level_1_annotation", levels1, color_key1)

#plotting SSP enrichment on level 3 annotations
# getting the color scheme
DAP3_color_key_l3_table <- read.csv("inputs/DAP3_level_3_plotting_color_scheme.csv")
DAP3_color_key_l3_vector <- DAP3_color_key_l3_table$color
names(DAP3_color_key_l3_vector) <- DAP3_color_key_l3_table$tissue

DAP5_color_key_l3_table <- read.csv("inputs/DAP5_level_3_plotting_color_scheme.csv")
DAP5_color_key_l3_vector <- DAP5_color_key_l3_table$color
names(DAP5_color_key_l3_vector) <- DAP5_color_key_l3_table$tissue

DAP7_color_key_l3_table <- read.csv("inputs/DAP7_level_3_plotting_color_scheme.csv")
DAP7_color_key_l3_vector <- DAP7_color_key_l3_table$color
names(DAP7_color_key_l3_vector) <- DAP7_color_key_l3_table$tissue

#need to do level 3 annotation separately
plotting_pep_mods_L3 <- function(seurat_object, timepoint, color_key, wi, he){
  
  color_order <- color_key[order(match(names(color_key), unique(Idents(seurat_object))))]
  
  # Extract numeric part of the cluster labels and reorder factors numerically
  numeric_order <- sapply(unique(seurat_object$level_3_annotation_full), function(x) {
    as.numeric(str_split(x, "\\.")[[1]][1])  # Extract and convert to numeric
  })
  
  # Reorder `cluster` factor levels based on numeric sorting, for both the full and plotting annotation
  ordered_full <- unique(seurat_object$level_3_annotation_full)[order(numeric_order)]
  ordered_plotting <- unique(seurat_object$level_3_plotting)[order(numeric_order)]
  
  seurat_object$level_3_annotation_full <- factor(seurat_object$level_3_annotation_full, levels = ordered_full)
  seurat_object$level_3_plotting <- factor(seurat_object$level_3_plotting, levels = ordered_plotting)
  
  seurat_object <- SetIdent(seurat_object, value = "level_3_plotting") 
  
  features_to_plot <- intersect(colnames(seurat_object[[]]), c("All_SSP_modsc1","Asn_SSP_modsc1", "IDA_modsc1", "GLV_modsc1", 
                                                              "EPF_modsc1", "CEP_modsc1", "PSK_modsc1", "PSY_modsc1", 
                                                              "RALFL_modsc1", "GASA_modsc1", "CLE_modsc1", "SCRL_modsc1", 
                                                              "TPD_modsc1", "LTP_modsc1", "LCR_modsc1", "DEFL_modsc1", 
                                                              "Thionin_modsc1", "PMEI_modsc1", "CX8_modsc1", "DIR_modsc1", "SI_modsc1"))
  
  for(y in features_to_plot){
    png(paste0("outputs/", timepoint, "/figures/pngs/",y, "_", timepoint, "_level_3_annotation_pep_mod_scores_vlnplots.png"),
        width     = wi,
        height    = he,
        units     = "in",
        res       = 480)
    p1 <- VlnPlot(object = seurat_object, features = y, alpha = 0.02, cols = color_order) +
      geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
      ggtitle(paste0(y))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 10), 
            legend.position = "none")
    plot(p1 + theme(legend.position = "none"))
    dev.off()
    
  }

  pdf(paste0("outputs/",timepoint, "/figures/",timepoint, "_level_3_annotation_mod_scores_vlnplots.pdf"), width = wi, height = he)
 
   for(y in features_to_plot){
    p1 <- VlnPlot(object = seurat_object, features = y, alpha = 0.02, cols = color_order) +
      geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
      ggtitle(paste0("Cell module score for genes with ",y))+
      theme_classic() +
      theme(axis.text.y = element_text(size = 18),
            axis.text.x = element_text(size = 14),  # Adjust size for visibility
            axis.ticks.x = element_blank())
    plot(p1 + theme(legend.position = "none"))
  }
  dev.off()
}

#3DAP
plotting_pep_mods_L3(dap3, "3DAP", DAP3_color_key_l3_vector, 13, 3)
#5DAP
plotting_pep_mods_L3(dap5, "5DAP", DAP5_color_key_l3_vector, 13, 3)
#7DAP
plotting_pep_mods_L3(dap7, "7DAP", DAP7_color_key_l3_vector, 13, 3)

#printing the level 3 key, 3DAP
dap3 <- SetIdent(dap3, value = "level_3_annotation_full")
color_order <- DAP3_color_key_l3_vector[order(match(names(DAP3_color_key_l3_vector), unique(Idents(dap3))))]

names(color_order) <- unique(Idents(dap3))
  
png(paste0("outputs/3DAP/figures/DAP3_level_3_annotation_pep_mod_scores_key.png"),
    width     = 12,
    height    = 6,
    units     = "in",
    res       = 1200,
    pointsize = 4)
p1 <- VlnPlot(object = dap3, features = "Asn_SSP_modsc1", alpha = 0.02, cols = color_order) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10))

key <- get_legend(p1)

plot(plot_grid(key, ncol = 2))

dev.off()

#5DAP
dap5 <- SetIdent(dap5, value = "level_3_annotation_full")
color_order <- DAP5_color_key_l3_vector[order(match(names(DAP5_color_key_l3_vector), unique(Idents(dap5))))]

names(color_order) <- unique(Idents(dap5))

png(paste0("outputs/5DAP/figures/5DAP_level_3_annotation_pep_mod_scores_key.png"),
    width     = 12,
    height    = 6,
    units     = "in",
    res       = 1200,
    pointsize = 4)
p1 <- VlnPlot(object = dap5, features = "Asn_SSP_modsc1", alpha = 0.02, cols = color_order) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10))

key <- get_legend(p1)

plot(plot_grid(key, ncol = 2))

dev.off()

#7DAP
dap7 <- SetIdent(dap7, value = "level_3_annotation_full")
color_order <- DAP7_color_key_l3_vector[order(match(names(DAP7_color_key_l3_vector), unique(Idents(dap7))))]

names(color_order) <- unique(Idents(dap7))

png(paste0("outputs/7DAP/figures/7DAP_level_3_annotation_pep_mod_scores_key.png"),
    width     = 15,
    height    = 6,
    units     = "in",
    res       = 1200,
    pointsize = 4)
p1 <- VlnPlot(object = dap7, features = "Asn_SSP_modsc1", alpha = 0.02, cols = color_order) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10))

key <- get_legend(p1)

plot(plot_grid(key, ncol = 2))

dev.off()




