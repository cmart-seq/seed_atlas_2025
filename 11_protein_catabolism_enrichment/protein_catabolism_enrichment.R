#performing an enrichment analysis for protein catabolism GO terms and gene sets, then performing the statistical analysis for all module scores in the timepoint and atlas objects

library(dplyr)
library(Seurat)
library(readr)
library(pheatmap)
library(RColorBrewer)


setwd("11_protein_catabolism_enrichment")

#reading in the atlas and individual timepoints
atlas <- readRDS("../10_peptide_enrichment/outputs/ATLAS/ATLAS_peps.rds")
DefaultAssay(atlas) <- "RNA"
atlas <- JoinLayers(atlas)

######## reading in data ######## 

#loading all of the timepoint annotated seurats
dap3 <- readRDS("../10_peptide_enrichment/outputs/3DAP/3DAP_peps.rds") 
DefaultAssay(dap3) <- "RNA"
dap3 <- JoinLayers(dap3)
dap5 <- readRDS("../10_peptide_enrichment/outputs/5DAP/5DAP_peps.rds")
DefaultAssay(dap5) <- "RNA"
dap5 <- JoinLayers(dap5)
dap7 <- readRDS("../10_peptide_enrichment/outputs/7DAP/7DAP_peps.rds")
DefaultAssay(dap7) <- "RNA"
dap7 <- JoinLayers(dap7)

#making the output directory if it doesnt exist
if (!dir.exists("outputs")) {
  dir.create("outputs", recursive = TRUE)
}

#Ub ligase keyword search from TAIR
ub_ligase <- read_delim("inputs/uniprot/ub_ligase_TAIRgene_results_2025-03-12.tsv", delim = "\t")

#GO term search results from uniprot. some genes are concatenated, need to split
read_uniprot <- function(path){
  genes_list <- read_delim(path, delim = "\t")%>% 
    mutate(TAIR = gsub(';', '', TAIR)) %>% 
    pull(TAIR) %>%
    na.omit()
  genes_list <- unname(sapply(genes_list, function(x){ifelse(nchar(x) == 18, c(substr(x,1,9), substr(x,10,18)), x)}))
}

#GO term search results from uniprot: 
#GO:0006473
acetylation <- read_uniprot("inputs/uniprot/acetylation_uniprotkb_GO_0006473_AND_taxonomy_id_37_2025_03_09.tsv")
#GO:0043543
acylation <- read_uniprot("inputs/uniprot/acylation_uniprotkb_GO_0043543_AND_taxonomy_id_37_2025_03_09.tsv")
#GO:0008213
alkylation <- read_uniprot("inputs/uniprot/alkylation_uniprotkb_GO_0008213_AND_taxonomy_id_37_2025_03_09.tsv")
#GO:0001666
hypox <- read_uniprot("inputs/uniprot/hypox_uniprotkb_go_0001666_AND_taxonomy_id_37_2025_03_09.tsv")
#GO:0000502
proteasome <- read_uniprot("inputs/uniprot/proteasome_uniprotkb_GO_0000502_AND_taxonomy_id_37_2025_03_11.tsv")
#GO:0030163
protein_catabol <- read_uniprot("inputs/uniprot/protein_catabol_uniprotkb_GO_0030163_AND_taxonomy_id_37_2025_03_09.tsv")
#GO:0031399
reg_of_prot_mod <- read_uniprot("inputs/uniprot/reg_of_prot_mod_uniprotkb_GO_0031399_AND_taxonomy_id_37_2025_03_09.tsv")
#GO:0016567
ubiquitin <- read_uniprot("inputs/uniprot/ubiquitin_uniprotkb_GO_0016567_AND_taxonomy_id_37_2025_03_11.tsv")

prot_mods <- list(ub_ligase = ub_ligase$Locus,
                  acetylation = acetylation,
                  acylation = acylation,
                  alkylation = alkylation,
                  hypox = hypox,
                  proteasome = proteasome,
                  protein_catabol = protein_catabol,
                  reg_of_prot_mod = reg_of_prot_mod,
                  ubiquitin = ubiquitin)

seu_list <- list(dap3, dap5, dap7, atlas)
seu_names <- list("DAP3", "DAP5", "DAP7", "ATLAS")

######## enrichment analysis ######## 

for(i in 1:length(seu_list)){
  #now adding mod scores for the complete gene sets
  seu <- AddModuleScore(
    object = seu_list[[i]],
    features = prot_mods,
    name = 'prot_mods'
  )
  
  #renaming the metadata
  new_cols <- c(grep("prot_mods",colnames(seu[[]]), invert = TRUE, value = TRUE), 
                "Ub_ligase_TAIR_modsc", 
                "Prot_acetylat_GO0006473_modsc", 
                "Prot_acylat_GO0043543_modsc",
                "Prot_alkylat_GO0008213_modsc", 
                "Resp_hypox_GO0001666_modsc",
                "Proteas_compl_GO0000502_modsc",
                "Protein_catabol_GO0030163_modsc", 
                "Reg_of_prot_mod_GO0031399_modsc", 
                "Prot_ubiq_GO0016567_modsc")
  
  length(new_cols) == length(colnames(seu[[]])) 
  
  colnames(seu[[]]) <- new_cols
  #this move just copies the modscores, getting rid of the old columns
  seu@meta.data <- seu@meta.data[, !grepl("prot_mods", colnames(seu@meta.data))]
  
  saveRDS(seu, paste0("outputs/", seu_names[[i]], "_merged_annotated_all_mods.rds"))
  
}

#since this is the last modification script, performing a correlation analysis between all mods
all_scores <- c(#signalling and transport modscores
  "Nucleoside_met_GO0009116_modsc", 
  "Cyto_resp_GO0009735_modsc", 
  "Cyto_bios_GO0009691_modsc", 
  "Isopren_GO0019840_modsc",
  "Lipid_drop_GO0005811_modsc",
  "Mono_storage_GO0012511_modsc",
  "Lipid_storage_GO0019915_modsc",
  "Resp_alc_GO0097306_modsc",
  "Alc_biosyn_GO0046165_modsc",
  "Cell_resp_ABA_GO0071215_modsc", 
  "Resp_ABA_GO0009737_modsc", 
  "ABA_bio_GO0009688_modsc", 
  "Nitrate_TAIR_transp_modsc", 
  "SWEET_TAIR_modsc", 
  "UMAMIT_TAIR_modsc", 
  "Autophagy_GO0006914_modsc", 
  "Brassino_homeo_GO0010268_modsc", 
  "Cadmium_resp_GO0046686_modsc", 
  "Cell_death_GO0008219_modsc", 
  "Metal_ion_GO0010038_modsc", 
  "PCD_GO0012501_modsc", 
  "Senescence_G00090693_modsc", 
  "Trans_metal_transp_GO0000041_modsc", 
  "Zinc_homeo_GO0006882_modsc_modsc", 
  "Brassino_biosyn_GO0016132_modsc", 
  "Photosynth_GO0015979_modsc", 
  "Hypoxia_GO0001666_modsc", 
  "Cell_wall_GO0042546_modsc", 
  "Mech_stim_GO0009612_modsc", 
  "Sec_vesc_GO0099503_modsc", 
  "Chitin_GO0010200_modsc", 
  "Hydroly_GO0004553_modsc", 
  "Cell_cell_sig_GO0007267_modsc", 
  "Polysac_binding_GO0030247_modsc", 
  "Xyloglucan_met_GO0010411_modsc", 
  "Xylogluc_Xyloglucos_GO0016762_modsc", 
  "Brassino_resp_GO0009741_modsc",
  "Callose_synth_0003843_modsc",
  "Auxin_biosynth_0009851_modsc",
  "Auxin_response_0009851_modsc",
  
  #peptide modscores
  "Asn_SSP_modsc1", "All_SSP_modsc1","IDA_modsc1", "GLV_modsc1", 
  "EPF_modsc1", "CEP_modsc1", "PSK_modsc1", "PSY_modsc1", 
  "RALFL_modsc1", "GASA_modsc1", "CLE_modsc1", "SCRL_modsc1", 
  "TPD_modsc1", "LTP_modsc1", "LCR_modsc1", "DEFL_modsc1", 
  "Thionin_modsc1", "PMEI_modsc1", "CX8_modsc1", "DIR_modsc1", "SI_modsc1",
  
  #protein catabolism scores
  "Ub_ligase_TAIR_modsc", 
  "Prot_acetylat_GO0006473_modsc", 
  "Prot_acylat_GO0043543_modsc",
  "Prot_alkylat_GO0008213_modsc", 
  "Resp_hypox_GO0001666_modsc",
  "Proteas_compl_GO0000502_modsc",
  "Protein_catabol_GO0030163_modsc", 
  "Reg_of_prot_mod_GO0031399_modsc", 
  "Prot_ubiq_GO0016567_modsc")


######## statistical analysis: cluster vs. all ######## 
#cluster vs. all other aggregated clusters with BH correction
#performing on the biological replicate means for each module

#performing statistical analysis
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
assign_scode <- function(pval) {
  case_when(
    pval == 0 ~ "*****",
    pval <= 0.0001 ~ "****",
    pval <= 0.001 ~ "***",
    pval <= 0.01 ~ "**",
    pval <= 0.05 ~ "*",
    pval <= 0.1 ~ "n.s. 0.1",
    pval <= 1 ~ "n.s.", 
    TRUE ~ "Invalid"
  )
}

#making a DAP3, DAP5 subset from the ATLAS for testing
dap3_5 <- subset(atlas, timepoint != "7DAP")

saveRDS(dap3_5, "outputs/DAP3_5_ATLAS.rds")

#making all L3 annotations distinct
dap3$level_3_annotation_full_timed <- paste(dap3$timepoint, dap3$level_3_annotation_full)
dap5$level_3_annotation_full_timed <- paste(dap5$timepoint, dap5$level_3_annotation_full)
dap7$level_3_annotation_full_timed <- paste(dap7$timepoint, dap7$level_3_annotation_full)
atlas$level_3_annotation_full_timed <- paste(atlas$timepoint, atlas$level_3_annotation_full)
dap3_5$level_3_annotation_full_timed <- paste(dap3_5$timepoint, dap3_5$level_3_annotation_full)

seu_list <- c(dap3, dap5, dap7, atlas, dap3_5)
names(seu_list) <- c("3DAP", "5DAP", "7DAP", "ATLAS", "3_5DAP")
level_list <- c("level_1_annotation", "level_2_annotation", "level_3_annotation_full_timed")

#function for module score statistics
mod_ttest <- function(seu, mod, cluster, level){
  mdata <-seu[[]] %>% dplyr::select("bio_rep", !!sym(level), !!sym(mod))
  mdata$group <- ifelse(mdata[[level]] == cluster, "group1", "group2")
  mdata_summary <- mdata %>% group_by(group, bio_rep) %>% summarise(avg = mean(!!sym(mod)))
  results<- t.test(avg ~ group, data = mdata_summary, alternative = "two.sided",var.equal = TRUE)  #student's t-test
  return(results)
}

t_test_table <- data.frame("cluster" = character(),
                          "comp" = character(),
                           "DAP" = double(),
                           "level" = character(),
                           "module" = character(),
                           "pval" = character(),
                           "mean_group1" = double(),
                           "mean_group2" = double())

for(i in 1:length(seu_list)){
  
  test_table <- seu_list[[i]][[]]
  
  timepoint <- names(seu_list)[i]
  
  for(l in level_list){
    
    features_to_test <- intersect(colnames(test_table), all_scores)
    for(f in features_to_test){
      for(c in unique(test_table[[l]])){
        print(paste0(c, ", ",f,", ",l))
        mod_results <- mod_ttest(seu_list[[i]], f, c, l)
        results_table_mods <- data.frame("cluster" = c, 
                                         "comp" = paste0(c, " vs. all"),
                                         "DAP" = timepoint,
                                         "level" = l,
                                         "module" = f,
                                         "pval" = mod_results$p.value, 
                                         "mean_group1" = mod_results$estimate[[1]],
                                         "mean_group2" = mod_results$estimate[[2]])
        
        # Append to main result table
        t_test_table <- rbind(t_test_table, results_table_mods)
      }
    }
  }
}


t_test_table <- arrange(t_test_table, pval)

write.csv(t_test_table, "t_test_table.csv")

#performing BH correction and assigning codes
t_test_table <- read_csv("t_test_table.csv")

t_test_table_corrected <- data.frame(cluster = character(),
                                     comp = character(),
                                     DAP = character(),
                                     level = character(),
                                     module = character(),
                                     pval = numeric(),
                                     mean_group1= numeric(),
                                     mean_group2= numeric(), 
                                     pval_BH_adj= numeric(), 
                                     pos_neg = character())

for(d in unique(t_test_table$DAP)){
  for(l in unique(t_test_table$level)){
    for(m in unique(t_test_table$module)){
      subgroup <- filter(t_test_table, DAP == d) %>% 
        filter(level == l) %>%
        filter(module == m)
      #correcting
      subgroup$pval_BH_adj <- p.adjust(subgroup$pval, method = 'BH')
      #noting whether the difference is positive or negative
      subgroup <- mutate(subgroup, pos_neg = ifelse(mean_group1 > mean_group2, "pos", "neg")) 
      subgroup <- dplyr::select(subgroup, cluster,comp, DAP, level, module, pval, mean_group1, mean_group2, pval_BH_adj, pos_neg)
      t_test_table_corrected <-rbind(t_test_table_corrected, subgroup)
    }
  }
}
  
  
  
t_test_table_corrected$code <- assign_scode(t_test_table_corrected$pval_BH_adj)

#saving
t_test_table_corrected <- t_test_table_corrected %>% dplyr::select(cluster, comp, DAP, level, module, pval,mean_group1,
                                               mean_group2, pval_BH_adj, code, pos_neg)

write_csv(t_test_table_corrected, "outputs/t_test_table_corrected.csv")
t_test_table_corrected <- read_csv("outputs/t_test_table_corrected.csv")

#performing another round to find timed L2 vs. atlas differences
seu_list <- c(atlas)
names(seu_list) <- c("ATLAS")
level_list <- c("level_2_annotation_timed")

t_test_table2 <- data.frame("cluster" = character(),
                           "comp" = character(),
                           "DAP" = double(),
                           "level" = character(),
                           "module" = character(),
                           "pval" = character(),
                           "mean_group1" = double(),
                           "mean_group2" = double())

for(i in 1:length(seu_list)){
  
  test_table <- seu_list[[i]][[]]
  
  timepoint <- names(seu_list)[i]
  
  for(l in level_list){
    
    features_to_test <- intersect(colnames(test_table), all_scores)
    for(f in features_to_test){
      for(c in unique(test_table[[l]])){
        mod_results <- mod_ttest(seu_list[[i]], f, c, l)
        results_table_mods <- data.frame("cluster" = c, 
                                         "comp" = paste0(c, " vs. all"),
                                         "DAP" = timepoint,
                                         "level" = l,
                                         "module" = f,
                                         "pval" = mod_results$p.value, 
                                         "mean_group1" = mod_results$estimate[[1]],
                                         "mean_group2" = mod_results$estimate[[2]])
        
        # Append to main result table
        t_test_table2 <- rbind(t_test_table2, results_table_mods)
      }
    }
  }
}



#performing BH correction and assigning codes
t_test_table_corrected2 <- data.frame(cluster = character(),
                                     comp = character(),
                                     DAP = character(),
                                     level = character(),
                                     module = character(),
                                     pval = numeric(),
                                     mean_group1= numeric(),
                                     mean_group2= numeric(), 
                                     pval_BH_adj= numeric(), 
                                     pos_neg = character())

for(d in unique(t_test_table2$DAP)){
  for(l in unique(t_test_table2$level)){
    for(m in unique(t_test_table2$module)){
      subgroup <- filter(t_test_table2, DAP == d) %>% 
        filter(level == l) %>%
        filter(module == m)
      
      subgroup$pval_BH_adj <- p.adjust(subgroup$pval, method = 'BH')
      #noting whether the difference is positive or negative
      subgroup <- mutate(subgroup, pos_neg = ifelse(mean_group1 > mean_group2, "pos", "neg")) 
      subgroup <- dplyr::select(subgroup, cluster,comp, DAP, level, module, pval, mean_group1, mean_group2, pval_BH_adj, pos_neg)
      t_test_table_corrected2 <-rbind(t_test_table_corrected2, subgroup)
    }
  }
}

t_test_table_corrected2$code <- assign_scode(t_test_table_corrected2$pval_BH_adj)

#saving
t_test_table_corrected2 <- t_test_table_corrected2 %>% dplyr::select(cluster, comp, DAP, level, module, pval,mean_group1,
                                                                   mean_group2, pval_BH_adj, code, pos_neg)

t_test_table_corrected <- read_csv("outputs/t_test_table_corrected.csv")
t_test_table_corrected <- rbind(t_test_table_corrected, t_test_table_corrected2)
write_csv(t_test_table_corrected, "outputs/t_test_table_corrected.csv")
