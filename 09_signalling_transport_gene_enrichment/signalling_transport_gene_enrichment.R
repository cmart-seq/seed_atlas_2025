#performing an enrichment analysis for signalling and transport gene sets in all timepoint objects and the atlas

library(dplyr)
library(Seurat)
library(readr)

setwd("09_signalling_transport_gene_enrichment")

#reading in the atlas and individual timepoints
atlas <- readRDS("../06_atlas_merging/outputs/ATLAS_integrated_annotated_rpca.rds")
DefaultAssay(atlas) <- "RNA"
atlas <- JoinLayers(atlas)

#loading all of the timepoint annotated seurats
dap3 <- readRDS("../04_manual_annotation/outputs/DAP3_wcze_subs/DAP3_wcze_subs_annotated.rds") 
DefaultAssay(dap3) <- "RNA"
dap3 <- JoinLayers(dap3)
dap5 <- readRDS("../04_manual_annotation/outputs/DAP5_wcze_subs/DAP5_wcze_subs_annotated.rds")
DefaultAssay(dap5) <- "RNA"
dap5 <- JoinLayers(dap5)
dap7 <- readRDS("../04_manual_annotation/outputs/DAP7/DAP7_annotated.rds")
DefaultAssay(dap7) <- "RNA"
dap7 <- JoinLayers(dap7)


#making the output directory if it doesnt exist
if (!dir.exists("outputs")) {
  dir.create("outputs", recursive = TRUE)
}

#reading in the input gene lists
#keyword search results from TAIR
nitrate_trans <- read_delim("inputs/tair/nitrate_transporters_gene_results_2025-03-1.tsv", delim = "\t")
sweet_trans <- read_delim("inputs/tair/SWEET_gene_results_2025-03-14.tsv", delim = "\t")
umami_trans <- read_delim("inputs/tair/UMAMIT_gene_results_2025-03-14.tsv", delim = "\t")

#GO term search results from uniprot. some genes are concatenated, need to split
read_uniprot <- function(path){
  genes_list <- read_delim(path, delim = "\t")%>% 
    mutate(TAIR = gsub(';', '', TAIR)) %>% 
    pull(TAIR) %>%
    na.omit()
  genes_list <- unname(sapply(genes_list, function(x){ifelse(nchar(x) == 18, c(substr(x,1,9), substr(x,10,18)), x)}))
}

# nucleoside metabolism, GO_0009116
nuc_met <- read_uniprot("inputs/uniprot/nucleoside_met_uniprotkb_GO_0009116_AND_taxonomy_id_37_2025_04_02.tsv")
# cytokinin response
cyto_resp <- read_uniprot("inputs/uniprot/cytokinin_resp_uniprotkb_go_0009735_AND_taxonomy_id_37_2025_04_02.tsv")
# cytokinin biosynthesis
cyto_bios <- read_uniprot("inputs/uniprot/cytokinin_biosynthetic_uniprotkb_go_0009691_AND_taxonomy_id_37_2025_04_02.tsv")
#GO:0019840, isoprenoid
isopren <- read_uniprot("inputs/uniprot/isopren_uniprotkb_go_0019840_AND_taxonomy_id_37_2025_04_02.tsv")
#GO:0005811, "lipid droplet" 
lipid_drop <- read_uniprot("inputs/uniprot/lipid_drop_uniprotkb_go_0005811_AND_taxonomy_id_37_2025_04_02.tsv")
#GO:0012511, "monolayer-surrounded lipid storage body" 
mono_storage <- read_uniprot("inputs/uniprot/mono_storage_uniprotkb_go_0012511_AND_taxonomy_id_37_2025_04_02.tsv")
#GO:0019915, "lipid storage"
lipid_storage <- read_uniprot("inputs/uniprot/lipid_storage_uniprotkb_go_0019915_AND_taxonomy_id_37_2025_04_02.tsv")
#GO:0097306, "cellular response to alcohol"
resp_alc <- read_uniprot("inputs/uniprot/resp_alc_uniprotkb_go_0097306_AND_taxonomy_id_37_2025_04_02.tsv")
#GO:0046165, "alcohol biosynthetic process"  
alc_bios <- read_uniprot("inputs/uniprot/alc_biosuniprotkb_go_0046165_AND_taxonomy_id_37_2025_04_02.tsv")
#GO:0071215 "cellular response to abscisic acid stimulus"
cell_resp_aba <- read_uniprot("inputs/uniprot/cell_resp_aba_uniprotkb_go_0071215_AND_taxonomy_id_37_2025_04_02.tsv")
#GO:0009737    response to abscisic acid
resp_aba <- read_uniprot("inputs/uniprot/resp_aba_uniprotkb_go_0009737_AND_taxonomy_id_37_2025_04_02.tsv")
#GO:0009688    abscisic acid biosynthetic process
aba_bio <- read_uniprot("inputs/uniprot/aba_bio_uniprotkb_go_0009688_AND_taxonomy_id_37_2025_04_02.tsv")
#GO:0006914
autophagy <- read_uniprot("inputs/uniprot/autophagy_uniprotkb_go_0006914_AND_taxonomy_id_37_2025_03_08.tsv")
#GO:0010268
br_homeo <- read_uniprot("inputs/uniprot/BR_homeo_uniprotkb_go_0010268_AND_taxonomy_id_37_2025_03_08.tsv")
#GO:0046686
cadmium <- read_uniprot("inputs/uniprot/cad_uniprotkb_go_0046686_AND_taxonomy_id_37_2025_03_08.tsv")
#GO:0008219
cell_death <- read_uniprot("inputs/uniprot/cell_death_uniprotkb_go_0008219_AND_taxonomy_id_37_2025_03_08.tsv")
#GO:0010038
metal_ion_resp <- read_uniprot("inputs/uniprot/metal_ion_response_uniprotkb_go_0010038_AND_taxonomy_id_37_2025_03_08.tsv")
#GO:0012501
pcd <- read_uniprot("inputs/uniprot/PCD_uniprotkb_go_0012501_AND_taxonomy_id_37_2025_03_08.tsv")
#GO:0090693
senes <- read_uniprot("inputs/uniprot/senescence_uniprotkb_go_0090693_AND_taxonomy_id_37_2025_03_08.tsv")
#GO:0000041
trans_metal <- read_uniprot("inputs/uniprot/trans_metal_uniprotkb_go_0000041_AND_taxonomy_id_37_2025_03_08.tsv")
#GO:0006882
zinc <- read_uniprot("inputs/uniprot/zinc_ion_uniprotkb_go_0006882_AND_taxonomy_id_37_2025_03_08.tsv")
#GO:0016132
brassino_bios <- read_uniprot("inputs/uniprot/brassino_biosynth_uniprotkb_GO_0016132_AND_taxonomy_id_37_2025_03_10.tsv")
#GO:0015979
photosynth <- read_uniprot("inputs/uniprot/photosynth_uniprotkb_GO_0015979_AND_taxonomy_id_37_2025_03_10.tsv")
#GO:0001666
hypoxia <- read_uniprot("inputs/uniprot/hypoxia_uniprotkb_GO_0001666_AND_taxonomy_id_37_2025_03_10.tsv")
#GO:0042546
cell_wall <- read_uniprot("inputs/uniprot/cell_wall_uniprotkb_GO_0042546_AND_taxonomy_id_37_2025_03_10.tsv")
#GO:0009612
mech_stim <- read_uniprot("inputs/uniprot/mech_stim_uniprotkb_GO_0009612_AND_taxonomy_id_37_2025_03_10.tsv")
#GO:0099503
sec_vesc <- read_uniprot("inputs/uniprot/sec_vesc_uniprotkb_GO_0099503_AND_taxonomy_id_37_2025_03_10.tsv")
#GO:0010200
chitin <- read_uniprot("inputs/uniprot/chitin_uniprotkb_GO_0010200_AND_taxonomy_id_37_2025_03_10.tsv")
#GO:0004553
hydroly <- read_uniprot("inputs/uniprot/hydroly_uniprotkb_GO_0004553_AND_taxonomy_id_37_2025_03_10.tsv")
#GO:0007267
cell_cell <- read_uniprot("inputs/uniprot/cell_cell_uniprotkb_GO_0007267_AND_taxonomy_id_37_2025_03_10.tsv")
#GO:0030247
poly_binding <- read_uniprot("inputs/uniprot/poly_binding_uniprotkb_GO_0030247_AND_taxonomy_id_37_2025_03_10.tsv")
#GO:0010411
xyl_met <- read_uniprot("inputs/uniprot/xyl_met_uniprotkb_GO_0010411_AND_taxonomy_id_37_2025_03_10.tsv")
#GO:0016762
xyl_trans <- read_uniprot("inputs/uniprot/xyl_trans_uniprotkb_GO_0016762_AND_taxonomy_id_37_2025_03_10.tsv")
#GO:0009741
br_response <- read_uniprot("inputs/uniprot/br_uniprotkb_GO_0009741_AND_taxonomy_id_37_2025_03_10.tsv")
#GO:0003843
cal_synth <- read_uniprot("inputs/uniprot/cal_synth_uniprotkb_taxonomy_id_3702_AND_GO_00038_2025_03_25.tsv")
#GO:0009851
aux_synth <- read_uniprot("inputs/uniprot/aux_biosynth_uniprotkb_GO_0009851_AND_model_organism_2025_03_26.tsv")
#GO:0009733
aux_resp <- read_uniprot("inputs/uniprot/resp_aux_uniprotkb_GO_0009733_AND_taxonomy_id_37_2025_04_20.tsv")

sig_mods <- list(nuc_met= nuc_met, 
                 cyto_resp = cyto_resp,
                 cyto_bios = cyto_bios,
                 isopren = isopren, 
                 lipid_drop = lipid_drop, 
                 mono_storage = mono_storage, 
                 lipid_storage = lipid_storage, 
                 resp_alc = resp_alc, 
                 alc_bios = alc_bios, 
                 cell_resp_aba = cell_resp_aba, 
                 resp_aba = resp_aba, 
                 aba_bio = aba_bio, 
                 nitrate_trans = nitrate_trans$Locus,
                 sweet_trans = sweet_trans$Locus,
                 umami_trans = umami_trans$Locus,
                 autophagy = autophagy,
                 br_homeo = br_homeo,
                 cadmium = cadmium,
                 cell_death = cell_death,
                 metal_ion_resp = metal_ion_resp,
                 pcd = pcd,
                 senes = senes,
                 trans_metal = trans_metal,
                 zinc = zinc, 
                 brassino_bios = brassino_bios, 
                 photosynth = photosynth, 
                 hypoxia = hypoxia, 
                 cell_wall = cell_wall, 
                 mech_stim = mech_stim, 
                 sec_vesc = sec_vesc, 
                 chitin = chitin, 
                 hydroly = hydroly,
                 cell_cell = cell_cell, 
                 poly_binding = poly_binding, 
                 xyl_met = xyl_met, 
                 xyl_trans = xyl_trans, 
                 br_response = br_response, 
                 cal_synth = cal_synth,
                 aux_synth = aux_synth, 
                 aux_resp = aux_resp)

seu_list <- list(dap3, dap5, dap7, atlas)
seu_names <- list("DAP3", "DAP5", "DAP7", "ATLAS")

for(i in 1:length(seu_list)){
  #now adding mod scores for the complete gene sets
  seu <- AddModuleScore(
    object = seu_list[[i]],
    features = sig_mods,
    name = 'sig_mods'
  )
  
  #renaming the metadata

  new_cols <- c(grep("sig_mods",colnames(seu[[]]), invert = TRUE, value = TRUE), 
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
                "Auxin_response_0009851_modsc")
  
  length(new_cols) == length(colnames(seu[[]])) 
  
  colnames(seu[[]]) <- new_cols
  #this move just copies the modscores, getting rid of the old columns
  seu@meta.data <- seu@meta.data[, !grepl("sig_mods", colnames(seu@meta.data))]
  
  saveRDS(seu, paste0("outputs/", seu_names[[i]], "_merged_annotated_sigmods.rds"))
  
}




