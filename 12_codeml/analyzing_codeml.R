library(dplyr)
library(readr)
library(clusterProfiler)
library(org.At.tair.db)
library(Seurat)
library(seqinr)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(viridis)
library(tidyr)

#function for codeml output
parse_positive_selection <- function(filename, prob_threshold = 0.95) {
  basename <- str_split(filename, "\\/|_dir")[[1]][5]
  file_path <- paste0(grep(filename,top_dirs, value = T),"/",basename, ".codeml.out"  )
  lines <- readLines(file_path)
  
  #locate where BEB results start
  start_line <- grep("(BEB)", lines)
  if (length(start_line) == 0) {
    stop("No BEB section found in the file.")
  }
  
  end_line_m2a <- grep("(NEB)", lines)
  
  # Extract lines after BEB header that look like selected site entries
  site_line_m2a <- lines[(start_line[1] + 1):end_line_m2a[2]]
  site_line_m8 <- lines[(start_line[2] + 1):length(lines)]
  
  site_lines_m2a <- site_line_m2a[grep("^\\s*\\d+\\s+", site_line_m2a)]
  site_lines_m8 <- site_line_m8[grep("^\\s*\\d+\\s+", site_line_m8)]
  
  # Parse lines
  site_data_m2a <- do.call(rbind, lapply(site_lines_m2a, function(line) {
    fields <- unlist(strsplit(trimws(line), "\\s+"))
    site_pos <- as.numeric(fields[1])
    aa <- fields[2]
    prob <- as.numeric(gsub("[^0-9.]", "", fields[3]))
    mark <- ifelse(grepl("\\*", fields[3]), TRUE, FALSE)
    return(data.frame(site = site_pos, amino_acid = aa, prob = prob, selected = mark))
  }))
  
  site_data_m8 <- do.call(rbind, lapply(site_lines_m8, function(line) {
    fields <- unlist(strsplit(trimws(line), "\\s+"))
    site_pos <- as.numeric(fields[1])
    aa <- fields[2]
    prob <- as.numeric(gsub("[^0-9.]", "", fields[3]))
    mark <- ifelse(grepl("\\*", fields[3]), TRUE, FALSE)
    return(data.frame(site = site_pos, amino_acid = aa, prob = prob, selected = mark))
  }))
  
  #label models
  site_data_m2a$model <- "M2a"
  site_data_m8$model <- "M8"
  
  site_data_all <- rbind(site_data_m2a, site_data_m8)
  
  # Filter by probability threshold
  positively_selected <- subset(site_data_all, prob >= prob_threshold)
  return(positively_selected)
}


######## processing codeml outputs ######## 
#reading in parsed codeml outputs, adding orthogroup information
ortho_table <- read_delim("orthofinder_inputs4/OrthoFinder/Results_Oct31/Orthogroups/Orthogroups.tsv", "\t") %>% dplyr::select(Orthogroup, Cleanthaliana_Proteins)
codeml <- read_csv("/lab/solexa_gehring/carly/snRNA_seed_dev/clean_analysis/orthofinder/clean/orthofinder_inputs4/codeml_LRT_summary_full.csv")
codeml$Orthogroup <- sapply(codeml$Gene, function(x){str_split(x, "_")[[1]][1]})
codeml <- left_join(codeml, ortho_table, by ="Orthogroup")
codeml$tair <- sapply(codeml$Cleanthaliana_Proteins, function(x){str_split(x, "_|[.]")[[1]][1]})

#splitting the protein names to add isoform
codeml$isoform <- sapply(codeml$Cleanthaliana_Proteins, function(x){str_split(x, ".Araport11.447")[[1]][1]})

#getting the lengths for all isoforms/genes so that I can filter for the longest isoform
#getting CDSs
cds_thaliana <- read.fasta("../../12_div_peptides_old/inputs/divergence/Athaliana_447_Araport11.cds.fa")
#getting the length
tair_lengths <- data.frame(tair= names(cds_thaliana), 
                           length = getLength(cds_thaliana))
#adding the isoform to the codeml table
codeml <- left_join(codeml, tair_lengths, by = c("isoform" = "tair"))

#saving a table with just the key info
codeml_simplified <- dplyr::select(codeml, Gene, tair, isoform, M2a_vs_M1a_LRT, M2a_vs_M1a_pval, M2a_vs_M1a_padj)
codeml_simplified <- codeml_simplified[,c("tair", "isoform",  "M2a_vs_M1a_LRT", "M2a_vs_M1a_pval", "M2a_vs_M1a_padj")]
write_csv(codeml_simplified, "outputs/codeml_simplified_final.csv")
#str_split("OG0000012_Proteins_AT1G49240.1_dir", "_")[[1]][3]

#beginning filtering, first by the longest isoform: 
#filtering for longest isoform
codeml_ordered <- group_by(codeml, tair) %>% arrange(desc(length)) %>% slice_head(n=1)
#saving all, selecting columns for supp table
codeml_ordered_write <- dplyr::select(codeml_ordered, tair, isoform, length, M2a_vs_M1a_LRT, M2a_vs_M1a_pval, M2a_vs_M1a_padj)
write_csv(codeml_ordered_write, "all_sc_orthos_M2a_vs_M1a.csv")

#now filtering for significant genes
codeml_sig <- filter(codeml_ordered, M2a_vs_M1a_padj < 0.05)
codeml_write <- dplyr::select(codeml_sig, tair, isoform, M2a_vs_M1a_LRT, M2a_vs_M1a_padj)
write_csv(codeml_write, "all_sig_sc_orthos_M2a_vs_M1a.csv")

####### looking at the intersection of interpro domain outs for all DE genes and those with evidence of positive selection ########

#interproscan results
all_interpro <- read_csv("all_m2am1a_gene_interpro_domains.csv", col_names = FALSE)
all_interpro$tair <- sapply(all_interpro$X1, function(x){strsplit(x, "[.]")[[1]][1]})

#getting the files for site parsing
top_dirs <- list.dirs(path = "orthofinder_inputs4/OrthoFinder/Results_Oct31/Gene_Trees", full.names = TRUE, recursive = TRUE)
codeml_sig_filenames <- grep(paste(codeml_sig$Gene, collapse="|"), top_dirs, value = TRUE)

#parsing positive selection, keeping threshhold at 0.95
pos_residues <- lapply(codeml_sig_filenames, function(x){parse_positive_selection(x)})
#naming the residues
names(pos_residues)<-codeml_sig_filenames
#making a table
pos_residues_df <- do.call(rbind, pos_residues)
#extracting the orthogroups
pos_residues_df$orthogroup <- sapply(rownames(pos_residues_df), function(x){strsplit(x, "orthofinder_inputs4/OrthoFinder/Results_Oct31/Gene_Trees/|[.]")[[1]][2]})
#mapping orthogroups to TAIRs
og_tair_map <- dplyr::select(codeml_sig, Gene, tair)
pos_residues_df <- left_join(pos_residues_df, og_tair_map, by = c("orthogroup" = "Gene"))
#focusing on M2a
pos_residues_df <- filter(pos_residues_df, model == "M2a")
#filtering 
pos_residues_df_seed <- filter(pos_residues_df, tair %in% atlas_marks_timed$gene)
#getting # unique sites by adding a gene-position column
pos_residues_df_seed$gene_pos <- paste(pos_residues_df_seed$tair, pos_residues_df_seed$site)
length(unique(pos_residues_df_seed$gene_pos))
#[1] 359

#looking at domain intersections
selected_domains <- pos_residues_df_seed %>%
  # join by tair first
  inner_join(all_interpro, by = "tair") %>%
  # filter on interval overlap
  filter(!is.na(X7), !is.na(X8)) %>%
  filter(site >= X7 & site <= X8)

selected_domains_write <- dplyr::select(selected_domains, site, amino_acid, prob, model, tair, X7, X8, X4,X5, X6)

#saving
write_csv(selected_domains_write, "selected_domain_families.csv")



  
