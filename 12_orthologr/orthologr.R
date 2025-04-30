library(orthologr)
library(biomartr)

setwd("12_orthologr")

lyrata_thaliana <- dNdS(query_file= "inputs/Athaliana_447_Araport11.cds.fa",
                         subject_file    = "inputs/Alyrata_384_v2.1.cds.fa",
                         delete_corrupt_cds = TRUE,
                         ortho_detection = "RBH",
                         aa_aln_type     = "pairwise",
                         aa_aln_tool     = "NW",
                         codon_aln_tool  = "pal2nal",
                         dnds_est.method = "NG",
                         store_locally = TRUE, 
                         comp_cores      = 8) 

write.csv(lyrata_thaliana, "outputs/lyrata_thaliana_NG.csv")




