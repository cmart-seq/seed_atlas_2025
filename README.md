## snRNA-seq and evolutionary analysis code for Martin et al. 2025, _"A transcriptional atlas of early Arabidopsis seed development suggests mechanisms for inter-tissue coordination"_

Below are the command line arguments for running all of the scripts implemented in Martin et al 2025.
Steps 1-11 are for snRNA-seq data preprocessing analysis, while step 12 and 13 are for genome-wide dN/dS and codeml analysis. 
Steps 1-11 have to be run sequentially, and all of the scripts assume a file structure like the one in this repository. 

### Step 1a: perform initial SoupX background correction, low quality cell removal, and clustering
```
cd 01_per_library_QCs
Rscript per_library_QCs.R --wd . --sample DAP7_2 --tenXlibname DAP7_colcol_4 --ccGenes inputs/Menges2005_Mphase_Sphase.csv
Rscript per_library_QCs.R --wd . --sample DAP7_3 --tenXlibname DAP7_colcol_5 --ccGenes inputs/Menges2005_Mphase_Sphase.csv
Rscript per_library_QCs.R --wd . --sample DAP5_2 --tenXlibname 231018_whole_seed_DAP5v2_colcol --ccGenes inputs/Menges2005_Mphase_Sphase.csv
Rscript per_library_QCs.R --wd . --sample DAP5_3 --tenXlibname DAP5_colcol_5 --ccGenes inputs/Menges2005_Mphase_Sphase.csv
Rscript per_library_QCs.R --wd . --sample DAP3_1a --tenXlibname DAP3_colcol_1 --ccGenes inputs/Menges2005_Mphase_Sphase.csv
Rscript per_library_QCs.R --wd . --sample DAP3_1b --tenXlibname DAP3_colcol_2 --ccGenes inputs/Menges2005_Mphase_Sphase.csv
Rscript per_library_QCs.R --wd . --sample DAP3_2 --tenXlibname DAP3_colcol_3 --ccGenes inputs/Menges2005_Mphase_Sphase.csv
```
### Step 1b: after inspecting clusters, remove low quality cells and detect doublets for each library
```
Rscript score_doublets.R --wd . --seu_path outputs/DAP3_1a/DAP3_1asoupx.rds --sample DAP3_1a --rm_clusters 16
Rscript score_doublets.R --wd . --seu_path outputs/DAP3_1b/DAP3_1bsoupx.rds --sample DAP3_1b --rm_clusters 0
Rscript score_doublets.R --wd . --seu_path outputs/DAP3_2/DAP3_2soupx.rds --sample DAP3_2
Rscript score_doublets.R --wd . --seu_path outputs/DAP5_2/DAP5_2soupx.rds --sample DAP5_2
Rscript score_doublets.R --wd . --seu_path outputs/DAP5_3/DAP5_3soupx.rds --sample DAP5_3
Rscript score_doublets.R --wd . --seu_path outputs/DAP7_2/DAP7_2soupx.rds --sample DAP7_2 --rm_clusters 16
Rscript score_doublets.R --wd . --seu_path outputs/DAP7_3/DAP7_3soupx.rds --sample DAP7_3
```
### Step 2: merge libraries into timepoint datasets and determine an appropriate # PCs and whether they should be integrated by timepoint 
```
cd ../02_merge_libraries_filtering_and_batch_effects

#ROUND1: 
#DAP3, no integration
Rscript merge_libraries_and_batch_effects_harmony.R --wd . --round 1 --merge_only TRUE --dataset DAP3 --libraries DAP3_1a DAP3_1b DAP3_2 --n_VarGenes 3000

#DAP5, YES integration
Rscript merge_libraries_and_batch_effects_harmony.R --wd . --round 1 --merge_only FALSE  --int_variable bio_rep --dataset DAP5 --libraries DAP5_2 DAP5_3 --n_VarGenes 3000

#DAP7, no integration
Rscript merge_libraries_and_batch_effects_harmony.R --wd . --round 1 --merge_only TRUE --dataset DAP7 --libraries DAP7_2 DAP7_3 --n_VarGenes 3000

#ROUND2: 
#DAP3, no integration
Rscript merge_libraries_and_batch_effects_harmony.R --wd . --round 2 --merge_only TRUE --dataset DAP3 --n_VarGenes 3000 --ndims 60

#DAP5, YES integration
Rscript merge_libraries_and_batch_effects_harmony.R --wd . --round 2 --merge_only FALSE --int_variable bio_rep --dataset DAP5 --n_VarGenes 3000 --ndims 59

#DAP7, no integration
Rscript merge_libraries_and_batch_effects_harmony.R --wd . --round 2 --merge_only TRUE --dataset DAP7 --n_VarGenes 3000 --ndims 57
```
### Step 3: clustering parameter sweep
```
cd ../03_clustering
Rscript clustering.R --wd . --integrated FALSE --timepoint DAP3 --dataname DAP3 --low_res 1 --high_res 2 --inc 0.1 --cze_marks inputs
Rscript clustering.R --wd . --integrated TRUE --timepoint DAP5 --dataname DAP5 --low_res 1 --high_res 2 --inc 0.1 --cze_marks inputs
Rscript clustering.R --wd . --integrated FALSE --timepoint DAP7 --dataname DAP7 --low_res 1 --high_res 2 --inc 0.1 --cze_marks inputs
```
### Step 4a: apply manual annotation to the de novo clusters to generate the L3 annotations
```
cd  ../04_manual_annotation
Rscript manual_annotation.R  --wd . --dataset DAP3  --harmony FALSE --seupath ../03_clustering/outputs/DAP3/DAP3_clustered.rds --mpath inputs/DAP3_level_3_manual_metadata.csv
Rscript manual_annotation.R  --wd . --dataset DAP5  --harmony TRUE --seupath ../03_clustering/outputs/DAP5/DAP5_clustered.rds --mpath inputs/DAP5_level_3_manual_metadata.csv
Rscript manual_annotation.R  --wd . --dataset DAP7  --harmony FALSE --seupath ../03_clustering/outputs/DAP7/DAP7_clustered.rds --mpath inputs/DAP7_level_3_manual_metadata.csv
```
### Step 4b: subclustering to ID endosperm subtupes
```
cd  subclustering_for_CZE_subtypes
#IDing the nodule, nodule-like, apical cyst, and basal cyst
04_manual_annotation/subclustering_for_CZE_subtypes/ID_RALFL3_and_nodule.R

#re-run the annotation script which appends the endosperm clusters ID'ed in the subclustering analysis 
cd  ../
Rscript manual_annotation.R  --wd . --dataset DAP3_wcze_subs  --harmony FALSE --seupath subclustering_for_CZE_subtypes/DAP3_clustered_wcze_subs.rds --mpath inputs/DAP3_level_3_wcze_subs.csv
Rscript manual_annotation.R  --wd . --dataset DAP5_wcze_subs  --harmony TRUE --seupath subclustering_for_CZE_subtypes/DAP5_clustered_wcze_subs.rds --mpath inputs/DAP5_level_3_wcze_subs.csv
```
### Step 5a: merge tissues across timepoints to enable pseudotime analysis
```
#peforming integration by timepoint and regressing out cell cycle genes
cd ../05_across_timepoints/01_merging
sbatch -p 20 --mem=180gb --mail-type=ALL --job-name l2merging --wrap "Rscript level_2_merging_harmony.R"
```
### Step 5b: Pseudotime, identify root after merging (step 5a)
```
#for MCE, CZE, EMB
cd ../05_across_timepoints/02_pseudotime
Rscript level_2_pseudotime_merged_timepoints.R

#The chalazal cyst and PEN 5 DAP pseudotime analysis had to be performed separetely, since it involved two L2 annotations
#relies only on the 5 DAP object
cd ../../02_pseudotime/czen_spectrum_pseudo
Rscript harmony_chalazal_endosperm_trajectory.R
```

### Step 6: merge all annotated timepoints into one atlas dataset
```
cd ../../06_atlas_merging
Rscript atlas_merging_rpca.R --wd . --dataset ATLAS --n_VarGenes 4000
```

### Step 7: differential expression for all annotated timepoints 
```
cd ../07_differential_expression
Rscript differential_expression.R --wd .  --sample DAP3
Rscript differential_expression.R --wd .  --sample DAP5
Rscript differential_expression.R --wd .  --sample DAP7
Rscript differential_expression.R --wd .  --sample ATLAS
```
### Step 8: GO analysis for all DE genes for all clusters
```
cd ../08_GO_analysis
#intermediate filtering
Rscript clusterprofiler_intermediate.R --mks ../07_differential_expression/outputs/DAP3/DAP3_level_3_annotation_full_markers.csv --sample DAP3_level_3_intermediate --wd .
Rscript clusterprofiler_intermediate.R --mks ../07_differential_expression/outputs/DAP5/DAP5_level_3_annotation_full_markers.csv --sample DAP5_level_3_intermediate --wd .
Rscript clusterprofiler_intermediate.R --mks ../07_differential_expression/outputs/DAP7/DAP7_level_3_annotation_full_markers.csv --sample DAP7_level_3_intermediate --wd .

Rscript clusterprofiler_intermediate.R --mks ../07_differential_expression/outputs/DAP3/DAP3_level_2_annotation_markers.csv --sample DAP3_level_2_intermediate --wd .
Rscript clusterprofiler_intermediate.R --mks ../07_differential_expression/outputs/DAP5/DAP5_level_2_annotation_markers.csv --sample DAP5_level_2_intermediate --wd .
Rscript clusterprofiler_intermediate.R --mks ../07_differential_expression/outputs/DAP7/DAP7_level_2_annotation_markers.csv --sample DAP7_level_2_intermediate --wd .

Rscript clusterprofiler_intermediate.R --mks ../07_differential_expression/outputs/DAP3/DAP3_level_1_annotation_markers.csv --sample DAP3_level_1_intermediate --wd .
Rscript clusterprofiler_intermediate.R --mks ../07_differential_expression/outputs/DAP5/DAP5_level_1_annotation_markers.csv --sample DAP5_level_1_intermediate --wd .
Rscript clusterprofiler_intermediate.R --mks ../07_differential_expression/outputs/DAP7/DAP7_level_1_annotation_markers.csv --sample DAP7_level_1_intermediate --wd .

Rscript clusterprofiler_intermediate.R --mks ../07_differential_expression/outputs/ATLAS/ATLAS_level_2_annotation_timed_final_markers.csv --sample ATLAS_level_2_timed_intermediate --wd .
```
### Step 9-11: Enrichment analysis 
```
#Signalling GO term enrichment
cd ../../09_signalling_transport_gene_enrichment
sbatch -p 20 --mem=64gb --mail-type=ALL --job-name signalling_enrichment --wrap "Rscript signalling_transport_gene_enrichment.R"

#Peptide enrichment
cd ../10_peptide_enrichment
sbatch -p 20  --dependency=afterok:7361923 --mem=64gb --mail-type=ALL --job-name proten --wrap "Rscript peptide_enrichment.R"

#Protein catabolism enrichment
cd ../11_protein_catabolism_enrichment
sbatch -p 20 --dependency=afterok:7361924 --mem=64gb --mail-type=ALL --job-name procat --wrap "Rscript protein_catabolism_enrichment.R"
```
### dN/dS analysis ###

### genome-wide codeml ###
see the ```orthofinder_to_codeml.sh``` script for generating alignments and gene trees for codeml analysis







