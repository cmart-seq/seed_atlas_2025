#signalp6 
cd 10_peptide_enrichment/inputs/signalP
#first, filtering all seqs to 250 or less
seqkit seq -M 250 ../../../00_references/Araport11_pep_20220914.fa > Araport11_pep_20220914_max250.fa

#shortening the FASTA headers
sed -E 's/(^>[^|]*)\|.*/\1/' Araport11_pep_20220914_max250.fa > Araport11_pep_20220914_max250_abbr_names.fa

#running signalp6 on the resulting sequences
conda activate signal6
sbatch -p 20 --mem=64gb --mail-type=ALL --job-name signal6 --wrap "signalp6 --fastafile Araport11_pep_20220914_max250_abbr_names.fa --output_dir signalP_results2 --model_dir models --format txt --organism eukarya --mode slow"
