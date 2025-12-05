#I performed an initial run of dnds analysis using the proteomes and transcriptomes of 
# Arabidopsis thaliana and lyrata. That did not work because the 
#transcripts from lyrata were not complete. I am trying again using this workflow, 
#which uses only the genome assembly and the annotation gff, so there should
#not be any issues between the loci and the protein coding sequences

#https://bioinformaticsworkbook.org/dataAnalysis/ComparativeGenomics/Finding_Positive_Selection_With_Codeml.html#gsc.tab=0

#get the lyrata, arabidopsis, brassica rapa, and Capsella grandiflora proteins
cd orthofinder

#have to un-mask b rapa
seqkit seq -u GCF_000309985.2_CAAS_Brap_v3.01_genomic.fna > GCF_000309985.2_CAAS_Brap_v3.01_genomic_upper.fna

gffread  -g Athaliana_447_TAIR10.fa Athaliana_447_Araport11.gene.gff3 -x clean/cds/thaliana_CDSs.fasta -y clean/proteins/thaliana_Proteins.fasta
gffread  -g Cgrandiflora_266_v1.fa Cgrandiflora_266_v1.1.gene.gff3 -x clean/cds/capsella_CDSs.fasta -y clean/proteins/capsella_Proteins.fasta
gffread  -g Alyrata_384_v1.fa Alyrata_384_v2.1.gene.gff3 -x clean/cds/lyrata_CDSs.fasta -y clean/proteins/lyrata_Proteins.fasta
gffread  -g GCA_905216605.1_AARE701a_genomic.fna GCA_905216605.1_genomic.gff -x clean/cds/a_arenosa_CDSs.fasta -y clean/proteins/a_arenosa_Proteins.fasta


cd clean/proteins
for f in *Proteins.fasta; do sed  '/^[^>]/s/\.//g' $f |awk '{print $1}' >Clean${f};done

#now setting up an input folder
cd ../
mkdir orthofinder_inputs4
mv proteins/Cleancapsella_Proteins.fasta orthofinder_inputs4
mv proteins/Cleanlyrata_Proteins.fasta orthofinder_inputs4
mv proteins/Cleanthaliana_Proteins.fasta orthofinder_inputs4
mv proteins/Cleana_arenosa_Proteins.fasta orthofinder_inputs4

sbatch -p 20 --mail-type=ALL --mem=64gb --cpus-per-task=8 --wrap "./orthofinder -a 32 -t 32 -f clean/orthofinder_inputs4"

#all proteins
cd  clean/orthofinder_inputs4
cat *fasta > AllCleanProteins.fasta
cdbfasta AllCleanProteins.fasta

cd  ../cds
for f in *_CDSs.fasta; do awk '{print $1}' $f >Clean${f};done
cat  Clean* > AllCleanCDSs.fasta
cdbfasta AllCleanCDSs.fasta

#select single copy orthologs
cd ../orthofinder_inputs4/OrthoFinder/Results_Oct31/Single_Copy_Orthologue_Sequences
ls | wc -l
ls OG* > sc_orthos.txt
sed 's/\.fa$//' sc_orthos.txt > sc_orthos_nofa.txt
cat sc_orthos_nofa.txt | wc -l

#Show which tree files would be deleted
find . -maxdepth 1 -type f -name "*_tree.txt" -printf "%f\n" \
  | grep -vFf <(awk '{sub(/\r$/,""); if (NF) print $0 "_tree.txt"}' ../Single_Copy_Orthologue_Sequences/sc_orthos_nofa.txt) \
  | wc -l
 #11900
  
#then delete
find . -maxdepth 1 -type f -name "*_tree.txt" -printf "%f\n" \
  | grep -vFf <(awk '{sub(/\r$/,""); if (NF) print $0 "_tree.txt"}' ../Single_Copy_Orthologue_Sequences/sc_orthos_nofa.txt) \
  | xargs -r rm

ls | wc -l
#7187

#prep for alignment
cd orthofinder/
for f in clean/orthofinder_inputs4/OrthoFinder/Results_Oct31/Gene_Trees/*_tree.txt; do sed 's/(//g' $f |sed 's/,/\n/g' |sed 's/:.*//g' |sort|uniq|sed 's/Proteins_/\t/1'|cut -f 2  |cdbyank clean/orthofinder_inputs4/AllCleanProteins.fasta.cidx >${f}_Proteins.fasta;done
for f in clean/orthofinder_inputs4/OrthoFinder/Results_Oct31/Gene_Trees/*_tree.txt; do sed 's/(//g' $f |sed 's/,/\n/g' |sed 's/:.*//g' |sort|uniq|sed 's/Proteins_/\t/1'|cut -f 2  |cdbyank clean/cds/AllCleanCDSs.fasta.cidx >${f}_CDSs.fasta;done

#running clustalo 

cd clean/orthofinder_inputs4/OrthoFinder/Results_Oct31/Gene_Trees

for f in *_tree.txt_Proteins.fasta; do
  base=${f%_tree.txt_Proteins.fasta}
  echo "clustalo -i $f -o ${base}_alignment.fasta" 
done > generateAlignment.sh

# Split the script into chunks of 100 lines each
split -l 100 generateAlignment.sh alignment_batch_

#wrap each in a SLURM script
for f in alignment_batch_*; do
  cat <<EOF > ${f}.slurm
#!/bin/bash
#SBATCH --job-name=codeml_${f}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=01:00:00

bash $f
EOF
done

for f in alignment_batch_*.slurm; do
  sbatch -p 20 $f
done

#generating codon alignments for all trees and protein alignments with pal2nal
bash generate_codon_alignment.sh

#cleaning tree names
for t in *_tree.txt; do
  sed -E 's/Clean[a-zA-Z_]+_Proteins_//g' "$t" > "${t%.txt}_stripped.txt"
done

#!/bin/bash
#OG0012453_alignment.fasta_pal2nal

#making codeml control files for all orthogroup alignments
for f in *_alignment.fasta_pal2nal; do 
    base="${f%_alignment.fasta_pal2nal}"
    tree="../${base}_tree_stripped.txt"
    outfile="${base}.codeml.out"
    ctlfile="${base}.codeml.ctl"

cat > "$ctlfile" <<EOF
      seqfile = ../${base}_alignment.fasta_pal2nal
      treefile = $tree
      outfile = $outfile

      verbose = 1
      runmode = 0
      seqtype = 1
      CodonFreq = 2

      model = 0
      NSsites = 0 1 2 7 8

      icode = 0
      fix_kappa = 0
      kappa = 2

      fix_omega = 0
      omega = 0.4

      cleandata = 1
EOF
done

#moving all codeml inputs to individual directories
for f in OG*.codeml.ctl; do 
    base=${f%.codeml.ctl}
    mkdir -p "${base}_dir"
    mv "$f" "${base}_dir/codeml.ctl"
done


for d in *_dir; do
   # Check if it's a directory
    if [ -d "$d" ] && [ -f "$d/codeml.ctl" ]; then
        # Append the command to the SLURM script for each subdirectory
        echo "cd $PWD/$d && codeml codeml.ctl" >> IterateCodeml.sh
    fi
done

#then batch the codeml script
split -l 100 IterateCodeml.sh codeml_batch_

for f in codeml_batch_*; do
  cat <<EOF > ${f}.slurm
#!/bin/bash
#SBATCH --job-name=codeml_${f}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G

bash $f
EOF
done

for f in codeml_batch_*.slurm; do
  sbatch -p 20 $f
done

#next, parse to get stats with parse_codeml_lrt.R


