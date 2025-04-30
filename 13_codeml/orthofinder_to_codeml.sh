#this bash script outlines the steps to perform genome-wide codeml in four Arabidopsis species

#get the arabidopsis lyrata, arabidopsis thaliana, brassica rapa, and Capsella grandiflora proteins
cd 13_codeml
mkdir cds
mkdir proteins
gffread  -g Athaliana_447_TAIR10.fa Athaliana_447_Araport11.gene.gff3 -x cds/thaliana_CDSs.fasta -y proteins/thaliana_Proteins.fasta
gffread  -g GCF_000309985.2_CAAS_Brap_v3.01_genomic.fna GCF_000309985.2_CAAS_Brap_v3.01_genomic.gff -x cds/brassica_CDSs.fasta -y proteins/brassica_Proteins.fasta
gffread  -g Cgrandiflora_266_v1.fa Cgrandiflora_266_v1.1.gene.gff3 -x cds/capsella_CDSs.fasta -y proteins/capsella_Proteins.fasta
gffread  -g Alyrata_384_v1.fa Alyrata_384_v2.1.gene.gff3 -x cds/lyrata_CDSs.fasta -y proteins/lyrata_Proteins.fasta

#cleaning up the FASTA headers, if necessary 
cd proteins
for f in *Proteins.fasta; do sed  '/^[^>]/s/\.//g' $f |awk '{print $1}' >Clean${f};done

#running orthofinder
sbatch -p 20 --mail-type=ALL --mem=64gb --cpus-per-task=8 --wrap "./orthofinder -a 32 -t 32 -f proteins"

#making an all-protein FASTA with an index
cd  proteins
cat *fasta >AllCleanProteins.fasta
cdbfasta AllCleanProteins.fasta

#making an all-CDS FASTA with an index
cd  ../cds
for f in *_CDSs.fasta; do awk '{print $1}' $f >Clean${f};done
cat  Clean* > AllCleanCDSs.fasta
cdbfasta AllCleanCDSs.fasta

#once orthofinder completes, prune the trees to generate trees with only one Arabidopsis thaliana gene and its best ortholog in 2-3 other species
cd  OrthoFinder/Results_Apr17/Gene_Trees 
Rscript ../../../prune_trees_best.R

#then, fix names, so that tree names and FASTA names match 
python ../../../fix_names.py

#prep for alignment, extracting proteins and CDSs from the gene trees
cd ../../../
for f in OrthoFinder/Results_Apr17/Gene_Trees/*_pruned.tree; do sed 's/(//g' $f |sed 's/,/\n/g' |sed 's/:.*//g' |sort|uniq|sed 's/Proteins_/\t/1'|cut -f 2  |cdbyank proteins/AllCleanProteins.fasta.cidx >${f}_Proteins.fasta;done
for f in OrthoFinder/Results_Apr17/Gene_Trees/*_pruned.tree; do sed 's/(//g' $f |sed 's/,/\n/g' |sed 's/:.*//g' |sort|uniq|sed 's/Proteins_/\t/1'|cut -f 2  |cdbyank cds/AllCleanCDSs.fasta.cidx >${f}_CDSs.fasta;done

#running clustalo for proteins in gene trees
cd OrthoFinder/Results_Apr17/Gene_Trees

for f in *Proteins.fasta; do
  base=${f%Araport11.447_pruned.tree_Proteins.fasta}
  echo "clustalo -i $f -o ${base}_alignment.fasta" 
done > generateAlignment.sh

#running alignment
bash generateAlignment.sh

#now creating codon alignments
for f in *alignment.fasta; do 
  base=${f%._alignment.fasta}
  echo "perl ../../../pal2nal.pl $f ${base}.Araport11.447_pruned.tree_CDSs.fasta -output paml -nogap > ${f}_pal2nal" 
done > pal2nal.sh

# Split the script into chunks of 100 lines each for parallelization
split -l 100 pal2nal.sh pal2nal_batch_

#wrap each in a SLURM script
for f in pal2nal_batch_*; do
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

for f in pal2nal_batch_*.slurm; do
  sbatch $f
done

#the number of files gets unwieldy, so I batched the pal2nal outputs into folders of 1k with move_to_dirs.sh script, then work within "dir_0001" etc
bash ../../../move_to_dirs.sh

#now, generating codeml scripts for all pal2nal files: 

for dir in dir_*; do
  echo "Processing $dir"
  cd "$dir" || continue
  
  for f in OG*_alignment.fasta_pal2nal; do 
    base=${f%._alignment.fasta_pal2nal}
    tree="../../${base}.Araport11.447_pruned.tree"
    outfile="${base}.codeml.out"
    ctlfile="${base}.codeml.ctl"

    cat <<EOF > $ctlfile
        seqfile = ../$f
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

  cd ..
done

#moving all control files into gene-specific subdirectories

for dir in dir_*; do
  echo "Organizing codeml.ctl files in $dir"
  cd "$dir" || continue

  for f in OG*.codeml.ctl; do 
    base=${f%.codeml.ctl}
    mkdir -p "${base}_dir"
    mv "$f" "${base}_dir/codeml.ctl"
  done

  cd ..
done

for d in dir_*; do
    if [ -d "$d" ]; then
        #looping through the subdirectories within each "dir_XXXX"
        for subd in "$d"/*/; do
            if [ -d "$subd" ] && [ -f "$subd/codeml.ctl" ]; then
                # generating a codeml command for each directory, too
                echo "cd $PWD/$subd && codeml codeml.ctl" >> IterateCodeml.sh
            fi
        done
    fi
done

#then batch the codeml script
split -l 100 IterateCodeml.sh codeml_batch_

#making a slurm submission script for codeml batches
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

chmod +x codeml_batch_*.slurm

for f in codeml_batch_*.slurm; do
  sbatch -p 20 $f
done

#get all stats with parse_codeml_lrt.R


