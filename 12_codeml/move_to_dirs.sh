#breaking the pal2nal outputs into groups of 1000 for parallelization

#batching into 1000-file directories
batch_size=1000
batch_num=1
file_count=0

mkdir -p dir_$(printf "%04d" $batch_num)

# looping through all pal2nal outputs
for file in OG*_alignment.fasta_pal2nal; do
  if (( file_count >= batch_size )); then
    ((batch_num++))
    file_count=0
    mkdir -p dir_$(printf "%04d" $batch_num)
  fi

  #moving pal2nal outputs into the new directory
  mv "$file" dir_$(printf "%04d" $batch_num)/
  ((file_count++))
done
