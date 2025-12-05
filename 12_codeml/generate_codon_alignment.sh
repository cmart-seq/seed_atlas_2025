#!/bin/bash
logfile="pal2nal_run.log"

echo "=== Starting PAL2NAL batch run: $(date) ===" > "$logfile"

for file in *_alignment.fasta; do
  base=${file%_alignment.fasta}  # get orthogroup base name
  echo "Processing $file" | tee -a "$logfile"
  
  if ! perl pal2nal.pl "$file" "${base}_tree.txt_CDSs.fasta" -output paml -nogap > "${file}_pal2nal" 2>>"$logfile"; then
    echo "failed on $file" | tee -a "$logfile"
    echo "  Check: $file and ${base}_tree.txt_CDSs.fasta" | tee -a "$logfile"
    echo "=== Error encountered, aborting at $(date) ===" | tee -a "$logfile"
    exit 1
  fi

  echo "completed $file" | tee -a "$logfile"
done

echo "=== All done: $(date) ===" | tee -a "$logfile"