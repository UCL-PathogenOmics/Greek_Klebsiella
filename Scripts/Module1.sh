#!/bin/bash

## Module 1 - Assemble genomes, detect plasmids, and run Kleborate

set -euo pipefail


INPUT="/FASTQ"
THREADS=16
OUTDIR="Assemblies"

mkdir -p "$OUTDIR/flye" "$OUTDIR/Medaka" "$OUTDIR/Contigs/" "$OUTDIR/snps" "$OUTDIR/matrix"

for R1 in "$INPUT"/*.fastq.gz; do
    base=$(basename "$R1" | sed 's/.fastq.gz')

	# Assemblies
	flye --nano-hq "$R1" --out-dir "$OUTDIR/flye/${base}" --threads "$THREADS"
	medaka_consensus -i "$R1" -d "$OUTDIR/flye/${base}/assembly.fasta" -o "$OUTDIR/Medaka/${base}" -t "$THREADS" -m r1041_e82_400bps_hac_v4.3.0 

	# MOB-SUITE
	
		
	# Separate Contigs into chromosome and plasmids
	python Separate_contigs.py "$OUTDIR/Medaka/${base}/consensus.fasta" "$OUTDIR/Contigs/${base}" 

	gtdbtk classify_wf \
 	 --genome_dir "$OUTDIR/Contigs/${base}" \
  	--out_dir GTDB_OUTPUT \
  	--extension fasta \
  	--cpus "$THREADS"


	checkm2 predict \
  	--input-directory INPUT_FASTA_DIR \
 	 --output-directory CHECKM2_OUTPUT \
 	 --threads THREADS \
 	 --force \
 	 --extension fasta

done

kleborate \
  --input INPUT_FASTA \
  --output KLEBORATE_OUTPUT
