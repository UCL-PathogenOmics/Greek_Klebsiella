#!/usr/bin/env bash

## Module 3 - per CG pangenome and core alignment tree

set -euo pipefail

DIR="Cluster"
THREADS=16

# PROKKA 

KINGDOM="Bacteria"
GENUS="Klebsiella"
SPECIES="pneumoniae"
mkdir -p "${DIR}/Prokka"

for f in "${DIR}/Chromosome_assemblies/"*.fasta; do
    base=$(basename "$f" .fasta)
    prokka \
    --outdir  "${DIR}/Prokka/${base}" \
    --prefix  "${base}" \
    --force \
    --cpus    "${THREADS}" \
    --kingdom "$KINGDOM" \
    --genus   "$GENUS" \
    --species "$SPECIES" \
    "$f"
done

# PANAROO
mkdir -p "${DIR}/Panaroo"
panaroo \
    -i "${DIR}/Prokka/"*/*.gff \
    -o "${DIR}/Panaroo" \
    -t "${THREADS}" \
    --clean-mode strict \
    --remove-invalid-genes \
    --alignment pan

# GUBBINS
mkdir -p "${DIR}/Gubbins"
ALIGN="${DIR}/Panaroo/core_gene_alignment.aln"
run_gubbins.py \
    --threads "${THREADS}" \
    --prefix "${DIR}/Gubbins/CG" \
    "$ALIGN"

# MASK ALIGNMENT
python mask_sites.py "$ALIGN" "${DIR}/Gubbins/CG.recombination_predictions.gff" "${DIR}/Panaroo/masked.fasta"

# SNP-SITES
mkdir -p "${DIR}/SNPs"
snp-sites -c -o "${DIR}/SNPs/SNPs.fasta" "${DIR}/Panaroo/masked.fasta"

# IQTREE3 tree build
SEED=${SEED:-12345}
mkdir -p "${DIR}/Tree"

iqtree3 \
    -s "${DIR}/SNPs/SNPs.fasta" \
    -m MFP+ASC \
    -nt "${THREADS}" \
    -B 1000 -bnni \
    -safe \
    -keep-ident \
    -seed "$SEED" \
    -pre "${DIR}/Tree/CG" \
    -redo

# IQTREE3 recalculate branch length
iqtree3 \
    -s  "${DIR}/SNPs/SNPs.fasta" \
    -te "${DIR}/Tree/CG.treefile" \
    -m  MFP+ASC \
    -nt "${THREADS}" \
    -seed "$SEED" \
    -keep-ident \
    -pre "${DIR}/Tree/CG_relength"

# Run BactDate in R to date tree (need CG_dates.txt file in Tree folder in "%Y/%m/%d" format)
Rscript run_bactdate.R "${DIR}/SNPs/SNPs.fasta" "${DIR}/Tree/CG_relength.treefile" "${DIR}/Tree/CG_dates.txt" "${DIR}/Tree/CG_timedtree.nwk"
