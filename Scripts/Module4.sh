#!/usr/bin/env bash

## Module 4 - Reference mapping, variant calling and recombination analysis for CG transmission clusters 

set -euo pipefail

READS="/FASTQ/Cluster1"
OUTDIR="Cluster1"
THREADS=16
REF="Klebsiella_reference.fasta"

mkdir -p "$OUTDIR/minimap" "$OUTDIR/vcf" "$OUTDIR/gubbins" "$OUTDIR/snps" "$OUTDIR/matrix" "$OUTDIR/consensus"

# Minimap2 + BAM processing
for R1 in "$READS"/*.fastq.gz; do
    base=$(basename "$R1" | sed 's/.fastq.*//')

    minimap2 -t "$THREADS" -a -x map-ont "$REF" "$R1" | \
    samtools sort -@ "$THREADS" -o "$OUTDIR/minimap/${base}.sorted.bam"
    
    samtools index "$OUTDIR/minimap/${base}.sorted.bam"
done

# Bcftools variant calling
for bam in "$OUTDIR/minimap"/*.sorted.bam; do
    sample=$(basename "$bam" .sorted.bam)

    bcftools mpileup -Ou -f "$REF" "$bam" | \
    bcftools call -mv -Oz -o "$OUTDIR/vcf/${sample}.vcf.gz"

    bcftools index "$OUTDIR/vcf/${sample}.vcf.gz"
done

for vcf in "$OUTDIR"/vcf/*.vcf.gz; do
    sample=$(basename "$vcf" .vcf.gz)

    bcftools consensus -f "$REF" -H 1 -M N "$vcf" > "$OUTDIR/consensus/${sample}.fasta"
done

# Combine consensus fastas for gubbins
cat "$OUTDIR"/consensus/*.fasta > "$OUTDIR/Cluster1_multi.fasta"

# Gubbins
run_gubbins.py --prefix "$OUTDIR/gubbins/Cluster1" --threads "$THREADS" --tree-builder iqtree "$OUTDIR/Cluster1_multi.fasta"

# SNP-dists
snp-dists -b "$OUTDIR/gubbins/Cluster1.filtered_polymorphic_sites.fasta" > "$OUTDIR/Cluster1_SNPdist.txt"

# Make transmission clusters
python Transmission_clustering.py --input "$OUTDIR/Cluster1_SNPdist.txt" --thresholds 21 --output "$OUTDIR/Cluster1_clusters" --symmetric
