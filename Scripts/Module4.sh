#!/usr/bin/env bash

## Module 4 - Reference mapping, variant calling and recombination analysis for CG transmission clusters 

set -euo pipefail

READS="/FASTQ/Cluster1"
OUTDIR="Cluster1"
THREADS=16
REF="Klebsiella_reference.fasta"

mkdir -p "$OUTDIR/minimap" "$OUTDIR/vcf" "$OUTDIR/gubbins" "$OUTDIR/snps" "$OUTDIR/matrix"

# Minimap2 + BAM processing
for R1 in "$READS"/*.fastq.gz; do
    base=$(basename "$R1" | sed 's/_R1.*//')

    minimap2 -t "$THREADS" -ax sr "$REF" "$R1" "$R2" | \
        samtools view -b | \
        samtools sort -o "$OUTDIR/minimap/${base}.sorted.bam"

    samtools index "$OUTDIR/minimap/${base}.sorted.bam"
done

# Bcftools variant calling
for bam in "$OUTDIR/minimap"/*.bam; do
    sample=$(basename "$bam" .sorted.bam)

    bcftools mpileup -Ou -f "$REF" "$bam" | \
        bcftools call -mv -Oz -o "$OUTDIR/vcf/${sample}.vcf.gz"

    bcftools index "$OUTDIR/vcf/${sample}.vcf.gz"
done

bcftools merge -Oz -o "$OUTDIR/merged.vcf.gz" "$OUTDIR/vcf"/*.vcf.gz
bcftools index "$OUTDIR/merged.vcf.gz"
bcftools consensus -f "$REF" -H 1 -M N \
    "$OUTDIR/merged.vcf.gz" \
    > "$OUTDIR/consensus.fasta"
    
# Gubbins
run_gubbins.py --prefix Cluster1 --threads "$THREADS" --tree-builder iqtree "$OUTDIR/consensus.fasta"

# SNP-dists
snp-dists -b Cluster1.filtered_polymorphic_sites.fasta > Cluster1_SNPdist.txt

# Make transmission clusters
python Transmission_clustering.py --input Cluster1_SNPdist.txt --thresholds 21 --output Cluster1_clusters
