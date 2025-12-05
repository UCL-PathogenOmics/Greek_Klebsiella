#!/usr/bin/env bash
set -euo pipefail


##### Script to make detect novel AMR genes.

## Variables (per sample per gene analysis)
SAMPLE="samplename"
GENE="lpoA"
FASTA="/kleborate_output/${SAMPLE}.fasta"
GENE_REF="${GENE}/reference.fasta"
OUTPUT="Output/${GENE}_${SAMPLE}"
mkdir -p "$OUTPUT"

# Prodigal gene prediction
prodigal -i "$FASTA" -o "$OUTPUT/${SAMPLE}.gff" -f gff -p meta -q
awk '$3=="CDS"{print $1"\t"$4"\t"$5}' "$OUTPUT/${SAMPLE}.gff" \
  > "$OUTPUT/${SAMPLE}_genes.tsv"

# Create 10× repeated gene
SEQ=$(awk 'BEGIN{ORS=""} /^>/ {next} {print toupper($0)} END{print "\n"}' "$GENE_REF")
FA10="$OUTPUT/${GENE}_10x.fasta"
: > "$FA10"
for i in {1..10}; do printf ">${GENE}_copy$i\n%s\n" "$SEQ" >> "$FA10"; done

# Map gene→assembly and extract first mapping position
minimap2 -a "$FASTA" "$FA10" > "$OUTPUT/map.sam"

samtools view -bS "$OUTPUT/map.sam" | samtools sort -o "$OUTPUT/map.bam"
samtools index "$OUTPUT/map.bam"
rm -f "$OUTPUT/map.sam"

samtools view "$OUTPUT/map.bam" \
  | awk '($3!="*"){print $3,$4,$4+length($10)-1; exit}' \
  > "$OUTPUT/mapped_pos.tsv"

# Get closest Prodigal CDS to mapping
read CONTIG START END < "$OUTPUT/mapped_pos.tsv"
MID=$(( (START+END)/2 ))

BEST=$(awk -v C="$CONTIG" -v M="$MID" '
  BEGIN{bestD=1e18}
  {
    if($1==C){
      mid=($2+$3)/2; d=M-mid; if(d<0)d=-d
      if(d<bestD){bestD=d; best=$1":"$2"-"$3}
    }
  }
  END{print best}
' "$OUTPUT/${SAMPLE}_genes.tsv")

samtools faidx "$FASTA" "$BEST" > "$OUTPUT/${GENE}_${SAMPLE}_extracted.fasta"

# Add reference gene and align
REFSEQ=$(awk 'BEGIN{ORS=""} /^>/ {next} {print toupper($0)} END{print "\n"}' "$GENE_REF")
printf ">ref_%s\n%s\n" "$GENE" "$REFSEQ" >> "$OUTPUT/${GENE}_${SAMPLE}_extracted.fasta"

mafft --auto "$OUTPUT/${GENE}_${SAMPLE}_extracted.fasta" \
  > "$OUTPUT/${GENE}_${SAMPLE}_aligned.fasta"
