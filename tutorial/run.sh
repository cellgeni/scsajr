#¡/bin/bash
# Script to run nf-scsajr with a sample dataset

# Prepare input files
mkdir -p input

url="https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-vdj/4.0.0/Parent_SC5v1_Human_Glioblastoma/Parent_SC5v1_Human_Glioblastoma"
curl "${url}_raw_feature_bc_matrix.h5" -o input/raw_feature_bc_matrix.h5
curl "${url}_filtered_feature_bc_matrix.h5" -o input/filtered_feature_bc_matrix.h5
curl "${url}_possorted_genome_bam.bam" -o input/alignment.bam
curl "${url}_possorted_genome_bam.bam.bai" -o input/alignment.bam.bai
curl "${url}_analysis.tar.gz" -o input/clusters.tar.gz

tar -xzf input/clusters.tar.gz 
echo -e "sample\t$(pwd)/input/alignment.bam" > input/samples.tsv
awk -F, 'BEGIN { OFS="\t" }; NR>1 {print  "sample",$1,$2}' analysis/clustering/graphclust/clusters.csv > input/barcodes.tsv

# Get nf-scsajr from GitHub
git clone https://github.com/cellgeni/nf-scsajr

# Run nf-scsajr
nextflow run ./nf-scsajr/main.nf \
  --SAMPLEFILE ./input/samples.tsv \
  --BARCODEFILE ./input/barcodes.tsv \
  --minsamples 1 \
  --ref ./nf-scsajr/ref/human_2020A_chr \
  --outdir ./nf-scsajr_output
