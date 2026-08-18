#!/bin/bash
#SBATCH --job-name=phoenix
#SBATCH --account=bphl-umbrella
#SBATCH --qos=bphl-umbrella
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tatyana.kiryutina@dph.ga.gov
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=phoenix.%j.out
#SBATCH --error=phoenix.%j.err

## Indicate project directory
PROJECT_DIR="/blue/bphl-georgia/tkiryutina/projects/260707-UGA-Training-CODE"

## Derive directory variables
FASTQ_DIR="${PROJECT_DIR}/fastqs"; mkdir -p "$FASTQ_DIR"
RUN_DIR="${PROJECT_DIR}/phoenix"; mkdir -p "$RUN_DIR"
SAMPLESHEET="${RUN_DIR}/samplesheet.csv"
KRAKEN2DB="/blue/bphl-georgia/tkiryutina/databases/k2_standard_08_GB_20260226.tar.gz"
export NXF_SINGULARITY_CACHEDIR="/blue/bphl-georgia/tkiryutina/databases/singularity_cache/phoenix-v2.3.2"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

## Move into run directory
cd "$RUN_DIR"

## Create samplesheet header
echo "sample,fastq_1,fastq_2" > "$SAMPLESHEET"

## Populate samplesheet with each pair of FASTQs
for r1 in "$FASTQ_DIR"/*_1.fastq.gz; do
    sample=$(basename "$r1" _1.fastq.gz)
    r2="${r1/_1.fastq.gz/_2.fastq.gz}"
    echo "${sample},${r1},${r2}" >> "$SAMPLESHEET"
done

## Load required modules
module load nextflow/25.10.4
module load singularity

## Run PHoeNIx
nextflow run CDCgov/phoenix \
    -r v2.3.2 \
    -profile singularity \
    --mode PHOENIX \
    --input "$SAMPLESHEET" \
    --kraken2db "$KRAKEN2DB" \
    --outdir "${RUN_DIR}/output" \
    -work-dir "${RUN_DIR}/work" \
    -resume



