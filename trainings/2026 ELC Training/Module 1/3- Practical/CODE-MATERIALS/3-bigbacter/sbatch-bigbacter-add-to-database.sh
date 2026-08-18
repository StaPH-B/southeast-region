#!/bin/bash
#SBATCH --job-name=BigBacter
#SBATCH --account=bphl-umbrella
#SBATCH --qos=bphl-umbrella
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tatyana.kiryutina@dph.ga.gov
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=bigbacter.%j.out
#SBATCH --error=bigbacter.%j.err

## Indicate project directory
PROJECT_DIR="/blue/bphl-georgia/tkiryutina/projects/260707-UGA-Training-CODE"

## Derived directories
FASTQ_DIR="${PROJECT_DIR}/fastqs"
PHOENIX_DIR="${PROJECT_DIR}/phoenix"
RUN_DIR="${PROJECT_DIR}/bigbacter"

## BigBacter variables
SAMPLESHEET="${RUN_DIR}/samplesheet.csv"
BB_DB="${RUN_DIR}/db"
BB_OVERRIDE_CONFIG="${RUN_DIR}/bb_override.config"

## Create BigBacter run directory
mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

## Load modules
module purge
module load nextflow/24.04.2

## Set cache/temp dirs
export NXF_SINGULARITY_CACHEDIR="/blue/bphl-georgia/tkiryutina/singularity/bigbacter-containers"
export TMPDIR="/tmp/"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

## Run BigBacter again to push results to the database
## Notice the "--push true" and "-resume" flags were added
nextflow run DOH-JDJ0303/bigbacter-nf \
    -r v1.0.0 \
    -profile singularity \
    --input "$SAMPLESHEET" \
    --db "$BB_DB" \
    --outdir "${RUN_DIR}/results" \
    --max_cpus 16 \
    -c "$BB_OVERRIDE_CONFIG" \
    --push true \
    -resume



