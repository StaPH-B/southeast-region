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
PROJECT_DIR="/blue/bphl-georgia/tkiryutina/projects/260707-UGA-Training-CODE/ADD-ISOLATES"

## Derived directories
FASTQ_DIR="${PROJECT_DIR}/fastqs"
PHOENIX_DIR="${PROJECT_DIR}/phoenix"
RUN_DIR="${PROJECT_DIR}/bigbacter"; mkdir -p "$RUN_DIR"

## PHoeNIx summary containing Final_Taxa_ID
PHOENIX_SUMMARY="${PHOENIX_DIR}/output/Phoenix_Summary.tsv"

## BigBacter variables
SAMPLESHEET="${RUN_DIR}/samplesheet.csv"
BB_DB="/blue/bphl-georgia/tkiryutina/projects/260707-UGA-Training-CODE/bigbacter/db"
BB_OVERRIDE_CONFIG="/blue/bphl-georgia/tkiryutina/projects/260707-UGA-Training-CODE/bigbacter/bb_override.config"

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

## Create BigBacter samplesheet
echo "sample,taxa,assembly,fastq_1,fastq_2" > "$SAMPLESHEET"

## Loop through the phoenix summary IDs that are only A. baumannii
grep "Acinetobacter baumannii" "$PHOENIX_SUMMARY" | cut -f1 | while read -r sample; do
    assembly="${PHOENIX_DIR}/output/${sample}/assembly/${sample}.filtered.scaffolds.fa.gz"
    fastq1=$(ls "${FASTQ_DIR}/${sample}"*_1.fastq.gz)
    fastq2=$(ls "${FASTQ_DIR}/${sample}"*_2.fastq.gz)
    echo "${sample},Acinetobacter_baumannii,${assembly},${fastq1},${fastq2}" >> "$SAMPLESHEET"
done

## Run BigBacter
nextflow run DOH-JDJ0303/bigbacter-nf \
    -r v1.0.0 \
    -profile singularity \
    --input "$SAMPLESHEET" \
    --db "$BB_DB" \
    --outdir "${RUN_DIR}/results" \
    --max_cpus 16 \
    -c "$BB_OVERRIDE_CONFIG"

## Optional: If you have more isolates to add later, either run or re-run the previous command with:
#   --push true \
#   -resume
