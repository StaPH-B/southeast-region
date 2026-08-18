#!/bin/bash
#SBATCH --job-name=srr-download
#SBATCH --account=bphl-umbrella
#SBATCH --qos=bphl-umbrella
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tatyana.kiryutina@dph.ga.gov
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --time=12:00:00
#SBATCH --output=srr-download.%j.log
#SBATCH --error=srr-download.%j.log

###_______________________________________________________
### Download FASTQs from SRR IDs (phase two: adding additional isolates to the outbreak)

## Run from a new work directory. For example: a directory "ADD-ISOLATES" from the previous run dir.
## "/blue/bphl-georgia/tkiryutina/projects/260707-UGA-Training-CODE/ADD-ISOLATES"

## Create directory for downloaded FASTQs
mkdir -p fastqs

## Loop through SRR ID file and download
while read -r srr; do
    echo "Downloading $srr"

    fasterq-dump "$srr" \
        --split-files \
        --threads 8 \
        --outdir fastqs

    pigz -p 8 fastqs/${srr}_*.fastq
done < srr_ids.txt


