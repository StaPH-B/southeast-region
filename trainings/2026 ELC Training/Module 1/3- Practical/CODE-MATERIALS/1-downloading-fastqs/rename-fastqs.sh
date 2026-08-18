#!/bin/bash

### Rename SRR IDs to Sample (WGS) IDs in FASTQs.

## If needed, copy fastqs to preserve original downloads.
cp -r fastqs/ fastqs-srr-original

## Rename downloaded FASTQs from SRR IDs to sample names
tail -n +2 sample_ids.txt | while IFS=$'\t' read -r sample accession; do
    mv -v "fastqs/${accession}_1.fastq.gz" "fastqs/${sample}_1.fastq.gz"
    mv -v "fastqs/${accession}_2.fastq.gz" "fastqs/${sample}_2.fastq.gz"
done

## Example outputs commented out using ctrl+/

# renamed 'fastqs/SRR33240013_1.fastq.gz' -> 'fastqs/2025LY00031R_1.fastq.gz'
# renamed 'fastqs/SRR33240013_2.fastq.gz' -> 'fastqs/2025LY00031R_2.fastq.gz'
# renamed 'fastqs/SRR33240012_1.fastq.gz' -> 'fastqs/2025LY00033_1.fastq.gz'
# renamed 'fastqs/SRR33240012_2.fastq.gz' -> 'fastqs/2025LY00033_2.fastq.gz'
# renamed 'fastqs/SRR33240001_1.fastq.gz' -> 'fastqs/2025LY00034_1.fastq.gz'
# renamed 'fastqs/SRR33240001_2.fastq.gz' -> 'fastqs/2025LY00034_2.fastq.gz'
# renamed 'fastqs/SRR33239997_1.fastq.gz' -> 'fastqs/2025LY00035_1.fastq.gz'
# renamed 'fastqs/SRR33239997_2.fastq.gz' -> 'fastqs/2025LY00035_2.fastq.gz'
# renamed 'fastqs/SRR33239996_1.fastq.gz' -> 'fastqs/2025LY00036_1.fastq.gz'
# renamed 'fastqs/SRR33239996_2.fastq.gz' -> 'fastqs/2025LY00036_2.fastq.gz'
# renamed 'fastqs/SRR33239995_1.fastq.gz' -> 'fastqs/2025LY00037_1.fastq.gz'
# renamed 'fastqs/SRR33239995_2.fastq.gz' -> 'fastqs/2025LY00037_2.fastq.gz'
# renamed 'fastqs/SRR33239994_1.fastq.gz' -> 'fastqs/2025LY00038_1.fastq.gz'
# renamed 'fastqs/SRR33239994_2.fastq.gz' -> 'fastqs/2025LY00038_2.fastq.gz'
# renamed 'fastqs/SRR33239993_1.fastq.gz' -> 'fastqs/2025LY00039_1.fastq.gz'
# renamed 'fastqs/SRR33239993_2.fastq.gz' -> 'fastqs/2025LY00039_2.fastq.gz'
# renamed 'fastqs/SRR33239992_1.fastq.gz' -> 'fastqs/2025LY00040_1.fastq.gz'
# renamed 'fastqs/SRR33239992_2.fastq.gz' -> 'fastqs/2025LY00040_2.fastq.gz'
# renamed 'fastqs/SRR33239991_1.fastq.gz' -> 'fastqs/2025LY00041_1.fastq.gz'
# renamed 'fastqs/SRR33239991_2.fastq.gz' -> 'fastqs/2025LY00041_2.fastq.gz'
# renamed 'fastqs/SRR33240011_1.fastq.gz' -> 'fastqs/2025LY00042_1.fastq.gz'
# renamed 'fastqs/SRR33240011_2.fastq.gz' -> 'fastqs/2025LY00042_2.fastq.gz'
# renamed 'fastqs/SRR33240010_1.fastq.gz' -> 'fastqs/2025LY00043_1.fastq.gz'
# renamed 'fastqs/SRR33240010_2.fastq.gz' -> 'fastqs/2025LY00043_2.fastq.gz'
# renamed 'fastqs/SRR33240009_1.fastq.gz' -> 'fastqs/2025LY00044_1.fastq.gz'
# renamed 'fastqs/SRR33240009_2.fastq.gz' -> 'fastqs/2025LY00044_2.fastq.gz'
# renamed 'fastqs/SRR33240008_1.fastq.gz' -> 'fastqs/2025LY00045_1.fastq.gz'
# renamed 'fastqs/SRR33240008_2.fastq.gz' -> 'fastqs/2025LY00045_2.fastq.gz'
# renamed 'fastqs/SRR33240007_1.fastq.gz' -> 'fastqs/2025LY00046_1.fastq.gz'
# renamed 'fastqs/SRR33240007_2.fastq.gz' -> 'fastqs/2025LY00046_2.fastq.gz'
# renamed 'fastqs/SRR33240006_1.fastq.gz' -> 'fastqs/2025LY00047_1.fastq.gz'
# renamed 'fastqs/SRR33240006_2.fastq.gz' -> 'fastqs/2025LY00047_2.fastq.gz'
# renamed 'fastqs/SRR33240005_1.fastq.gz' -> 'fastqs/2025LY00048_1.fastq.gz'
# renamed 'fastqs/SRR33240005_2.fastq.gz' -> 'fastqs/2025LY00048_2.fastq.gz'
# renamed 'fastqs/SRR33240004_1.fastq.gz' -> 'fastqs/2025LY00049_1.fastq.gz'
# renamed 'fastqs/SRR33240004_2.fastq.gz' -> 'fastqs/2025LY00049_2.fastq.gz'
# renamed 'fastqs/SRR33240003_1.fastq.gz' -> 'fastqs/2025LY00050_1.fastq.gz'
# renamed 'fastqs/SRR33240003_2.fastq.gz' -> 'fastqs/2025LY00050_2.fastq.gz'
# renamed 'fastqs/SRR33240002_1.fastq.gz' -> 'fastqs/2025LY00051_1.fastq.gz'
# renamed 'fastqs/SRR33240002_2.fastq.gz' -> 'fastqs/2025LY00051_2.fastq.gz'
# renamed 'fastqs/SRR33138450_1.fastq.gz' -> 'fastqs/2025LY00052R_1.fastq.gz'
# renamed 'fastqs/SRR33138450_2.fastq.gz' -> 'fastqs/2025LY00052R_2.fastq.gz'
# renamed 'fastqs/SRR33240000_1.fastq.gz' -> 'fastqs/2025LY00053_1.fastq.gz'
# renamed 'fastqs/SRR33240000_2.fastq.gz' -> 'fastqs/2025LY00053_2.fastq.gz'
# renamed 'fastqs/SRR33239999_1.fastq.gz' -> 'fastqs/2025LY00054_1.fastq.gz'
# renamed 'fastqs/SRR33239999_2.fastq.gz' -> 'fastqs/2025LY00054_2.fastq.gz'
# renamed 'fastqs/SRR33239998_1.fastq.gz' -> 'fastqs/2025LY00055_1.fastq.gz'
# renamed 'fastqs/SRR33239998_2.fastq.gz' -> 'fastqs/2025LY00055_2.fastq.gz'

