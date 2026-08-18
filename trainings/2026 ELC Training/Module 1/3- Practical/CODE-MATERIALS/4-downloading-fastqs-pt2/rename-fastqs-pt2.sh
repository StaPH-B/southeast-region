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

# renamed 'fastqs/SRR36255135_1.fastq.gz' -> 'fastqs/2025LY-00187_1.fastq.gz'
# renamed 'fastqs/SRR36255135_2.fastq.gz' -> 'fastqs/2025LY-00187_2.fastq.gz'
# renamed 'fastqs/SRR36255124_1.fastq.gz' -> 'fastqs/2025LY-00188_1.fastq.gz'
# renamed 'fastqs/SRR36255124_2.fastq.gz' -> 'fastqs/2025LY-00188_2.fastq.gz'
# renamed 'fastqs/SRR36255106_1.fastq.gz' -> 'fastqs/2025LY-00194_1.fastq.gz'
# renamed 'fastqs/SRR36255106_2.fastq.gz' -> 'fastqs/2025LY-00194_2.fastq.gz'
# renamed 'fastqs/SRR36255126_1.fastq.gz' -> 'fastqs/2025LY-00204_1.fastq.gz'
# renamed 'fastqs/SRR36255126_2.fastq.gz' -> 'fastqs/2025LY-00204_2.fastq.gz'
# renamed 'fastqs/SRR36255125_1.fastq.gz' -> 'fastqs/2025LY-00205_1.fastq.gz'
# renamed 'fastqs/SRR36255125_2.fastq.gz' -> 'fastqs/2025LY-00205_2.fastq.gz'
# renamed 'fastqs/SRR36788088_1.fastq.gz' -> 'fastqs/2025LY-00206_1.fastq.gz'
# renamed 'fastqs/SRR36788088_2.fastq.gz' -> 'fastqs/2025LY-00206_2.fastq.gz'
# renamed 'fastqs/SRR36255118_1.fastq.gz' -> 'fastqs/2025LY-00212_1.fastq.gz'
# renamed 'fastqs/SRR36255118_2.fastq.gz' -> 'fastqs/2025LY-00212_2.fastq.gz'
# renamed 'fastqs/SRR33984279_1.fastq.gz' -> 'fastqs/2025LY00093_1.fastq.gz'
# renamed 'fastqs/SRR33984279_2.fastq.gz' -> 'fastqs/2025LY00093_2.fastq.gz'
# renamed 'fastqs/SRR34157520_1.fastq.gz' -> 'fastqs/2025LY00106_1.fastq.gz'
# renamed 'fastqs/SRR34157520_2.fastq.gz' -> 'fastqs/2025LY00106_2.fastq.gz'
# renamed 'fastqs/SRR34157519_1.fastq.gz' -> 'fastqs/2025LY00107_1.fastq.gz'
# renamed 'fastqs/SRR34157519_2.fastq.gz' -> 'fastqs/2025LY00107_2.fastq.gz'
# renamed 'fastqs/SRR34157518_1.fastq.gz' -> 'fastqs/2025LY00108_1.fastq.gz'
# renamed 'fastqs/SRR34157518_2.fastq.gz' -> 'fastqs/2025LY00108_2.fastq.gz'
# renamed 'fastqs/SRR34157517_1.fastq.gz' -> 'fastqs/2025LY00109_1.fastq.gz'
# renamed 'fastqs/SRR34157517_2.fastq.gz' -> 'fastqs/2025LY00109_2.fastq.gz'
# renamed 'fastqs/SRR34323655_1.fastq.gz' -> 'fastqs/2025LY00111_1.fastq.gz'
# renamed 'fastqs/SRR34323655_2.fastq.gz' -> 'fastqs/2025LY00111_2.fastq.gz'
# renamed 'fastqs/SRR34559556_1.fastq.gz' -> 'fastqs/2025LY-00121_1.fastq.gz'
# renamed 'fastqs/SRR34559556_2.fastq.gz' -> 'fastqs/2025LY-00121_2.fastq.gz'
# renamed 'fastqs/SRR34657797_1.fastq.gz' -> 'fastqs/2025LY-00129_1.fastq.gz'
# renamed 'fastqs/SRR34657797_2.fastq.gz' -> 'fastqs/2025LY-00129_2.fastq.gz'
# renamed 'fastqs/SRR34559555_1.fastq.gz' -> 'fastqs/2025LY-00122_1.fastq.gz'
# renamed 'fastqs/SRR34559555_2.fastq.gz' -> 'fastqs/2025LY-00122_2.fastq.gz'
# renamed 'fastqs/SRR33893821_1.fastq.gz' -> 'fastqs/2025LY00094_1.fastq.gz'
# renamed 'fastqs/SRR33893821_2.fastq.gz' -> 'fastqs/2025LY00094_2.fastq.gz'
# renamed 'fastqs/SRR37210572_1.fastq.gz' -> 'fastqs/2026LY-00018_1.fastq.gz'
# renamed 'fastqs/SRR37210572_2.fastq.gz' -> 'fastqs/2026LY-00018_2.fastq.gz'
# renamed 'fastqs/SRR37210564_1.fastq.gz' -> 'fastqs/2026LY-00023_1.fastq.gz'
# renamed 'fastqs/SRR37210564_2.fastq.gz' -> 'fastqs/2026LY-00023_2.fastq.gz'
# renamed 'fastqs/SRR37211280_1.fastq.gz' -> 'fastqs/2026LY-00034_1.fastq.gz'
# renamed 'fastqs/SRR37211280_2.fastq.gz' -> 'fastqs/2026LY-00034_2.fastq.gz'
# renamed 'fastqs/SRR37211273_1.fastq.gz' -> 'fastqs/2026LY-00036_1.fastq.gz'
# renamed 'fastqs/SRR37211273_2.fastq.gz' -> 'fastqs/2026LY-00036_2.fastq.gz'
# renamed 'fastqs/SRR37398810_1.fastq.gz' -> 'fastqs/2026LY-00070_1.fastq.gz'
# renamed 'fastqs/SRR37398810_2.fastq.gz' -> 'fastqs/2026LY-00070_2.fastq.gz'
# renamed 'fastqs/SRR37398817_1.fastq.gz' -> 'fastqs/2026LY-00082_1.fastq.gz'
# renamed 'fastqs/SRR37398817_2.fastq.gz' -> 'fastqs/2026LY-00082_2.fastq.gz'
# renamed 'fastqs/SRR37468327_1.fastq.gz' -> 'fastqs/2026LY-00129_1.fastq.gz'
# renamed 'fastqs/SRR37468327_2.fastq.gz' -> 'fastqs/2026LY-00129_2.fastq.gz'
# renamed 'fastqs/SRR37468326_1.fastq.gz' -> 'fastqs/2026LY-00130_1.fastq.gz'
# renamed 'fastqs/SRR37468326_2.fastq.gz' -> 'fastqs/2026LY-00130_2.fastq.gz'
# renamed 'fastqs/SRR37210571_1.fastq.gz' -> 'fastqs/2026LY-00019_1.fastq.gz'
# renamed 'fastqs/SRR37210571_2.fastq.gz' -> 'fastqs/2026LY-00019_2.fastq.gz'
# renamed 'fastqs/SRR37398811_1.fastq.gz' -> 'fastqs/2026LY-00069_1.fastq.gz'
# renamed 'fastqs/SRR37398811_2.fastq.gz' -> 'fastqs/2026LY-00069_2.fastq.gz'
# renamed 'fastqs/SRR30224866_1.fastq.gz' -> 'fastqs/2024LY00067_1.fastq.gz'
# renamed 'fastqs/SRR30224866_2.fastq.gz' -> 'fastqs/2024LY00067_2.fastq.gz'
# renamed 'fastqs/SRR34547107_1.fastq.gz' -> 'fastqs/2025LY-00117_1.fastq.gz'
# renamed 'fastqs/SRR34547107_2.fastq.gz' -> 'fastqs/2025LY-00117_2.fastq.gz'
# renamed 'fastqs/SRR29437671_1.fastq.gz' -> 'fastqs/2024LY00046_1.fastq.gz'
# renamed 'fastqs/SRR29437671_2.fastq.gz' -> 'fastqs/2024LY00046_2.fastq.gz'
# renamed 'fastqs/SRR29517191_1.fastq.gz' -> 'fastqs/2024LY00054_1.fastq.gz'
# renamed 'fastqs/SRR29517191_2.fastq.gz' -> 'fastqs/2024LY00054_2.fastq.gz'
# renamed 'fastqs/SRR30224867_1.fastq.gz' -> 'fastqs/2024LY00066_1.fastq.gz'
# renamed 'fastqs/SRR30224867_2.fastq.gz' -> 'fastqs/2024LY00066_2.fastq.gz'
# renamed 'fastqs/SRR37211287_1.fastq.gz' -> 'fastqs/2026LY-00047_1.fastq.gz'
# renamed 'fastqs/SRR37211287_2.fastq.gz' -> 'fastqs/2026LY-00047_2.fastq.gz'
# renamed 'fastqs/SRR37114319_1.fastq.gz' -> 'fastqs/2025LY-00225_1.fastq.gz'
# renamed 'fastqs/SRR37114319_2.fastq.gz' -> 'fastqs/2025LY-00225_2.fastq.gz'
# renamed 'fastqs/SRR37468329_1.fastq.gz' -> 'fastqs/2026LY-00127_1.fastq.gz'
# renamed 'fastqs/SRR37468329_2.fastq.gz' -> 'fastqs/2026LY-00127_2.fastq.gz'


