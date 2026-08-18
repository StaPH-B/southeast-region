#!/bin/bash

##_!!_## Run line by line, or by each small code block!

###_______________________________________________________
### Prerequisite: install SRA Toolkit to run fasterq-dump

## Download the Linux release of SRA Toolkit from:
## https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit

## Make a personal bin directory if needed
mkdir -p "$HOME/bin"

## Move the downloaded tar.gz file into $HOME/bin, then:
cd "$HOME/bin"

## Extract the archive
tar -xzf sratoolkit.3.4.1-ubuntu64.tar.gz

## Add SRA Toolkit to your PATH for this session
export PATH="$HOME/bin/sratoolkit.3.4.1-ubuntu64/bin:$PATH"

## Optional: permanently add it to ~/.bashrc
echo 'export PATH="$HOME/bin/sratoolkit.3.4.1-ubuntu64/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

## Verify installation
fasterq-dump --version

###_______________________________________________________
### Download FASTQs from SRR IDs:

## Go to your working directory (choose your own path and replace)
# cd /path/to/your/project
cd /blue/bphl-georgia/tkiryutina/projects/260707-UGA-Training-CODE

## Create directory for downloaded FASTQs
mkdir -p fastqs

## Loop through SRR ID file and download
while read -r srr; do
    echo "Downloading $srr"

    fasterq-dump "$srr" \
        --split-files \
        --threads 8 \
        --outdir fastqs

    gzip fastqs/${srr}_*.fastq
done < srr_ids.txt

## Expected output from this loop: (commented out using ctrl+/)

# Downloading SRR33240013
# spots read      : 807,899
# reads read      : 1,615,798
# reads written   : 1,615,798
# Downloading SRR33240012
# spots read      : 601,104
# reads read      : 1,202,208
# reads written   : 1,202,208
# Downloading SRR33240001
# spots read      : 529,197
# reads read      : 1,058,394
# reads written   : 1,058,394
# Downloading SRR33239997
# spots read      : 557,384
# reads read      : 1,114,768
# reads written   : 1,114,768
# Downloading SRR33239996
# spots read      : 484,104
# reads read      : 968,208
# reads written   : 968,208
# Downloading SRR33239995
# spots read      : 562,313
# reads read      : 1,124,626
# reads written   : 1,124,626
# Downloading SRR33239994
# spots read      : 520,232
# reads read      : 1,040,464
# reads written   : 1,040,464
# Downloading SRR33239993
# spots read      : 861,945
# reads read      : 1,723,890
# reads written   : 1,723,890
# Downloading SRR33239992
# spots read      : 450,278
# reads read      : 900,556
# reads written   : 900,556
# Downloading SRR33239991
# spots read      : 618,594
# reads read      : 1,237,188
# reads written   : 1,237,188
# Downloading SRR33240011
# spots read      : 543,757
# reads read      : 1,087,514
# reads written   : 1,087,514
# Downloading SRR33240010
# spots read      : 508,474
# reads read      : 1,016,948
# reads written   : 1,016,948
# Downloading SRR33240009
# spots read      : 336,592
# reads read      : 673,184
# reads written   : 673,184
# Downloading SRR33240008
# spots read      : 419,175
# reads read      : 838,350
# reads written   : 838,350
# Downloading SRR33240007
# spots read      : 581,082
# reads read      : 1,162,164
# reads written   : 1,162,164
# Downloading SRR33240006
# spots read      : 546,746
# reads read      : 1,093,492
# reads written   : 1,093,492
# Downloading SRR33240005
# spots read      : 582,624
# reads read      : 1,165,248
# reads written   : 1,165,248
# Downloading SRR33240004
# spots read      : 377,527
# reads read      : 755,054
# reads written   : 755,054
# Downloading SRR33240003
# spots read      : 418,959
# reads read      : 837,918
# reads written   : 837,918
# Downloading SRR33240002
# spots read      : 355,092
# reads read      : 710,184
# reads written   : 710,184
# Downloading SRR33138450
# spots read      : 1,021,087
# reads read      : 2,042,174
# reads written   : 2,042,174
# Downloading SRR33240000
# spots read      : 381,595
# reads read      : 763,190
# reads written   : 763,190
# Downloading SRR33239999
# spots read      : 476,831
# reads read      : 953,662
# reads written   : 953,662
# Downloading SRR33239998
# spots read      : 479,052
# reads read      : 958,104
# reads written   : 958,104

