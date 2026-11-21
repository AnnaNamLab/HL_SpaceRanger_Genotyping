## Spatial Genotyping

This repository contains the bash scripts used to generate and clean the Space Ranger genotyping results for Visium HD data. The initial Space Ranger output must be cleaned because some reads contain the reverse complement of the genotyping probe sequence, which leads to missed/incorrect genotyping calls.

## Requirements

The initial genotyping script (scripts/example_spaceranger_count.sh) requires the 10X Space Ranger tool (most versions should be compatible):\
spaceranger-3.0.0 (https://www.10xgenomics.com/support/software/space-ranger/3.0)

The cleaning script (scripts/hl_sr_gen_cleaning.sh) requires Bash (most versions should be compatible):\
GNU bash, version 5.1.8(1)-release (x86_64-redhat-linux-gnu)

## Running

The genotyping can be run like any other Space Ranger count run, just with the appropriate genotyping fastqs and probe barcodes.\
\
Once the Space Ranger output bam is generated, we can then use this as the input for the subsequent cleaning step (data/test_sr_gen_bam.bam). This step identifies the reads assigned to a valid genotyping probe based off of the bam tags, then extracts the relevant read sections to compare to a reference file (data/std_rev.txt). This file contains the forward and reverse complement sequences for a portion of each probe, as well as the expected MUT/WT alleles and their position (i.e. base 25 or 26), which varies between the two orientations. 
