## Spatial Genotyping

This repository contains the bash scripts used to clean the Space Ranger genotyping results for Visium HD data. The initial Space Ranger output must be cleaned because some reads contain the reverse complement of the genotyping probe sequence, which leads to missed/incorrect genotyping calls.

## Requirements

The cleaning script (scripts/hl_sr_gen_cleaning.sh) requires Bash (most versions should be compatible):\
GNU bash, version 5.1.8(1)-release (x86_64-redhat-linux-gnu)

## Running

Once the Space Ranger output bam is generated, we can then use this as the input for the subsequent cleaning step (data/test_sr_gen_bam.bam). This step identifies the reads assigned to a valid genotyping probe based off of the bam tags, then extracts the relevant read sections to compare to a reference file (data/std_rev.txt). This file contains the forward and reverse complement sequences for a portion of each probe, as well as the expected MUT/WT alleles and their position (i.e. base 25 or 26), which varies between the two orientations. Once orientation is determined, mismatched reads (i.e. assigned to MUT probe but is actually the reverse comp of the WT probe) are reassigned to their correct genotype, and the results are saved as .sam files.
\
\
On this subsetted bam file, this demo should take < 1 minute to run.
