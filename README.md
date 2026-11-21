## Spatial Genotyping

This repository contains the bash scripts used to generate and clean the Space Ranger genotyping results for Visium HD data. The initial Space Ranger output must be cleaned because some reads contain the reverse complement of the genotyping probe sequence, which leads to missed/incorrect genotyping calls.

## Requirements

The initial Space Ranger genotyping script (example_spaceranger_count.sh) requires the 10X Space Ranger tool (most versions should be compatible):\
spaceranger-3.0.0 (https://www.10xgenomics.com/support/software/space-ranger/3.0)

The cleaning script (hl_sr_gen_clean.sh) requires Bash (most versions should be compatible):\
GNU bash, version 5.1.8(1)-release (x86_64-redhat-linux-gnu)

