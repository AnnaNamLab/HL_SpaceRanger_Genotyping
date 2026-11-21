#!/bin/bash

#SBATCH --partition=scu-gpu   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=SR_GEN_Output
#SBATCH --time=48:00:00   # HH/MM/SS
#SBATCH --mem=64G   


spaceranger_path=/path/to/spaceranger
spatial_base_files=/path/to/input/files

spaceranger count --id="Sample_VHD_SR_Output"_$(date +%F) \
	  --transcriptome=$spaceranger_path/reference_genomes/refdata-gex-GRCh38-2020-A \
          --probe-set=$spatial_base_files/20240530_VisiumHD_HLGenotyping_probe_set.csv \
          --fastqs=$spatial_base_files/fastqs/GEN/fastqs \
          --image=$spatial_base_files/images/HE_images/Sample_20X-03.jpeg \
          --cytaimage=$spatial_base_files/images/CytAssist_images/Sample_run1_assay_CAVG10562_2024-05-01_16-38-17_H1-P6FCYFC_1714548378_CytAssist/CAVG10562_2024-05-01_16-56-18_2024-05-01_16-38-17_H1-P6FCYFC_A1_Sample.tif \
          --loupe-alignment=$spatial_base_files/images/Alignment_files/Sample_H1-P6FCYFC-A1-fiducials-image-registration.json \
	  --slide=H1-P6FCYFC \
          --area=A1 \
	  --create-bam=true \
	  --no-probe-filter  
exit



