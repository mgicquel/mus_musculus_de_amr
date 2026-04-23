#!/bin/bash

#SBATCH --array=1-875%50
#SBATCH --cpus-per-task=8
#SBATCH --job-name=rgi_annotation
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=10G

# Load necessary environment and activate conda environment
source /home/user/.bashrc
conda activate rgi

# Author: Víctor Hugo Jarquín Díaz
# Group: AG Forslund
# Aim: Script for ARG annotation of reads using RGIbwt
# Location/date: MDC - 01.2024

# Run script from root directory of all samples (e.g. ../data/HMHZ_Shotgun)
# Define variables

SAMPLE_LIST="sample_list.txt"
INPUT_DIR="$PWD/preprocessed"

# Main loop
cat "$SAMPLE_LIST" | while read -r sample; do
    echo "Processing sample: $sample"
    
    #Indicate the output directory
    annotation_dir="$PWD/output/${sample}"
    mkdir -p "$annotation_dir"
    
    # Check if output already exists
    if [[ -f "$annotation_dir/${sample}.gene_mapping_data.txt" ]]; then
        echo "Annotation output already exists for sample $sample. Skipping."
        continue
    fi
    
    # Check input files existence
    if [[ ! -f "$INPUT_DIR/${sample}_filtered_reads.pair.1.fq.gz" || ! -f "$INPUT_DIR/${sample}_filtered_reads.pair.2.fq.gz" ]]; then
        echo "Error: Input files not found for sample $sample"
        continue
    fi

    ## Run the RGI pipeline
    
    rgi bwt -1 "$INPUT_DIR/${sample}_filtered_reads.pair.1.fq.gz" -2 "$INPUT_DIR/${sample}_filtered_reads.pair.2.fq.gz" -a 'bowtie2' -n 8 -o  "$annotation_dir/${sample}" --clean
done
