#!/bin/bash

#!/bin/sh

#SBATCH --job-name=03_PICRUSt2
#SBATCH --error=data/logs/%x-%j.err
#SBATCH --output=data/logs/%x-%j.out

#SBATCH --partition=general # This is the default partition
#SBATCH --qos=regular
#SBATCH --cpus-per-task=48
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=10:00:00
#SBATCH --mem=24000

# Time and memory consumption estimates are orientative. Please adjust them according to you requirments.

# Loas modules
module load PICRUSt2/2.6.1-foss-2022b

# Define the threshold variable
Th=0.5  # Adjust as needed

# Define input and output file names
ASV_TABLE="ASV_table_${Th}.tsv"
FASTA_FILE="2023_16S_GorBEEa_prj_ASVs_filt_${Th}.fa"
OUTPUT_DIR="picrust2_output_${Th}"
CPU="48"

# Run the PICRUSt2 pipeline
picrust2_pipeline.py -s ${FASTA_FILE} -i ${ASV_TABLE} -o ${OUTPUT_DIR} -p ${CPU}
