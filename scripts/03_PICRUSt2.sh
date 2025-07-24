#!/bin/bash

# Usage: ./03_PICRUSt2.sh <Th>
# Example: ./03_PICRUSt2.sh 0.1

if [ $# -ne 1 ]; then
  echo "Usage: $0 <Th>"
  echo "Example: $0 0.1"
  exit 1
fi

Th=$1  # threshold from command line argument

# Define input/output paths based on Th
DATA_DIR="./data"
ASV_TABLE="${DATA_DIR}/ASV_table_${Th}.tsv"
FASTA_FILE="${DATA_DIR}/2023_16S_GorBEEa_prj_ASVs_filt_c.${Th}.fa"
OUTPUT_DIR="${DATA_DIR}/picrust2_output_${Th}"
CPU=4

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Define log and error files inside the output directory
LOG_FILE="${OUTPUT_DIR}/picrust2_${Th}.log"
ERR_FILE="${OUTPUT_DIR}/picrust2_${Th}.err"

# Run PICRUSt2 pipeline
picrust2_pipeline.py -s "${FASTA_FILE}" -i "${ASV_TABLE}" -o "${OUTPUT_DIR}" -p "${CPU}" > "${LOG_FILE}" 2> "${ERR_FILE}"

