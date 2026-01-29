#!/usr/bin/env bash
#
# Initial preprocessing pipeline for Aedes aegypti RNA-seq data
#
# Steps:
#   1) Update system and install dependencies (one-time setup)
#   2) Verify installed tool versions
#   3) Run quality control with FastQC
#   4) Prepare reference genome and annotation files
#   5) Build HISAT2 index
#   6) Align paired-end reads with HISAT2
#   7) Convert SAM → sorted BAM with samtools
#   8) Count reads per gene with featureCounts
#
# Notes:
#   - Adjust paths, file names, and number of threads (THREADS)
#     according to your project organization.
#   - This example assumes NCBI AaegL5 reference:
#       GCF_002204515.2_AaegL5.0_genomic.fna(.gz)
#       GCF_002204515.2_AaegL5.0_genomic.gtf(.gz)
#   - Commands with `sudo` must be run in a regular terminal
#     with appropriate privileges (not inside a restricted job).
#

###############################################################################
# 1. System update and tool installation (run once on a new machine)
###############################################################################

# (Optional) Update package index and upgrade existing packages
# sudo apt update
# sudo apt upgrade -y

# Install required tools
# sudo apt install -y fastqc samtools hisat2 subread

###############################################################################
# 2. Verify tool installation
###############################################################################

# These commands should return version information.
# If any of them fail, check the installation.

# fastqc --version
# hisat2 --version
# samtools --version
# featureCounts -v

###############################################################################
# 3. Quality control with FastQC
###############################################################################
# Run this inside the folder containing FASTQ files for a given batch/experiment.

# Example project layout:
#   project/
#     references/
#     batch_1/
#       sample1_R1.fastq.gz
#       sample1_R2.fastq.gz
#       sample2_R1.fastq.gz
#       sample2_R2.fastq.gz
#       ...
#     batch_2/
#       ...

# Go to the batch folder:
# cd /path/to/project/batch_1

# Create an output folder for FastQC reports
# mkdir -p fastqc_reports

# Run FastQC on all compressed FASTQ files in this folder
# fastqc *.fastq.gz \
#   -o fastqc_reports \
#   -t 8

###############################################################################
# 4. Reference genome and annotation preparation (NCBI AaegL5)
###############################################################################
# Run this in the references folder of the project.

# cd /path/to/project/references

# Decompress NCBI files if they are still in .gz format
# gunzip GCF_002204515.2_AaegL5.0_genomic.fna.gz
# gunzip GCF_002204515.2_AaegL5.0_genomic.gtf.gz

# Check that the decompressed files exist:
# ls GCF_002204515.2_AaegL5.0_genomic.fna
# ls GCF_002204515.2_AaegL5.0_genomic.gtf

###############################################################################
# 5. Build HISAT2 index
###############################################################################
# Still inside the references folder.

# hisat2-build GCF_002204515.2_AaegL5.0_genomic.fna AaegL5_index

# This will generate several index files:
#   AaegL5_index.1.ht2, AaegL5_index.2.ht2, ..., AaegL5_index.8.ht2

###############################################################################
# 6. Alignment with HISAT2 (example for a single sample)
###############################################################################
# Run this inside the batch folder where FASTQ files are stored.

# cd /path/to/project/batch_1

# Define the base sample name (without R1/R2 suffix or extension)
# For example:
# SAMPLE="control_rep1"
# Then the FASTQ files are:
#   control_rep1_R1.fastq.gz
#   control_rep1_R2.fastq.gz

# SAMPLE="control_rep1"
# R1="${SAMPLE}_R1.fastq.gz"
# R2="${SAMPLE}_R2.fastq.gz"

# Path to HISAT2 index (references folder)
# REF_DIR="/path/to/project/references"

# Number of threads to use
# THREADS=8

# Paired-end alignment
# hisat2 \
#   -p "${THREADS}" \
#   -x "${REF_DIR}/AaegL5_index" \
#   -1 "${R1}" \
#   -2 "${R2}" \
#   -S "${SAMPLE}.sam"

# Output: SAM file with alignments: ${SAMPLE}.sam

###############################################################################
# 7. SAM to sorted BAM conversion with samtools
###############################################################################
# Convert the SAM file generated above into a sorted BAM file.

# samtools view -bS "${SAMPLE}.sam" \
#   | samtools sort -o "${SAMPLE}_sorted.bam"

# (Optional but recommended) Index the sorted BAM file
# samtools index "${SAMPLE}_sorted.bam"

# Key alignment file for downstream counting:
#   ${SAMPLE}_sorted.bam

###############################################################################
# 8. Gene-level read counting with featureCounts
###############################################################################
# This step links sorted BAM files to the GTF annotation.

# featureCounts \
#   -T 4 \
#   -a "${REF_DIR}/GCF_002204515.2_AaegL5.0_genomic.gtf" \
#   -o "${SAMPLE}_counts.txt" \
#   "${SAMPLE}_sorted.bam"

# Parameters:
#   -T 4  : number of threads
#   -a    : GTF annotation file
#   -o    : output file with raw gene counts
#
# The file `${SAMPLE}_counts.txt` will be one of the inputs
# for building the count matrix in R (DESeq2).

###############################################################################
# 9. Example loop for multiple samples
###############################################################################
# You can automate steps 6–8 for multiple samples with a simple loop:
#
# cd /path/to/project/batch_1
#
# REF_DIR="/path/to/project/references"
# THREADS=8
#
# for SAMPLE in control_rep1 control_rep2 control_rep3 \
#                treatmentA_rep1 treatmentA_rep2 treatmentA_rep3; do
#
#   R1="${SAMPLE}_R1.fastq.gz"
#   R2="${SAMPLE}_R2.fastq.gz"
#
#   # Alignment
#   hisat2 -p "${THREADS}" \
#     -x "${REF_DIR}/AaegL5_index" \
#     -1 "${R1}" \
#     -2 "${R2}" \
#     -S "${SAMPLE}.sam"
#
#   # SAM → sorted BAM
#   samtools view -bS "${SAMPLE}.sam" \
#     | samtools sort -o "${SAMPLE}_sorted.bam"
#
#   samtools index "${SAMPLE}_sorted.bam"
#
#   # Gene-level counts
#   featureCounts -T 4 \
#     -a "${REF_DIR}/GCF_002204515.2_AaegL5.0_genomic.gtf" \
#     -o "${SAMPLE}_counts.txt" \
#     "${SAMPLE}_sorted.bam"
#
# done
#
# Make sure file names in this loop match your actual FASTQ and BAM names.
