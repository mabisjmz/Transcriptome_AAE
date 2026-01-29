# Aedes aegypti RNA-seq analysis pipeline

This repository contains a complete, reproducible pipeline for transcriptome analysis of *Aedes aegypti* mosquitoes, from raw FASTQ files to differential expression and Gene Ontology (GO) enrichment.

The workflow has two main components:

1. **`00_preprocessing_rnaseq_aedes.sh`**  
   Bash script that performs **quality control**, **read alignment** and **gene-level counting**.

2. **`Transcriptome_AAE.R`**  
   R script that performs **count matrix construction**, **DESeq2-based differential expression**, **QC visualizations**, **GO enrichment**, and **figure generation** for the manuscript.

Although the pipeline was developed for a specific experiment (two *A. aegypti* populations exposed to broflanilide and imidacloprid), it can be adapted to similar RNA-seq designs.

---

## 1. Repository contents

- `00_preprocessing_rnaseq_aedes.sh`  
  Shell script with all commands needed to:
  - run FastQC on raw FASTQ files,
  - prepare the reference genome and annotation,
  - build a HISAT2 index,
  - align paired-end reads,
  - convert SAM → sorted BAM with `samtools`,
  - generate gene-level count files with `featureCounts`.

- `Transcriptome_AAE.R`  
  R script that:
  - reads all `*_counts.txt` files,
  - builds and saves a unified count matrix and sample metadata,
  - runs DESeq2 using a **condition-based design**,
  - produces QC plots (PCA, sample–sample distance heatmap),
  - computes pairwise differential expression contrasts,
  - generates volcano plots,
  - performs GO over-representation analysis (BP, MF, CC),
  - produces dotplots and barplots for enriched GO terms,
  - generates supplementary tables (counts, DESeq2 results, DEG lists, GO enrichment).

---

## 2. Experimental design (example)

The current version of the R script assumes six **condition** codes, corresponding to population × treatment combinations:

- **New Orleans (lab strain)**
  - `CNO`: control  
  - `NOB`: broflanilide  
  - `NOI`: imidacloprid  

- **San Nicolás (field population)**
  - `CSN`: control  
  - `SNB`: broflanilide  
  - `SNI`: imidacloprid  

Each condition is expected to have **three biological replicates**, named like:

- `CNO1`, `CNO2`, `CNO3`
- `NOB1`, `NOB2`, `NOB3`
- `NOI1`, `NOI2`, `NOI3`
- `CSN1`, `CSN2`, `CSN3`
- `SNB1`, `SNB2`, `SNB3`
- `SNI1`, `SNI2`, `SNI3`

These sample names are used consistently throughout the pipeline and are mapped to more reader-friendly labels (e.g. `NOC1`, `SNC1`) in the final figures.

---

## 3. Software requirements

### System tools (used in `00_preprocessing_rnaseq_aedes.sh`)

- Linux with `apt` (or equivalent)
- `fastqc`
- `hisat2`
- `samtools`
- `subread` (for `featureCounts`)

The script includes example `apt` commands:

```bash
sudo apt update
sudo apt upgrade -y
sudo apt install -y fastqc samtools hisat2 subread

---

### R environment (used in `Transcriptome_AAE.R`)

- R (≥ 4.0 recommended)
- Packages (CRAN + Bioconductor), for example:

```r
install.packages(c(
  "tidyverse",
  "patchwork",
  "matrixStats",
  "writexl",
  "data.table"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2",
  "pheatmap",
  "EnhancedVolcano",
  "clusterProfiler"
))
```

The script itself loads the required packages at the top; make sure all of them install cleanly before running the full workflow.

---

## 4. Input data and directory structure

A typical project layout might look like:

```text
project_root/
├── references/
│   ├── GCF_002204515.2_AaegL5.0_genomic.fna(.gz)
│   ├── GCF_002204515.2_AaegL5.0_genomic.gtf(.gz)
│   └── AaegL5_index.*.ht2       # created by hisat2-build
├── batch_1/                      # raw FASTQ files (one or more batches)
│   ├── CNO1_R1.fastq.gz
│   ├── CNO1_R2.fastq.gz
│   ├── ...
│   └── SNI3_R2.fastq.gz
├── counts/                       # optional folder to store *_sorted.bam & *_counts.txt
├── 00_preprocessing_rnaseq_aedes.sh
└── Transcriptome_AAE.R
```

The R script expects all `*_counts.txt` files (featureCounts output) to be available in the working directory (or a specified folder) when it starts.

---

## 5. Usage

### Step 1 – Preprocessing: QC, alignment, and counting

1. **Edit paths** in `00_preprocessing_rnaseq_aedes.sh`:
   - `project_root` or batch folder paths,
   - `REF_DIR` (reference genome and GTF),
   - sample names in the loop section.

2. Make the script executable:

```bash
chmod +x 00_preprocessing_rnaseq_aedes.sh
```

3. Run it (or copy-paste the relevant sections into your terminal):

```bash
./00_preprocessing_rnaseq_aedes.sh
```

This will:

- run FastQC on all `.fastq.gz` files in each batch,
- build the HISAT2 index (if not already present),
- align each pair of FASTQ files to AaegL5,
- create sorted BAM files (`*_sorted.bam`),
- generate per-sample count files (`*_counts.txt`).

> **Important:** Make sure that the final `*_counts.txt` filenames follow the expected naming convention (e.g. `CNO1_counts.txt`, `NOB2_counts.txt`, etc.).

---

### Step 2 – R-based transcriptome analysis

1. Open `Transcriptome_AAE.R` in R or RStudio.

2. At the top of the script, adjust the working directory:

```r
setwd("/path/to/project_root")  # EDIT THIS
```

3. Confirm that the `*_counts.txt` files are present in the working directory (or in the directory used in the script).

4. Source or run the script:

```r
source("Transcriptome_AAE.R")
```

The script will:

- construct a unified count matrix (`matriz_counts_todos.txt`),
- create a `sample_info` table with population, treatment and condition,
- perform DESeq2 analysis with `design = ~ condition`,
- generate QC plots (PCA, sample–sample distance heatmap),
- produce:
  - multi-panel volcano plots for all pairwise contrasts,
  - a heatmap of the top 10 most variable genes (with gene ID + symbol labels),
  - GO enrichment dotplots (BP, MF, CC),
  - GO BP barplots per population and insecticide,
- export supplementary Excel files with:
  - raw and filtered counts + metadata,
  - full DESeq2 results per contrast,
  - DEG lists and GO enrichment tables.

Outputs are saved to folders such as:

- `counts_matriz/`
- `qc/`
- `heatmaps/` or `figures_final/`
- `degs/`
- `results/GO/`
- `supplementary/`

---

## 6. Reproducing manuscript figures

The R pipeline is organized so that each major figure corresponds to a specific code block, for example:

- **Figure 1** – sample-to-sample distance heatmap.  
- **Figure 2** – multi-panel volcano plots for all DESeq2 contrasts.  
- **Figure 3** – heatmap of the top 10 most variable genes, with gene ID + symbol labels.  
- **Figures 4–6** – GO BP/CC/MF dotplots across insecticide and population contrasts.  
- **Figures 7–8** – GO BP barplots stratified by population and insecticide.

All figure files are written as high-resolution `.png` files suitable for journal submission.

---

## 7. Adapting the pipeline to other experiments

To reuse this workflow for other RNA-seq experiments:

1. Adjust sample naming and metadata logic in `Transcriptome_AAE.R`:
   - edit how `condition`, `population`, and `treatment` are derived from sample names,
   - update factor levels and contrast definitions.

2. Update reference genome and annotation paths in:
   - `00_preprocessing_rnaseq_aedes.sh` (HISAT2 index, GTF),
   - `Transcriptome_AAE.R` if needed (e.g. for GO annotations).

3. Confirm that all required packages install correctly and that the directory structure matches the new project.

---

## 8. Contact / authorship

Please update this section with the appropriate project information:
- **Repository maintainer:** <Mariana Jiménez, https://github.com/mabisjmz>

If you use or adapt this pipeline in your own work, please acknowledge the original authors and cite the associated manuscript once published.


