# RNA-Seq Transcriptome Analysis of Tomato Response to Bacterial Wilt Infection

## Project Overview

This project presents a comprehensive RNA-Seq–based transcriptomic analysis to investigate gene expression changes in tomato following infection with *Ralstonia solanacearum*. The analysis is based on and inspired by the study:

**Transcriptome analysis reveals differential transcription in tomato (Solanum lycopersicum) following inoculation with Ralstonia solanacearum**

The goal of this project is to reproduce and explore transcriptomic patterns associated with host-pathogen interaction, focusing on identifying differentially expressed genes (DEGs) and key defense-related pathways.

---

## Objectives

* To analyze RNA-Seq data from infected vs control tomato samples
* To identify differentially expressed genes associated with bacterial wilt response
* To explore transcriptional reprogramming during infection
* To perform functional enrichment analysis of DEGs

---

## Biological Context

* Host organism: Solanum lycopersicum
* Pathogen: Ralstonia solanacearum
* Disease: Bacterial wilt
* Study focus: Host transcriptional response to infection

---

## Dataset

* Source: Public RNA-Seq dataset associated with the published study
* Repository: NCBI Sequence Read Archive (SRA)
* Accession ID: PRJNA787007
* Experimental design:

  * Control (uninfected) samples
  * Infected samples post-inoculation

---

## Workflow

### 1. Quality Control

* Tool: FastQC
* Purpose: Evaluate sequencing quality metrics

### 2. Preprocessing

* Adapter trimming and quality filtering

### 3. Alignment

* Tool: BWA MEM
* Reference genome: Tomato reference genome
* Output: Sorted BAM files


### 4. Differential Expression Analysis

* Tool: DESeq2 (R)
* Thresholds:

  * Adjusted p-value (FDR) < 0.05
  * |log2 Fold Change| ≥ 1

### 4. Quantification
  * Tool: featureCounts
  * Output: Gene expression count matrix
---

## Key Results

* Identification of significantly upregulated and downregulated genes following infection
* Distinct separation of control and infected samples in PCA analysis
* Enrichment of defense-related biological processes, including:

  * Response to biotic stimulus
  * Signal transduction pathways
  * Stress-responsive gene regulation
* Detection of candidate genes potentially involved in resistance mechanisms

---


## Tools and Software

* FastQC
* HISAT2 / STAR
* SAMtools
* featureCounts / HTSeq
* DESeq2 (R)
* Python (for preprocessing and visualization)

---

## Outputs

* Quality control reports
* Alignment files (BAM)
* Gene count matrix
* Differential expression results (CSV format)
* Visualization:

  * PCA plots
  * Heatmaps
  * Volcano plots

---

## Biological Insights

The transcriptomic analysis reveals extensive transcriptional reprogramming in tomato following *R. solanacearum* infection. Several genes associated with plant defense, signaling pathways, and stress response are significantly regulated, providing insights into molecular mechanisms underlying bacterial wilt resistance.

---

## Reproducibility

All scripts, parameters, and computational environments are documented to ensure reproducibility of the analysis.

---

## Disclaimer

This project is an independent academic reproduction and analysis based on the published study by Chen et al. (2022). It is intended for educational and research demonstration purposes.

---

## Author

Chandrika Ghoshal
SRF, Division of Vegetable Science
Indian Agricultural Research institute
Pusa, New Delhi, India

---

## Reference

Chen, N., Shao, Q., Lu, Q., Li, X., & Gao, Y. (2022).
**Transcriptome analysis reveals differential transcription in tomato (Solanum lycopersicum) following inoculation with Ralstonia solanacearum**
Published in Scientific Reports
