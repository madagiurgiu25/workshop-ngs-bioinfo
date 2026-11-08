# Snakemake Bioinformatics Workflow for NGS 

This repository contains a reproducible **bioinformatics workflow** implemented using **Snakemake** for NGS analysis explained in [tutorial_from_fastq_to_bam.md](tutorial_from_fastq_to_bam.md) and [tutorial_from_bam_to_vcf.md](tutorial_from_bam_to_vcf.md). This step by step tutorial does not require upfront knowledge about **Snakemake**.


## 1. Requirements

- [Snakemake](https://snakemake.readthedocs.io/)
- Python ≥ 3.8
- Conda or Mamba (optional, for environment management)

## 2. Setup

```
git clone https://github.com/madagiurgiu25/workshop-ngs-bioinfo.git

cd workshop-ngs-bioinfo
mkdir -p data/raw_fastq && mkdir -p data/index
```

### 2.1 Environment installation

Bioinformatics packages (see [tutorial_install.md](tutorial_install.md)).
After having the environment please do:

```bash
conda activate bioinfo2026
```

### 2.2 Download input data
Download the raw fastq from [zenodo](https://zenodo.org/records/17598455) and place them under `workshop-ngs-bioinfo/data/raw_fastq`.

### 2.3 Prepare reference genome

Check [Step 2](https://github.com/madagiurgiu25/workshop-ngs-bioinfo/blob/main/tutorials/tutorial_from_fastq_to_bam.md#2-prepare-reference-genome) in [tutorial_from_fastq_to_bam.md](tutorial_from_fastq_to_bam.md).

And move this trimmed reference genome under `workshop-ngs-bioinfo/data/index`.

### 2.4 Project structure

```
workshop-ngs-bioinfo/
├── Snakefile
├── config.yaml
├── data/
│   ├── raw_fastq/
│   │   ├── HG00152_male_SRR769545_1.cut.fastq.gz
│   │   ├── HG00152_male_SRR769545_2.cut.fastq.gz
│   │   ├── HG00171_female_ERR034564_1.cut.fastq.gz
│   │   └── HG00171_female_ERR034564_2.cut.fastq.gz
│   └── index/GRCh38.primary_assembly.genome.cut.fa
```

### 2.5 Configuration

The `config.yaml` file contains parameters such as sample name, reference genome, output directory etc., and can be configured from the commandline.

Example:

```yaml
# ========== CONFIG PARAMETERS ==============

# Fastq 1
R1: "data/raw_fastq/HG00171_female_ERR034564_1.cut.fastq.gz"

# Fastq 2
R2: "data/raw_fastq/HG00171_female_ERR034564_2.cut.fastq.gz"

# Sample name (without _1/_2 suffix)
sample: "HG00171_female_ERR034564"

# Reference genome FASTA file
ref: "data/index/GRCh38.primary_assembly.genome.cut.fa"

# Output directory for all results
outdir: "results/HG00171_female_ERR034564"

# Number of threads for multithreaded tools
threads: 4

# Minimum read length for trimming
minlen: 50

# Quality threshold for trimming
q: 20

# Minimal depth of coverage for variant calling
DP: 20

```

## 3. Running the Pipeline

Run Snakemake with 4 cores:

```bash
snakemake --cores 4 \
  --config \
  R1=data/raw_fastq/HG00171_female_ERR034564_1.cut.fastq.gz \
  R2=data/raw_fastq/HG00171_female_ERR034564_2.cut.fastq.gz \
  ref=data/index/GRCh38.primary_assembly.genome.cut.fa \
  outdir=results/HG00171_female_ERR034564
```

## 4. Workflow structure

![pipeline](../dag.png)

## 5. Homework

Exercise 1: Try to change the pipeline such that you include FastQC and MultiQC step after `mark_deduplicates` step. Hint: This is similar to `fastq_reports` and `multiqc_report` rules.

Exercise 2: Make the pipeline more configurable by adding as parameter the `SOR` score. This means the `yaml` file and `variant_filtration` rule need adjustment. Hint: This is similar to `DP` parameter.

