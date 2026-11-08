

#  Bioinformatics Environment Setup Tutorial

Before starting, ensure you have on your laptop/work station:

* **At least 10 GB** of free disk space
* **At least 8 GB RAM**
* **A stable internet connection** for downloading reference and sample data

To download the code from this repository you need to install Git. If you are not familiar with Git check for example this [video tutorial](https://www.youtube.com/watch?v=QzvA7r-WndM).

```
# MacOS
brew install git

# Linux
sudo apt get install git

# Windows
# install from https://git-scm.com/install/windows
```

## **1. Prepare a Linux-like Terminal**

Most bioinformatics tools are designed for Linux.

* **If you already use macOS or Linux:**
  You’re good to go — just open your default terminal.

* **If you’re on Windows:**

  1. Download and install **Git Bash** from [https://gitforwindows.org/](https://gitforwindows.org/)
  2. Once installed, open **Git Bash**. It provides a Unix-like shell similar to Linux.
  

## **2. Install Conda (via Miniforge)**

If you already have Conda installed please ignore this step.
Conda is a package and environment manager that simplifies installation of scientific tools. 

### Steps:

1. Visit the **Miniforge**
   [https://conda-forge.org/download/](https://conda-forge.org/download/)

2. Download the installer for your operating system.

3. Run the installer from your terminal, for example:

   ```bash
   # Linux or macOS
   bash <name of the downloaded script>.sh

   # Windows - run the .exe file
   ```

4. Restart your terminal and test if Conda works:

   ```bash
   conda --version
   ```

You should see something like `conda 24.x.x`.


## **3. Create the Bioinformatics Environment**

Create a Conda environment containing all the bioinformatics tools necessary for the next-generation sequencing processing workshop.

### Run the following command:

```bash
conda create -n bioinfo2025 -c conda-forge -c bioconda \
  python==3.10 \
  fastqc==0.11.8 \
  trim-galore==0.6.10 \
  bwa==0.7.19 \
  picard==3.4.0 \
  samtools==1.22.1 \
  bedtools==2.31.1 \
  deeptools==3.5.6 \
  multiqc==1.31 \
  gatk4==4.6.2.0 \
  bcftools==1.22 \
  snakemake
```

Once installed, activate the environment:

```bash
conda activate bioinfo2025
```

Check that tools are accessible, e.g.:

```bash
fastqc --version
bwa
samtools --version
```

## **4. Download the Human Reference Genome**

We will use as human reference genome the GRCh38 primary assembly from **GENCODE Release 49**.

### Option A: Download via browser

Visit GENCODE website [https://www.gencodegenes.org/human/](https://www.gencodegenes.org/human/) to **Release 49** → Fasta files → **Genome sequence primary assembly (GRCh38)**:

![Gencode](./img/refgenome.png)

### Option B: Download via command line

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
```

## **5. Download Raw Data for Two Coriell Samples**

We’ll use two example datasets from the **Coriell Institute** [https://www.coriell.org/], `HG00152` and `HG00171`. The raw data is Illumina short-read, paired-end whole-exome sequencing (WES) and you can find this on <b>The International Genome Sample Resource (IGSR)</b> platform [https://www.internationalgenome.org/data-portal/sample]. 

For practicality, in the workshop we will work only with regions from a subset of genes: APC, BRCA1, BRCA2, PMS2, DMD, spanning chr5,chr7,chr13,chr17,chrX. The modified FASTQ files can be downloaded from [Raw_FastQ_OnTarget_Coriell](https://zenodo.org/records/17598455)

```
wget https://zenodo.org/records/17598455/files/HG00152_male_SRR769545_1.cut.fastq.gz
wget https://zenodo.org/records/17598455/files/HG00152_male_SRR769545_2.cut.fastq.gz
wget https://zenodo.org/records/17598455/files/HG00171_female_ERR034564_1.cut.fastq.gz
wget https://zenodo.org/records/17598455/files/HG00171_female_ERR034564_2.cut.fastq.gz
```

## You’re Ready!

Once all steps are completed, you’ll have a fully functional bioinformatics environment ready for:

* Quality control (`FastQC`, `MultiQC`)
* Adapter trimming (`Trim Galore`)
* Read alignment (`BWA`)
* Post-processing (`Picard`, `Samtools`)
* Coverage analysis (`deeptools`, `bedtools`)
* Variant calling (`GATK`, `bcftools`)



