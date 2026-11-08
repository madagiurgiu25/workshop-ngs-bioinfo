#############################################
# Run with:
# snakemake --cores 4 --config R1=data/raw_fastq/HG00171_female_ERR034564_1.cut.fastq.gz R2=data/raw_fastq/HG00171_female_ERR034564_2.cut.fastq.gz ref=data/index/GRCh38.primary_assembly.genome.cut.fa outdir=results/HG00171_female_ERR034564 sample=HG00171_female_ERR034564
# snakemake --cores 4 --config R1=data/raw_fastq/HG00152_male_SRR769545_1.cut.fastq.gz R2=data/raw_fastq/HG00152_male_SRR769545_2.cut.fastq.gz ref=data/index/GRCh38.primary_assembly.genome.cut.fa outdir=results/HG00152_male_SRR769545 sample=HG00152_male_SRR769545 
#############################################

configfile: "config.yaml"

R1 = config["R1"]
R2 = config["R2"]
SAMPLE = config["sample"]
REF = config["ref"]
OUTDIR = config["outdir"]
THREADS = config["threads"]
MINLEN = config["minlen"]
Q = config["q"] 
DP = config["DP"]

# ========== CREATE OUTPUT DIRECTORY ==========
os.makedirs(OUTDIR, exist_ok=True)

# ========== TARGET RULE ===================

rule all:
    input:
        f"{OUTDIR}/fastqc_reports/{SAMPLE}_1.cut_fastqc.html",
        f"{OUTDIR}/fastqc_reports/{SAMPLE}_2.cut_fastqc.html",
        f"{OUTDIR}/multiqc_report/multiqc_report.html",
        f"{OUTDIR}/fastqc_reports_trimmed/{SAMPLE}_1.cut_val_1_fastqc.html",
        f"{OUTDIR}/fastqc_reports_trimmed/{SAMPLE}_2.cut_val_2_fastqc.html",
        f"{OUTDIR}/multiqc_report_trimmed/multiqc_report.html",
        # f"{OUTDIR}/mapping/{SAMPLE}.sorted.bam.bai",
        # f"{OUTDIR}/dedup/{SAMPLE}.dedup.bam.bai",
        # f"{OUTDIR}/ppaired/{SAMPLE}.proper.sorted.bam.bai",
        f"{OUTDIR}/fastqc_reports_clean/{SAMPLE}.proper_R1_fastqc.html",
        f"{OUTDIR}/fastqc_reports_clean/{SAMPLE}.proper_R2_fastqc.html",
        f"{OUTDIR}/multiqc_report_clean/multiqc_report.html",
        f"{OUTDIR}/variants/{SAMPLE}_filtered.vcf.gz"


# ------------------------------------------------------------
# 1. Run FastQC on raw FASTQs
# ------------------------------------------------------------
rule fastqc_reports:
    input:
        R1, R2
    output:
        f"{OUTDIR}/fastqc_reports/{SAMPLE}_1.cut_fastqc.html",
        f"{OUTDIR}/fastqc_reports/{SAMPLE}_2.cut_fastqc.html"
    threads: THREADS
    shell:
        """
        mkdir -p {OUTDIR}/fastqc_reports
        fastqc -t {threads} -o {OUTDIR}/fastqc_reports {input}
        """

# ------------------------------------------------------------
# 2. Run MultiQC on raw FastQC reports
# ------------------------------------------------------------
rule multiqc_report:
    input:
        f"{OUTDIR}/fastqc_reports/{SAMPLE}_1.cut_fastqc.html",
        f"{OUTDIR}/fastqc_reports/{SAMPLE}_2.cut_fastqc.html",
    output:
        f"{OUTDIR}/multiqc_report/multiqc_report.html"
    shell:
        "mkdir -p {OUTDIR}/multiqc_report && multiqc {OUTDIR}/fastqc_reports -o {OUTDIR}/multiqc_report"

# ------------------------------------------------------------
# 3. Trim with Trim Galore
# ------------------------------------------------------------
rule trim_galore:
    input:
        R1, R2
    output:
        R1t=f"{OUTDIR}/trimmed/{SAMPLE}_1.cut_val_1.fq.gz",
        R2t=f"{OUTDIR}/trimmed/{SAMPLE}_2.cut_val_2.fq.gz"
    params:
        outdir=f"{OUTDIR}/trimmed",
        q=Q,
        minlen=MINLEN
    shell:
        """
        mkdir -p {params.outdir}
        trim_galore -q {params.q} --phred33 --illumina --paired -j 1 --length {params.minlen} --fastqc -o {params.outdir} {input}
        """

# ------------------------------------------------------------
# 4. QC after trimming
# ------------------------------------------------------------
rule fastqc_reports_trimmed:
    input:
        R1t=f"{OUTDIR}/trimmed/{SAMPLE}_1.cut_val_1.fq.gz",
        R2t=f"{OUTDIR}/trimmed/{SAMPLE}_2.cut_val_2.fq.gz"
    output:
        f"{OUTDIR}/fastqc_reports_trimmed/{SAMPLE}_1.cut_val_1_fastqc.html",
        f"{OUTDIR}/fastqc_reports_trimmed/{SAMPLE}_2.cut_val_2_fastqc.html"
    threads: THREADS
    shell:
        """
        mkdir -p {OUTDIR}/fastqc_reports_trimmed
        fastqc -t {threads} -o {OUTDIR}/fastqc_reports_trimmed {input}
        """

rule multiqc_report_trimmed:
    input:
        f"{OUTDIR}/fastqc_reports_trimmed/{SAMPLE}_1.cut_val_1_fastqc.html",
        f"{OUTDIR}/fastqc_reports_trimmed/{SAMPLE}_2.cut_val_2_fastqc.html",
    output:
        f"{OUTDIR}/multiqc_report_trimmed/multiqc_report.html"
    shell:
        "mkdir -p {OUTDIR}/multiqc_report_trimmed && multiqc {OUTDIR}/fastqc_reports_trimmed -o {OUTDIR}/multiqc_report_trimmed"

# ------------------------------------------------------------
# 5. Index reference (BWA)
# ------------------------------------------------------------
rule bwa_index:
    input:
        REF
    output:
        expand("{ref}.{ext}", ref=REF, ext=["amb", "ann", "bwt", "pac", "sa"])
    shell:
        "bwa index {input}"

# ------------------------------------------------------------
# 6. Align reads using BWA and sort
# ------------------------------------------------------------
rule bwa_mem:
    input:
        ref=REF,
        indexbwa=f"{REF}.sa",
        R1=f"{OUTDIR}/trimmed/{SAMPLE}_1.cut_val_1.fq.gz",
        R2=f"{OUTDIR}/trimmed/{SAMPLE}_2.cut_val_2.fq.gz"
    output:
        bam=f"{OUTDIR}/mapping/{SAMPLE}.sorted.bam",
        bai=f"{OUTDIR}/mapping/{SAMPLE}.sorted.bam.bai",
    threads: THREADS
    shell:
        """
        mkdir -p {OUTDIR}/mapping
        bwa mem -t {threads} -R \"@RG\\tID:{SAMPLE}\\tSM:{SAMPLE}\\tPL:ILLUMINA\\tLB:lib1\\tPU:unit1\" {input.ref} {input.R1} {input.R2} |
        samtools view -bS - | samtools sort -o {output.bam}
        samtools index {output.bam}
        """

# ------------------------------------------------------------
# 7. Deduplicate with Picard
# ------------------------------------------------------------
rule mark_duplicates:
    input:
        bam=f"{OUTDIR}/mapping/{SAMPLE}.sorted.bam"
    output:
        dedup=f"{OUTDIR}/dedup/{SAMPLE}.dedup.bam",
        dedupbai=f"{OUTDIR}/dedup/{SAMPLE}.dedup.bam.bai",
        metrics=f"{OUTDIR}/dedup/{SAMPLE}.dedup.metrics.txt"
    shell:
        """
        mkdir -p {OUTDIR}/dedup
        picard MarkDuplicates I={input.bam} O={output.dedup} M={output.metrics} REMOVE_DUPLICATES=true
        samtools index {output.dedup}
        """

# ------------------------------------------------------------
# 8. QC after deduplication
# ------------------------------------------------------------
#
#   HOMEWORK to add the FastQC and MultiQC analysis
#
#

# ------------------------------------------------------------
# 9. Keep only properly paired reads
# ------------------------------------------------------------
rule filter_properly_paired:
    input:
        dedup=f"{OUTDIR}/dedup/{SAMPLE}.dedup.bam"
    output:
        bam=f"{OUTDIR}/ppaired/{SAMPLE}.proper.sorted.bam",
        bai=f"{OUTDIR}/ppaired/{SAMPLE}.proper.sorted.bam.bai"
    shell:
        """
        mkdir -p {OUTDIR}/ppaired
        samtools view -b -f 2 -F 4 -F 8 {input.dedup} > {OUTDIR}/ppaired/{SAMPLE}.proper.bam
        samtools sort -o {output.bam} {OUTDIR}/ppaired/{SAMPLE}.proper.bam
        samtools index {output.bam}
        """

# ------------------------------------------------------------
# 10. QC after properly paired filtering
# ------------------------------------------------------------
rule convert_bam_to_fastq:
    input:
        bam=f"{OUTDIR}/ppaired/{SAMPLE}.proper.sorted.bam",
    output:
        R1p=f"{OUTDIR}/ppaired/{SAMPLE}.proper_R1.fastq",
        R2p=f"{OUTDIR}/ppaired/{SAMPLE}.proper_R2.fastq"
    shell:
        """
        samtools fastq {input.bam} -1 {output.R1p} -2 {output.R2p}
        """

rule fastqc_reports_clean:
    input:
        R1t=f"{OUTDIR}/ppaired/{SAMPLE}.proper_R1.fastq",
        R2t=f"{OUTDIR}/ppaired/{SAMPLE}.proper_R2.fastq"
    output:
        f"{OUTDIR}/fastqc_reports_clean/{SAMPLE}.proper_R1_fastqc.html",
        f"{OUTDIR}/fastqc_reports_clean/{SAMPLE}.proper_R2_fastqc.html"
    threads: THREADS
    shell:
        """
        mkdir -p {OUTDIR}/fastqc_reports_clean
        fastqc -t {threads} -o {OUTDIR}/fastqc_reports_clean {input}
        """

rule multiqc_report_clean:
    input:
        f"{OUTDIR}/fastqc_reports_clean/{SAMPLE}.proper_R1_fastqc.html",
        f"{OUTDIR}/fastqc_reports_clean/{SAMPLE}.proper_R2_fastqc.html"
    output:
        f"{OUTDIR}/multiqc_report_clean/multiqc_report.html"
    shell:
        "mkdir -p {OUTDIR}/multiqc_report_clean && multiqc {OUTDIR}/fastqc_reports_clean -o {OUTDIR}/multiqc_report_clean"

# ------------------------------------------------------------
# 11. Variant Calling (GATK HaplotypeCaller)
# ------------------------------------------------------------
rule haplotype_caller:
    input:
        bam=f"{OUTDIR}/ppaired/{SAMPLE}.proper.sorted.bam",
        ref=REF
    output:
        vcf=f"{OUTDIR}/variants/{SAMPLE}.vcf.gz"
    shell:
        """
        mkdir -p {OUTDIR}/variants
        gatk HaplotypeCaller -R {input.ref} -I {input.bam} -O {output.vcf} -A AlleleFraction
        """

# ------------------------------------------------------------
# 12. Variant Filtering
# ------------------------------------------------------------
rule variant_filtration:
    input:
        vcf=f"{OUTDIR}/variants/{SAMPLE}.vcf.gz",
        ref=REF
    output:
        flagged=f"{OUTDIR}/variants/{SAMPLE}_flagged.vcf.gz"
    params:
        dp=DP
    shell:
        """
        gatk VariantFiltration \
          -R {input.ref} \
          -V {input.vcf} \
          -O {output.flagged} \
          --filter-name \"LowDP\" --filter-expression \"DP < {params.dp}\" \
          --filter-name \"StrandBiasFS\" --filter-expression \"FS > 60.0\" \
          --filter-name \"StrandBiasSOR\" --filter-expression \"SOR > 3.0" \
          --filter-name \"TailDistBias\" --filter-expression \"ReadPosRankSum < -8.0\" \
          --filter-name \"NonGermline\" --filter-expression \"(AD[1]/(AD[0]+AD[1]) < 0.4)\" \
          --filter-name \"LowQual\" --filter-expression \"QUAL < 10\"
        """

# ------------------------------------------------------------
# 13. Select PASS Variants
# ------------------------------------------------------------
rule select_variants:
    input:
        f"{OUTDIR}/variants/{SAMPLE}_flagged.vcf.gz"
    output:
        f"{OUTDIR}/variants/{SAMPLE}_filtered.vcf.gz"
    shell:
        "gatk SelectVariants -V {input} --exclude-filtered -O {output}"
