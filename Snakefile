"""
Author: Arlet Culhaci
Affiliation: UvA
Aim: A simple Snakemake workflow to process RAD-seq data and output a VCF (Variant Calling Format) matrix. This matrix serves as the input for QTL analysis.
Date of creation: Sun Feb 24, 2019
Dry run: snakemake -np
Run: snakemake --use-conda --cores N
Latest modification:
  - todo
"""
from snakemake.utils import min_version
min_version("5.4.3")
import pandas as pd # to read barcode to sample correspondence

####################
# Configuration file
####################
configfile: "configure.yaml"

TEMP_DIR =   config["temp_dir"]    # will be removed upon pipeline completion
RESULT_DIR = config["result_dir"]  # will be kept upon pipeline completion

############################################
# Get the sample names from the barcode file
############################################
barcodes = pd.read_csv("barcodes.tsv",delimiter="\t",header=None,names=["barcode","sample"]) # imports the barcodes as a Pandas dataframe
SAMPLES = barcodes["sample"].to_list()


################
# Desired output
################
rule all:
    input:
        expand(RESULT_DIR + "{sample}.bam", sample=SAMPLES)
    message:"All done! Removing {TEMP_DIR} directory"
    shell:
        "rm -r {TEMP_DIR}"

#######
# Rules
#######
rule demultiplex_filter_pe:
    input:
        forward = config["fastq"]["forward"],
        reverse = config["fastq"]["reverse"]
    output:
        expand(TEMP_DIR + "fastq/{sample}.{mate}.fq.gz",sample=SAMPLES,mate=["1","2"])
    conda:
        "envs/stacks.yaml"
    params:
        barcodes = config["barcodes"],
        output_dir = TEMP_DIR + "fastq/",
        restriction_enzyme = config["restriction_enzyme"]
    log:
        RESULT_DIR + "logs/demultiplex.log"
    message: "demultiplexing {input.forward} and {input.reverse} reads"
    shell:
        "process_radtags -1 {input.forward} -2 {input.reverse} "
        "-b {params.barcodes} "
        "-o {params.output_dir} "
        "-e {params.restriction_enzyme} "
        "--rescue " # rescue barcodes and RAD-Tags
        "--clean "  # clean data, remove read with an uncalled base
        "--quality " # discard reads with low quality scores
        "2>{log}"

rule calc_stat:
    input:
        directory("results/")
    output:
        "calc/initial_reads.txt"
    params:
        num_ind = config["num_ind"]
    shell:
        "cat {input}process_radtags* | tail -n +12 | awk '{{print $2, $3, $6}}' |"
        " head -{params.num_ind} | tr ' ' ';' > {output}"

rule vis_reads:
    input: "calc/initial_reads.txt"
    output: "calc/initial_reads.pdf"
    conda:
        "envs/ggplot2.yaml"
    shell:
        "Rscript scripts/initial_reads_IBEDs_pipeline.R {input} {output}"

#glob_wildcards('results/{sample}.fq.gz')
rule unzip_fastq:
    input:
        forward = TEMP_DIR + "fastq/{sample}.1.fq.gz",
        reverse = TEMP_DIR + "fastq/{sample}.2.fq.gz"
    output:
        forward = TEMP_DIR + "fastq/{sample}.1.fq",
        reverse = TEMP_DIR + "fastq/{sample}.2.fq"
    message:"unzipping {wildcards.sample} fastq files"
    shell:
        "gunzip {input}"

rule bowtie2_index:
    input:
        config["genome"]
    output:
        [TEMP_DIR + "genome." + str(i) + ".bt2" for i in range(1,5)]
    message:"indexing genome reference: {input}"
    conda:
        "envs/bowtie2.yaml"
    params:
        index_name = TEMP_DIR + "genome"
    threads: 10
    shell:
        "bowtie2-build --threads {threads} --quiet {input} {params}"

rule bowtie2_align_pe:
    input:
        forward = TEMP_DIR + "fastq/{sample}.1.fq",
        reverse = TEMP_DIR + "fastq/{sample}.2.fq",
        index = [TEMP_DIR + "genome." + str(i) + ".bt2" for i in range(1,5)]
    output:
        RESULT_DIR + "{sample}.bam"
    conda:
        "envs/bowtie2.yaml"
    threads: 10
    params:
        index_name = TEMP_DIR + "genome"
    log:
        RESULT_DIR + "logs/{sample}.bowtie2.log"
    message:"aligning {wildcards.sample} reads using Bowtie2"
    shell:
        "bowtie2 -p {threads} "
        "-x  {params.index_name} "
        "-1 {input.forward} -2 {input.reverse} "
        "-S - | samtools view -Sb -F 4 -o {output} " # takes standard input and converts on the fly to .bam format (removes unmapped)
