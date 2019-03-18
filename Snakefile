"""
Author: A. Culhaci
Affiliation: UvA
Aim: A simple Snakemake workflow to process RAD-seq data.
Date: Sun Feb  24  2019
Run: snakemake   -s Snakefile
Latest modification:
  - todo
"""
from snakemake.utils import min_version
min_version("5.4.3")


configfile: "configure.yaml"

################
# Desired output
################
rule all:
    input: directory("calc/initial_reads.pdf")
#######
# Rules
#######
if config["paired"]:
    rule demultiplex_filter_pe:
        input:
                forward = "data/sub_forward.fq.gz",
                reverse = "data/sub_reverse.fq.gz"
        output:
            directory("results/")
        #threads: CLUSTER["align"]["cpu"]
        conda:
            "envs/stacks.yaml"
        params:
            barcodes = config["barcodes"]
        log:
            "results/"
        shell:
            "process_radtags -1 {input.forward} -2 {input.reverse} "
            "-b {params.barcodes} "
            "-o {output} "
            "-e sbfI "
            "--rescue " # rescue barcodes and RAD-Tags
            "--clean "  # clean data, remove read with an uncalled base
            "--quality" # discard reads with low quality scores
else:
    rule demultiplex_filter_se:
        input:
                forward = "data/sub_forward.fq.gz",
        output:
            directory("results/")
        #threads: CLUSTER["align"]["cpu"]
        conda:
            "envs/stacks.yaml"
        params:
            barcodes = config["barcodes"]
        log:
            "results/"
        shell:
            "process_radtags -f {input.forward} "
            "-b {params.barcodes} "
            "-o {output} "
            "-e sbfI "
            "--rescue " # rescue barcodes and RAD-Tags
            "--clean "  # clean data, remove read with an uncalled base
            "--quality" # discard reads with low quality scores

rule calc_stat:
    input:
        directory("results/")
    output: "calc/initial_reads.txt"
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
        "Rscript initial_reads_IBEDs_pipeline.R {input} {output}"
