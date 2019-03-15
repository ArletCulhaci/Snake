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
    input: directory("results/")

#######
# Rules
#######

rule demultiplex_filter:
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
