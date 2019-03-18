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

SAMPLES = []

f = open( config["barcodes"], 'rU' ) #open the file in read universal mode
for line in f:
    cells = line.split( "\t" )
    SAMPLES.append(cells[1].rstrip("\n")) #since we want the first, second and third column
f.close()

################
# Desired output
################
rule all:
    input: expand("results/{sample}.bam", sample=SAMPLES)
    #input: directory("results/{sample}.bam")
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

glob_wildcards('results/{sample}.fq.gz')

if config["paired"]:
    rule bowtie2_align_pe:
        input:
            glob_wildcards('results/{sample}.fq.gz').example
            #"results/{}"
        output:
            directory("results/{sample}.bam")
        #threads: CLUSTER["align"]["cpu"]
        #conda:
        #    "envs/stacks.yaml"
        params:
            reference = config["reference"]
        log:
            "results/"
        shell:
          "bowtie2  -p 20 -x  {params.reference} -1 results/{sample}.1.fq.gz -2 results/{samples}.2.fq.gz -S {output}"
          #samtools view -Sb ${5}/"STACKS_$DATE-$N/alignment"/${names}.sam > ${5}/"STACKS_$DATE-$N/alignment"/${names}.bam
          #samtools sort ${5}/"STACKS_$DATE-$N/alignment"/${names}.bam ${5}/"STACKS_$DATE-$N/alignment"/${names}
          #samtools index ${5}/"STACKS_$DATE-$N/alignment"/${names}.bam

else:
    rule bowtie2_align_se:
