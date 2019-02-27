"""
Author: A. Culhaci
Affiliation: UvA
Aim: A simple Snakemake workflow to process RAD-seq data.
Date: Sun Feb  24  2019
Run: snakemake   -s Snakefile
Latest modification:
  - todo
"""
configfile: "configure.yaml"

rule dempulitplex_filter:
    input:  "ata/sub_forward.fq.gz", "ata/sub_reverse.fq.gz"
    output: directory("/snake")
    #threads: CLUSTER["align"]["cpu"]
    params:
            #jobname = "{sample}",
            ## add read group for bwa mem mapping, change accordingly if you know PL:ILLUMINA, LB:library1 PI:200 etc...
            #rg = "@RG\tID:{sample}\tSM:{sample}"
    #message: "aligning bwa {input}: {threads} threads"
    #log:
        #bwa = "00log/{sample}.align",
        #markdup = "00log/{sample}.markdup"
    shell:
        "process_radtags -1 {input[0]} -2 {input[1]} -b {configure[barcodes]} -o {output} -e sbfI -r -c -q"
    #run:
        ## paired end reads
    #    if config["paired_end"]:
            ## short reads < 70bp
            ## Probably one of the most important is how many mismatches you will allow between a read and a potential mapping location for that location to be considered a match.
            ## The default is 4% of the read length, but you can set this to be either another proportion of the read length, or a fixed integer
    #        shell(
    #            r"""
    #            process_radtags -1 {input[0]} -2 {input[1]} -b {configure[barcodes]} -o {output} -e sbfI -r -c -q
    #            """)
        # single end  reads
    #    else:
            ## short reads < 70bp
    #        shell(
    #            r"""
    #            process_radtags -f {input[0]} -b {configure[barcodes]} -o {output} -e sbfI -r -c -q
    #            """)
rule all:
    input: directory("/snake")
