"""Computation Hub omic data processing pipeline"""

import datetime
import sys
import os
import pandas as pd
import json

timestamp = ('{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now()))

configfile:"config.yaml"
project_id = config["project_id"]

SAMPLES, = glob_wildcards("samples/raw/{sample}.fastq")

filt_type = ['out','filtered']

with open('cluster.json') as json_file:
    json_dict = json.load(json_file)

rule_dirs = list(json_dict.keys())

for rule in rule_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'logs',rule)):
        log_out = os.path.join(os.getcwd(), 'logs', rule)
        os.makedirs(log_out)
        print(log_out)

def message(mes):
    sys.stderr.write("|--- " + mes + "\n")

for sample in SAMPLES:
    message("Sample " + sample + " will be processed")

rule all:
    input:
        expand("samples/readLength/{sample}_unprocessed_readLength.txt", sample = SAMPLES),
        "results/multiqc/{project_id}_multiqc.html".format(project_id=config["project_id"]),
#        expand("results/tables/hsa_miRNA_count_from_hg38_{filt_type}.txt", filt_type = filt_type),
        expand("samples/readLength/{sample}_readLength.txt", sample = SAMPLES),
#        expand("samples/star_miRNA/{sample}_bam/{sample}.filtered.sorted.bw", sample = SAMPLES),
#        "results/tables/miR_counts_bowtie2.txt",
#        "results/tables/miR_counts_bowtie1.txt",
#        expand("samples/bowtie2/{sample}.bw", sample = SAMPLES),
        "results/tables/preprocessing_read_summary.txt",
#        "results/tables/alignment_read_summary.txt", 
        "results/tables/qiagen_counts.txt"

include: "rules/preprocess.smk"
#include: "rules/align_star.smk"
#include: "rules/align_bowtie.smk"
include: "rules/summary.smk"
#include: "rules/mirpro.smk"
include: "rules/align_qiagen.smk"
