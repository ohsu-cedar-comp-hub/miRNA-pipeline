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

with open('cluster.json') as json_file:
    json_dict = json.load(json_file)

rule_dirs = list(json_dict.keys())
rule_dirs.pop(rule_dirs.index('__default__'))

for rule in rule_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'logs',rule)):
        log_out = os.path.join(os.getcwd(), 'logs', rule)
        os.makedirs(log_out)
        print(log_out)

def message(mes):
    sys.stderr.write("|--- " + mes + "\n")

def get_deseq2_threads(wildcards=None):
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(config["omic_meta_data"]) < 100 or few_coeffs else 6

def get_contrast(wildcards):
    """Return each contrast provided in the configuration file"""
    return config["diffexp"]["contrasts"][wildcards.contrast]

for sample in SAMPLES:
    message("Sample " + sample + " will be processed")

rule all:
    input:
        expand("samples/ampumi/{sample}_ampumi.fastq",sample = SAMPLES),
        "results/mirge/miR.Counts.csv",
        "results/diffexp/pca.pdf",
        expand("results/diffexp/{contrast}.diffexp.01.VolcanoPlot.pdf", contrast = config["diffexp"]["contrasts"]),

include: "rules/align.smk"
include: "rules/deseq.smk"
