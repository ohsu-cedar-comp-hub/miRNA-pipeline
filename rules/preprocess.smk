rule get_initial_length_distr:
    input:
        "samples/raw/{sample}.fastq"
    output:
        "samples/readLength/{sample}_unprocessed_readLength.txt"
    shell:
        """cat {input} | awk '{{if(NR%4==2) print length($1)}}' | sort -n | uniq -c > {output}"""

#rule ampUMI:
#    input: 
#        "samples/raw/{sample}.fastq"
#    output: 
#        "samples/ampumi/{sample}_ampumi.fastq"
#    params:
#        ampUMI = config["ampUMI_tool"]
#    conda:
#        "../envs/ampumi.yaml"
#    shell:
#        """python {params.ampUMI} Process \
#                   --fastq {input} --fastq_out {output} \
#                   --umi_regex "AACTGTAGGCACCATCAATIIIIIIIIIIII.*" --write_UMI_counts"""

rule fastqc:
    input:
        "samples/raw/{sample}.fastq"
    output:
        "samples/fastqc/{sample}/{sample}_fastqc.zip",
    conda:
        "../envs/fastqc.yaml"
    shell:
        """fastqc --outdir  samples/fastqc/{wildcards.sample} --extract  -f fastq {input}"""

rule multiqc:
    input:
        expand("samples/fastqc/{sample}/{sample}_fastqc.zip", sample = SAMPLES),
    output:
        "results/multiqc/{project_id}_multiqc.html".format(project_id=config["project_id"])
    params:
        project_id=config["project_id"]
    conda:
        "../envs/multiqc.yaml"
    shell:
        """multiqc --outdir results/multiqc --filename {params.project_id}_multiqc samples/fastqc/"""

rule cutadapt:
    input:
        "samples/raw/{sample}.fastq"
    output:
        "samples/trimmed/{sample}_trimmed.fastq"
    conda:
        "../envs/cutadapt.yaml"
    shell:
        """cutadapt -a AACTGTAGGCACCATCAAT \
                    -o {output} {input}"""

rule filter_short_reads:
    input:
        "samples/trimmed/{sample}_trimmed.fastq"
    output:
        "samples/trimmed/{sample}_trimmed_filtered.fastq"
    shell:
        "scripts/filter_short_reads.sh -i {input} -o {output}"

rule get_length_distr:
    input:
        "samples/trimmed/{sample}_trimmed_filtered.fastq"
    output:
        "samples/readLength/{sample}_readLength.txt"
    shell:
        """cat {input} | awk '{{if(NR%4==2) print length($1)}}' | sort -n | uniq -c > {output}"""

