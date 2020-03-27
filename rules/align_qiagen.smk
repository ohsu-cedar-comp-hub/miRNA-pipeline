rule qiagen_bowtie:
    input:
        "samples/trimmed/{sample}_trimmed_filtered.fastq"
    output:
        sam = "samples/qiagen_mirna/{sample}.sam",
        log = "samples/qiagen_mirna/{sample}.log",
        unaligned = "samples/qiagen_mirna/{sample}.unaligned.fastq"
    params:
       index = "/home/exacloud/lustre1/CEDAR/anurpa/mirna/qiagen/qiaseq_ref/bowtie_index/Human/mirna"
    conda:
        "../envs/bowtie.yaml" 
    shell:
        """bowtie {params.index} -n 0 --un {output.unaligned} -q {input} 1>{output.sam} 2>{output.log}"""

rule secondary_mapping:
    input:
        "samples/qiagen_mirna/{sample}.unaligned.fastq"
    output:
        sam = "samples/qiagen_mirna_secondPass/{sample}.sam",
        log = "samples/qiagen_mirna_secondPass/{sample}.log"
    params:
       index = "/home/exacloud/lustre1/CEDAR/anurpa/mirna/qiagen/qiaseq_ref/bowtie_index/Human/mirna"
    conda:
        "../envs/bowtie.yaml"
    shell:
        """bowtie {params.index} -q -l 15 -5 1 -3 2 -n 2 {input} 1>{output.sam} 2>{output.log}"""

rule merge_sams:
    input:
        firstPass = "samples/qiagen_mirna/{sample}.sam",
        secondPass = "samples/qiagen_mirna_secondPass/{sample}.sam"
    output:
        "samples/qiagen_mirna/{sample}.final.merged.sam"
    shell:
        """cat {input.firstPass} {input.secondPass} > {output}"""

rule count_qiagen:
    input:
        "samples/qiagen_mirna/{sample}.final.merged.sam"
    output:
        "samples/count_qiagen/{sample}.counts.txt"
    conda:
        "../envs/bowtie.yaml"
    shell:
        "cat {input} | cut -f 3 | sort | uniq -c > {output}"

rule combine_qiagen_counts:
    input:
        expand("samples/count_qiagen/{sample}.counts.txt", sample = SAMPLES)
    output:
        "results/tables/qiagen_counts.txt"
    script:
        "../scripts/combine_bowtie_counts.R"

