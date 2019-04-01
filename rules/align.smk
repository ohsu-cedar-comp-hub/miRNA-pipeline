rule ampUMI:
    input: 
        "samples/raw/{sample}.fastq"
    output: 
        "samples/ampumi/{sample}_ampumi.fastq"
    params:
        ampUMI = config["ampUMI_tool"]
    conda:
        "../envs/ampumi.yaml"
    shell:
        """python {params.ampUMI} Process \
                   --fastq {input} --fastq_out {output} \
                   --umi_regex "AACTGTAGGCACCATCAATIIIIIIIIIIII.*" --write_UMI_counts > logs/ampUMI/{wildcards.sample}_ampumi.log 2>&1"""

rule count:
    input:
        expand("samples/ampumi/{sample}_ampumi.fastq", sample = SAMPLES)
    output:
        "results/mirge/miR.Counts.csv"
    params:
        bowtie = config["bowtie_bin"],
        refs = config["mirge_refs"],
    conda:
        "../envs/mirge.yaml"
    shell:
        """
        miRge2.0 annotate -s {input} -o results/mirge \
                          -d miRBase -pb {params.bowtie} \
                          -lib {params.refs} -sp human -gff -trf -cpu 4
        mv results/miRge.*/* results/
        rm -r results/miRge.*
        """

