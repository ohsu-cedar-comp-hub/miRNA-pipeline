rule summarize_preprocess:
    input:
        expand("samples/raw/{sample}.fastq", sample = SAMPLES),
        #expand("samples/ampumi/{sample}_ampumi.fastq", sample = SAMPLES),
        expand("samples/trimmed/{sample}_trimmed.fastq", sample = SAMPLES),
        expand("samples/trimmed/{sample}_trimmed_filtered.fastq", sample = SAMPLES)
    output:
        "results/tables/preprocessing_read_summary.txt"
    shell:
        """
        ls -1 samples/raw/*fastq | while read line; do samp=$( basename $line | sed 's/.fastq//g' ); reads=$( grep "^@" $line | wc -l) ; echo $samp $reads | tr ' ' '\t'; done > results/summary_stats/initial_reads.txt
        #ls -1 samples/ampumi/*ampumi.fastq | while read line; do samp=$( basename $line | sed 's/_ampumi.fastq//g' ); reads=$( grep "^@" $line | wc -l) ; echo $samp $reads | tr ' ' '\t'; done > results/summary_stats/UMI_collapsed_reads.txt
        ls -1 samples/trimmed/*trimmed.fastq | while read line; do samp=$( basename $line | sed 's/_trimmed.fastq//g' ); reads=$( grep "^@" $line | wc -l) ; echo $samp $reads | tr ' ' '\t'; done > results/summary_stats/trimmed_reads.txt
        ls -1 samples/trimmed/*trimmed_filtered.fastq | while read line; do samp=$( basename $line | sed 's/_trimmed_filtered.fastq//g' ); reads=$( grep "^@" $line | wc -l) ; echo $samp $reads | tr ' ' '\t'; done > results/summary_stats/filtered_by_length_reads.txt
        
        #paste results/summary_stats/{{initial,UMI_collapsed,trimmed,filtered_by_length}}_reads.txt | awk '{{print $1,$2,$4,$6,$8}}' | tr ' ' '\t' > {output}
        paste results/summary_stats/{{initial,trimmed,filtered_by_length}}_reads.txt | awk '{{print $1,$2,$4,$6}}' | tr ' ' '\t' > {output}
 
        echo -e "sample\tinitial_reads\ttrimmed_reads\tfiltered_reads" | cat - {output} > /tmp/out && mv /tmp/out {output}
        """

#rule summarize_alignment:
#    input:
#        expand("samples/star_miRNA/{sample}_bam/Aligned.out.sam", sample = SAMPLES),
#        expand("samples/star_miRNA/{sample}_bam/Aligned.filtered.sam", sample = SAMPLES),
#        expand("samples/bowtie/{sample}.bam", sample = SAMPLES),
#        expand("samples/bowtie2/{sample}.bam", sample = SAMPLES)
#    output:
#        "results/tables/alignment_read_summary.txt"
#    conda:
#        "../envs/bowtie2.yaml"
#    shell:
#        """
#        ls -1 samples/star_miRNA/*/Aligned.out.sam | while read line; do samp=$( echo $line | cut -d / -f 3 | sed 's/_bam//g' ); numreads=$( samtools view -q255 $line | wc -l ); echo $samp $numreads | tr ' ' '\t'; done > results/summary_stats/star_miRNA_unique_reads.txt
#        ls -1 samples/star_miRNA/*/Aligned.filtered.sam | while read line; do samp=$( echo $line | cut -d / -f 3 | sed 's/_bam//g' ); numreads=$( samtools view -q255 $line | wc -l ); echo $samp $numreads | tr ' ' '\t'; done > results/summary_stats/star_miRNA_filtered_unique_reads.txt
#        ls -1 samples/bowtie/*bam | while read line; do samp=$( basename $line | sed 's/.bam//g' ); numreads=$( samtools view -q255 $line | wc -l ); echo $samp $numreads | tr ' ' '\t'; done > results/summary_stats/bowtie_miRNA_unique_reads.txt
#        ls -1 samples/bowtie2/*bam | while read line; do samp=$( basename $line | sed 's/.bam//g' ); numreads=$( samtools view -q255 $line | wc -l ); echo $samp $numreads | tr ' ' '\t'; done > results/summary_stats/bowtie2_miRNA_unique_reads.txt
#        
#        paste results/summary_stats/{{star_miRNA,star_miRNA_filtered,bowtie_miRNA,bowtie2_miRNA}}_unique_reads.txt | awk '{{print $1,$2,$4,$6,$8}}' | tr ' ' '\t' > {output}
#
#        echo -e "sample\tstar_unique_reads\tstar_filtered_reads\tbowtie_reads\tbowtie2_reads" | cat - {output} > /tmp/out && mv /tmp/out {output}
#        """

