contrast = get_contrast

rule deseq2_init:
    input:
        counts = "results/mirge/miR.Counts.csv"
    output:
        rds="results/diffexp/{project_id}_all.rds".format(project_id=project_id),
        normed_counts="results/tables/{project_id}_normed_counts.txt".format(project_id = project_id),
        rld_out = "results/diffexp/{project_id}_rlog_dds.rds".format(project_id = project_id),
    params:
        samples=config["omic_meta_data"],
        design=config["pca"]["labels"],
        row_names=config["sample_id"],
    conda:
        "../envs/deseq2.yaml",
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2-init.R"

rule deseq2_plots:
    input:
        rds_object = "results/diffexp/{project_id}_rlog_dds.rds".format(project_id = project_id),
    output:
        pca="results/diffexp/pca.pdf",
        sd_mean_plot="results/diffexp/sd_mean_plot.pdf",
        distance_plot = "results/diffexp/distance_plot.pdf",
        ggplot_pca_factor = "results/diffexp/ggplot_factor_pca.pdf",
    params:
        pca_labels=config["pca"]["labels"]
    conda:
        "../envs/deseq2_plots.yaml"
    script:
        "../scripts/plot-pca.R"


rule deseq2:
    input:
        rds="results/diffexp/{project_id}_all.rds".format(project_id=project_id),
        rld = "results/diffexp/{project_id}_rlog_dds.rds".format(project_id = project_id)
    output:
        table="results/diffexp/{contrast}.diffexp.tsv",
        ma_plot="results/diffexp/{contrast}.ma_plot.pdf",
        p_hist="results/diffexp/{contrast}.phist_plot.pdf",
        heatmap_plot = "results/diffexp/{contrast}.heatmap_plot.pdf",
        panel_ma = "results/diffexp/{contrast}.panel_ma.pdf",
        var_heat = "results/diffexp/{contrast}.variance_heatmap.pdf"
    params:
        contrast=get_contrast,
        condition = config["linear_model"]
    conda:
        "../envs/deseq2_plots.yaml"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2.R"

rule volcano:
    input:
        "results/diffexp/{contrast}.diffexp.tsv"
    output:
        "results/diffexp/{contrast}.diffexp.01.VolcanoPlot.pdf"
    params:
       contrast = get_contrast
    shell: """Rscript scripts/RNAseq_makeVolcano.R --degFile={input} --adjp=0.01 --FC=2"""

