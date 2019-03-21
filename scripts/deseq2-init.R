library("DESeq2")

counts = snakemake@input[['counts']]

params = snakemake@params[['samples']]

output = snakemake@output[['rds']]

dds_design = snakemake@params[['design']]

row_names = snakemake@params[['row_names']]

out_table = snakemake@output[['normed_counts']]

rld_out = snakemake@output[['rld_out']]

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

cts <- read.csv(counts, row.names=1)
# We remove the first row "miRNAtotal" as it is not necessary for deseq
cts = cts[-1,]
coldata <- read.delim(params, header=TRUE, row.names=1)

# Extract just the sampleID from the column name (default is {sampleID}_ampumi.fastq) and rename the columns of counts
samps <- c()
for (i in 1:(ncol(cts))) {
  name <- colnames(cts)[i]
  samp <- strsplit(name, "_")[[1]][1]
  samps[i] <- samp
}
names(cts) <- samps

coldata <- coldata[colnames(cts),]
cts <- cts[, rownames(coldata)]
# Check
stopifnot(rownames(coldata)==colnames(cts))

dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design= as.formula(paste('~',dds_design)))

# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 1, ]
# normalization and preprocessing
dds <- DESeq(dds, parallel=parallel)

saveRDS(dds, file=output)

normed_counts <-counts(dds,normalized=TRUE)
write.table(normed_counts,quote=F,sep='\t',file=out_table)

# obtain normalized counts
rld <- rlog(dds, blind=FALSE)
saveRDS(rld, file=rld_out)
