files <- snakemake@input
files <- unlist(files)

for (i in 1:length(files)) {
    file <- files[i]
    samp <- unlist(strsplit(unlist(strsplit(file, "/"))[3], "[.]"))[1]

    df <- read.table(file, header = TRUE)
    df <- df[,c(2,1)]
    colnames(df) <- c("miRNA", samp)

    if (i==1) {
        all <- df
    } else {
        all <- merge(all, df, by = "miRNA", all = TRUE)
    }
}

write.table(all, file = snakemake@output[[1]], row.names = F, quote = F, sep = "\t")

