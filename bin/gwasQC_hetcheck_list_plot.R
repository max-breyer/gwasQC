#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Exactly two arguments must be supplied. Arg1 = Input file, Arg2 = Output file prefix")
}
infile <- args[1]
outfile <- args[2]

het <- read.table(infile, header = T)

het$het_dist_std <- (het$Het_rate - mean(het$Het_rate))/(sd(het$Het_rate))
het_fail = subset(het, (het$het_dist_std < -4 |  het$het_dist_std > 4))
het_fail_ids <- het_fail[c("FID", "IID")]
write.table(het_fail, paste0(outfile, "excess_het_outliers.txt"), row.names = FALSE, quote = FALSE)
write.table(het_fail_ids, paste0(outfile, "excess_het_outliers_ids_only.txt"), row.names = FALSE, quote = FALSE) 

pdf(paste0(outfile, "heterozygosity.pdf"))
par(mfrow=c(1,2))
hist(het$Het_rate, xlab="Heterozygosity Rate", ylab="Frequency")
hist(het$het_dist_std, xlab="Standardized Heterozygosity Rate", ylab="Frequency")
dev.off()




