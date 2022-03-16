#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("Exactly one arguments must be supplied. Arg1 = Output file")
}
outfile<- args[1]
outfile <- "HapMap_3_r3_1_dropSummary.txt"

#bashcmd <- "grep \"Before .*pruning, there are\" *.log | cut -d: -f1" 
bashcmd <- 'grep "people (" *.log | cut -d: -f1'
fname <- system(bashcmd, intern=T)

bashcmd <- 'grep "people (" *.log | cut -d: -f2 | cut -d" " -f1'
before_indiv <- system(bashcmd, intern=T)
bashcmd <- 'grep "pass filters" *.log | cut -d: -f2 | cut -d" " -f4'
after_indiv <- system(bashcmd, intern=T)

#bashcmd <- "grep \"Before .*pruning, there are\" *.log | cut -d : -f2 | cut -d' ' -f8"
bashcmd <- 'grep "variants loaded" *.log | cut -d: -f2 | cut -d" " -f1'
snps_start <- system(bashcmd, intern=T)
#bashcmd <- "grep \"After .*pruning, there are\" *.log | cut -d : -f2 | cut -d' ' -f8"
bashcmd <- 'grep "pass filters" *.log | cut -d: -f2 | cut -d" " -f1'
snps_end <- system(bashcmd, intern=T)

snps_dropped <- as.numeric(snps_start) - as.numeric(snps_end)
indiv_dropped <- as.numeric(before_indiv) - as.numeric(after_indiv)

df <- data.frame(fname, snps_dropped, indiv_dropped)

df$fname[grepl("P2.log", df$fname)] <- "MAF-0.01"
df$fname[grepl("P3.log", df$fname)] <- "QC-subjects"
df$fname[grepl("P3.1.log", df$fname)] <- "MIND-0.05"
df$fname[grepl("P4.log", df$fname)] <- "excludeNonAutosomalSNP"
df$fname[grepl("P5.log", df$fname)] <- "GENO-0.05"
df$fname[grepl("P6.log", df$fname)] <- "excludeSex"
df$fname[grepl("P7.log", df$fname)] <- "excludeIBD"
df$fname[grepl("P8.log", df$fname)] <- "HWE-1E-6"
df$fname[grepl("IBD.het.log", df$fname)] <- "excludeHeterozygosity"
df$fname[grepl("QC_subjects.log", df$fname)] <- "dropQCsubjects"
df$fname[grepl("P6.IBD.log", df$fname)] <- "excludeRelatedness"
df$fname[grepl("het_check.log", df$fname)] <- "checkHeterozygosity"

write.table(df, file=outfile, row.name=F, quote=F)
