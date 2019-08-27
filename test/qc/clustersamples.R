#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("VariantAnnotation"))

vcf_file = "metrics/summary/snps-filtered.vcf"
vcf = readVcf(vcf_file, "b37")
gt = geno(vcf)$GT
ad = geno(vcf)$AD
af = structure(sapply(ad, function(x) x[2] / sum(x)), dim = dim(ad))
X = matrix(0, nrow = nrow(gt), ncol = ncol(gt), dimnames = list(rownames(gt), colnames(gt)))
X[is.na(af)] = NA
X[af > 0.15 & af < 0.95] = 1
X[af >= 0.95] = 2
X[!gt %in% c("0/0", "0/1", "1/1")] = NA
gt = matrix(as.integer(factor(X)), nrow = nrow(gt), ncol = ncol(gt), dimnames = list(rownames(gt), colnames(gt)))
write.table(gt, file="metrics/summary/snps-filtered.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
