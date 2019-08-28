#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("VariantAnnotation"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--library", default = NA, type = 'character', help = "libray type"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

if (as.character(opt$library)=="STANDARD") {

	vcf_file = "metrics/summary/snps_filtered-standard.vcf"
	vcf = readVcf(vcf_file, "b37")
	gt = geno(vcf)$GT
	ad = geno(vcf)$AD
	af = structure(sapply(ad, function(x) x[2] / sum(x)), dim = dim(ad))
	X = matrix(0, nrow = nrow(gt), ncol = ncol(gt), dimnames = list(rownames(gt), colnames(gt)))
	X[is.na(af)] = NA
	X[af > 0.3 & af < 0.7] = 1
	X[af >= 0.7] = 2
	X[!gt %in% c("0/0", "0/1", "1/1")] = NA
	gt = matrix(as.integer(factor(X)), nrow = nrow(gt), ncol = ncol(gt), dimnames = list(rownames(gt), colnames(gt)))
	write.table(gt, file="metrics/summary/snps_filtered-standard.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

} else if (as.character(opt$library)=="UNFILTERED") {

	vcf_file = "metrics/summary/snps_filtered-unfiltered.vcf"
	vcf = readVcf(vcf_file, "b37")
	gt = geno(vcf)$GT
	ad = geno(vcf)$AD
	af = structure(sapply(ad, function(x) x[2] / sum(x)), dim = dim(ad))
	X = matrix(0, nrow = nrow(gt), ncol = ncol(gt), dimnames = list(rownames(gt), colnames(gt)))
	X[is.na(af)] = NA
	X[af > 0.15 & af < 0.95] = 1
	X[af >= 0.7] = 2
	X[!gt %in% c("0/0", "0/1", "1/1")] = NA
	gt = matrix(as.integer(factor(X)), nrow = nrow(gt), ncol = ncol(gt), dimnames = list(rownames(gt), colnames(gt)))
	write.table(gt, file="metrics/summary/snps_filtered-unfiltered.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

} else if (as.character(opt$library)=="SIMPLEX") {

	vcf_file = "metrics/summary/snps_filtered-simplex.vcf"
	vcf = readVcf(vcf_file, "b37")
	gt = geno(vcf)$GT
	ad = geno(vcf)$AD
	af = structure(sapply(ad, function(x) x[2] / sum(x)), dim = dim(ad))
	X = matrix(0, nrow = nrow(gt), ncol = ncol(gt), dimnames = list(rownames(gt), colnames(gt)))
	X[is.na(af)] = NA
	X[af > 0.3 & af < 0.7] = 1
	X[af >= 0.7] = 2
	X[!gt %in% c("0/0", "0/1", "1/1")] = NA
	gt = matrix(as.integer(factor(X)), nrow = nrow(gt), ncol = ncol(gt), dimnames = list(rownames(gt), colnames(gt)))
	write.table(gt, file="metrics/summary/snps_filtered-simplex.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

} else if (as.character(opt$library)=="DUPLEX") {

	vcf_file = "metrics/summary/snps_filtered-duplex.vcf"
	vcf = readVcf(vcf_file, "b37")
	gt = geno(vcf)$GT
	ad = geno(vcf)$AD
	af = structure(sapply(ad, function(x) x[2] / sum(x)), dim = dim(ad))
	X = matrix(0, nrow = nrow(gt), ncol = ncol(gt), dimnames = list(rownames(gt), colnames(gt)))
	X[is.na(af)] = NA
	X[af > 0.3 & af < 0.7] = 1
	X[af >= 0.7] = 2
	X[!gt %in% c("0/0", "0/1", "1/1")] = NA
	gt = matrix(as.integer(factor(X)), nrow = nrow(gt), ncol = ncol(gt), dimnames = list(rownames(gt), colnames(gt)))
	write.table(gt, file="metrics/summary/snps_filtered-duplex.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

}
