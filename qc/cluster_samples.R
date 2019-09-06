#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("VariantAnnotation"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--switch", default = NA, type = 'character', help = "libray type"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

if (as.numeric(opt$switch)==1) {

	vcf_file = "metrics/summary/snps_filtered.vcf"
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
	write.table(gt, file="metrics/summary/snps_filtered.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

} else if (as.numeric(opt$switch)==2) {

	suppressPackageStartupMessages(library("superheat"))
	suppressPackageStartupMessages(library("viridis"))

	data = read.csv(file="metrics/summary/snps_filtered.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
	data = data	 %>%
		   rename_all(funs(gsub(pattern=".", replacement="-", x=make.names(names(data)), fixed=TRUE))) %>%
		   type_convert()
	data = data[apply(data, 1, function(x) {sum(is.na(x))})!=ncol(data),,drop=FALSE]
	data[is.na(data)] = 1
	data[data==3] = 1
	for (i in 1:2) {
		data = data[apply(data, 1, function(x) {sum(x==i, na.rm=TRUE)})!=ncol(data),,drop=FALSE]
	}
	dm = as.matrix(dist(t(data), method="manhattan", diag=TRUE, upper=TRUE))
	dm = 1-((max(dm)-dm)/(max(dm) - min(dm)))
	pdf(file="metrics/report/snps_clustering.pdf", width=14, height=14)
	superheat(X = dm, smooth.heat = TRUE, scale = FALSE, legend = TRUE, grid.hline = TRUE, grid.vline = TRUE,
			  row.dendrogram = TRUE, col.dendrogram = TRUE, force.grid.hline = TRUE, force.grid.vline = TRUE,
			  bottom.label.text.angle = 90, bottom.label.text.size = 3.5, bottom.label.size = .15,
			  left.label.size = .15, left.label.text.size = 3.5, grid.hline.col = "grey90",
			  grid.vline.col = "grey90", heat.pal = viridis(n=100), heat.pal.values = seq(from=0, to=1, length=100),
			  print.plot = TRUE)
	dev.off()

}
