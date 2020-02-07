#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--switch", default = NA, type = 'character', help = "libray type"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

'read_yaml' <- function(file)
{
	yaml = read.csv(file=file, header=FALSE, sep="\t", stringsAsFactors=FALSE)
	index = grep("- name: ", yaml[,1])
	x = gsub(pattern="- name: ", replacement="", yaml[index,1])
	y = list()
	for (i in 1:length(index)) {
		ii = index[i] + 1
		if (i == length(index)) {
			jj = nrow(yaml)
		} else {
			jj = index[i+1] - 1
		}
		string = yaml[ii:jj,1]
		marker = c(" ", ":", "[", "]", "normal", "tumor")
		for (j in 1:length(marker)) {
			string = gsub(pattern=marker[j], replacement="", string, fixed=TRUE)
		}
		string = unlist(lapply(string, function(x) { unlist(strsplit(x, ",", fixed=TRUE)) }))
		y[[i]] = string
	}
	names(y) = x
	return(invisible(y))
}

if (as.numeric(opt$switch)==1) {

	suppressPackageStartupMessages(library("VariantAnnotation"))

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

	suppressPackageStartupMessages(library("readr"))
	suppressPackageStartupMessages(library("dplyr"))
	suppressPackageStartupMessages(library("magrittr"))
	suppressPackageStartupMessages(library("superheat"))
	suppressPackageStartupMessages(library("viridis"))

	yaml = read_yaml("samples.yaml")
	data = read.csv(file="metrics/summary/snps_filtered.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
	data = data	 %>%
		   rename_all(funs(gsub(pattern=".", replacement="-", x=make.names(names(data)), fixed=TRUE))) %>%
		   type_convert()
	data = data[apply(data, 1, function(x) {sum(is.na(x))})<=(.5*ncol(data)),,drop=FALSE]
	index = is.na(data)
	data[index] = 1
	data[data==3] = 1
	data[index] = NA
	for (i in 1:2) {
		data = data[apply(data, 1, function(x) {sum(x==i, na.rm=TRUE)})!=ncol(data),,drop=FALSE]
	}
	dm = as.matrix(dist(t(data), method="manhattan", diag=TRUE, upper=TRUE))
	dm = 1-((max(dm)-dm)/(max(dm) - min(dm)))
	col_groups = list()
	for (i in 1:length(yaml)) {
		col_groups[[i]] = rep(names(yaml)[i], length(yaml[[i]]))
	}
	col_groups = unlist(col_groups)
	dm = dm[as.vector(unlist(yaml)),as.vector(unlist(yaml)),drop=FALSE]
	
	d = dist(dm, method = "euclidean")
	h = hclust(d, method = "ward.D")
	dm = dm[,h$order,drop=FALSE]
	col_groups = col_groups[h$order]
	
	pdf(file="metrics/report/snps_clustering.pdf", width=14, height=14)
	superheat(X = dm,
			  smooth.heat = FALSE,
			  scale = FALSE,
			  legend = TRUE,
			  grid.hline = FALSE,
			  grid.vline = FALSE,
			  force.grid.hline = FALSE,
			  force.grid.vline = FALSE,
			  row.dendrogram = TRUE,
			  col.dendrogram = FALSE,
			  membership.cols = col_groups,
			  clustering.method = "hierarchical",
			  dist.method = "euclidean",
			  linkage.method = "ward.D",
			  bottom.label.text.angle = 90,
			  bottom.label.text.size = 3.5,
			  bottom.label.size = .1,
			  left.label.text.size = 3.5,
			  left.label.size = .1,
			  left.label.text.alignment = "center",
			  heat.pal = viridis(n=100),
			  heat.pal.values = c(seq(from=0, to=.4, length=70), seq(from=.41, to=1, length=30)),
			  print.plot = TRUE)
	dev.off()

}
