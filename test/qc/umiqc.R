#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--type", default = NA, type = 'character', help = "analysis type"),
				 make_option("--sample_names", default = NA, type = 'character', help = "sample names"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

if (as.numeric(opt$type)==1) {

	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	umi_frequencies = list()
	umi_list = list()
	for (i in 1:length(sample_names)) {
		umi_frequencies[[i]] = read_tsv(file=paste0("marianas/", sample_names[i], "/umi-frequencies.txt"), col_names=FALSE, col_types = cols(.default = col_character())) %>%
		   					   type_convert()
		umi_list[[i]] = umi_frequencies %>%
						.[["X1"]]
	}
	umi_list = unique(unlist(umi_list))
	umis = data.frame(matrix(0, nrow=length(umi_list), ncol=length(sample_names), dimnames=list(umi_list, sample_names)))
	for (i in 1:length(sample_names)) {
		ix = umi_frequencies[[i]] %>%
			 .[["X1"]]
		iy = umi_frequencies[[i]] %>%
			 .[["X2"]]
		umis[ix,i] = iy
	}
	colnames(umis) = gsub(pattern=".", replacement="-", x=colnames(umis), fixed=TRUE)
	write.table(umis, file="metrics/summary/umi-frequencies.tsv", sep="\t", col.names=TRUE, row.names=TRUE, append=FALSE, quote=FALSE)

}
	
