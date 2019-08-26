#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("ggplot2"))


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--type", default = NA, type = 'character', help = "analysis type"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

if (as.numeric(opt$type)==1) {

	data = read.csv(file="metrics/summary/umi_frequencies.tsv", sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE) %>%
		   rename_all(funs(gsub(pattern=".", replacement="-", x=make.names(names(data)), fixed=TRUE))) %>%
		   type_convert() %>%
		   replace(is.na(.), 0)
	pdf(file="metrics/report/umi_frequencies.pdf", width=7, height=14)
	heatmap(data, Rowv=NA, Colv=NA, scale="column")
	dev.off()
	
}
