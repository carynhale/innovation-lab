#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("readr"))


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--metric", default = NA, type = 'character', help = "metric type"),
				 make_option("--samples", default = NA, type = 'character', help = "sample names"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

if (as.numeric(opt$metric)==1) {

	sample_names = unlist(strsplit(x=as.character(opt$samples), split=" ", fixed=TRUE))
	idx_metrics = list()
	for (i in 1:length(sample_names)) {
		idx_metrics[[i]] = read_tsv(file=paste0("metrics/picard/", sample_names[i], "-idx_stats.txt"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   				   type_convert() %>%
		   				   select(X3, X4) %>%
		   				   filter(!(is.na(X3) | is.na(X4))) %>%
		   				   mutate(X3 = as.numeric(gsub("Aligned= ", "", X3, fixed=TRUE))) %>%
		   				   mutate(X4 = as.numeric(gsub("Unaligned= ", "", X4, fixed=TRUE))) %>%
		   				   summarize(N_ALIGNED = sum(X3),
		   				   			 N_UNALIGNED = sum(X4)
		   				   ) %>%
		   				   mutate(N_TOTAL = N_ALIGNED + N_UNALIGNED) %>%
		   				   mutate(SAMPLE = sample_names[i])
	}
	idx_metrics = do.call(rbind, idx_metrics)
	write_tsv(x=idx_metrics, path="metrics/summary/metrics_idx.tsv", na = "NA", append = FALSE, col_names = TRUE)
	
}
