#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("readr"))


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}


args_list <- list(
                	make_option("--samples", default = NA, type = 'character', help = "sample names")
                 )
                 
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

sample_names = unlist(strsplit(x=opt$samples, split=" ", fixed=TRUE))
fusion_summary = list()
for (i in 1:length(sample_names)) {
	fusion_summary[[i]] = read_tsv(file=paste0("integrate_rnaseq/oncofuse/", sample_names[i], ".oncofuse.txt"), col_types = cols(.default = col_character())) %>%
		   			 	  type_convert() %>%
		   			 	  mutate(sample_id = sample_names[i])
}
fusion_summary = do.call(rbind, fusion_summary)
write_tsv(fusion_summary, path="integrate_rnaseq/summary.tsv")
