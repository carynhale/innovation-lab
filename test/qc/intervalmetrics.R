#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("readr"))


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--metric_type", default = NA, type = 'character', help = "metric type"),
				 make_option("--sample_names", default = NA, type = 'character', help = "sample names"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

if (as.numeric(opt$metric_type)==1) {

	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	idx_metrics = list()
	for (i in 1:length(sample_names)) {
		idx_metrics[[i]] = read_tsv(file=paste0("metrics/standard/", sample_names[i], ".idx_stats.txt"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   				   type_convert() %>%
		   				   select(X3, X4) %>%
		   				   filter(!is.na(X3)) %>%
		   				   filter(!is.na(X4)) %>%
		   				   mutate(X3 = as.numeric(gsub("Aligned= ", "", X3, fixed=TRUE))) %>%
		   				   mutate(X4 = as.numeric(gsub("Unaligned= ", "", X4, fixed=TRUE))) %>%
		   				   summarize(N_ALIGNED = sum(X3),
		   				   			 N_UNALIGNED = sum(X4)
		   				   ) %>%
		   				   mutate(N_TOTAL = N_ALIGNED + N_UNALIGNED) %>%
		   				   mutate(SAMPLE = sample_names[i])
	}
	idx_metrics = do.call(rbind, idx_metrics)
	write_tsv(x=idx_metrics, path="metrics/standard/idx_metrics.tsv", na = "NA", append = FALSE, col_names = TRUE)
	
} else if (as.numeric(opt$metric_type)==2) {

	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	aln_metrics = list()
	for (i in 1:length(sample_names)) {
		aln_metrics[[i]] = read_tsv(file=paste0("metrics/standard/", sample_names[i], ".aln_metrics.txt"), comment="#", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		   				   type_convert() %>%
		   				   slice(2:n()) %>%
		   				   mutate(TOTAL_READS = as.numeric(TOTAL_READS)) %>%
		   				   mutate(PF_READS = as.numeric(PF_READS)) %>%
		   				   mutate(PCT_PF_READS = 100*as.numeric(PCT_PF_READS)) %>%
		   				   mutate(PF_NOISE_READS = as.numeric(PF_NOISE_READS)) %>%
		   				   mutate(PF_READS_ALIGNED = as.numeric(PF_READS_ALIGNED)) %>%
		   				   mutate(PCT_PF_READS_ALIGNED = 100*as.numeric(PCT_PF_READS_ALIGNED)) %>%
		   				   mutate(PF_ALIGNED_BASES = as.numeric(PF_ALIGNED_BASES)) %>%
		   				   mutate(PF_HQ_ALIGNED_READS = as.numeric(PF_HQ_ALIGNED_READS)) %>%
		   				   mutate(PF_HQ_ALIGNED_BASES = as.numeric(PF_HQ_ALIGNED_BASES)) %>%
		   				   mutate(PF_HQ_ALIGNED_Q20_BASES = as.numeric(PF_HQ_ALIGNED_Q20_BASES)) %>%
		   				   mutate(PF_HQ_MEDIAN_MISMATCHES = as.numeric(PF_HQ_MEDIAN_MISMATCHES)) %>%
		   				   mutate(PF_MISMATCH_RATE = 100*as.numeric(PF_MISMATCH_RATE)) %>%
		   				   mutate(PF_HQ_ERROR_RATE = 100*as.numeric(PF_HQ_ERROR_RATE)) %>%
		   				   mutate(PF_INDEL_RATE = 100*as.numeric(PF_INDEL_RATE)) %>%
		   				   mutate(MEAN_READ_LENGTH = as.numeric(MEAN_READ_LENGTH)) %>%
		   				   mutate(READS_ALIGNED_IN_PAIRS = as.numeric(READS_ALIGNED_IN_PAIRS)) %>%
		   				   mutate(PCT_READS_ALIGNED_IN_PAIRS = 100*as.numeric(PCT_READS_ALIGNED_IN_PAIRS)) %>%
		   				   mutate(BAD_CYCLES = as.numeric(BAD_CYCLES)) %>%
		   				   mutate(STRAND_BALANCE = as.numeric(STRAND_BALANCE)) %>%
		   				   mutate(PCT_CHIMERAS = 100*as.numeric(PCT_CHIMERAS)) %>%
		   				   mutate(PCT_ADAPTER = as.numeric(PCT_ADAPTER)) %>%
		   				   mutate(SAMPLE = sample_names[i])
		   				   
	}
	aln_metrics = do.call(rbind, aln_metrics)
	write_tsv(x=aln_metrics, path="metrics/standard/aln_metrics.tsv", na = "NA", append = FALSE, col_names = TRUE)

}
