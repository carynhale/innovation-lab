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
	write_tsv(x=idx_metrics, path="metrics/standard/metrics_idx.tsv", na = "NA", append = FALSE, col_names = TRUE)
	
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
	write_tsv(x=aln_metrics, path="metrics/standard/metrics_aln.tsv", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$metric_type)==3) {

	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	insert_metrics = list()
	for (i in 1:length(sample_names)) {
		insert_metrics[[i]] = read_tsv(file=paste0("metrics/standard/", sample_names[i], ".insert_metrics.txt"), comment="#", col_names = TRUE, n_max = 2, col_types = cols(.default = col_character())) %>%
		   				      slice(2:n()) %>%
		   				      type_convert() %>%
		   				      mutate(SAMPLE = sample_names[i])
		   				   
	}
	insert_metrics = do.call(rbind, insert_metrics)
	write_tsv(x=insert_metrics, path="metrics/standard/metrics_insert.tsv", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$metric_type)==4) {

	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	insert_metrics = list()
	for (i in 1:length(sample_names)) {
		insert_metrics[[i]] = read_tsv(file=paste0("metrics/standard/", sample_names[i], ".insert_metrics.txt"), comment="#", col_names = TRUE, skip = 3, col_types = cols(.default = col_character())) %>%
		   				      slice(2:n()) %>%
		   				      type_convert()
	}
	n = max(unlist(lapply(insert_metrics, function(x) {max(x$`insert_size`)})))
	insert_distribution = matrix(NA, nrow=n, ncol=length(insert_metrics))
	for (i in 1:length(insert_metrics)) {
		ix = as.numeric(as.character(insert_metrics[[i]]$`insert_size`))
		iy = as.numeric(as.character(insert_metrics[[i]]$`All_Reads.fr_count`))
		insert_distribution[ix,i] = iy
	}
	colnames(insert_distribution) = sample_names
	insert_distribution = as.data.frame(insert_distribution)
	write_tsv(x=insert_distribution, path="metrics/standard/metrics_insert_distribution.tsv", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$metric_type)==5) {

	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	oxog_metrics = list()
	for (i in 1:length(sample_names)) {
		oxog_metrics[[i]] = read_tsv(file=paste0("metrics/standard/", sample_names[i], ".oxog_metrics.txt"), comment="#", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		   				    slice(2:n()) %>%
		   				    type_convert()
	}
	oxog_metrics = do.call(rbind, oxog_metrics)
	write_tsv(x=oxog_metrics, path="metrics/standard/metrics_oxog.tsv", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$metric_type)==6) {

	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	hs_metrics_a = hs_metrics_b = list()
	for (i in 1:length(sample_names)) {
		hs_metrics_a[[i]] = read_tsv(file=paste0("metrics/standard/", sample_names[i], ".probe-A.hs_metrics.txt"), comment="#", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		   				    slice(2:n()) %>%
		   				    type_convert() %>%
		   				    mutate(SAMPLE = sample_names[i])
		hs_metrics_b[[i]] = read_tsv(file=paste0("metrics/standard/", sample_names[i], ".probe-B.hs_metrics.txt"), comment="#", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		   				    slice(2:n()) %>%
		   				    type_convert() %>%
		   				    mutate(SAMPLE = sample_names[i])
		   				   
	}
	hs_metrics_a = do.call(rbind, hs_metrics_a)
	hs_metrics_b = do.call(rbind, hs_metrics_b)
	hs_metrics = rbind(hs_metrics_a, hs_metrics_b) %>%
				 arrange(SAMPLE) %>%
				 dplyr::select(GENOME_SIZE,
				 			   BAIT_SET,
				 			   TARGET_TERRITORY,
				 			   ON_TARGET_BASES,
				 			   MEAN_TARGET_COVERAGE,
				 			   PCT_TARGET_BASES_2X,
				 			   PCT_TARGET_BASES_10X,
				 			   PCT_TARGET_BASES_20X,
				 			   PCT_TARGET_BASES_30X,
				 			   PCT_TARGET_BASES_40X,
				 			   PCT_TARGET_BASES_50X,
				 			   PCT_TARGET_BASES_100X,
				 			   AT_DROPOUT,
				 			   GC_DROPOUT,
				 			   SAMPLE) %>%
				 mutate(PCT_TARGET_BASES_2X = 100*PCT_TARGET_BASES_2X,
				 		PCT_TARGET_BASES_10X = 100*PCT_TARGET_BASES_10X,
				 		PCT_TARGET_BASES_20X = 100*PCT_TARGET_BASES_20X,
				 		PCT_TARGET_BASES_30X = 100*PCT_TARGET_BASES_30X,
				 		PCT_TARGET_BASES_40X = 100*PCT_TARGET_BASES_40X,
				 		PCT_TARGET_BASES_50X = 100*PCT_TARGET_BASES_50X,
				 		PCT_TARGET_BASES_100X = 100*PCT_TARGET_BASES_100X)
	write_tsv(x=hs_metrics, path="metrics/standard/metrics_hs.tsv", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$metric_type)==7) {

	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	idx_metrics = list()
	for (i in 1:length(sample_names)) {
		idx_metrics[[i]] = read_tsv(file=paste0("metrics/unfiltered/", sample_names[i], ".idx_stats.txt"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
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
	write_tsv(x=idx_metrics, path="metrics/unfiltered/metrics_idx.tsv", na = "NA", append = FALSE, col_names = TRUE)
	
} else if (as.numeric(opt$metric_type)==8) {

	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	aln_metrics = list()
	for (i in 1:length(sample_names)) {
		aln_metrics[[i]] = read_tsv(file=paste0("metrics/unfiltered/", sample_names[i], ".aln_metrics.txt"), comment="#", col_names = TRUE, col_types = cols(.default = col_character())) %>%
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
	write_tsv(x=aln_metrics, path="metrics/unfiltered/metrics_aln.tsv", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$metric_type)==9) {

	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	insert_metrics = list()
	for (i in 1:length(sample_names)) {
		insert_metrics[[i]] = read_tsv(file=paste0("metrics/unfiltered/", sample_names[i], ".insert_metrics.txt"), comment="#", col_names = TRUE, n_max = 2, col_types = cols(.default = col_character())) %>%
		   				      slice(2:n()) %>%
		   				      type_convert() %>%
		   				      mutate(SAMPLE = sample_names[i])
		   				   
	}
	insert_metrics = do.call(rbind, insert_metrics)
	write_tsv(x=insert_metrics, path="metrics/unfiltered/metrics_insert.tsv", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$metric_type)==10) {

	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	insert_metrics = list()
	for (i in 1:length(sample_names)) {
		insert_metrics[[i]] = read_tsv(file=paste0("metrics/unfiltered/", sample_names[i], ".insert_metrics.txt"), comment="#", col_names = TRUE, skip = 3, col_types = cols(.default = col_character())) %>%
		   				      slice(2:n()) %>%
		   				      type_convert()
	}
	n = max(unlist(lapply(insert_metrics, function(x) {max(x$`insert_size`)})))
	insert_distribution = matrix(NA, nrow=n, ncol=length(insert_metrics))
	for (i in 1:length(insert_metrics)) {
		ix = as.numeric(as.character(insert_metrics[[i]]$`insert_size`))
		iy = as.numeric(as.character(insert_metrics[[i]]$`All_Reads.fr_count`))
		insert_distribution[ix,i] = iy
	}
	colnames(insert_distribution) = sample_names
	insert_distribution = as.data.frame(insert_distribution)
	write_tsv(x=insert_distribution, path="metrics/unfiltered/metrics_insert_distribution.tsv", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$metric_type)==11) {

	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	hs_metrics_a = hs_metrics_b = list()
	for (i in 1:length(sample_names)) {
		hs_metrics_a[[i]] = read_tsv(file=paste0("metrics/unfiltered/", sample_names[i], ".probe-A.hs_metrics.txt"), comment="#", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		   				    slice(2:n()) %>%
		   				    type_convert() %>%
		   				    mutate(SAMPLE = sample_names[i])
		hs_metrics_b[[i]] = read_tsv(file=paste0("metrics/unfiltered/", sample_names[i], ".probe-B.hs_metrics.txt"), comment="#", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		   				    slice(2:n()) %>%
		   				    type_convert() %>%
		   				    mutate(SAMPLE = sample_names[i])
		   				   
	}
	hs_metrics_a = do.call(rbind, hs_metrics_a)
	hs_metrics_b = do.call(rbind, hs_metrics_b)
	hs_metrics = rbind(hs_metrics_a, hs_metrics_b) %>%
				 arrange(SAMPLE) %>%
				 dplyr::select(GENOME_SIZE,
				 			   BAIT_SET,
				 			   TARGET_TERRITORY,
				 			   ON_TARGET_BASES,
				 			   MEAN_TARGET_COVERAGE,
				 			   PCT_TARGET_BASES_2X,
				 			   PCT_TARGET_BASES_10X,
				 			   PCT_TARGET_BASES_20X,
				 			   PCT_TARGET_BASES_30X,
				 			   PCT_TARGET_BASES_40X,
				 			   PCT_TARGET_BASES_50X,
				 			   PCT_TARGET_BASES_100X,
				 			   AT_DROPOUT,
				 			   GC_DROPOUT,
				 			   SAMPLE) %>%
				 mutate(PCT_TARGET_BASES_2X = 100*PCT_TARGET_BASES_2X,
				 		PCT_TARGET_BASES_10X = 100*PCT_TARGET_BASES_10X,
				 		PCT_TARGET_BASES_20X = 100*PCT_TARGET_BASES_20X,
				 		PCT_TARGET_BASES_30X = 100*PCT_TARGET_BASES_30X,
				 		PCT_TARGET_BASES_40X = 100*PCT_TARGET_BASES_40X,
				 		PCT_TARGET_BASES_50X = 100*PCT_TARGET_BASES_50X,
				 		PCT_TARGET_BASES_100X = 100*PCT_TARGET_BASES_100X)
	write_tsv(x=hs_metrics, path="metrics/unfiltered/metrics_hs.tsv", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$metric_type)==12) {

	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	idx_metrics = list()
	for (i in 1:length(sample_names)) {
		idx_metrics[[i]] = read_tsv(file=paste0("metrics/duplex/", sample_names[i], ".idx_stats.txt"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
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
	write_tsv(x=idx_metrics, path="metrics/duplex/metrics_idx.tsv", na = "NA", append = FALSE, col_names = TRUE)
	
} else if (as.numeric(opt$metric_type)==13) {

	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	aln_metrics = list()
	for (i in 1:length(sample_names)) {
		aln_metrics[[i]] = read_tsv(file=paste0("metrics/duplex/", sample_names[i], ".aln_metrics.txt"), comment="#", col_names = TRUE, col_types = cols(.default = col_character())) %>%
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
	write_tsv(x=aln_metrics, path="metrics/duplex/metrics_aln.tsv", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$metric_type)==14) {

	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	insert_metrics = list()
	for (i in 1:length(sample_names)) {
		insert_metrics[[i]] = read_tsv(file=paste0("metrics/duplex/", sample_names[i], ".insert_metrics.txt"), comment="#", col_names = TRUE, n_max = 2, col_types = cols(.default = col_character())) %>%
		   				      slice(2:n()) %>%
		   				      type_convert() %>%
		   				      mutate(SAMPLE = sample_names[i])
		   				   
	}
	insert_metrics = do.call(rbind, insert_metrics)
	write_tsv(x=insert_metrics, path="metrics/duplex/metrics_insert.tsv", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$metric_type)==15) {

	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	insert_metrics = list()
	for (i in 1:length(sample_names)) {
		insert_metrics[[i]] = read_tsv(file=paste0("metrics/duplex/", sample_names[i], ".insert_metrics.txt"), comment="#", col_names = TRUE, skip = 3, col_types = cols(.default = col_character())) %>%
		   				      slice(2:n()) %>%
		   				      type_convert()
	}
	n = max(unlist(lapply(insert_metrics, function(x) {max(x$`insert_size`)})))
	insert_distribution = matrix(NA, nrow=n, ncol=length(insert_metrics))
	for (i in 1:length(insert_metrics)) {
		ix = as.numeric(as.character(insert_metrics[[i]]$`insert_size`))
		iy = as.numeric(as.character(insert_metrics[[i]]$`All_Reads.fr_count`))
		insert_distribution[ix,i] = iy
	}
	colnames(insert_distribution) = sample_names
	insert_distribution = as.data.frame(insert_distribution)
	write_tsv(x=insert_distribution, path="metrics/duplex/metrics_insert_distribution.tsv", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$metric_type)==16) {

	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	hs_metrics_a = hs_metrics_b = list()
	for (i in 1:length(sample_names)) {
		hs_metrics_a[[i]] = read_tsv(file=paste0("metrics/duplex/", sample_names[i], ".probe-A.hs_metrics.txt"), comment="#", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		   				    slice(2:n()) %>%
		   				    type_convert() %>%
		   				    mutate(SAMPLE = sample_names[i])
		hs_metrics_b[[i]] = read_tsv(file=paste0("metrics/duplex/", sample_names[i], ".probe-B.hs_metrics.txt"), comment="#", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		   				    slice(2:n()) %>%
		   				    type_convert() %>%
		   				    mutate(SAMPLE = sample_names[i])
		   				   
	}
	hs_metrics_a = do.call(rbind, hs_metrics_a)
	hs_metrics_b = do.call(rbind, hs_metrics_b)
	hs_metrics = rbind(hs_metrics_a, hs_metrics_b) %>%
				 arrange(SAMPLE) %>%
				 dplyr::select(GENOME_SIZE,
				 			   BAIT_SET,
				 			   TARGET_TERRITORY,
				 			   ON_TARGET_BASES,
				 			   MEAN_TARGET_COVERAGE,
				 			   PCT_TARGET_BASES_2X,
				 			   PCT_TARGET_BASES_10X,
				 			   PCT_TARGET_BASES_20X,
				 			   PCT_TARGET_BASES_30X,
				 			   PCT_TARGET_BASES_40X,
				 			   PCT_TARGET_BASES_50X,
				 			   PCT_TARGET_BASES_100X,
				 			   AT_DROPOUT,
				 			   GC_DROPOUT,
				 			   SAMPLE) %>%
				 mutate(PCT_TARGET_BASES_2X = 100*PCT_TARGET_BASES_2X,
				 		PCT_TARGET_BASES_10X = 100*PCT_TARGET_BASES_10X,
				 		PCT_TARGET_BASES_20X = 100*PCT_TARGET_BASES_20X,
				 		PCT_TARGET_BASES_30X = 100*PCT_TARGET_BASES_30X,
				 		PCT_TARGET_BASES_40X = 100*PCT_TARGET_BASES_40X,
				 		PCT_TARGET_BASES_50X = 100*PCT_TARGET_BASES_50X,
				 		PCT_TARGET_BASES_100X = 100*PCT_TARGET_BASES_100X)
	write_tsv(x=hs_metrics, path="metrics/duplex/metrics_hs.tsv", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$metric_type)==17) {

	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	idx_metrics = list()
	for (i in 1:length(sample_names)) {
		idx_metrics[[i]] = read_tsv(file=paste0("metrics/simplex/", sample_names[i], ".idx_stats.txt"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
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
	write_tsv(x=idx_metrics, path="metrics/simplex/metrics_idx.tsv", na = "NA", append = FALSE, col_names = TRUE)
	
} else if (as.numeric(opt$metric_type)==18) {

	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	aln_metrics = list()
	for (i in 1:length(sample_names)) {
		aln_metrics[[i]] = read_tsv(file=paste0("metrics/simplex/", sample_names[i], ".aln_metrics.txt"), comment="#", col_names = TRUE, col_types = cols(.default = col_character())) %>%
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
	write_tsv(x=aln_metrics, path="metrics/simplex/metrics_aln.tsv", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$metric_type)==19) {

	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	insert_metrics = list()
	for (i in 1:length(sample_names)) {
		insert_metrics[[i]] = read_tsv(file=paste0("metrics/simplex/", sample_names[i], ".insert_metrics.txt"), comment="#", col_names = TRUE, n_max = 2, col_types = cols(.default = col_character())) %>%
		   				      slice(2:n()) %>%
		   				      type_convert() %>%
		   				      mutate(SAMPLE = sample_names[i])
		   				   
	}
	insert_metrics = do.call(rbind, insert_metrics)
	write_tsv(x=insert_metrics, path="metrics/simplex/metrics_insert.tsv", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$metric_type)==20) {

	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	insert_metrics = list()
	for (i in 1:length(sample_names)) {
		insert_metrics[[i]] = read_tsv(file=paste0("metrics/simplex/", sample_names[i], ".insert_metrics.txt"), comment="#", col_names = TRUE, skip = 3, col_types = cols(.default = col_character())) %>%
		   				      slice(2:n()) %>%
		   				      type_convert()
	}
	n = max(unlist(lapply(insert_metrics, function(x) {max(x$`insert_size`)})))
	insert_distribution = matrix(NA, nrow=n, ncol=length(insert_metrics))
	for (i in 1:length(insert_metrics)) {
		ix = as.numeric(as.character(insert_metrics[[i]]$`insert_size`))
		iy = as.numeric(as.character(insert_metrics[[i]]$`All_Reads.fr_count`))
		insert_distribution[ix,i] = iy
	}
	colnames(insert_distribution) = sample_names
	insert_distribution = as.data.frame(insert_distribution)
	write_tsv(x=insert_distribution, path="metrics/simplex/metrics_insert_distribution.tsv", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$metric_type)==21) {

	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	hs_metrics_a = hs_metrics_b = list()
	for (i in 1:length(sample_names)) {
		hs_metrics_a[[i]] = read_tsv(file=paste0("metrics/simplex/", sample_names[i], ".probe-A.hs_metrics.txt"), comment="#", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		   				    slice(2:n()) %>%
		   				    type_convert() %>%
		   				    mutate(SAMPLE = sample_names[i])
		hs_metrics_b[[i]] = read_tsv(file=paste0("metrics/simplex/", sample_names[i], ".probe-B.hs_metrics.txt"), comment="#", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		   				    slice(2:n()) %>%
		   				    type_convert() %>%
		   				    mutate(SAMPLE = sample_names[i])
		   				   
	}
	hs_metrics_a = do.call(rbind, hs_metrics_a)
	hs_metrics_b = do.call(rbind, hs_metrics_b)
	hs_metrics = rbind(hs_metrics_a, hs_metrics_b) %>%
				 arrange(SAMPLE) %>%
				 dplyr::select(GENOME_SIZE,
				 			   BAIT_SET,
				 			   TARGET_TERRITORY,
				 			   ON_TARGET_BASES,
				 			   MEAN_TARGET_COVERAGE,
				 			   PCT_TARGET_BASES_2X,
				 			   PCT_TARGET_BASES_10X,
				 			   PCT_TARGET_BASES_20X,
				 			   PCT_TARGET_BASES_30X,
				 			   PCT_TARGET_BASES_40X,
				 			   PCT_TARGET_BASES_50X,
				 			   PCT_TARGET_BASES_100X,
				 			   AT_DROPOUT,
				 			   GC_DROPOUT,
				 			   SAMPLE) %>%
				 mutate(PCT_TARGET_BASES_2X = 100*PCT_TARGET_BASES_2X,
				 		PCT_TARGET_BASES_10X = 100*PCT_TARGET_BASES_10X,
				 		PCT_TARGET_BASES_20X = 100*PCT_TARGET_BASES_20X,
				 		PCT_TARGET_BASES_30X = 100*PCT_TARGET_BASES_30X,
				 		PCT_TARGET_BASES_40X = 100*PCT_TARGET_BASES_40X,
				 		PCT_TARGET_BASES_50X = 100*PCT_TARGET_BASES_50X,
				 		PCT_TARGET_BASES_100X = 100*PCT_TARGET_BASES_100X)
	write_tsv(x=hs_metrics, path="metrics/simplex/metrics_hs.tsv", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$metric_type)==22) {

	standard = read_tsv(file="metrics/standard/metrics_idx.tsv", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		   	   type_convert() %>%
		   	   mutate(TYPE = "STANDARD")
	unfiltered = read_tsv(file="metrics/unfiltered/metrics_idx.tsv", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		   	     type_convert() %>%
		   	     mutate(TYPE = "UNFILTERED")
	simplex = read_tsv(file="metrics/simplex/metrics_idx.tsv", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		   	  type_convert() %>%
		   	  mutate(TYPE = "SIMPLEX")
	duplex = read_tsv(file="metrics/duplex/metrics_idx.tsv", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		   	 type_convert() %>%
		   	 mutate(TYPE = "DUPLEX")
	metrics_idx = bind_rows(standard, unfiltered, simplex, duplex)
	write_tsv(x=metrics_idx, path="metrics/summary/metrics_idx.tsv", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$metric_type)==23) {

	standard = read_tsv(file="metrics/standard/metrics_aln.tsv", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		   	   type_convert() %>%
		   	   mutate(TYPE = "STANDARD")
	unfiltered = read_tsv(file="metrics/unfiltered/metrics_aln.tsv", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		   	     type_convert() %>%
		   	     mutate(TYPE = "UNFILTERED")
	simplex = read_tsv(file="metrics/simplex/metrics_aln.tsv", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		   	  type_convert() %>%
		   	  mutate(TYPE = "SIMPLEX")
	duplex = read_tsv(file="metrics/duplex/metrics_aln.tsv", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		   	 type_convert() %>%
		   	 mutate(TYPE = "DUPLEX")
	metrics_aln = bind_rows(standard, unfiltered, simplex, duplex)
	write_tsv(x=metrics_aln, path="metrics/summary/metrics_aln.tsv", na = "NA", append = FALSE, col_names = TRUE)

}
