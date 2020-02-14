#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--type", default = NA, type = 'character', help = "analysis type"),
				 make_option("--target_file", default = NA, type= 'character', help = 'target file'),
				 make_option("--sample_names", default = NA, type = 'character', help = "sample names"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

cutoffAF = 0.02
chr = 21

if (as.numeric(opt$type)==1) {
	target_positions = read_tsv(file=as.character(opt$target_file), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   			   type_convert() %>%
		   			   rename(chrom = X1, pos = X2) %>%
		   			   mutate(uuid = paste0(chrom, ":", pos))
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	x1 = x2 = x3 = x4 = list()
	for (i in 1:length(sample_names)) {
		df = read_tsv(file=paste0("waltz/", sample_names[i], "-pileup.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   	 readr::type_convert() %>%
		   	 dplyr::select(`chrom`	= X1,
		   	 			   `pos`	= X2,
		   	 			   `ref`	= X3,
		   	 			   `total`	= X4,
		   	 			   `a`		= X5,
		   	 			   `g`		= X6,
		   	 			   `c`		= X7,
		   	 			   `t`		= X8,
		   	 			   `ins`	= X9,
		   	 			   `del`	= X10) %>%
		   	 dplyr::mutate(`total_n` = a+g+c+t,
		   	 			   `uuid` = paste0(chrom, ":", pos)) %>%
		   	 dplyr::filter(uuid %in% target_positions$uuid) %>%
		   	 dplyr::mutate(`counts` = apply(dplyr::tibble(.$ref, .$a, .$g, .$c, .$t, .$total_n), 1, function(x) {
		   	 						ref = as.character(x[1])
		   	 						alt = as.numeric(x[2:6])
		   	 						if (ref=="A") {
		   	 							index = c(1:4)[-1]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						} else if (ref=="G") {
										index = c(1:4)[-2]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						} else if (ref=="C") {
		   	 							index = c(1:4)[-3]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						} else if (ref=="T") {
		   	 							index = c(1:4)[-4]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						}
		   	 						return(y)
		   	 					}))	%>%
		   	 dplyr::filter(!is.na(counts)) %>%
		   	 dplyr::summarize(total_bases = sum(total_n),
		   	 				  alt_counts = sum(counts),
		   	 				  ref_counts = total_bases - alt_counts,
		   	 				  contributing_sites = n())
		   	 x1[[i]] = df
	}
	x1 = do.call(rbind, x1) %>%
		 mutate(sample_names = sample_names, bam_file = "standard")
	for (i in 1:length(sample_names)) {
		df = read_tsv(file=paste0("waltz/", sample_names[i], "__aln_srt_IR_FX-pileup.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   	 readr::type_convert() %>%
		   	 dplyr::select(`chrom`	= X1,
		   	 			   `pos`	= X2,
		   	 			   `ref`	= X3,
		   	 			   `total`	= X4,
		   	 			   `a`		= X5,
		   	 			   `g`		= X6,
		   	 			   `c`		= X7,
		   	 			   `t`		= X8,
		   	 			   `ins`	= X9,
		   	 			   `del`	= X10) %>%
		   	 dplyr::mutate(`total_n` = a+g+c+t,
		   	 			   `uuid` = paste0(chrom, ":", pos)) %>%
		   	 dplyr::filter(uuid %in% target_positions$uuid) %>%
		   	 dplyr::mutate(`counts` = apply(dplyr::tibble(.$ref, .$a, .$g, .$c, .$t, .$total_n), 1, function(x) {
		   	 						ref = as.character(x[1])
		   	 						alt = as.numeric(x[2:6])
		   	 						if (ref=="A") {
		   	 							index = c(1:4)[-1]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						} else if (ref=="G") {
										index = c(1:4)[-2]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						} else if (ref=="C") {
		   	 							index = c(1:4)[-3]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						} else if (ref=="T") {
		   	 							index = c(1:4)[-4]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						}
		   	 						return(y)
		   	 					}))	%>%
		   	 dplyr::filter(!is.na(counts)) %>%
		   	 dplyr::summarize(total_bases = sum(total_n),
		   	 				  alt_counts = sum(counts),
		   	 				  ref_counts = total_bases - alt_counts,
		   	 				  contributing_sites = n())
		   	 x2[[i]] = df
	}
	x2 = do.call(rbind, x2) %>%
		 mutate(sample_names = sample_names, bam_file = "collapsed")
	for (i in 1:length(sample_names)) {
		df = read_tsv(file=paste0("waltz/", sample_names[i], "__aln_srt_IR_FX-simplex-pileup.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   	 readr::type_convert() %>%
		   	 dplyr::select(`chrom`	= X1,
		   	 			   `pos`	= X2,
		   	 			   `ref`	= X3,
		   	 			   `total`	= X4,
		   	 			   `a`		= X5,
		   	 			   `g`		= X6,
		   	 			   `c`		= X7,
		   	 			   `t`		= X8,
		   	 			   `ins`	= X9,
		   	 			   `del`	= X10) %>%
		   	 dplyr::mutate(`total_n` = a+g+c+t,
		   	 			   `uuid` = paste0(chrom, ":", pos)) %>%
		   	 dplyr::filter(uuid %in% target_positions$uuid) %>%
		   	 dplyr::mutate(`counts` = apply(dplyr::tibble(.$ref, .$a, .$g, .$c, .$t, .$total_n), 1, function(x) {
		   	 						ref = as.character(x[1])
		   	 						alt = as.numeric(x[2:6])
		   	 						if (ref=="A") {
		   	 							index = c(1:4)[-1]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						} else if (ref=="G") {
										index = c(1:4)[-2]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						} else if (ref=="C") {
		   	 							index = c(1:4)[-3]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						} else if (ref=="T") {
		   	 							index = c(1:4)[-4]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						}
		   	 						return(y)
		   	 					}))	%>%
		   	 dplyr::filter(!is.na(counts)) %>%
		   	 dplyr::summarize(total_bases = sum(total_n),
		   	 				  alt_counts = sum(counts),
		   	 				  ref_counts = total_bases - alt_counts,
		   	 				  contributing_sites = n())
		   	 x3[[i]] = df
	}
	x3 = do.call(rbind, x3) %>%
		 mutate(sample_names = sample_names, bam_file = "simplex")
	for (i in 1:length(sample_names)) {
		df = read_tsv(file=paste0("waltz/", sample_names[i], "__aln_srt_IR_FX-duplex-pileup.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   	 readr::type_convert() %>%
		   	 dplyr::select(`chrom`	= X1,
		   	 			   `pos`	= X2,
		   	 			   `ref`	= X3,
		   	 			   `total`	= X4,
		   	 			   `a`		= X5,
		   	 			   `g`		= X6,
		   	 			   `c`		= X7,
		   	 			   `t`		= X8,
		   	 			   `ins`	= X9,
		   	 			   `del`	= X10) %>%
		   	 dplyr::mutate(`total_n` = a+g+c+t,
		   	 			   `uuid` = paste0(chrom, ":", pos)) %>%
		   	 dplyr::filter(uuid %in% target_positions$uuid) %>%
		   	 dplyr::mutate(`counts` = apply(dplyr::tibble(.$ref, .$a, .$g, .$c, .$t, .$total_n), 1, function(x) {
		   	 						ref = as.character(x[1])
		   	 						alt = as.numeric(x[2:6])
		   	 						if (ref=="A") {
		   	 							index = c(1:4)[-1]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						} else if (ref=="G") {
										index = c(1:4)[-2]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						} else if (ref=="C") {
		   	 							index = c(1:4)[-3]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						} else if (ref=="T") {
		   	 							index = c(1:4)[-4]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						}
		   	 						return(y)
		   	 					}))	%>%
		   	 dplyr::filter(!is.na(counts)) %>%
		   	 dplyr::summarize(total_bases = sum(total_n),
		   	 				  alt_counts = sum(counts),
		   	 				  ref_counts = total_bases - alt_counts,
		   	 				  contributing_sites = n())
		   	 x4[[i]] = df
	}
	x4 = do.call(rbind, x4) %>%
		 mutate(sample_names = sample_names, bam_file = "duplex")
	x = rbind(x1, x2, x3, x4)
	write_tsv(x, path="waltz/noise_metrics_with_duplicates.txt", na = "NA", append = FALSE, col_names = TRUE)
} else if (as.numeric(opt$type)==2) {
	target_positions = read_tsv(file=as.character(opt$target_file), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   			   type_convert() %>%
		   			   rename(chrom = X1, pos = X2) %>%
		   			   mutate(uuid = paste0(chrom, ":", pos))
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	x1 = x2 = x3 = x4 = list()
	for (i in 1:length(sample_names)) {
		df = read_tsv(file=paste0("waltz/", sample_names[i], "-pileup-without-duplicates.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   	 readr::type_convert() %>%
		   	 dplyr::select(`chrom`	= X1,
		   	 			   `pos`	= X2,
		   	 			   `ref`	= X3,
		   	 			   `total`	= X4,
		   	 			   `a`		= X5,
		   	 			   `g`		= X6,
		   	 			   `c`		= X7,
		   	 			   `t`		= X8,
		   	 			   `ins`	= X9,
		   	 			   `del`	= X10) %>%
		   	 dplyr::mutate(`total_n` = a+g+c+t,
		   	 			   `uuid` = paste0(chrom, ":", pos)) %>%
		   	 dplyr::filter(uuid %in% target_positions$uuid) %>%
		   	 dplyr::mutate(`counts` = apply(dplyr::tibble(.$ref, .$a, .$g, .$c, .$t, .$total_n), 1, function(x) {
		   	 						ref = as.character(x[1])
		   	 						alt = as.numeric(x[2:6])
		   	 						if (ref=="A") {
		   	 							index = c(1:4)[-1]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						} else if (ref=="G") {
										index = c(1:4)[-2]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						} else if (ref=="C") {
		   	 							index = c(1:4)[-3]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						} else if (ref=="T") {
		   	 							index = c(1:4)[-4]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						}
		   	 						return(y)
		   	 					}))	%>%
		   	 dplyr::filter(!is.na(counts)) %>%
		   	 dplyr::summarize(total_bases = sum(total_n),
		   	 				  alt_counts = sum(counts),
		   	 				  ref_counts = total_bases - alt_counts,
		   	 				  contributing_sites = n())
		   	 x1[[i]] = df
	}
	x1 = do.call(rbind, x1) %>%
		 mutate(sample_names = sample_names, bam_file = "standard")
	for (i in 1:length(sample_names)) {
		df = read_tsv(file=paste0("waltz/", sample_names[i], "__aln_srt_IR_FX-pileup-without-duplicates.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   	 readr::type_convert() %>%
		   	 dplyr::select(`chrom`	= X1,
		   	 			   `pos`	= X2,
		   	 			   `ref`	= X3,
		   	 			   `total`	= X4,
		   	 			   `a`		= X5,
		   	 			   `g`		= X6,
		   	 			   `c`		= X7,
		   	 			   `t`		= X8,
		   	 			   `ins`	= X9,
		   	 			   `del`	= X10) %>%
		   	 dplyr::mutate(`total_n` = a+g+c+t,
		   	 			   `uuid` = paste0(chrom, ":", pos)) %>%
		   	 dplyr::filter(uuid %in% target_positions$uuid) %>%
		   	 dplyr::mutate(`counts` = apply(dplyr::tibble(.$ref, .$a, .$g, .$c, .$t, .$total_n), 1, function(x) {
		   	 						ref = as.character(x[1])
		   	 						alt = as.numeric(x[2:6])
		   	 						if (ref=="A") {
		   	 							index = c(1:4)[-1]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						} else if (ref=="G") {
										index = c(1:4)[-2]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						} else if (ref=="C") {
		   	 							index = c(1:4)[-3]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						} else if (ref=="T") {
		   	 							index = c(1:4)[-4]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						}
		   	 						return(y)
		   	 					}))	%>%
		   	 dplyr::filter(!is.na(counts)) %>%
		   	 dplyr::summarize(total_bases = sum(total_n),
		   	 				  alt_counts = sum(counts),
		   	 				  ref_counts = total_bases - alt_counts,
		   	 				  contributing_sites = n())
		   	 x2[[i]] = df
	}
	x2 = do.call(rbind, x2) %>%
		 mutate(sample_names = sample_names, bam_file = "collapsed")
	for (i in 1:length(sample_names)) {
		df = read_tsv(file=paste0("waltz/", sample_names[i], "__aln_srt_IR_FX-simplex-pileup-without-duplicates.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   	 readr::type_convert() %>%
		   	 dplyr::select(`chrom`	= X1,
		   	 			   `pos`	= X2,
		   	 			   `ref`	= X3,
		   	 			   `total`	= X4,
		   	 			   `a`		= X5,
		   	 			   `g`		= X6,
		   	 			   `c`		= X7,
		   	 			   `t`		= X8,
		   	 			   `ins`	= X9,
		   	 			   `del`	= X10) %>%
		   	 dplyr::mutate(`total_n` = a+g+c+t,
		   	 			   `uuid` = paste0(chrom, ":", pos)) %>%
		   	 dplyr::filter(uuid %in% target_positions$uuid) %>%
		   	 dplyr::mutate(`counts` = apply(dplyr::tibble(.$ref, .$a, .$g, .$c, .$t, .$total_n), 1, function(x) {
		   	 						ref = as.character(x[1])
		   	 						alt = as.numeric(x[2:6])
		   	 						if (ref=="A") {
		   	 							index = c(1:4)[-1]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						} else if (ref=="G") {
										index = c(1:4)[-2]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						} else if (ref=="C") {
		   	 							index = c(1:4)[-3]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						} else if (ref=="T") {
		   	 							index = c(1:4)[-4]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						}
		   	 						return(y)
		   	 					}))	%>%
		   	 dplyr::filter(!is.na(counts)) %>%
		   	 dplyr::summarize(total_bases = sum(total_n),
		   	 				  alt_counts = sum(counts),
		   	 				  ref_counts = total_bases - alt_counts,
		   	 				  contributing_sites = n())
		   	 x3[[i]] = df
	}
	x3 = do.call(rbind, x3) %>%
		 mutate(sample_names = sample_names, bam_file = "simplex")
	for (i in 1:length(sample_names)) {
		df = read_tsv(file=paste0("waltz/", sample_names[i], "__aln_srt_IR_FX-duplex-pileup-without-duplicates.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   	 readr::type_convert() %>%
		   	 dplyr::select(`chrom`	= X1,
		   	 			   `pos`	= X2,
		   	 			   `ref`	= X3,
		   	 			   `total`	= X4,
		   	 			   `a`		= X5,
		   	 			   `g`		= X6,
		   	 			   `c`		= X7,
		   	 			   `t`		= X8,
		   	 			   `ins`	= X9,
		   	 			   `del`	= X10) %>%
		   	 dplyr::mutate(`total_n` = a+g+c+t,
		   	 			   `uuid` = paste0(chrom, ":", pos)) %>%
		   	 dplyr::filter(uuid %in% target_positions$uuid) %>%
		   	 dplyr::mutate(`counts` = apply(dplyr::tibble(.$ref, .$a, .$g, .$c, .$t, .$total_n), 1, function(x) {
		   	 						ref = as.character(x[1])
		   	 						alt = as.numeric(x[2:6])
		   	 						if (ref=="A") {
		   	 							index = c(1:4)[-1]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						} else if (ref=="G") {
										index = c(1:4)[-2]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						} else if (ref=="C") {
		   	 							index = c(1:4)[-3]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						} else if (ref=="T") {
		   	 							index = c(1:4)[-4]
		   	 							if (any(alt[index]>cutoffAF*alt[5])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(alt[index])
		   	 							}
		   	 						}
		   	 						return(y)
		   	 					}))	%>%
		   	 dplyr::filter(!is.na(counts)) %>%
		   	 dplyr::summarize(total_bases = sum(total_n),
		   	 				  alt_counts = sum(counts),
		   	 				  ref_counts = total_bases - alt_counts,
		   	 				  contributing_sites = n())
		   	 x4[[i]] = df
	}
	x4 = do.call(rbind, x4) %>%
		 mutate(sample_names = sample_names, bam_file = "duplex")
	x = rbind(x1, x2, x3, x4)
	write_tsv(x, path="waltz/noise_metrics_without_duplicates.txt", na = "NA", append = FALSE, col_names = TRUE)
} else if (as.numeric(opt$type)==3) {
	target_positions = read_tsv(file=as.character(opt$target_file), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   			   type_convert() %>%
		   			   rename(chrom = X1, pos = X2) %>%
		   			   mutate(uuid = paste0(chrom, ":", pos))
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	x1 = list()
	for (i in 1:length(sample_names)) {
		df = read_tsv(file=paste0("waltz/", sample_names[i], "-pileup.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
			 readr::type_convert() %>%
			 dplyr::select(chrom = X1,
									   Position = X2,
									   Reference_Allele = X3,
									   Total_Depth = X4,
									   A = X5,
									   C = X6,
									   G = X7,
									   T = X8) %>%
						 mutate(AF_A = 100*A/Total_Depth,
								AF_C = 100*C/Total_Depth,
								AF_G = 100*G/Total_Depth,
								AF_T = 100*T/Total_Depth)
		nuc_metrics[[i]] = tibble(Chromosome = rep(pileup_metrics$Chromosome, 4),
								  Position = rep(pileup_metrics$Position, 4),
								  Reference_Allele = rep(pileup_metrics$Reference_Allele, 4),
								  Alternate_Allele = c(rep("A", nrow(pileup_metrics)),
													   rep("C", nrow(pileup_metrics)),
													   rep("G", nrow(pileup_metrics)),
													   rep("T", nrow(pileup_metrics))),
								  Allele_Frequency = c(pileup_metrics$AF_A,
													   pileup_metrics$AF_C,
													   pileup_metrics$AF_G,
													   pileup_metrics$AF_T)) %>%
								  filter(Reference_Allele!=Alternate_Allele) %>%
								  filter(Allele_Frequency<AF) %>%
								  filter(Chromosome==CHR) %>%
								  arrange(Position)
	}
	standard_bam = nuc_metrics[[1]][,c("Chromosome", "Position", "Reference_Allele", "Alternate_Allele"),drop=FALSE]
	for (i in 1:length(sample_names)) {
		cat(i, "\n")
		standard_bam = left_join(standard_bam, nuc_metrics[[i]], by=c("Chromosome", "Position", "Reference_Allele", "Alternate_Allele"))
	}
	colnames(standard_bam)[5:ncol(standard_bam)] = sample_names

	nuc_metrics = list()
	for (i in 1:length(sample_names)) {
		pileup_metrics = read_tsv(file=paste0("metrics/standard/", sample_names[i], "-pileup-without-duplicates.txt"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
						 type_convert() %>%
						 dplyr::select(Chromosome = X1,
									   Position = X2,
									   Reference_Allele = X3,
									   Total_Depth = X4,
									   A = X5,
									   C = X6,
									   G = X7,
									   T = X8) %>%
						 mutate(AF_A = 100*A/Total_Depth,
								AF_C = 100*C/Total_Depth,
								AF_G = 100*G/Total_Depth,
								AF_T = 100*T/Total_Depth)
		nuc_metrics[[i]] = tibble(Chromosome = rep(pileup_metrics$Chromosome, 4),
								  Position = rep(pileup_metrics$Position, 4),
								  Reference_Allele = rep(pileup_metrics$Reference_Allele, 4),
								  Alternate_Allele = c(rep("A", nrow(pileup_metrics)),
													   rep("C", nrow(pileup_metrics)),
													   rep("G", nrow(pileup_metrics)),
													   rep("T", nrow(pileup_metrics))),
								  Allele_Frequency = c(pileup_metrics$AF_A,
													   pileup_metrics$AF_C,
													   pileup_metrics$AF_G,
													   pileup_metrics$AF_T)) %>%
								  filter(Reference_Allele!=Alternate_Allele) %>%
								  filter(Allele_Frequency<AF) %>%
								  filter(Chromosome==CHR) %>%
								  arrange(Position)
	}
	standard_bam_dedup = nuc_metrics[[1]][,c("Chromosome", "Position", "Reference_Allele", "Alternate_Allele"),drop=FALSE]
	for (i in 1:length(sample_names)) {
		cat(i, "\n")
		standard_bam_dedup = left_join(standard_bam_dedup, nuc_metrics[[i]], by=c("Chromosome", "Position", "Reference_Allele", "Alternate_Allele"))
	}
	colnames(standard_bam_dedup)[5:ncol(standard_bam_dedup)] = sample_names

	nuc_metrics = list()
	for (i in 1:length(sample_names)) {
		pileup_metrics = read_tsv(file=paste0("metrics/simplex/", sample_names[i], "-pileup.txt"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
						 type_convert() %>%
						 dplyr::select(Chromosome = X1,
									   Position = X2,
									   Reference_Allele = X3,
									   Total_Depth = X4,
									   A = X5,
									   C = X6,
									   G = X7,
									   T = X8) %>%
						 mutate(AF_A = 100*A/Total_Depth,
								AF_C = 100*C/Total_Depth,
								AF_G = 100*G/Total_Depth,
								AF_T = 100*T/Total_Depth)
		nuc_metrics[[i]] = tibble(Chromosome = rep(pileup_metrics$Chromosome, 4),
								  Position = rep(pileup_metrics$Position, 4),
								  Reference_Allele = rep(pileup_metrics$Reference_Allele, 4),
								  Alternate_Allele = c(rep("A", nrow(pileup_metrics)),
													   rep("C", nrow(pileup_metrics)),
													   rep("G", nrow(pileup_metrics)),
													   rep("T", nrow(pileup_metrics))),
								  Allele_Frequency = c(pileup_metrics$AF_A,
													   pileup_metrics$AF_C,
													   pileup_metrics$AF_G,
													   pileup_metrics$AF_T)) %>%
								  filter(Reference_Allele!=Alternate_Allele) %>%
								  filter(Allele_Frequency<AF) %>%
								  filter(Chromosome==CHR) %>%
								  arrange(Position)
	}
	simplex_bam = nuc_metrics[[1]][,c("Chromosome", "Position", "Reference_Allele", "Alternate_Allele"),drop=FALSE]
	for (i in 1:length(sample_names)) {
		cat(i, "\n")
		simplex_bam = left_join(simplex_bam, nuc_metrics[[i]], by=c("Chromosome", "Position", "Reference_Allele", "Alternate_Allele"))
	}
	colnames(simplex_bam)[5:ncol(simplex_bam)] = sample_names

	nuc_metrics = list()
	for (i in 1:length(sample_names)) {
		pileup_metrics = read_tsv(file=paste0("metrics/duplex/", sample_names[i], "-pileup.txt"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
						 type_convert() %>%
						 dplyr::select(Chromosome = X1,
									   Position = X2,
									   Reference_Allele = X3,
									   Total_Depth = X4,
									   A = X5,
									   C = X6,
									   G = X7,
									   T = X8) %>%
						 mutate(AF_A = 100*A/Total_Depth,
								AF_C = 100*C/Total_Depth,
								AF_G = 100*G/Total_Depth,
								AF_T = 100*T/Total_Depth)
		nuc_metrics[[i]] = tibble(Chromosome = rep(pileup_metrics$Chromosome, 4),
								  Position = rep(pileup_metrics$Position, 4),
								  Reference_Allele = rep(pileup_metrics$Reference_Allele, 4),
								  Alternate_Allele = c(rep("A", nrow(pileup_metrics)),
													   rep("C", nrow(pileup_metrics)),
													   rep("G", nrow(pileup_metrics)),
													   rep("T", nrow(pileup_metrics))),
								  Allele_Frequency = c(pileup_metrics$AF_A,
													   pileup_metrics$AF_C,
													   pileup_metrics$AF_G,
													   pileup_metrics$AF_T)) %>%
								  filter(Reference_Allele!=Alternate_Allele) %>%
								  filter(Allele_Frequency<AF) %>%
								  filter(Chromosome==CHR) %>%
								  arrange(Position)
	}
	duplex_bam = nuc_metrics[[1]][,c("Chromosome", "Position", "Reference_Allele", "Alternate_Allele"),drop=FALSE]
	for (i in 1:length(sample_names)) {
		cat(i, "\n")
		duplex_bam = left_join(duplex_bam, nuc_metrics[[i]], by=c("Chromosome", "Position", "Reference_Allele", "Alternate_Allele"))
	}
	colnames(duplex_bam)[5:ncol(duplex_bam)] = sample_names

	nuc_pileup = left_join(standard_bam,
						   standard_bam_dedup,
						   by = c("Chromosome", "Position", "Reference_Allele", "Alternate_Allele")) %>%
				 left_join(simplex_bam,
						   by = c("Chromosome", "Position", "Reference_Allele", "Alternate_Allele")) %>%
				 left_join(duplex_bam,
						   by = c("Chromosome", "Position", "Reference_Allele", "Alternate_Allele"))
	index = order(apply(nuc_pileup[,5:ncol(nuc_pileup),drop=FALSE], 1, mean, na.rm=TRUE), decreasing=FALSE)
	nuc_pileup = nuc_pileup[index,,drop=FALSE] %>%
				 arrange(Alternate_Allele) %>%
				 arrange(Reference_Allele)
			 
			 
	col_groups = rep(c("STANDARD\nWITH DUPLICATES", "STANDARD\nDEDUPLICATED", "COLLAPSED\nSIMPLEX", "COLLAPSED\nDUPLEX"), each=length(sample_names))
	row_groups = paste0(nuc_pileup$Reference_Allele, " > ", nuc_pileup$Alternate_Allele, "         ")
	
	write_tsv(nuc_pileup, path="metrics/summary/snp_pileup_non_ref.tsv", na = "NA", append = FALSE, col_names = TRUE)

	nuc_pileup = nuc_pileup %>%
				 dplyr::select(-Chromosome, -Position, -Reference_Allele, -Alternate_Allele)

	index = apply(nuc_pileup, 1, function(x) {sum(is.na(x))})==0
	pdf(file="metrics/report/non_reference_calls.pdf", height=14, width=14)
	superheat(X = as.matrix(nuc_pileup[index,,drop=FALSE]),
			  smooth.heat = FALSE,
			  scale = FALSE,
			  legend = TRUE,
			  grid.hline = FALSE,
			  grid.vline = FALSE,
			  membership.cols=col_groups,
			  membership.rows=row_groups[index],
			  row.dendrogram = FALSE,
			  col.dendrogram = FALSE,
			  force.grid.hline = FALSE,
			  force.grid.vline = FALSE,
			  bottom.label.text.angle = 0,
			  bottom.label.text.size = 3.5,
			  bottom.label.size = .15,
			  left.label.size = .15,
			  left.label.text.size = 3.5,
			  print.plot = TRUE,
			  heat.pal = viridis(n=10),
			  heat.pal.values = c(seq(0,.4,l=9), 1))
	dev.off()
	

}
