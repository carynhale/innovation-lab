#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--type", default = NA, type = 'character', help = "analysis type"),
				 make_option("--target_file", default = NA, type= 'character', help = 'target file')
				 make_option("--sample_names", default = NA, type = 'character', help = "sample names"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

if (as.numeric(opt$type)==1) {
	target_positions = read_tsv(file=as.character(opt$target_file), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   			   type_convert() %>%
		   			   rename(chrom = X1, pos = X2) %>%
		   			   mutate(uuid = paste0(chrom, ":", pos))
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	x1 = x2 = x3 = x4 = list()
	for (i in 1:length(sample_names)) {
		df = read_tsv(file=paste0("waltz/", sample_names[i], "-STANDARD-pileup.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
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
		   	 						REF = x[1]
		   	 						x[1] = 0
		   	 						x = as.numeric(x)
		   	 						if (REF=="A") {
		   	 							index = c(2:5)[-1]
		   	 							if (any(x[index]>.02*x[6])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(x[index])
		   	 							}
		   	 						} else if (REF=="G") {
										index = c(2:5)[-2]
		   	 							if (any(x[index]>.02*x[6])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(x[index])
		   	 							}
		   	 						} else if (REF=="C") {
		   	 							index = c(2:5)[-3]
		   	 							if (any(x[index]>.02*x[6])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(x[index])
		   	 							}
		   	 						} else if (REF=="T") {
		   	 							index = c(2:5)[-4]
		   	 							if (any(x[index]>.02*x[6])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(x[index])
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
		df = read_tsv(file=paste0("waltz/", sample_names[i], "-COLLAPSED-pileup.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
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
		   	 						REF = x[1]
		   	 						x[1] = 0
		   	 						x = as.numeric(x)
		   	 						if (REF=="A") {
		   	 							index = c(2:5)[-1]
		   	 							if (any(x[index]>.02*x[6])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(x[index])
		   	 							}
		   	 						} else if (REF=="G") {
										index = c(2:5)[-2]
		   	 							if (any(x[index]>.02*x[6])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(x[index])
		   	 							}
		   	 						} else if (REF=="C") {
		   	 							index = c(2:5)[-3]
		   	 							if (any(x[index]>.02*x[6])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(x[index])
		   	 							}
		   	 						} else if (REF=="T") {
		   	 							index = c(2:5)[-4]
		   	 							if (any(x[index]>.02*x[6])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(x[index])
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
		df = read_tsv(file=paste0("waltz/", sample_names[i], "-SIMPLEX-pileup.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
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
		   	 						REF = x[1]
		   	 						x[1] = 0
		   	 						x = as.numeric(x)
		   	 						if (REF=="A") {
		   	 							index = c(2:5)[-1]
		   	 							if (any(x[index]>.02*x[6])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(x[index])
		   	 							}
		   	 						} else if (REF=="G") {
										index = c(2:5)[-2]
		   	 							if (any(x[index]>.02*x[6])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(x[index])
		   	 							}
		   	 						} else if (REF=="C") {
		   	 							index = c(2:5)[-3]
		   	 							if (any(x[index]>.02*x[6])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(x[index])
		   	 							}
		   	 						} else if (REF=="T") {
		   	 							index = c(2:5)[-4]
		   	 							if (any(x[index]>.02*x[6])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(x[index])
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
		df = read_tsv(file=paste0("waltz/", sample_names[i], "-DUPLEX-pileup.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
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
		   	 						REF = x[1]
		   	 						x[1] = 0
		   	 						x = as.numeric(x)
		   	 						if (REF=="A") {
		   	 							index = c(2:5)[-1]
		   	 							if (any(x[index]>.02*x[6])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(x[index])
		   	 							}
		   	 						} else if (REF=="G") {
										index = c(2:5)[-2]
		   	 							if (any(x[index]>.02*x[6])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(x[index])
		   	 							}
		   	 						} else if (REF=="C") {
		   	 							index = c(2:5)[-3]
		   	 							if (any(x[index]>.02*x[6])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(x[index])
		   	 							}
		   	 						} else if (REF=="T") {
		   	 							index = c(2:5)[-4]
		   	 							if (any(x[index]>.02*x[6])) {
		   	 								y = NA
		   	 							} else {
		   	 								y = sum(x[index])
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
	write_tsv(x, path="waltz/noise_metrics.txt", na = "NA", append = FALSE, col_names = !append)
}
