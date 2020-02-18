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
chr = 1:21
epsilon = 0.00001

if (as.numeric(opt$type)==1) {

	target_positions = read_tsv(file=as.character(opt$target_file), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   			   readr::type_convert() %>%
		   			   dplyr::rename(chrom = X1, pos = X2) %>%
		   			   dplyr::mutate(uuid = paste0(chrom, ":", pos))
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	x1 = x2 = x3 = x4 = list()
	for (i in 1:length(sample_names)) {
		x = read_tsv(file=paste0("waltz/", sample_names[i], "-pileup.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   	readr::type_convert() %>%
		   	dplyr::select(`chrom`	= X1,
		   				  `pos`	= X2,
		   	 			  `ref`	= X3,
		   	 			  `total`	= X4,
		   	 			  `a`		= X5,
		   	 			  `g`		= X6,
		   	 			  `c`		= X7,
		   	 			  `t`		= X8) %>%
		   	dplyr::mutate(`total_n` = a+g+c+t,
		   	 			  `n` = total - total_n,
		   	 			  `uuid` = paste0(chrom, ":", pos)) %>%
		   	dplyr::filter(uuid %in% target_positions$uuid) %>%
		   	dplyr::mutate(af_a = 100*a/(total_n+epsilon),
						  af_c = 100*c/(total_n+epsilon),
						  af_g = 100*g/(total_n+epsilon),
						  af_t = 100*t/(total_n+epsilon))
		y = dplyr::tibble(`ref` = rep(x$ref, 4),
						  `alt` = rep(c("A", "C", "G", "T"), each=nrow(x)),
						  `freq` = c(x$af_a, x$af_c, x$af_g,x$af_t),
						  `count` = c(x$a, x$c, x$g,x$t),
						  `total_count` = rep(x$total_n, 4),
						  `n_count` =  rep(x$n, 4),
						  `uuid` = rep(x$uuid, 4)) %>%
			dplyr::filter(alt != ref,
						  freq < cutoffAF*100) %>%
			dplyr::group_by(uuid) %>%
			dplyr::summarize(total_count = mean(total_count),
		   	 				 alt_count = sum(count),
		   	 				 ref_count = total_count - alt_count,
		   	 				 n_count = mean(n_count)) %>%
		   	dplyr::summarize(total_counts = sum(total_count),
		   	 				 alt_counts = sum(alt_count),
		   	 				 ref_counts = sum(ref_count),
		   	 				 n_counts = sum(n_count),
		   	 				 contributing_sites = n())
		x1[[i]] = y
	}
	x1 = do.call(rbind, x1) %>%
		 mutate(sample_names = sample_names, bam_file = "standard")
	for (i in 1:length(sample_names)) {
		x = read_tsv(file=paste0("waltz/", sample_names[i], "__aln_srt_IR_FX-pileup.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   	readr::type_convert() %>%
		   	dplyr::select(`chrom`	= X1,
		   				  `pos`	= X2,
		   	 			  `ref`	= X3,
		   	 			  `total`	= X4,
		   	 			  `a`		= X5,
		   	 			  `g`		= X6,
		   	 			  `c`		= X7,
		   	 			  `t`		= X8) %>%
		   	dplyr::mutate(`total_n` = a+g+c+t,
		   	 			  `n` = total - total_n,
		   	 			  `uuid` = paste0(chrom, ":", pos)) %>%
		   	dplyr::filter(uuid %in% target_positions$uuid) %>%
		   	dplyr::mutate(af_a = 100*a/(total_n+epsilon),
						  af_c = 100*c/(total_n+epsilon),
						  af_g = 100*g/(total_n+epsilon),
						  af_t = 100*t/(total_n+epsilon))
		y = dplyr::tibble(`ref` = rep(x$ref, 4),
						  `alt` = rep(c("A", "C", "G", "T"), each=nrow(x)),
						  `freq` = c(x$af_a, x$af_c, x$af_g,x$af_t),
						  `count` = c(x$a, x$c, x$g,x$t),
						  `total_count` = rep(x$total_n, 4),
						  `n_count` =  rep(x$n, 4),
						  `uuid` = rep(x$uuid, 4)) %>%
			dplyr::filter(alt != ref,
						  freq < cutoffAF*100) %>%
			dplyr::group_by(uuid) %>%
			dplyr::summarize(total_count = mean(total_count),
		   	 				 alt_count = sum(count),
		   	 				 ref_count = total_count - alt_count,
		   	 				 n_count = mean(n_count)) %>%
		   	dplyr::summarize(total_counts = sum(total_count),
		   	 				 alt_counts = sum(alt_count),
		   	 				 ref_counts = sum(ref_count),
		   	 				 n_counts = sum(n_count),
		   	 				 contributing_sites = n())
		   	 x2[[i]] = y
	}
	x2 = do.call(rbind, x2) %>%
		 mutate(sample_names = sample_names, bam_file = "collapsed")
	for (i in 1:length(sample_names)) {
		x = read_tsv(file=paste0("waltz/", sample_names[i], "__aln_srt_IR_FX-simplex-pileup.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   	readr::type_convert() %>%
		   	dplyr::select(`chrom`	= X1,
		   				  `pos`	= X2,
		   	 			  `ref`	= X3,
		   	 			  `total`	= X4,
		   	 			  `a`		= X5,
		   	 			  `g`		= X6,
		   	 			  `c`		= X7,
		   	 			  `t`		= X8) %>%
		   	dplyr::mutate(`total_n` = a+g+c+t,
		   	 			  `n` = total - total_n,
		   	 			  `uuid` = paste0(chrom, ":", pos)) %>%
		   	dplyr::filter(uuid %in% target_positions$uuid) %>%
		   	dplyr::mutate(af_a = 100*a/(total_n+epsilon),
						  af_c = 100*c/(total_n+epsilon),
						  af_g = 100*g/(total_n+epsilon),
						  af_t = 100*t/(total_n+epsilon))
		y = dplyr::tibble(`ref` = rep(x$ref, 4),
						  `alt` = rep(c("A", "C", "G", "T"), each=nrow(x)),
						  `freq` = c(x$af_a, x$af_c, x$af_g,x$af_t),
						  `count` = c(x$a, x$c, x$g,x$t),
						  `total_count` = rep(x$total_n, 4),
						  `n_count` =  rep(x$n, 4),
						  `uuid` = rep(x$uuid, 4)) %>%
			dplyr::filter(alt != ref,
						  freq < cutoffAF*100) %>%
			dplyr::group_by(uuid) %>%
			dplyr::summarize(total_count = mean(total_count),
		   	 				 alt_count = sum(count),
		   	 				 ref_count = total_count - alt_count,
		   	 				 n_count = mean(n_count)) %>%
		   	dplyr::summarize(total_counts = sum(total_count),
		   	 				 alt_counts = sum(alt_count),
		   	 				 ref_counts = sum(ref_count),
		   	 				 n_counts = sum(n_count),
		   	 				 contributing_sites = n())
		   	 x3[[i]] = y
	}
	x3 = do.call(rbind, x3) %>%
		 mutate(sample_names = sample_names, bam_file = "simplex")
	for (i in 1:length(sample_names)) {
		x = read_tsv(file=paste0("waltz/", sample_names[i], "__aln_srt_IR_FX-duplex-pileup.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   	readr::type_convert() %>%
		   	dplyr::select(`chrom`	= X1,
		   				  `pos`	= X2,
		   	 			  `ref`	= X3,
		   	 			  `total`	= X4,
		   	 			  `a`		= X5,
		   	 			  `g`		= X6,
		   	 			  `c`		= X7,
		   	 			  `t`		= X8) %>%
		   	dplyr::mutate(`total_n` = a+g+c+t,
		   	 			  `n` = total - total_n,
		   	 			  `uuid` = paste0(chrom, ":", pos)) %>%
		   	dplyr::filter(uuid %in% target_positions$uuid) %>%
		   	dplyr::mutate(af_a = 100*a/(total_n+epsilon),
						  af_c = 100*c/(total_n+epsilon),
						  af_g = 100*g/(total_n+epsilon),
						  af_t = 100*t/(total_n+epsilon))
		y = dplyr::tibble(`ref` = rep(x$ref, 4),
						  `alt` = rep(c("A", "C", "G", "T"), each=nrow(x)),
						  `freq` = c(x$af_a, x$af_c, x$af_g,x$af_t),
						  `count` = c(x$a, x$c, x$g,x$t),
						  `total_count` = rep(x$total_n, 4),
						  `n_count` =  rep(x$n, 4),
						  `uuid` = rep(x$uuid, 4)) %>%
			dplyr::filter(alt != ref,
						  freq < cutoffAF*100) %>%
			dplyr::group_by(uuid) %>%
			dplyr::summarize(total_count = mean(total_count),
		   	 				 alt_count = sum(count),
		   	 				 ref_count = total_count - alt_count,
		   	 				 n_count = mean(n_count)) %>%
		   	dplyr::summarize(total_counts = sum(total_count),
		   	 				 alt_counts = sum(alt_count),
		   	 				 ref_counts = sum(ref_count),
		   	 				 n_counts = sum(n_count),
		   	 				 contributing_sites = n())
		   	 x4[[i]] = y
	}
	x4 = do.call(rbind, x4) %>%
		 mutate(sample_names = sample_names, bam_file = "duplex")
	x = rbind(x1, x2, x3, x4)
	write_tsv(x, path="waltz/noise_metrics_with_duplicates.txt", na = "NA", append = FALSE, col_names = TRUE)
	
} else if (as.numeric(opt$type)==2) {

	target_positions = read_tsv(file=as.character(opt$target_file), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   			   readr::type_convert() %>%
		   			   dplyr::rename(chrom = X1, pos = X2) %>%
		   			   dplyr::mutate(uuid = paste0(chrom, ":", pos))
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	x1 = x2 = x3 = x4 = list()
	for (i in 1:length(sample_names)) {
		x = read_tsv(file=paste0("waltz/", sample_names[i], "-pileup-without-duplicates.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   	readr::type_convert() %>%
		   	dplyr::select(`chrom`	= X1,
		   				  `pos`	= X2,
		   	 			  `ref`	= X3,
		   	 			  `total`	= X4,
		   	 			  `a`		= X5,
		   	 			  `g`		= X6,
		   	 			  `c`		= X7,
		   	 			  `t`		= X8) %>%
		   	dplyr::mutate(`total_n` = a+g+c+t,
		   	 			  `n` = total - total_n,
		   	 			  `uuid` = paste0(chrom, ":", pos)) %>%
		   	dplyr::filter(uuid %in% target_positions$uuid) %>%
		   	dplyr::mutate(af_a = 100*a/(total_n+epsilon),
						  af_c = 100*c/(total_n+epsilon),
						  af_g = 100*g/(total_n+epsilon),
						  af_t = 100*t/(total_n+epsilon))
		y = dplyr::tibble(`ref` = rep(x$ref, 4),
						  `alt` = rep(c("A", "C", "G", "T"), each=nrow(x)),
						  `freq` = c(x$af_a, x$af_c, x$af_g,x$af_t),
						  `count` = c(x$a, x$c, x$g,x$t),
						  `total_count` = rep(x$total_n, 4),
						  `n_count` =  rep(x$n, 4),
						  `uuid` = rep(x$uuid, 4)) %>%
			dplyr::filter(alt != ref,
						  freq < cutoffAF*100) %>%
			dplyr::group_by(uuid) %>%
			dplyr::summarize(total_count = mean(total_count),
		   	 				 alt_count = sum(count),
		   	 				 ref_count = total_count - alt_count,
		   	 				 n_count = mean(n_count)) %>%
		   	dplyr::summarize(total_counts = sum(total_count),
		   	 				 alt_counts = sum(alt_count),
		   	 				 ref_counts = sum(ref_count),
		   	 				 n_counts = sum(n_count),
		   	 				 contributing_sites = n())
		   	 x1[[i]] = y
	}
	x1 = do.call(rbind, x1) %>%
		 mutate(sample_names = sample_names, bam_file = "standard")
	for (i in 1:length(sample_names)) {
		x = read_tsv(file=paste0("waltz/", sample_names[i], "__aln_srt_IR_FX-pileup-without-duplicates.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   	readr::type_convert() %>%
		   	dplyr::select(`chrom`	= X1,
		   				  `pos`	= X2,
		   	 			  `ref`	= X3,
		   	 			  `total`	= X4,
		   	 			  `a`		= X5,
		   	 			  `g`		= X6,
		   	 			  `c`		= X7,
		   	 			  `t`		= X8) %>%
		   	dplyr::mutate(`total_n` = a+g+c+t,
		   	 			  `n` = total - total_n,
		   	 			  `uuid` = paste0(chrom, ":", pos)) %>%
		   	dplyr::filter(uuid %in% target_positions$uuid) %>%
		   	dplyr::mutate(af_a = 100*a/(total_n+epsilon),
						  af_c = 100*c/(total_n+epsilon),
						  af_g = 100*g/(total_n+epsilon),
						  af_t = 100*t/(total_n+epsilon))
		y = dplyr::tibble(`ref` = rep(x$ref, 4),
						  `alt` = rep(c("A", "C", "G", "T"), each=nrow(x)),
						  `freq` = c(x$af_a, x$af_c, x$af_g,x$af_t),
						  `count` = c(x$a, x$c, x$g,x$t),
						  `total_count` = rep(x$total_n, 4),
						  `n_count` =  rep(x$n, 4),
						  `uuid` = rep(x$uuid, 4)) %>%
			dplyr::filter(alt != ref,
						  freq < cutoffAF*100) %>%
			dplyr::group_by(uuid) %>%
			dplyr::summarize(total_count = mean(total_count),
		   	 				 alt_count = sum(count),
		   	 				 ref_count = total_count - alt_count,
		   	 				 n_count = mean(n_count)) %>%
		   	dplyr::summarize(total_counts = sum(total_count),
		   	 				 alt_counts = sum(alt_count),
		   	 				 ref_counts = sum(ref_count),
		   	 				 n_counts = sum(n_count),
		   	 				 contributing_sites = n())
		   	 x2[[i]] = y
	}
	x2 = do.call(rbind, x2) %>%
		 mutate(sample_names = sample_names, bam_file = "collapsed")
	for (i in 1:length(sample_names)) {
		x = read_tsv(file=paste0("waltz/", sample_names[i], "__aln_srt_IR_FX-simplex-pileup-without-duplicates.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   	readr::type_convert() %>%
		   	dplyr::select(`chrom`	= X1,
		   				  `pos`	= X2,
		   	 			  `ref`	= X3,
		   	 			  `total`	= X4,
		   	 			  `a`		= X5,
		   	 			  `g`		= X6,
		   	 			  `c`		= X7,
		   	 			  `t`		= X8) %>%
		   	dplyr::mutate(`total_n` = a+g+c+t,
		   	 			  `n` = total - total_n,
		   	 			  `uuid` = paste0(chrom, ":", pos)) %>%
		   	dplyr::filter(uuid %in% target_positions$uuid) %>%
		   	dplyr::mutate(af_a = 100*a/(total_n+epsilon),
						  af_c = 100*c/(total_n+epsilon),
						  af_g = 100*g/(total_n+epsilon),
						  af_t = 100*t/(total_n+epsilon))
		y = dplyr::tibble(`ref` = rep(x$ref, 4),
						  `alt` = rep(c("A", "C", "G", "T"), each=nrow(x)),
						  `freq` = c(x$af_a, x$af_c, x$af_g,x$af_t),
						  `count` = c(x$a, x$c, x$g,x$t),
						  `total_count` = rep(x$total_n, 4),
						  `n_count` =  rep(x$n, 4),
						  `uuid` = rep(x$uuid, 4)) %>%
			dplyr::filter(alt != ref,
						  freq < cutoffAF*100) %>%
			dplyr::group_by(uuid) %>%
			dplyr::summarize(total_count = mean(total_count),
		   	 				 alt_count = sum(count),
		   	 				 ref_count = total_count - alt_count,
		   	 				 n_count = mean(n_count)) %>%
		   	dplyr::summarize(total_counts = sum(total_count),
		   	 				 alt_counts = sum(alt_count),
		   	 				 ref_counts = sum(ref_count),
		   	 				 n_counts = sum(n_count),
		   	 				 contributing_sites = n())
		   	 x3[[i]] = y
	}
	x3 = do.call(rbind, x3) %>%
		 mutate(sample_names = sample_names, bam_file = "simplex")
	for (i in 1:length(sample_names)) {
		x = read_tsv(file=paste0("waltz/", sample_names[i], "__aln_srt_IR_FX-duplex-pileup-without-duplicates.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   	readr::type_convert() %>%
		   	dplyr::select(`chrom`	= X1,
		   				  `pos`	= X2,
		   	 			  `ref`	= X3,
		   	 			  `total`	= X4,
		   	 			  `a`		= X5,
		   	 			  `g`		= X6,
		   	 			  `c`		= X7,
		   	 			  `t`		= X8) %>%
		   	dplyr::mutate(`total_n` = a+g+c+t,
		   	 			  `n` = total - total_n,
		   	 			  `uuid` = paste0(chrom, ":", pos)) %>%
		   	dplyr::filter(uuid %in% target_positions$uuid) %>%
		   	dplyr::mutate(af_a = 100*a/(total_n+epsilon),
						  af_c = 100*c/(total_n+epsilon),
						  af_g = 100*g/(total_n+epsilon),
						  af_t = 100*t/(total_n+epsilon))
		y = dplyr::tibble(`ref` = rep(x$ref, 4),
						  `alt` = rep(c("A", "C", "G", "T"), each=nrow(x)),
						  `freq` = c(x$af_a, x$af_c, x$af_g,x$af_t),
						  `count` = c(x$a, x$c, x$g,x$t),
						  `total_count` = rep(x$total_n, 4),
						  `n_count` =  rep(x$n, 4),
						  `uuid` = rep(x$uuid, 4)) %>%
			dplyr::filter(alt != ref,
						  freq < cutoffAF*100) %>%
			dplyr::group_by(uuid) %>%
			dplyr::summarize(total_count = mean(total_count),
		   	 				 alt_count = sum(count),
		   	 				 ref_count = total_count - alt_count,
		   	 				 n_count = mean(n_count)) %>%
		   	dplyr::summarize(total_counts = sum(total_count),
		   	 				 alt_counts = sum(alt_count),
		   	 				 ref_counts = sum(ref_count),
		   	 				 n_counts = sum(n_count),
		   	 				 contributing_sites = n())
		   	 x4[[i]] = y
	}
	x4 = do.call(rbind, x4) %>%
		 mutate(sample_names = sample_names, bam_file = "duplex")
	x = rbind(x1, x2, x3, x4)
	write_tsv(x, path="waltz/noise_metrics_without_duplicates.txt", na = "NA", append = FALSE, col_names = TRUE)
	
} else if (as.numeric(opt$type)==3) {

	target_positions = read_tsv(file=as.character(opt$target_file), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   			   readr::type_convert() %>%
		   			   dplyr::rename(chrom = X1, pos = X2) %>%
		   			   dplyr::mutate(uuid = paste0(chrom, ":", pos))
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	x = list()
	for (i in 1:length(sample_names)) {
		y = read_tsv(file=paste0("waltz/", sample_names[i], "-pileup.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
			readr::type_convert() %>%
			dplyr::select(`chrom`	= X1,
						  `pos`		= X2,
						  `ref`		= X3,
						  `total`	= X4,
						  `a`		= X5,
						  `c`		= X6,
						  `g`		= X7,
						  `t`		= X8) %>%
			dplyr::mutate(`total_n` = a+c+g+t,
			 			  `uuid` = paste0(chrom, ":", pos)) %>%
			dplyr::filter(uuid %in% target_positions$uuid) %>%
			dplyr::mutate(af_a = 100*a/(total_n+epsilon),
						  af_c = 100*c/(total_n+epsilon),
						  af_g = 100*g/(total_n+epsilon),
						  af_t = 100*t/(total_n+epsilon))
		z = dplyr::tibble(`chrom` = rep(df$chrom, 4),
						  `pos` = rep(df$pos, 4),
						  `ref` = rep(df$ref, 4),
						  `alt` = rep(c("A", "C", "G", "T"), each=nrow(df)),
						  `freq` = c(df$af_a, df$af_c, df$af_g,df$af_t)) %>%
			dplyr::filter(chrom %in% chr,
						  alt != ref)
			dplyr::mutate(freq = ifelse(freq>(cutoffAF*100), 0, freq)) %>%
			dplyr::arrange(chrom, pos, ref, alt) %>%
			dplyr::rename_at("freq", funs(sample_names[i]))
		if (i!=1) {
			z = z %>% dplyr::select(-chrom, -pos, -ref, -alt)
		}
		x[[i]] = z
	}
	x = do.call(cbind, x)
	write_tsv(x, path="waltz/noise_by_position_standard_with_duplicates.txt", na = "NA", append = FALSE, col_names = TRUE)
	
} else if (as.numeric(opt$type)==4) {

	target_positions = read_tsv(file=as.character(opt$target_file), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   			   readr::type_convert() %>%
		   			   dplyr::rename(chrom = X1, pos = X2) %>%
		   			   dplyr::mutate(uuid = paste0(chrom, ":", pos))
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	x = list()
	for (i in 1:length(sample_names)) {
		y = read_tsv(file=paste0("waltz/", sample_names[i], "-pileup-without-duplicates.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
			 readr::type_convert() %>%
			 dplyr::select(`chrom`	= X1,
						   `pos`	= X2,
						   `ref`	= X3,
						   `total`	= X4,
						   `a`		= X5,
						   `c`		= X6,
						   `g`		= X7,
						   `t`		= X8) %>%
			 dplyr::mutate(`total_n` = a+c+g+t,
			 			   `uuid` = paste0(chrom, ":", pos)) %>%
			 dplyr::filter(uuid %in% target_positions$uuid) %>%
			 dplyr::mutate(af_a = 100*a/(total_n+epsilon),
						   af_c = 100*c/(total_n+epsilon),
						   af_g = 100*g/(total_n+epsilon),
						   af_t = 100*t/(total_n+epsilon))
		x[[i]] = dplyr::tibble(`chrom` = rep(df$chrom, 4),
							   `pos` = rep(df$pos, 4),
							   `ref` = rep(df$ref, 4),
							   `alt` = rep(c("A", "C", "G", "T"), each=nrow(df)),
							   `af` = c(df$af_a, df$af_c, df$af_g,df$af_t)) %>%
				 dplyr::mutate(af = ifelse(af>cutoffAF*100, 0, af)) %>%
				 dplyr::filter(chrom %in% chr,
				 			   ref != alt) %>%
				 dplyr::arrange(chrom, pos, ref, alt) %>%
				 dplyr::rename_at("af", funs(sample_names[i]))
		if (i!=1) {
			x[[i]] = x[[i]] %>%
					 dplyr::select(-chrom, -pos, -ref, -alt)
		}
	}
	x = do.call(cbind, x)
	write_tsv(x, path="waltz/noise_by_position_standard_without_duplicates.txt", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$type)==5) {
	target_positions = read_tsv(file=as.character(opt$target_file), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   			   readr::type_convert() %>%
		   			   dplyr::rename(chrom = X1, pos = X2) %>%
		   			   dplyr::mutate(uuid = paste0(chrom, ":", pos))
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	x = list()
	for (i in 1:length(sample_names)) {
		y = read_tsv(file=paste0("waltz/", sample_names[i], "__aln_srt_IR_FX-simplex-pileup-without-duplicates.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
			readr::type_convert() %>%
			dplyr::select(`chrom`	= X1,
						  `pos`		= X2,
						  `ref`		= X3,
						  `total`	= X4,
						  `a`		= X5,
						  `c`		= X6,
						  `g`		= X7,
						  `t`		= X8) %>%
			dplyr::mutate(`total_n` = a+c+g+t,
			 			  `uuid` = paste0(chrom, ":", pos)) %>%
			dplyr::filter(uuid %in% target_positions$uuid) %>%
			dplyr::mutate(af_a = 100*a/(total_n+epsilon),
						  af_c = 100*c/(total_n+epsilon),
						  af_g = 100*g/(total_n+epsilon),
						  af_t = 100*t/(total_n+epsilon))
		z = dplyr::tibble(`chrom` = rep(df$chrom, 4),
						  `pos` = rep(df$pos, 4),
						  `ref` = rep(df$ref, 4),
						  `alt` = rep(c("A", "C", "G", "T"), each=nrow(df)),
						  `freq` = c(df$af_a, df$af_c, df$af_g,df$af_t)) %>%
			dplyr::filter(chrom %in% chr,
						  alt != ref)
			dplyr::mutate(freq = ifelse(freq>(cutoffAF*100), 0, freq)) %>%
			dplyr::arrange(chrom, pos, ref, alt) %>%
			dplyr::rename_at("freq", funs(sample_names[i]))
		if (i!=1) {
			z = z %>% dplyr::select(-chrom, -pos, -ref, -alt)
		}
		x[[i]] = z
	}
	x = do.call(cbind, x)
	write_tsv(x, path="waltz/noise_by_position_simplex_without_duplicates.txt", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$type)==6) {

	target_positions = read_tsv(file=as.character(opt$target_file), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   			   readr::type_convert() %>%
		   			   dplyr::rename(chrom = X1, pos = X2) %>%
		   			   dplyr::mutate(uuid = paste0(chrom, ":", pos))
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	x = list()
	for (i in 1:length(sample_names)) {
		y = read_tsv(file=paste0("waltz/", sample_names[i], "__aln_srt_IR_FX-duplex-pileup-without-duplicates.txt.gz"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
			readr::type_convert() %>%
			dplyr::select(`chrom`	= X1,
						  `pos`		= X2,
						  `ref`		= X3,
						  `total`	= X4,
						  `a`		= X5,
						  `c`		= X6,
						  `g`		= X7,
						  `t`		= X8) %>%
			dplyr::mutate(`total_n` = a+c+g+t,
			 			  `uuid` = paste0(chrom, ":", pos)) %>%
			dplyr::filter(uuid %in% target_positions$uuid) %>%
			dplyr::mutate(af_a = 100*a/(total_n+epsilon),
						  af_c = 100*c/(total_n+epsilon),
						  af_g = 100*g/(total_n+epsilon),
						  af_t = 100*t/(total_n+epsilon))
		z = dplyr::tibble(`chrom` = rep(df$chrom, 4),
						  `pos` = rep(df$pos, 4),
						  `ref` = rep(df$ref, 4),
						  `alt` = rep(c("A", "C", "G", "T"), each=nrow(df)),
						  `freq` = c(df$af_a, df$af_c, df$af_g,df$af_t)) %>%
			dplyr::filter(chrom %in% chr,
						  alt != ref)
			dplyr::mutate(freq = ifelse(freq>(cutoffAF*100), 0, freq)) %>%
			dplyr::arrange(chrom, pos, ref, alt) %>%
			dplyr::rename_at("freq", funs(sample_names[i]))
		if (i!=1) {
			z = z %>% dplyr::select(-chrom, -pos, -ref, -alt)
		}
		x[[i]] = z
	}
	x = do.call(cbind, x)
	write_tsv(x, path="waltz/noise_by_position_duplex_without_duplicates.txt", na = "NA", append = FALSE, col_names = TRUE)

} else if (as.numeric(opt$type)==7) {

	suppressPackageStartupMessages(library("superheat"))
	suppressPackageStartupMessages(library("viridis"))
	
	target_positions = read_tsv(file=as.character(opt$target_file), col_names = FALSE, col_types = cols(.default = col_character())) %>%
		   			   readr::type_convert() %>%
		   			   dplyr::rename(chrom = X1, pos = X2) %>%
		   			   dplyr::mutate(uuid = paste0(chrom, ":", pos))
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	standard_w_dup = read_tsv(file="waltz/noise_by_position_standard_with_duplicates.txt", col_names = TRUE, col_types = cols(.default = col_character())) %>%
					 readr::type_convert()
	standard_w0_dup = read_tsv(file="waltz/noise_by_position_standard_without_duplicates.txt", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		  			  readr::type_convert()
	simplex_w0_dup = read_tsv(file="waltz/noise_by_position_simplex_without_duplicates.txt", col_names = TRUE, col_types = cols(.default = col_character())) %>%
					 readr::type_convert()
	duplex_w0_dup = read_tsv(file="waltz/noise_by_position_duplex_without_duplicates.txt", col_names = TRUE, col_types = cols(.default = col_character())) %>%
					readr::type_convert()

	nc_pu = dplyr::full_join(standard_w_dup, standard_w0_dup, by = c("chrom", "pos", "ref", "alt")) %>%
			dplyr::full_join(simplex_w0_dup, by = c("chrom", "pos", "ref", "alt")) %>%
			dplyr::full_join(duplex_w0_dup, by = c("chrom", "pos", "ref", "alt"))
	index = order(apply(nc_pu[,5:ncol(nc_pu),drop=FALSE], 1, mean, na.rm=TRUE), decreasing=FALSE)
	nc_pu = nc_pu[index,,drop=FALSE] %>%
			dplyr::arrange(ref, alt)
	col_groups = factor(rep(c("STANDARD\nWITH DUPLICATES", "STANDARD\nDEDUPLICATED", "COLLAPSED\nSIMPLEX", "COLLAPSED\nDUPLEX"), each=length(sample_names)),
						levels = c("STANDARD\nWITH DUPLICATES", "STANDARD\nDEDUPLICATED", "COLLAPSED\nSIMPLEX", "COLLAPSED\nDUPLEX"), ordered = TRUE)
	row_groups = paste0(nc_pu$ref, " > ", nc_pu$alt, "         ")
	nc_pu = nc_pu %>%
			dplyr::select(-chrom, -pos, -ref, -alt)
	index = apply(nc_pu, 1, function(x) {sum(is.na(x), na.rm=TRUE)})==0 |
			apply(nc_pu, 1, function(x) {sum(x==0, na.rm=TRUE)})>1
	pdf(file="waltz/noise_by_position.pdf", height=14, width=14)
	superheat(X = as.matrix(nc_pu[index,,drop=FALSE]),
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
			  heat.pal = viridis(n=15),
			  heat.pal.values = c(seq(0,.4,l=9), 1))
	dev.off()

}
