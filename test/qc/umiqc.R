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
	write.table(umis, file="metrics/summary/umi_frequencies.tsv", sep="\t", col.names=TRUE, row.names=TRUE, append=FALSE, quote=FALSE)

} else if (as.numeric(opt$type)==2) {

	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	umi_frequencies = list()
	umi_list = list()
	for (i in 1:length(sample_names)) {
		umi_frequencies[[i]] = read_tsv(file=paste0("marianas/", sample_names[i], "/composite-umi-frequencies.txt"), col_names=FALSE, col_types = cols(.default = col_character())) %>%
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
	write.table(umis, file="metrics/summary/umi_composite.tsv", sep="\t", col.names=TRUE, row.names=TRUE, append=FALSE, quote=FALSE)

} else if (as.numeric(opt$type)==3) {

	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))

	umi_families = list()
	for (i in 1:length(sample_names)) {
		umi_families[[i]] = read_tsv(file=paste0("marianas/", sample_names[i], "/family-sizes.txt"), col_names=TRUE, col_types = cols(.default = col_character())) %>%
		   					type_convert()
	}
	umi_families = do.call(rbind, umi_families)
	write_tsv(x=umi_families, path="metrics/summary/umi_families.tsv", na = "NA", append = FALSE, col_names = TRUE)
	
	family_types_a = list()
	for (i in 1:length(sample_names)) {
		family_types_a[[i]] = read_tsv(file=paste0("marianas/", sample_names[i], "/family-types-A.txt"), col_names=TRUE, col_types = cols(.default = col_character())) %>%
		   					  type_convert()
	}
	family_types_a = do.call(rbind, family_types_a) %>%
					 mutate(BAIT_SET = "MSK-ACCESS_v1.0_probe-A")
	family_types_b = list()
	for (i in 1:length(sample_names)) {
		family_types_b[[i]] = read_tsv(file=paste0("marianas/", sample_names[i], "/family-types-B.txt"), col_names=TRUE, col_types = cols(.default = col_character())) %>%
		   					  type_convert()
	}
	family_types_b = do.call(rbind, family_types_b) %>%
					 mutate(BAIT_SET = "MSK-ACCESS_v1.0_probe-B")
					 
	family_types = bind_rows(family_types_a, family_types_b)
	write_tsv(x=family_types, path="metrics/summary/umi_family_types.tsv", na = "NA", append = FALSE, col_names = TRUE)

}
