#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("ExomeDepth"))


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--sample_name", default = NA, type = 'character', help = "sample name"),
				 make_option("--probe", default = NA, type = 'character', help = "probe"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

if (as.character(opt$probe)=="A") {

	bed = read_tsv(file="~/share/reference/cnvkit_reference/MSK-ACCESS-v1_0-probe-A.bigWig", col_types = cols(.default = col_character())) %>%
		  type_convert() %>%
		  dplyr::select(chromosome=chr,
		  				start,
		  				end) %>%
		  mutate(name = paste0(chromosome, ":", start, "-", end))
	bed = as.data.frame(bed)
	counts = getBamCounts(bed.frame = bed,
                          bam.files = paste0("bam/", as.character(opt$sample_name), "-standard.bam"),
                          include.chr = FALSE)
    counts = as.data.frame(counts)
    colnames(counts)[6] = "counts"
    counts = counts %>%
    		 dplyr::select(
    		 	chr = space,
    		 	start = start,
    		 	end = end,
    		 	counts = counts)
    bed = read_tsv(file="~/share/reference/cnvkit_reference/MSK-ACCESS-v1_0-probe-A.bigWig", col_types = cols(.default = col_character())) %>%
		  type_convert() %>%
		  dplyr::select(gc, mappability)
	counts = bind_cols(counts, bed)
    write_tsv(counts, path=paste0("cnvaccess/cov/", as.character(opt$sample_name), ".probe-A.txt"), na = "NA", append = FALSE, col_names = TRUE)

} else if (as.character(opt$probe)=="B") {
	
	bed = read_tsv(file="~/share/reference/cnvkit_reference/MSK-ACCESS-v1_0-probe-B.bigWig", col_types = cols(.default = col_character())) %>%
		  type_convert() %>%
		  dplyr::select(chromosome=chr,
		  				start,
		  				end) %>%
		  mutate(name = paste0(chromosome, ":", start, "-", end))
	bed = as.data.frame(bed)
	counts = getBamCounts(bed.frame = bed,
                          bam.files = paste0("bam/", as.character(opt$sample_name), "-standard.bam"),
                          include.chr = FALSE)
    counts = as.data.frame(counts)
    colnames(counts)[6] = "counts"
    counts = counts %>%
    		 dplyr::select(
    		 	chr = space,
    		 	start = start,
    		 	end = end,
    		 	counts = counts)
    bed = read_tsv(file="~/share/reference/cnvkit_reference/MSK-ACCESS-v1_0-probe-B.bigWig", col_types = cols(.default = col_character())) %>%
		  type_convert() %>%
		  dplyr::select(gc, mappability)
	counts = bind_cols(counts, bed)
    write_tsv(counts, path=paste0("cnvaccess/cov/", as.character(opt$sample_name), ".probe-B.txt"), na = "NA", append = FALSE, col_names = TRUE)

} else if (as.character(opt$probe)=="NA") {
	
	bed = read_tsv(file="~/share/reference/cnvkit_reference/MSK-ACCESS-v1_0-noprobe.bigWig", col_types = cols(.default = col_character())) %>%
		  type_convert() %>%
		  dplyr::select(chromosome=chr,
		  				start,
		  				end) %>%
		  mutate(name = paste0(chromosome, ":", start, "-", end))
	bed = as.data.frame(bed)
	counts = getBamCounts(bed.frame = bed,
                          bam.files = paste0("bam/", as.character(opt$sample_name), "-standard.bam"),
                          include.chr = FALSE)
    counts = as.data.frame(counts)
    colnames(counts)[6] = "counts"
    counts = counts %>%
    		 dplyr::select(
    		 	chr = space,
    		 	start = start,
    		 	end = end,
    		 	counts = counts)
    bed = read_tsv(file="~/share/reference/cnvkit_reference/MSK-ACCESS-v1_0-noprobe.bigWig", col_types = cols(.default = col_character())) %>%
		  type_convert() %>%
		  dplyr::select(gc, mappability)
	counts = bind_cols(counts, bed)
    write_tsv(counts, path=paste0("cnvaccess/cov/", as.character(opt$sample_name), ".probe-AB.txt"), na = "NA", append = FALSE, col_names = TRUE)

}

