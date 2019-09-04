#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("CNaccess"))


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--sample_name", default = NA, type = 'character', help = "sample name"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

data(MSK_ACCESS_A_coord_)
data(MSK_ACCESS_A_counts_)
	
test_ = read_tsv(file=paste0("cnvaccess/cov/", opt$sample_name, ".probe-A.txt"), col_names=TRUE, col_types = cols(.default = col_character())) %>%
		type_convert() %>%
		dplyr::select(chr, start, end, counts)
indx = findmatch_(test_counts=test_$counts, panel_label = "A")
index = index_[, indx] & test_$counts>0
test_$counts[index] = test_$counts[index]/median(test_$counts[index])
t1 = tryspan_(counts = test_$counts[index],
			  gc = coord$gc[index],
			  mappability = coord$mappability[index])
optspan = t1$span[which.min(t1$ssr)]
l1 = loessfit_(counts = test_$counts[index],
			   gc = coord$gc[index],
			   mappability = coord$mappability[index],
			   span = optspan)
counts_[index,indx] = counts_[index,indx]/mean(counts_[index,indx])
t2 = tryspan_(counts = counts_[index,indx],
			  gc = coord$gc[index],
			  mappability = coord$mappability[index])
optspan = t2$span[which.min(t2$ssr)]
l2 = loessfit_(counts = counts_[index,indx],
			   gc = coord$gc[index],
			   mappability = coord$mappability[index],
			   span = optspan)
				   
A_ = test_ %>%
	 mutate(log2 = NA) %>%
	 mutate(panel = "A")
		
x = MASS::lqs(l1~l2)$residuals

A_[index,"log2"] = x

data(MSK_ACCESS_B_coord_)
data(MSK_ACCESS_B_counts_)
	
test_ = read_tsv(file=paste0("cnvaccess/cov/", opt$sample_name, ".probe-B.txt"), col_names=TRUE, col_types = cols(.default = col_character())) %>%
		type_convert() %>%
		dplyr::select(chr, start, end, counts)
indx = findmatch_(test_counts=test_$counts, panel_label = "B")
index = index_[, indx] & test_$counts>0
test_$counts[index] = test_$counts[index]/median(test_$counts[index])
t1 = tryspan_(counts = test_$counts[index],
			  gc = coord$gc[index],
			  mappability = coord$mappability[index])
optspan = t1$span[which.min(t1$ssr)]
l1 = loessfit_(counts = test_$counts[index],
			   gc = coord$gc[index],
			   mappability = coord$mappability[index],
			   span = optspan)
counts_[index,indx] = counts_[index,indx]/mean(counts_[index,indx])
t2 = tryspan_(counts = counts_[index,indx],
			  gc = coord$gc[index],
			  mappability = coord$mappability[index])
optspan = t2$span[which.min(t2$ssr)]
l2 = loessfit_(counts = counts_[index,indx],
			   gc = coord$gc[index],
			   mappability = coord$mappability[index],
			   span = optspan)
				   
B_ = test_ %>%
	 mutate(log2 = NA) %>%
	 mutate(panel = "B")
		
x = MASS::lqs(l1~l2)$residuals

B_[index,"log2"] = x

AB_ = rbind(A_, B_) %>%
	  mutate(chr = ifelse(chr=="X", "23", chr)) %>%
	  mutate(chr = ifelse(chr=="Y", "24", chr)) %>%
	  mutate(chr = as.numeric(chr)) %>%
	  arrange(chr, start)
	  
write_table(AB_, path=paste0("cnvaccess/log2/", opt$sample_name, ".txt"), na = "NA", append = FALSE, col_names = TRUE)
