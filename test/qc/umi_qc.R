#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--type", default = NA, type = 'character', help = "analysis type"),
				 make_option("--sample_name", default = NA, type = 'character', help = "sample name"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

if (as.numeric(opt$type)==0) {

	umi_frequencies = read_tsv(file=paste0("marianas/", opt$sample_name, "/umi-frequencies.txt"), col_names=FALSE, col_types = cols(.default = col_character())) %>%
		   			  type_convert()
	umi_w_n = umi_frequencies %>%
			  filter(grepl("N", X1))
	umi_w0_n = umi_frequencies %>%
			   filter(!grepl("N", X1))
	
	# % UMI with 'N'
	res_0 = 100 * sum(umi_w_n$X2) / sum(umi_w0_n$X2)
	
	# % UMI with 'N' at 1st or 2nd position
	umi_0n = do.call(rbind, strsplit(umi_w0_n$X1, split="", fixed=TRUE))
	res_1 = list()
	for (i in 1:nrow(umi_w_n)) {
		umi_1n = unlist(strsplit(umi_w_n$X1[i], split="", fixed=TRUE))
		index = which(umi_1n=="N")
		if (index==1) {
			inx = umi_0n[,2]==umi_1n[2] & umi_0n[,3]==umi_1n[3]
		} else if (index==2) {
			inx = umi_0n[,1]==umi_1n[1] & umi_0n[,3]==umi_1n[3]
		}
		res_1[[i]] = 100 * umi_w_n$X2[i] / sum(umi_w0_n$X2[inx])
	}
	res_1 = unlist(res_1)
	names(res_1) = umi_w_n$X1
	
	umi_composite = read_tsv(file=paste0("marianas/", opt$sample_name, "/composite-umi-frequencies.txt"), col_names=FALSE, col_types = cols(.default = col_character())) %>%
		   			type_convert()
		   			
	# % UMI with x 'N'
	index_n = list()
	for (i in 1:6) {
		index_n[[i]] = apply(do.call(rbind, strsplit(umi_composite$X1, "")), 1, function(x) {sum(x=="N")==i})
	}	
	index_n = do.call(cbind, index_n)
	apply(index_n, 2, sum)
	res_2 = 100 * sum(umi_composite$X2[index_n[,1]]) / sum(umi_composite$X2[apply(index_n, 1, sum)==0])
	res_3 = 100 * sum(umi_composite$X2[index_n[,2]]) / sum(umi_composite$X2[apply(index_n, 1, sum)==0])
	res_4 = 100 * sum(umi_composite$X2[index_n[,3]]) / sum(umi_composite$X2[apply(index_n, 1, sum)==0])
	res_5 = 100 * sum(umi_composite$X2[index_n[,4]]) / sum(umi_composite$X2[apply(index_n, 1, sum)==0])
	res_6 = 100 * sum(umi_composite$X2[index_n[,5]]) / sum(umi_composite$X2[apply(index_n, 1, sum)==0])
	res_7 = 100 * sum(umi_composite$X2[index_n[,6]]) / sum(umi_composite$X2[apply(index_n, 1, sum)==0])
	
	# Matrix of UMI combinations
	umi_composite = umi_composite %>%
					filter(!grepl("N", X1))
	umi_1 = unlist(lapply(strsplit(umi_composite$X1, split="+", fixed=TRUE), function(x) {x[1]}))
	umi_2 = unlist(lapply(strsplit(umi_composite$X1, split="+", fixed=TRUE), function(x) {x[2]}))
	u_umi_1 = sort(unique(umi_1))
	u_umi_2 = sort(unique(umi_2))
	res_8 = matrix(0, nrow=length(u_umi_1), ncol=length(u_umi_2))
	for (i in 1:length(u_umi_1)) {
		for (j in 1:length(u_umi_2)) {
			index = which(umi_1 == u_umi_1[i] & umi_2 == u_umi_2[j])
			res_8[i,j] = umi_composite$X2[index]
		}
	}
	rownames(res_8) = u_umi_1
	colnames(res_8) = u_umi_2
	
	
	# % Reads with UMI
	umi_info = read_tsv(file=paste0("marianas/", opt$sample_name, "/info.txt"), col_names=FALSE, col_types = cols(.default = col_character())) %>%
		   	   type_convert()
	x = as.numeric(gsub("Total read pairs: ", "", x=as.character(umi_info[1,1]), fixed=TRUE))
	y = as.numeric(strsplit(gsub("Read pairs with UMIs: ", "", x=as.character(umi_info[2,1]), fixed=TRUE), "/", fixed=TRUE)[[1]][1])
	
	
	
	# Save results
	res = list(percent_read_umi = c("Total N read pairs" = x,
									"N read pairs with UMI" = y),
			   percent_umi_any_n = res_0,
			   percent_umi_adj = res_1,
			   percent_umi_x_n = c("1" = res_2,
			   					   "2" = res_3,
			   					   "3" = res_4,
			   					   "4" = res_5,
			   					   "5" = res_6,
			   					   "6" = res_7),
			   umi_combination = res_8)
	save(res, file=paste0("marianas/", opt$sample_name, "/umi-info.RData"))
	
} else if (as.numeric(opt$type)==1) {
	
	load(paste0("marianas/", opt$sample_name, "/umi-info.RData"))

	# Heatmap of UMI combinations
	pdf(paste0("marianas/", opt$sample_name, "/umi-composite.pdf"))
	heatmap(x=res$umi_combination, scale="none", symm = TRUE, col=colorRampPalette(rev(c("white", "#deebf7", "#9ecae1", "#3182bd")))(12))
	dev.off()

}
