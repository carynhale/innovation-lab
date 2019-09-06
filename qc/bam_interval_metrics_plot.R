#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("scales"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--type", default = NA, type = 'character', help = "analysis type"),
				 make_option("--sample_names", default = NA, type = 'character', help = "sample names"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

AF = 5
CHR = "21"

if (as.numeric(opt$type)==1) {

	data = read_tsv(file="metrics/summary/metrics_hs.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   arrange(desc(MEAN_TARGET_COVERAGE)) %>%
		   mutate(`Sample ID` = factor(SAMPLE, levels=unique(SAMPLE), ordered=TRUE))
		   		   
	pdf(file="metrics/report/target_coverage.pdf", width=14)
	plot.0 = ggplot(data, aes(x=`Sample ID`, y=MEAN_TARGET_COVERAGE)) +
			 geom_bar(stat="identity", fill="#d7191c") +
			 theme_classic(base_size=15) +
			 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
			 labs(x="Sample ID", y="Depth\n", title="MEAN DEDUPLICATED COVERAGE") +
			 theme(plot.title = element_text(hjust = 0.5, size=16))
	print(plot.0)
	dev.off()
	
} else if (as.numeric(opt$type)==2) {

	data = read_tsv(file="metrics/summary/metrics_idx.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert()
		   
	x = data %>%
		dplyr::select(N = N_ALIGNED, SAMPLE) %>%
		mutate(TYPE = "Aligned")
	y = data %>%
		dplyr::select(N = N_UNALIGNED, SAMPLE) %>%
		mutate(TYPE = "Unaligned")
	data = bind_rows(x, y) %>%
		   arrange(desc(N)) %>%
		   dplyr::rename(`Sample ID` = SAMPLE, Type = TYPE) %>%
		   mutate(`Sample ID` = factor(`Sample ID`, levels=unique(`Sample ID`), ordered=TRUE))
		   
	pdf(file="metrics/report/alignment_summary.pdf", width=14)
	plot.0 = ggplot(data, aes(x=`Sample ID`, y=N, fill=Type)) +
			 geom_bar(stat="identity") +
			 scale_fill_manual(values=c("Aligned"="#d7191c", "Unaligned"="#2c7bb6")) +
		 	 theme_classic(base_size=15) +
		 	 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
		 	 labs(x="Sample ID", y="Number of read pairs\n", title="READ PAIR ALIGNMENT") +
		 	 theme(plot.title = element_text(hjust = 0.5, size=16))
	print(plot.0)
	dev.off()

} else if (as.numeric(opt$type)==3) {

	data = read_tsv(file="metrics/summary/metrics_insert.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   arrange(desc(MEDIAN_INSERT_SIZE)) %>%
		   dplyr::rename(`Sample ID` = SAMPLE) %>%
		   mutate(`Sample ID` = factor(`Sample ID`, levels=unique(`Sample ID`), ordered=TRUE))
		   
	pdf(file="metrics/report/insert_size_summary.pdf", width=14)
	plot.0 = ggplot(data, aes(x=`Sample ID`, y=MEDIAN_INSERT_SIZE)) +
			 geom_bar(stat="identity", fill="#d7191c") +
		 	 theme_classic(base_size=15) +
		 	 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
		 	 labs(x="Sample ID", y="Size (bp)\n", title="MEAN INSERT SIZE") +
		 	 theme(plot.title = element_text(hjust = 0.5, size=16)) +
		 	 ylim(0,250)
	print(plot.0)
	dev.off()

} else if (as.numeric(opt$type)==4) {

	data = read_tsv(file="metrics/summary/metrics_insert_distribution.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert()
		   
	pdf(file="metrics/report/insert_size_distribution.pdf", width=14)
	tmp = list()
	for (j in 1:ncol(data)) {
		x = 1:nrow(data)
		y = data[,j,drop=TRUE]
		x = x[!is.na(y)]
		y = y[!is.na(y)]
		z = smooth.spline(x=x, y=y, spar=0.3)
		tmp[[j]] = data.frame(x=z$x, y=z$y) %>%
			       mutate(`Sample ID` = colnames(data)[j])
	}
	tmp = do.call(rbind, tmp)
	plot.0 = ggplot(tmp, aes(x=x, y=y, color=`Sample ID`)) +
			 geom_point(size=.5) +
			 geom_line() +
			 theme_classic(base_size=15) +
			 labs(x="\nInsert size (bp)", y="Frequency\n", title="INSERT SIZE DISTRIBUTION") +
			 theme(plot.title = element_text(hjust = 0.5, size=16)) +
			 xlim(0, 400)
	print(plot.0)
	dev.off()
	
} else if (as.numeric(opt$type)==5) {

	suppressPackageStartupMessages(library("superheat"))
	suppressPackageStartupMessages(library("viridis"))
	
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	
	nuc_metrics = list()
	for (i in 1:length(sample_names)) {
		pileup_metrics = read_tsv(file=paste0("metrics/pileup/", sample_names[i], "-pileup.txt"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
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
	standard_bam = nuc_metrics[[1]][,c("Chromosome", "Position", "Reference_Allele", "Alternate_Allele"),drop=FALSE]
	for (i in 1:length(sample_names)) {
		cat(i, "\n")
		standard_bam = left_join(standard_bam, nuc_metrics[[i]], by=c("Chromosome", "Position", "Reference_Allele", "Alternate_Allele"))
	}
	colnames(standard_bam)[5:ncol(standard_bam)] = sample_names

	nuc_metrics = list()
	for (i in 1:length(sample_names)) {
		pileup_metrics = read_tsv(file=paste0("metrics/pileup/", sample_names[i], "-pileup-without-duplicates.txt"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
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

	nuc_pileup = left_join(standard_bam,
						   standard_bam_dedup,
						   by = c("Chromosome", "Position", "Reference_Allele", "Alternate_Allele"))
	index = order(apply(nuc_pileup[,5:ncol(nuc_pileup),drop=FALSE], 1, mean, na.rm=TRUE), decreasing=FALSE)
	nuc_pileup = nuc_pileup[index,,drop=FALSE] %>%
				 arrange(Alternate_Allele) %>%
				 arrange(Reference_Allele)
			 
	col_groups = rep(c("PILEUP\nWITH DUPLICATES", "DEDUPLICATED\nPILEUP"), each=length(sample_names))
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
			  heat.pal.values = c(seq(0,1,l=9), 5))
	dev.off()

} else if (as.numeric(opt$type)==6) {

	suppressPackageStartupMessages(library("superheat"))
	suppressPackageStartupMessages(library("viridis"))

	data = read_tsv(file="metrics/summary/metrics_oxog.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert()
	sample_id = unique(data$SAMPLE_ALIAS)
	nuc_context = unique(data$CONTEXT)
	oxog = matrix(0, nrow=length(nuc_context), ncol=length(sample_id), dimnames=list(nuc_context, sample_id))
	for (i in 1:length(sample_id)) {
		indx = data %>%
			   filter(SAMPLE_ALIAS == sample_id[i]) %>%
			   .[["CONTEXT"]]
		indy = sample_id[i]
		oxog[indx,indy] = as.numeric(data %>%
						  			 filter(SAMPLE_ALIAS == sample_id[i]) %>%
						  			 .[["OXIDATION_ERROR_RATE"]])
	}
	index = order(apply(oxog, 2, mean))
	oxog = oxog[,index,drop=FALSE]
	index = order(apply(oxog, 1, mean))
	oxog = oxog[index,,drop=FALSE]
	
	pdf(file="metrics/report/oxog_error_rate.pdf", height=14, width=14)
	superheat(X = t(oxog),
			  smooth.heat = TRUE,
			  scale = FALSE,
			  legend = TRUE,
			  grid.hline = FALSE,
			  grid.vline = FALSE,
			  row.dendrogram = FALSE,
			  col.dendrogram = FALSE,
			  force.grid.hline = FALSE,
			  force.grid.vline = FALSE,
			  bottom.label.text.angle = 0,
			  bottom.label.text.size = 3.5,
			  bottom.label.size = .15,
			  left.label.size = .15,
			  left.label.text.size = 3.5,
			  print.plot = TRUE)
	dev.off()	   		   
	
} else if (as.numeric(opt$type)==7) {

	data = read_tsv(file="metrics/summary/metrics_coverage.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert()
		   
	pdf(file="metrics/report/read_alignment_summary.pdf", width=14)
	tmp.1 = data %>%
			arrange(desc(N_TOTAL)) %>%
			filter(BAIT_SET=="On-target") %>%
			mutate(`Sample ID` = factor(SAMPLE, levels = unique(SAMPLE), ordered=TRUE)) %>%
			select(`Sample ID`, N = N_ALIGNED) %>%
			mutate(Type = "On-target")
			
	tmp.2 = data %>%
			arrange(desc(N_TOTAL)) %>%
			filter(BAIT_SET=="Off-target") %>%
			mutate(`Sample ID` = factor(SAMPLE, levels = unique(SAMPLE), ordered=TRUE)) %>%
			select(`Sample ID`, N = N_ALIGNED) %>%
			mutate(Type = "Off-target")
			
	tmp.0 = bind_rows(tmp.1, tmp.2)
			
			
	plot.0 = ggplot(tmp.0, aes(x=`Sample ID`, y=N, fill=Type)) +
			 geom_bar(stat="identity") +
			 scale_fill_manual(values=c("On-target"="#d7191c", "Off-target"="#2c7bb6")) +
			 theme_classic(base_size=15) +
			 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
			 labs(fill = "Type", x=" ", title="DISTRIBUTION OF READS", y="Number of read pairs\n") +
			 theme(plot.title = element_text(hjust = 0.5, size=16))
	print(plot.0)
	dev.off()

}
