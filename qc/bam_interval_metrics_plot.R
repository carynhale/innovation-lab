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

if (as.numeric(opt$type)==1) {

	suppressPackageStartupMessages(library("superheat"))

	data = read.csv(file="metrics/summary/umi_frequencies.tsv", sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)
	data = data	 %>%
		   rename_all(funs(gsub(pattern=".", replacement="-", x=make.names(names(data)), fixed=TRUE))) %>%
		   type_convert() %>%
		   replace(is.na(.), 0)
	index = order(apply(data, 1, sum))
	data = data[index,,drop=FALSE]
	index = order(apply(data, 2, sum))
	data = data[,index,drop=FALSE]
	pdf(file="metrics/report/umi_frequencies.pdf", width=14, height=14)
	superheat(X = t(data), smooth.heat = TRUE, scale = FALSE, legend = TRUE, grid.hline = TRUE, grid.vline = TRUE,
			  force.grid.hline = TRUE, force.grid.vline = TRUE, bottom.label.text.angle = 90, bottom.label.text.size = 3.5,
			  bottom.label.size = .05, left.label.size = .15, left.label.text.size = 3.5, grid.hline.col = "grey90",
			  grid.vline.col = "grey90", title="     ", print.plot = TRUE)
	dev.off()
	
} else if (as.numeric(opt$type)==2) {

	data = read_tsv(file="metrics/summary/umi_family_types.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert()

	pdf(file="metrics/report/umi_family_types_probe-A.pdf", width=14)
	tmp.0 = data %>%
			filter(BAIT_SET=="MSK-ACCESS_v1.0_probe-A") %>%
			arrange(Count) %>%
			mutate(`Sample ID` = factor(SAMPLE, levels=unique(SAMPLE), ordered=TRUE))
	plot.0 = ggplot(tmp.0, aes(x=`Sample ID`, y=Count, fill=Type)) +
			 geom_bar(stat="identity") +
			 scale_fill_manual(values=c("Duplex"      = "#d7191c",
			 							"Simplex"     = "#fdae61",
			 							"Singletons"  = "#abd9e9",
			 							"Sub-Simplex" = "#2c7bb6")) +
			 theme_classic(base_size=15) +
			 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
			 labs(fill = "Family type", x=" ", title="PROBE-A BAIT SET", y="Count\n") +
			 theme(plot.title = element_text(hjust = 0.5, size=16))
	print(plot.0)
	dev.off()

} else if (as.numeric(opt$type)==3) {

	data = read_tsv(file="metrics/summary/umi_family_types.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert()

	pdf(file="metrics/report/umi_family_types_probe-B.pdf", width=14)
	tmp.0 = data %>%
			filter(BAIT_SET=="MSK-ACCESS_v1.0_probe-B") %>%
			arrange(Count) %>%
			mutate(`Sample ID` = factor(SAMPLE, levels=unique(SAMPLE), ordered=TRUE))
	plot.0 = ggplot(tmp.0, aes(x=`Sample ID`, y=Count, fill=Type)) +
			 geom_bar(stat="identity") +
			 scale_fill_manual(values=c("Duplex"      = "#d7191c",
			 							"Simplex"     = "#fdae61",
			 							"Singletons"  = "#abd9e9",
			 							"Sub-Simplex" = "#2c7bb6")) +
			 theme_classic(base_size=15) +
			 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
			 labs(fill = "Family type", x=" ", title="PROBE-B BAIT SET", y="Count\n") +
			 theme(plot.title = element_text(hjust = 0.5, size=16))
	print(plot.0)
	dev.off()

} else if (as.numeric(opt$type)==4) {

	data = read_tsv(file="metrics/summary/umi_families.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   filter(FamilyType=="All") %>%
		   dplyr::rename(`Sample ID` = SAMPLE)
	
	pdf(file="metrics/report/umi_family_sizes_all.pdf", width=14)
	plot.0 = ggplot(data, aes(x=FamilySize, y=Frequency, color=`Sample ID`)) +
			 geom_point(size=.5) +
			 geom_line() +
			 scale_x_continuous(trans = 'log2',
			 					breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50),
    							labels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50),
    							limits = c(1, 50)) +
			 theme_classic(base_size=15) +
			 labs(x="\nFamily size", y="Frequency\n", title="ALL FAMILIY SIZES") +
			 theme(plot.title = element_text(hjust = 0.5, size=16))
	print(plot.0)
	dev.off()
	
} else if (as.numeric(opt$type)==5) {

	data = read_tsv(file="metrics/summary/umi_families.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   filter(FamilyType=="Duplex") %>%
		   dplyr::rename(`Sample ID` = SAMPLE)
	
	pdf(file="metrics/report/umi_family_sizes_duplex.pdf", width=14)
	plot.0 = ggplot(data, aes(x=FamilySize, y=Frequency, color=`Sample ID`)) +
			 geom_point(size=.5) +
			 geom_line() +
			 scale_x_continuous(trans = 'log2',
			 					breaks = c(2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50),
    							labels = c(2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50),
    							limits = c(2, 50)) +
			 theme_classic(base_size=15) +
			 labs(x="\nFamily size", y="Frequency\n", title="DUPLEX FAMILY SIZES") +
			 theme(plot.title = element_text(hjust = 0.5, size=16))
	print(plot.0)
	dev.off()
	
} else if (as.numeric(opt$type)==6) {

	data = read_tsv(file="metrics/summary/umi_families.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   filter(FamilyType=="Simplex") %>%
		   dplyr::rename(`Sample ID` = SAMPLE)
	
	pdf(file="metrics/report/umi_family_sizes_simplex.pdf", width=14)
	plot.0 = ggplot(data, aes(x=FamilySize, y=Frequency, color=`Sample ID`)) +
			 geom_point(size=.5) +
			 geom_line() +
			 scale_x_continuous(trans = 'log2',
			 					breaks = c(3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50),
    							labels = c(3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50),
    							limits = c(3, 50)) +
			 theme_classic(base_size=15) +
			 labs(x="\nFamily size", y="Frequency\n", title="SIMPLEX FAMILY SIZES") +
			 theme(plot.title = element_text(hjust = 0.5, size=16))
	print(plot.0)
	dev.off()
	
} else if (as.numeric(opt$type)==7) {

	data = read_tsv(file="metrics/summary/metrics_hs.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   filter(LIBRARY=="STANDARD") %>%
		   arrange(desc(MEAN_TARGET_COVERAGE)) %>%
		   mutate(`Sample ID` = factor(SAMPLE, levels=unique(SAMPLE), ordered=TRUE)) %>%
		   mutate(BAIT_SET = ifelse(BAIT_SET=="MSK-ACCESS-v1_0-probe-A", "Probe-A", "Probe-B"))
		   
	pdf(file="metrics/report/mean_standard_target_coverage-dedup.pdf", width=14)
	plot.0 = ggplot(data, aes(x=`Sample ID`, y=MEAN_TARGET_COVERAGE, fill=BAIT_SET)) +
			 geom_bar(stat="identity", position="dodge") +
			 scale_fill_manual(values=c("Probe-A"="#d7191c", "Probe-B"="#2c7bb6")) +
			 theme_classic(base_size=15) +
			 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
			 labs(x="Sample ID", y="Depth\n", fill="Bait set", title="MEAN DEDUPLICATED COVERAGE") +
			 theme(plot.title = element_text(hjust = 0.5, size=16))
	print(plot.0)
	dev.off()
	
} else if (as.numeric(opt$type)==8) {

	data = read_tsv(file="metrics/summary/metrics_hs.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   filter(LIBRARY=="STANDARD") %>%
		   arrange(desc(MEAN_TARGET_COVERAGE_NO_DEDUP)) %>%
		   mutate(`Sample ID` = factor(SAMPLE, levels=unique(SAMPLE), ordered=TRUE)) %>%
		   mutate(BAIT_SET = ifelse(BAIT_SET=="MSK-ACCESS-v1_0-probe-A", "Probe-A", "Probe-B"))
		   
	pdf(file="metrics/report/mean_standard_target_coverage-nodedup.pdf", width=14)
	plot.0 = ggplot(data, aes(x=`Sample ID`, y=MEAN_TARGET_COVERAGE_NO_DEDUP, fill=BAIT_SET)) +
			 geom_bar(stat="identity", position="dodge") +
			 scale_fill_manual(values=c("Probe-A"="#d7191c", "Probe-B"="#2c7bb6")) +
			 theme_classic(base_size=15) +
			 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
			 labs(x="Sample ID", y="Depth\n", fill="Bait set", title="MEAN COVERAGE WITH DUPLICATES") +
			 theme(plot.title = element_text(hjust = 0.5, size=16))
	print(plot.0)
	dev.off()
	
} else if (as.numeric(opt$type)==9) {

	data = read_tsv(file="metrics/summary/metrics_hs.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   filter(LIBRARY=="UNFILTERED") %>%
		   arrange(desc(MEAN_TARGET_COVERAGE)) %>%
		   mutate(`Sample ID` = factor(SAMPLE, levels=unique(SAMPLE), ordered=TRUE)) %>%
		   mutate(BAIT_SET = ifelse(BAIT_SET=="MSK-ACCESS-v1_0-probe-A", "Probe-A", "Probe-B"))
	
	pdf(file="metrics/report/mean_unfiltered_target_coverage.pdf", width=14)
	plot.0 = ggplot(data, aes(x=`Sample ID`, y=MEAN_TARGET_COVERAGE, fill=BAIT_SET)) +
			 geom_bar(stat="identity", position="dodge") +
			 scale_fill_manual(values=c("Probe-A"="#d7191c", "Probe-B"="#2c7bb6")) +
			 theme_classic(base_size=15) +
			 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
			 labs(x="Sample ID", y="Depth\n", fill="Bait set", title="MEAN COLLAPSED COVERAGE") +
			 theme(plot.title = element_text(hjust = 0.5, size=16))
	print(plot.0)
	dev.off()
	
} else if (as.numeric(opt$type)==10) {

	data = read_tsv(file="metrics/summary/metrics_hs.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   filter(LIBRARY=="DUPLEX") %>%
		   arrange(desc(MEAN_TARGET_COVERAGE)) %>%
		   mutate(`Sample ID` = factor(SAMPLE, levels=unique(SAMPLE), ordered=TRUE)) %>%
		   mutate(BAIT_SET = ifelse(BAIT_SET=="MSK-ACCESS-v1_0-probe-A", "Probe-A", "Probe-B"))
		   
	
	pdf(file="metrics/report/mean_duplex_target_coverage.pdf", width=14)
	plot.0 = ggplot(data, aes(x=`Sample ID`, y=MEAN_TARGET_COVERAGE, fill=BAIT_SET)) +
			 geom_bar(stat="identity", position="dodge") +
			 scale_fill_manual(values=c("Probe-A"="#d7191c", "Probe-B"="#2c7bb6")) +
			 theme_classic(base_size=15) +
			 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
			 labs(x="Sample ID", y="Depth\n", fill="Bait set", title="MEAN DUPLEX COVERAGE") +
			 theme(plot.title = element_text(hjust = 0.5, size=16))
	print(plot.0)
	dev.off()
	
} else if (as.numeric(opt$type)==11) {

	data = read_tsv(file="metrics/summary/metrics_hs.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   filter(LIBRARY=="SIMPLEX") %>%
		   arrange(desc(MEAN_TARGET_COVERAGE)) %>%
		   mutate(`Sample ID` = factor(SAMPLE, levels=unique(SAMPLE), ordered=TRUE)) %>%
		   mutate(BAIT_SET = ifelse(BAIT_SET=="MSK-ACCESS-v1_0-probe-A", "Probe-A", "Probe-B"))
		   
	
	pdf(file="metrics/report/mean_simplex_target_coverage.pdf", width=14)
	plot.0 = ggplot(data, aes(x=`Sample ID`, y=MEAN_TARGET_COVERAGE, fill=BAIT_SET)) +
			 geom_bar(stat="identity", position="dodge") +
			 scale_fill_manual(values=c("Probe-A"="#d7191c", "Probe-B"="#2c7bb6")) +
			 theme_classic(base_size=15) +
			 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
			 labs(x="Sample ID", y="Depth\n", fill="Bait set", title="MEAN SIMPLEX COVERAGE") +
			 theme(plot.title = element_text(hjust = 0.5, size=16))
	print(plot.0)
	dev.off()
	
} else if (as.numeric(opt$type)==12) {

	data = read_tsv(file="metrics/summary/metrics_idx.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert()
		   
	x = data %>%
		dplyr::select(N = N_ALIGNED, SAMPLE, LIBRARY) %>%
		mutate(TYPE = "Aligned")
	y = data %>%
		dplyr::select(N = N_UNALIGNED, SAMPLE, LIBRARY) %>%
		mutate(TYPE = "Unaligned")
	data = bind_rows(x, y) %>%
		   arrange(N) %>%
		   dplyr::rename(`Sample ID` = SAMPLE, Type = TYPE) %>%
		   mutate(`Sample ID` = factor(`Sample ID`, levels=unique(`Sample ID`), ordered=TRUE))
		   
	pdf(file="metrics/report/aligment_summary.pdf", width=14)
	libnames = c("STANDARD" = "STANDARD", "COLLAPSED" = "UNFILTERED", "DUPLEX" = "DUPLEX", "SIMPLEX" = "SIMPLEX")
	for (i in 1:length(libnames)) {
		plot.0 = ggplot(data %>% filter(LIBRARY==libnames[i]), aes(x=`Sample ID`, y=N, fill=Type)) +
				 geom_bar(stat="identity") +
				 scale_fill_manual(values=c("Aligned"="#d7191c", "Unaligned"="#2c7bb6")) +
		 		 theme_classic(base_size=15) +
		 		 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
		 		 labs(x="Sample ID", y="Number of read pairs\n", title=names(libnames)[i]) +
		 		 theme(plot.title = element_text(hjust = 0.5, size=16))
		print(plot.0)
	}
	dev.off()

} else if (as.numeric(opt$type)==13) {

	data = read_tsv(file="metrics/summary/metrics_insert.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   arrange(desc(MEDIAN_INSERT_SIZE)) %>%
		   dplyr::rename(`Sample ID` = SAMPLE) %>%
		   mutate(LIBRARY = case_when(
		   		LIBRARY == "STANDARD" ~ "Standard",
		   		LIBRARY == "UNFILTERED" ~ "Collapsed",
		   		LIBRARY == "DUPLEX" ~ "Duplex",
		   		LIBRARY == "SIMPLEX" ~ "Simplex"))
		   
	pdf(file="metrics/report/insert_size_summary.pdf", width=14)
	plot.0 = ggplot(data, aes(x=`Sample ID`, y=MEDIAN_INSERT_SIZE, fill=LIBRARY)) +
			 geom_bar(stat="identity", position="dodge") +
			 scale_fill_manual(values=c("Collapsed" = "#2c7bb6",
			 							"Duplex"    = "#d7191c",
			 							"Simplex"   = "#fdae61",
			 							"Standard"  = "#abd9e9")) +
		 	 theme_classic(base_size=15) +
		 	 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
		 	 labs(x="Sample ID", y="Size (bp)\n", fill="Library", title="MEAN INSERT SIZE") +
		 	 theme(plot.title = element_text(hjust = 0.5, size=16)) +
		 	 ylim(0,200)
	print(plot.0)
	dev.off()

} else if (as.numeric(opt$type)==14) {

	data = read_tsv(file="metrics/summary/metrics_insert_distribution.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert()
		   
	pdf(file="metrics/report/insert_size_distribution.pdf", width=14)
	libnames = c("STANDARD" = "STANDARD", "COLLAPSED" = "UNFILTERED", "DUPLEX" = "DUPLEX", "SIMPLEX" = "SIMPLEX")
	for (i in 1:length(libnames)) {
		tmp.0 = data %>%
				filter(LIBRARY==libnames[i]) %>%
				dplyr::select(-LIBRARY)
		tmp.1 = list()
		for (j in 1:ncol(tmp.0)) {
			x = 1:nrow(tmp.0)
			y = tmp.0[,j,drop=TRUE]
			x = x[!is.na(y)]
			y = y[!is.na(y)]
			z = smooth.spline(x=x, y=y, spar=0.3)
			tmp.1[[j]] = data.frame(x=z$x, y=z$y) %>%
				    	 mutate(`Sample ID` = colnames(tmp.0)[j])
		}
		tmp.1 = do.call(rbind, tmp.1)
		plot.0 = ggplot(tmp.1, aes(x=x, y=y, color=`Sample ID`)) +
				 geom_point(size=.5) +
				 geom_line() +
				 theme_classic(base_size=15) +
				 labs(x="\nInsert size (bp)", y="Frequency\n", title=names(libnames)[i]) +
				 theme(plot.title = element_text(hjust = 0.5, size=16)) +
				 xlim(0, 400)
		print(plot.0)
	}
	dev.off()
	
} else if (as.numeric(opt$type)==15) {

	suppressPackageStartupMessages(library("superheat"))
	suppressPackageStartupMessages(library("viridis"))

	data = read.csv(file="metrics/summary/snps_filtered-standard.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
	data = data	 %>%
		   rename_all(funs(gsub(pattern=".", replacement="-", x=make.names(names(data)), fixed=TRUE))) %>%
		   type_convert()
	data = data[apply(data, 1, function(x) {sum(is.na(x))})!=ncol(data),,drop=FALSE]
	data[is.na(data)] = 1
	data[data==3] = 1
	for (i in 1:2) {
		data = data[apply(data, 1, function(x) {sum(x==i, na.rm=TRUE)})!=ncol(data),,drop=FALSE]
	}
	dm = as.matrix(dist(t(data), method="manhattan", diag=TRUE, upper=TRUE))
	dm = 1-((max(dm)-dm)/(max(dm) - min(dm)))
	pdf(file="metrics/report/snps_clustering-standard.pdf", width=14, height=14)
	superheat(X = dm, smooth.heat = TRUE, scale = FALSE, legend = TRUE, grid.hline = TRUE, grid.vline = TRUE,
			  row.dendrogram = TRUE, col.dendrogram = TRUE, force.grid.hline = TRUE, force.grid.vline = TRUE,
			  bottom.label.text.angle = 90, bottom.label.text.size = 3.5, bottom.label.size = .15,
			  left.label.size = .15, left.label.text.size = 3.5, grid.hline.col = "grey90",
			  grid.vline.col = "grey90", heat.pal = viridis(n=100), heat.pal.values = seq(from=0, to=1, length=100),
			  print.plot = TRUE)
	dev.off()
	
} else if (as.numeric(opt$type)==16) {

	suppressPackageStartupMessages(library("superheat"))
	suppressPackageStartupMessages(library("viridis"))

	data = read.csv(file="metrics/summary/snps_filtered-unfiltered.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
	data = data	 %>%
		   rename_all(funs(gsub(pattern=".", replacement="-", x=make.names(names(data)), fixed=TRUE))) %>%
		   type_convert()
	data = data[apply(data, 1, function(x) {sum(is.na(x))})!=ncol(data),,drop=FALSE]
	data[is.na(data)] = 1
	data[data==3] = 1
	for (i in 1:2) {
		data = data[apply(data, 1, function(x) {sum(x==i, na.rm=TRUE)})!=ncol(data),,drop=FALSE]
	}
	dm = as.matrix(dist(t(data), method="manhattan", diag=TRUE, upper=TRUE))
	dm = 1-((max(dm)-dm)/(max(dm) - min(dm)))
	pdf(file="metrics/report/snps_clustering-unfiltered.pdf", width=14, height=14)
	superheat(X = dm, smooth.heat = TRUE, scale = FALSE, legend = TRUE, grid.hline = TRUE, grid.vline = TRUE,
			  row.dendrogram = TRUE, col.dendrogram = TRUE, force.grid.hline = TRUE, force.grid.vline = TRUE,
			  bottom.label.text.angle = 90, bottom.label.text.size = 3.5, bottom.label.size = .15,
			  left.label.size = .15, left.label.text.size = 3.5, grid.hline.col = "grey90",
			  grid.vline.col = "grey90", heat.pal = viridis(n=100), heat.pal.values = seq(from=0, to=1, length=100),
			  print.plot = TRUE)
	dev.off()
	
} else if (as.numeric(opt$type)==17) {

	suppressPackageStartupMessages(library("superheat"))
	suppressPackageStartupMessages(library("viridis"))

	data = read.csv(file="metrics/summary/snps_filtered-simplex.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
	data = data	 %>%
		   rename_all(funs(gsub(pattern=".", replacement="-", x=make.names(names(data)), fixed=TRUE))) %>%
		   type_convert()
	data = data[apply(data, 1, function(x) {sum(is.na(x))})!=ncol(data),,drop=FALSE]
	data[is.na(data)] = 1
	data[data==3] = 1
	for (i in 1:2) {
		data = data[apply(data, 1, function(x) {sum(x==i, na.rm=TRUE)})!=ncol(data),,drop=FALSE]
	}
	dm = as.matrix(dist(t(data), method="manhattan", diag=TRUE, upper=TRUE))
	dm = 1-((max(dm)-dm)/(max(dm) - min(dm)))
	pdf(file="metrics/report/snps_clustering-simplex.pdf", width=14, height=14)
	superheat(X = dm, smooth.heat = TRUE, scale = FALSE, legend = TRUE, grid.hline = TRUE, grid.vline = TRUE,
			  row.dendrogram = TRUE, col.dendrogram = TRUE, force.grid.hline = TRUE, force.grid.vline = TRUE,
			  bottom.label.text.angle = 90, bottom.label.text.size = 3.5, bottom.label.size = .15,
			  left.label.size = .15, left.label.text.size = 3.5, grid.hline.col = "grey90",
			  grid.vline.col = "grey90", heat.pal = viridis(n=100), heat.pal.values = seq(from=0, to=1, length=100),
			  print.plot = TRUE)
	dev.off()
	
} else if (as.numeric(opt$type)==18) {

	suppressPackageStartupMessages(library("superheat"))
	suppressPackageStartupMessages(library("viridis"))

	data = read.csv(file="metrics/summary/snps_filtered-duplex.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
	data = data	 %>%
		   rename_all(funs(gsub(pattern=".", replacement="-", x=make.names(names(data)), fixed=TRUE))) %>%
		   type_convert()
	data = data[apply(data, 1, function(x) {sum(is.na(x))})!=ncol(data),,drop=FALSE]
	data[is.na(data)] = 1
	data[data==3] = 1
	for (i in 1:2) {
		data = data[apply(data, 1, function(x) {sum(x==i, na.rm=TRUE)})!=ncol(data),,drop=FALSE]
	}
	dm = as.matrix(dist(t(data), method="manhattan", diag=TRUE, upper=TRUE))
	dm = 1-((max(dm)-dm)/(max(dm) - min(dm)))
	pdf(file="metrics/report/snps_clustering-duplex.pdf", width=14, height=14)
	superheat(X = dm, smooth.heat = TRUE, scale = FALSE, legend = TRUE, grid.hline = TRUE, grid.vline = TRUE,
			  row.dendrogram = TRUE, col.dendrogram = TRUE, force.grid.hline = TRUE, force.grid.vline = TRUE,
			  bottom.label.text.angle = 90, bottom.label.text.size = 3.5, bottom.label.size = .15,
			  left.label.size = .15, left.label.text.size = 3.5, grid.hline.col = "grey90",
			  grid.vline.col = "grey90", heat.pal = viridis(n=100), heat.pal.values = seq(from=0, to=1, length=100),
			  print.plot = TRUE)
	dev.off()
	
} else if (as.numeric(opt$type)==19) {

	data = read_tsv(file="metrics/summary/metrics_ts.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   dplyr::rename(`Sample ID` = SAMPLE)

 	pdf(file="metrics/report/read_alignment_summary.pdf", width=14)
	tmp.1 = data %>%
			arrange(desc(N_TOTAL)) %>%
			filter(BAIT_SET=="Probe-A") %>%
			mutate(`Sample ID` = factor(`Sample ID`, levels = `Sample ID`, ordered=TRUE)) %>%
			select(`Sample ID`, N = N_ALIGNED) %>%
			mutate(Type = "Probe-A")
			
	tmp.2 = data %>%
			arrange(desc(N_TOTAL)) %>%
			filter(BAIT_SET=="Probe-B") %>%
			mutate(`Sample ID` = factor(`Sample ID`, levels = `Sample ID`, ordered=TRUE)) %>%
			select(`Sample ID`, N = N_ALIGNED) %>%
			mutate(Type = "Probe-B")
			
	tmp.3 = data %>%
			arrange(desc(N_TOTAL)) %>%
			filter(BAIT_SET=="Probe-AB") %>%
			mutate(`Sample ID` = factor(`Sample ID`, levels = `Sample ID`, ordered=TRUE)) %>%
			select(`Sample ID`, N = N_ALIGNED) %>%
			mutate(Type = "Off-target")
			
	tmp.0 = bind_rows(tmp.1, tmp.2, tmp.3)
			
			
	plot.0 = ggplot(tmp.0, aes(x=`Sample ID`, y=N, fill=Type)) +
			 geom_bar(stat="identity") +
			 scale_fill_manual(values=c("Probe-A"="#d7191c", "Probe-B"="#2c7bb6", "Off-target"="#fdae61")) +
			 theme_classic(base_size=15) +
			 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
			 labs(fill = "Type", x=" ", title="DISTRIBUTION OF READS", y="Number of read pairs\n") +
			 theme(plot.title = element_text(hjust = 0.5, size=16))
	print(plot.0)
	dev.off()

} else if (as.numeric(opt$type)==20) {

	suppressPackageStartupMessages(library("superheat"))
	suppressPackageStartupMessages(library("viridis"))
	
	sample_names = unlist(strsplit(x=as.character(opt$sample_names), split=" ", fixed=TRUE))
	nuc_metrics = list()
	for (i in 1:length(sample_names)) {
		pileup_metrics = read_tsv(file=paste0("metrics/standard/", sample_names[i], "-pileup.txt"), col_names = FALSE, col_types = cols(.default = col_character())) %>%
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
								  filter(Allele_Frequency<10) %>%
								  filter(Chromosome=="21") %>%
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
								  filter(Allele_Frequency<10) %>%
								  filter(Chromosome=="21") %>%
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
								  filter(Allele_Frequency<10) %>%
								  filter(Chromosome=="21") %>%
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
								  filter(Allele_Frequency<10) %>%
								  filter(Chromosome=="21") %>%
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

	nuc_pileup = nuc_pileup %>%
				 dplyr::select(-Chromosome, -Position, -Reference_Allele, -Alternate_Allele)

	index = apply(nuc_pileup, 1, function(x) {sum(is.na(x))})==0
	pdf(file="metrics/report/non_reference_calls.pdf", height=14, width=14)
	superheat(X = as.matrix(nuc_pileup[index,,drop=FALSE]),
			  smooth.heat = FALSE,
			  scale = FALSE,
			  legend = FALSE,
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
			  heat.pal.values = c(seq(0,.1,l=9), 10))
	dev.off()

}
