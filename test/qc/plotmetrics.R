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

args_list = list(make_option("--type", default = NA, type = 'character', help = "analysis type"))
				  
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
	superheat(X = t(data), smooth.heat = FALSE, scale = FALSE, legend = TRUE, grid.hline = TRUE, grid.vline = TRUE,
			  force.grid.hline = TRUE, force.grid.vline = TRUE, bottom.label.text.angle = 90, bottom.label.text.size = 3.5,
			  bottom.label.size = .05, left.label.size = .15, left.label.text.size = 3.5, grid.hline.col = "grey90",
			  grid.vline.col = "grey90", print.plot = TRUE)
	dev.off()
	
} else if (as.numeric(opt$type)==2) {

	data = read_tsv(file="metrics/summary/umi_family_types.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert()

	pdf(file="metrics/report/umi_family_types_probe-A.pdf", width=14)
	tmp.0 = data %>%
			filter(BAIT_SET=="MSK-ACCESS_v1.0_probe-A") %>%
			arrange(Count) %>%
			mutate(`Sample ID` = factor(SAMPLE, levels=unique(SAMPLE), ordered=TRUE))
	plot.0 = ggplot(tmp.0, aes(x=`Sample ID`, y=Count, fill = Type)) +
			 geom_bar(stat="identity") +
			 theme_classic(base_size=15) +
			 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
			 labs(fill = "Family type", x=" ", title="Probe-A bait set", y="Count\n") +
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
	plot.0 = ggplot(tmp.0, aes(x=`Sample ID`, y=Count, fill = Type)) +
			 geom_bar(stat="identity") +
			 theme_classic(base_size=15) +
			 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
			 labs(fill = "Family type", x=" ", title="Probe-B bait set", y="Count\n") +
			 theme(plot.title = element_text(hjust = 0.5, size=16))
	print(plot.0)
	dev.off()

} else if (as.numeric(opt$type)==4) {

	data = read_tsv(file="metrics/summary/umi_family_types.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert()
	x = data %>%
		filter(Type=="Singletons") %>%
		filter(BAIT_SET=="MSK-ACCESS_v1.0_probe-A")
	y = data %>%
		filter(Type=="Singletons") %>%
		filter(BAIT_SET=="MSK-ACCESS_v1.0_probe-B")
	z.0 = left_join(x, y, by="SAMPLE")
	x = data %>%
		filter(Type=="Sub-Simplex") %>%
		filter(BAIT_SET=="MSK-ACCESS_v1.0_probe-A")
	y = data %>%
		filter(Type=="Sub-Simplex") %>%
		filter(BAIT_SET=="MSK-ACCESS_v1.0_probe-B")
	z.1 = left_join(x, y, by="SAMPLE")
	x = data %>%
		filter(Type=="Simplex") %>%
		filter(BAIT_SET=="MSK-ACCESS_v1.0_probe-A")
	y = data %>%
		filter(Type=="Simplex") %>%
		filter(BAIT_SET=="MSK-ACCESS_v1.0_probe-B")
	z.2 = left_join(x, y, by="SAMPLE")
	x = data %>%
		filter(Type=="Duplex") %>%
		filter(BAIT_SET=="MSK-ACCESS_v1.0_probe-A")
	y = data %>%
		filter(Type=="Duplex") %>%
		filter(BAIT_SET=="MSK-ACCESS_v1.0_probe-B")
	z.3 = left_join(x, y, by="SAMPLE")
	tmp.0 = bind_rows(z.0, z.1, z.2, z.3)
	pdf(file="metrics/report/umi_family_types_probe-AB.pdf")
	plot.0 = ggplot(tmp.0, aes(x=Count.x, y=Count.y)) +
			 geom_point(alpha=.75, size=1.95, shape=21, color="black", fill="red") +
			 geom_smooth(method="lm", color="goldenrod3", size=1, fullrange=TRUE) +
			 theme_bw(base_size=15) +
			 theme(axis.text=element_text(size=12)) +
			 labs(x="\nProbe-A", y="Probe-B\n") +
			 facet_wrap(~Type.x, scales="free")
	print(plot.0)
	dev.off()
	
} else if (as.numeric(opt$type)==5) {

	data = read_tsv(file="metrics/summary/umi_families.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   filter(FamilyType=="All") %>%
		   dplyr::rename(`Sample ID` = SAMPLE)
	
	pdf(file="metrics/report/umi_family_sizes_all.pdf", width=14)
	plot.0 = ggplot(data, aes(x=FamilySize, y=Frequency, color=`Sample ID`)) +
			 geom_point(size=.5) +
			 geom_line() +
			 scale_x_continuous(trans = 'log2',
			 					breaks = c(1, 10, 20, 30, 40, 50),
    							labels = c(1, 10, 20, 30, 40, 50),
    							limits = c(1, 50)) +
			 theme_classic(base_size=15) +
			 labs(x="\nFamily size", y="Frequency\n", title="All") +
			 theme(plot.title = element_text(hjust = 0.5, size=16))
	print(plot.0)
	dev.off()
	
} else if (as.numeric(opt$type)==6) {

	data = read_tsv(file="metrics/summary/umi_families.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   filter(FamilyType=="Duplex") %>%
		   dplyr::rename(`Sample ID` = SAMPLE)
	
	pdf(file="metrics/report/umi_family_sizes_duplex.pdf", width=14)
	plot.0 = ggplot(data, aes(x=FamilySize, y=Frequency, color=`Sample ID`)) +
			 geom_point(size=.5) +
			 geom_line() +
			 scale_x_continuous(trans = 'log2',
			 					breaks = c(1, 10, 20, 30, 40, 50),
    							labels = c(1, 10, 20, 30, 40, 50),
    							limits = c(1, 50)) +
			 theme_classic(base_size=15) +
			 labs(x="\nFamily size", y="Frequency\n", title="Duplex") +
			 theme(plot.title = element_text(hjust = 0.5, size=16))
	print(plot.0)
	dev.off()
	
} else if (as.numeric(opt$type)==7) {

	data = read_tsv(file="metrics/summary/umi_families.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   filter(FamilyType=="Simplex") %>%
		   dplyr::rename(`Sample ID` = SAMPLE)
	
	pdf(file="metrics/report/umi_family_sizes_simplex.pdf", width=14)
	plot.0 = ggplot(data, aes(x=FamilySize, y=Frequency, color=`Sample ID`)) +
			 geom_point(size=.5) +
			 geom_line() +
			 scale_x_continuous(trans = 'log2',
			 					breaks = c(1, 10, 20, 30, 40, 50),
    							labels = c(1, 10, 20, 30, 40, 50),
    							limits = c(1, 50)) +
			 theme_classic(base_size=15) +
			 labs(x="\nFamily size", y="Frequency\n", title="Simplex") +
			 theme(plot.title = element_text(hjust = 0.5, size=16))
	print(plot.0)
	dev.off()
	
} else if (as.numeric(opt$type)==8) {

	data = read_tsv(file="metrics/summary/metrics_hs.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   filter(LIBRARY=="STANDARD") %>%
		   arrange(desc(MEAN_TARGET_COVERAGE)) %>%
		   mutate(`Sample ID` = factor(SAMPLE, levels=unique(SAMPLE), ordered=TRUE)) %>%
		   mutate(BAIT_SET = ifelse(BAIT_SET=="MSK-ACCESS-v1_0-probe-A", "Probe-A", "Probe-B"))
		   
	pdf(file="metrics/report/mean_standard_target_coverage-dedup.pdf", width=14)
	plot.0 = ggplot(data, aes(x=`Sample ID`, y=MEAN_TARGET_COVERAGE, fill=BAIT_SET)) +
			 geom_bar(stat="identity", position="dodge") +
			 scale_fill_manual(values=c("Probe-A"="#619CFF", "Probe-B"="#00BA38")) +
			 theme_classic(base_size=15) +
			 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
			 labs(x="Sample ID", y="Depth\n", fill="Bait set", title="Deduplicated coverage") +
			 theme(plot.title = element_text(hjust = 0.5, size=16))
	print(plot.0)
	dev.off()
	
} else if (as.numeric(opt$type)==9) {

	data = read_tsv(file="metrics/summary/metrics_hs.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   filter(LIBRARY=="STANDARD") %>%
		   filter(!is.na(MEAN_TARGET_COVERAGE_NO_DEDUP)) %>%
		   arrange(desc(MEAN_TARGET_COVERAGE_NO_DEDUP)) %>%
		   mutate(`Sample ID` = factor(SAMPLE, levels=unique(SAMPLE), ordered=TRUE))
		   
	pdf(file="metrics/report/mean_standard_target_coverage-nodedup.pdf", width=14)
	plot.0 = ggplot(data, aes(x=`Sample ID`, y=MEAN_TARGET_COVERAGE_NO_DEDUP)) +
			 geom_bar(stat="identity", fill="#F8766D") +
			 theme_classic(base_size=15) +
			 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
			 labs(x="Sample ID", y="Depth\n", title="Duplicate coverage") +
			 theme(plot.title = element_text(hjust = 0.5, size=16))
	print(plot.0)
	dev.off()
	
} else if (as.numeric(opt$type)==10) {

	data = read_tsv(file="metrics/summary/metrics_hs.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   filter(LIBRARY=="UNFILTERED") %>%
		   arrange(desc(MEAN_TARGET_COVERAGE)) %>%
		   mutate(`Sample ID` = factor(SAMPLE, levels=unique(SAMPLE), ordered=TRUE)) %>%
		   mutate(BAIT_SET = ifelse(BAIT_SET=="MSK-ACCESS-v1_0-probe-A", "Probe-A", "Probe-B"))
	
	pdf(file="metrics/report/mean_unfiltered_target_coverage.pdf", width=14)
	plot.0 = ggplot(data, aes(x=`Sample ID`, y=MEAN_TARGET_COVERAGE, fill=BAIT_SET)) +
			 geom_bar(stat="identity", position="dodge") +
			 scale_fill_manual(values=c("Probe-A"="#619CFF", "Probe-B"="#00BA38")) +
			 theme_classic(base_size=15) +
			 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
			 labs(x="Sample ID", y="Depth", fill="Bait set", title="Collapsed coverage") +
			 theme(plot.title = element_text(hjust = 0.5, size=16))
	print(plot.0)
	dev.off()
	
} else if (as.numeric(opt$type)==11) {

	data = read_tsv(file="metrics/summary/metrics_hs.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   filter(LIBRARY=="DUPLEX") %>%
		   arrange(desc(MEAN_TARGET_COVERAGE)) %>%
		   mutate(`Sample ID` = factor(SAMPLE, levels=unique(SAMPLE), ordered=TRUE)) %>%
		   mutate(BAIT_SET = ifelse(BAIT_SET=="MSK-ACCESS-v1_0-probe-A", "Probe-A", "Probe-B"))
		   
	
	pdf(file="metrics/report/mean_duplex_target_coverage.pdf", width=14)
	plot.0 = ggplot(data, aes(x=`Sample ID`, y=MEAN_TARGET_COVERAGE, fill=BAIT_SET)) +
			 geom_bar(stat="identity", position="dodge") +
			 scale_fill_manual(values=c("Probe-A"="#619CFF", "Probe-B"="#00BA38")) +
			 theme_classic(base_size=15) +
			 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
			 labs(x="Sample ID", y="Depth", fill="Bait set", title="Duplex coverage") +
			 theme(plot.title = element_text(hjust = 0.5, size=16))
	print(plot.0)
	dev.off()
	
} else if (as.numeric(opt$type)==12) {

	data = read_tsv(file="metrics/summary/metrics_hs.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   filter(LIBRARY=="SIMPLEX") %>%
		   arrange(desc(MEAN_TARGET_COVERAGE)) %>%
		   mutate(`Sample ID` = factor(SAMPLE, levels=unique(SAMPLE), ordered=TRUE)) %>%
		   mutate(BAIT_SET = ifelse(BAIT_SET=="MSK-ACCESS-v1_0-probe-A", "Probe-A", "Probe-B"))
		   
	
	pdf(file="metrics/report/mean_simplex_target_coverage.pdf", width=14)
	plot.0 = ggplot(data, aes(x=`Sample ID`, y=MEAN_TARGET_COVERAGE, fill=BAIT_SET)) +
			 geom_bar(stat="identity", position="dodge") +
			 scale_fill_manual(values=c("Probe-A"="#619CFF", "Probe-B"="#00BA38")) +
			 theme_classic(base_size=15) +
			 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
			 labs(x="Sample ID", y="Depth", fill="Bait set", title="Simplex coverage") +
			 theme(plot.title = element_text(hjust = 0.5, size=16))
	print(plot.0)
	dev.off()
	
} else if (as.numeric(opt$type)==13) {

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
	libnames = c("Standard" = "STANDARD", "Collapsed" = "UNFILTERED", "Duplex" = "DUPLEX", "Simplex" = "SIMPLEX")
	for (i in 1:length(libnames)) {
		plot.0 = ggplot(data %>% filter(LIBRARY==libnames[i]), aes(x=`Sample ID`, y=N, fill=Type)) +
				 geom_bar(stat="identity") +
				 scale_fill_manual(values=c("Aligned"="#619CFF", "Unaligned"="#00BA38")) +
		 		 theme_classic(base_size=15) +
		 		 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
		 		 labs(x="Sample ID", y="Number of read pairs\n", title=names(libnames)[i]) +
		 		 theme(plot.title = element_text(hjust = 0.5, size=16))
		print(plot.0)
	}
	dev.off()

} else if (as.numeric(opt$type)==14) {

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
		 	 theme_classic(base_size=15) +
		 	 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
		 	 labs(x="Sample ID", y="Size (bp)\n", fill="Library", title="Mean insert size") +
		 	 theme(plot.title = element_text(hjust = 0.5, size=16)) +
		 	 ylim(0,200)
	print(plot.0)
	dev.off()

} else if (as.numeric(opt$type)==15) {

	data = read_tsv(file="metrics/summary/metrics_insert_distribution.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert()
		   
	pdf(file="metrics/report/insert_size_distribution.pdf", width=14)
	libnames = c("Standard" = "STANDARD", "Collapsed" = "UNFILTERED", "Duplex" = "DUPLEX", "Simplex" = "SIMPLEX")
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
	
} else if (as.numeric(opt$type)==16) {

	suppressPackageStartupMessages(library("superheat"))

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
	pdf(file="metrics/report/snps_clustering-standard.pdf", width=14, height=14)
	superheat(X = dm, smooth.heat = TRUE, scale = FALSE, legend = TRUE, grid.hline = TRUE, grid.vline = TRUE,
			  row.dendrogram = TRUE, col.dendrogram = TRUE, force.grid.hline = TRUE, force.grid.vline = TRUE,
			  bottom.label.text.angle = 90, bottom.label.text.size = 3.5, bottom.label.size = .15,
			  left.label.size = .15, left.label.text.size = 3.5, grid.hline.col = "grey90",
			  grid.vline.col = "grey90", print.plot = TRUE,
			  heat.lim = c(50, max(dm)), extreme.values.na = FALSE)
	dev.off()
	
} else if (as.numeric(opt$type)==17) {

	suppressPackageStartupMessages(library("superheat"))

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
	pdf(file="metrics/report/snps_clustering-unfiltered.pdf", width=14, height=14)
	superheat(X = dm, smooth.heat = TRUE, scale = FALSE, legend = TRUE, grid.hline = TRUE, grid.vline = TRUE,
			  row.dendrogram = TRUE, col.dendrogram = TRUE, force.grid.hline = TRUE, force.grid.vline = TRUE,
			  bottom.label.text.angle = 90, bottom.label.text.size = 3.5, bottom.label.size = .15,
			  left.label.size = .15, left.label.text.size = 3.5, grid.hline.col = "grey90",
			  grid.vline.col = "grey90", print.plot = TRUE)
	dev.off()
	
} else if (as.numeric(opt$type)==18) {

	suppressPackageStartupMessages(library("superheat"))

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
	pdf(file="metrics/report/snps_clustering-simplex.pdf", width=14, height=14)
	superheat(X = dm, smooth.heat = TRUE, scale = FALSE, legend = TRUE, grid.hline = TRUE, grid.vline = TRUE,
			  row.dendrogram = TRUE, col.dendrogram = TRUE, force.grid.hline = TRUE, force.grid.vline = TRUE,
			  bottom.label.text.angle = 90, bottom.label.text.size = 3.5, bottom.label.size = .15,
			  left.label.size = .15, left.label.text.size = 3.5, grid.hline.col = "grey90",
			  grid.vline.col = "grey90", print.plot = TRUE)
	dev.off()
	
} else if (as.numeric(opt$type)==19) {

	suppressPackageStartupMessages(library("superheat"))

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
	pdf(file="metrics/report/snps_clustering-duplex.pdf", width=14, height=14)
	superheat(X = dm, smooth.heat = TRUE, scale = FALSE, legend = TRUE, grid.hline = TRUE, grid.vline = TRUE,
			  row.dendrogram = TRUE, col.dendrogram = TRUE, force.grid.hline = TRUE, force.grid.vline = TRUE,
			  bottom.label.text.angle = 90, bottom.label.text.size = 3.5, bottom.label.size = .15,
			  left.label.size = .15, left.label.text.size = 3.5, grid.hline.col = "grey90",
			  grid.vline.col = "grey90", print.plot = TRUE)
	dev.off()
	
}
