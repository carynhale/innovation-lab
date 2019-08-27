#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("ggplot2"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--type", default = NA, type = 'character', help = "analysis type"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

if (as.numeric(opt$type)==1) {

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
	heatmap(t(as.matrix(data)), Rowv=NA, Colv=NA, scale="row", col=colorRampPalette(c("#f2f0f7", "#cbc9e2", "#9e9ac8", "#756bb1"))( 256 ))
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
			 labs(fill = "Family type")
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
			 labs(fill = "Family type")
	print(plot.0)
	dev.off()

} else if (as.numeric(opt$type)==4) {

	data = read_tsv(file="metrics/summary/umi_family_types.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert()
	
	pdf(file="metrics/report/umi_family_types_probe-AB.pdf")
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
	plot.0 = ggplot(tmp.0, aes(x=Count.x, y=Count.y)) +
			 geom_point() +
			 theme_bw(base_size=15) +
			 labs(x="Probe-A", y="Probe-B") +
			 facet_wrap(~Type.x)
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
			 theme_classic(base_size=15) +
			 labs(x="Family size", y="Frequency") +
			 xlim(0, 50)
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
			 theme_classic(base_size=15) +
			 labs(x="Family size", y="Frequency") +
			 xlim(0, 50)
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
			 theme_classic(base_size=15) +
			 labs(x="Family size", y="Frequency") +
			 xlim(0, 50)
	print(plot.0)
	dev.off()
	
} else if (as.numeric(opt$type)==8) {

	data = read_tsv(file="metrics/summary/metrics_hs.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   filter(LIBRARY=="STANDARD") %>%
		   arrange(desc(MEAN_TARGET_COVERAGE)) %>%
		   mutate(`Sample ID` = factor(SAMPLE, levels=unique(SAMPLE), ordered=TRUE)) %>%
		   mutate(BAIT_SET = ifelse(BAIT_SET=="MSK-ACCESS-v1_0-probe-A", "Probe-A", "Probe-B"))
		   
	
	pdf(file="metrics/report/mean_standard_target_coverage.pdf", width=14)
	plot.0 = ggplot(data, aes(x=`Sample ID`, y=MEAN_TARGET_COVERAGE, fill=BAIT_SET)) +
			 geom_bar(stat="identity", position="dodge") +
			 theme_classic(base_size=15) +
			 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
			 labs(x="Sample ID", y="Depth", fill="Bait set")
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
			 theme_classic(base_size=15) +
			 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
			 labs(x="Sample ID", y="Depth", fill="Bait set")
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
			 theme_classic(base_size=15) +
			 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
			 labs(x="Sample ID", y="Depth", fill="Bait set")
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
			 theme_classic(base_size=15) +
			 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
			 labs(x="Sample ID", y="Depth", fill="Bait set")
	print(plot.0)
	dev.off()
	
} else if (as.numeric(opt$type)==12) {

	data = read_tsv(file="metrics/summary/metrics_hs.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert() %>%
		   filter(LIBRARY=="STANDARD") %>%
		   filter(!is.na(MEAN_TARGET_COVERAGE_NO_DEDUP)) %>%
		   arrange(desc(MEAN_TARGET_COVERAGE_NO_DEDUP)) %>%
		   mutate(`Sample ID` = factor(SAMPLE, levels=unique(SAMPLE), ordered=TRUE))
		   
	pdf(file="metrics/report/mean_standard_target_coverage-nodedup.pdf", width=14)
	plot.0 = ggplot(data, aes(x=`Sample ID`, y=MEAN_TARGET_COVERAGE_NO_DEDUP)) +
			 geom_bar(stat="identity") +
			 theme_classic(base_size=15) +
			 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
			 labs(x="Sample ID", y="Depth")
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
		   
	
	LIBRARIES = c("STANDARD", "UNFILTERED", "DUPLEX", "SIMPLEX")
	for (i in LIBRARIES) {
		pdf(file=paste0("metrics/report/aligment_summary_", tolower(i), ".pdf"), width=14)
		plot.0 = ggplot(data %>% filter(LIBRARY==i), aes(x=`Sample ID`, y=N, fill=Type)) +
				 geom_bar(stat="identity") +
		 		 theme_classic(base_size=15) +
		 		 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
		 		 labs(x="Sample ID", y="Number of reads\n")
		print(plot.0)
		dev.off()
	}
	
}

