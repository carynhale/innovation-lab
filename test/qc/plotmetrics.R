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
			mutate(SAMPLE = factor(SAMPLE, levels=unique(SAMPLE), ordered=TRUE))
	plot.0 = ggplot(tmp.0, aes(x=SAMPLE, y=Count, fill = Type)) +
			 geom_bar(stat="identity") +
			 theme_classic(base_size=15) +
			 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
			 labs(fill = "Family type\n")
	print(plot.0)
	dev.off()

} else if (as.numeric(opt$type)==3) {

	data = read_tsv(file="metrics/summary/umi_family_types.tsv", col_types = cols(.default = col_character())) %>%
		   type_convert()

	pdf(file="metrics/report/umi_family_types_probe-B.pdf", width=14)
	tmp.0 = data %>%
			filter(BAIT_SET=="MSK-ACCESS_v1.0_probe-B") %>%
			arrange(Count) %>%
			mutate(SAMPLE = factor(SAMPLE, levels=unique(SAMPLE), ordered=TRUE))
	plot.0 = ggplot(tmp.0, aes(x=SAMPLE, y=Count, fill = Type)) +
			 geom_bar(stat="identity") +
			 theme_classic(base_size=15) +
			 theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_text(size=13)) +
			 labs(fill = "Family type\n")
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
	
}
