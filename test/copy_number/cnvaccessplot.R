#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("copynumber"))
suppressPackageStartupMessages(library("colorspace"))
suppressPackageStartupMessages(library("ASCAT"))
suppressPackageStartupMessages(library("GAP"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--type", default = NA, type = 'character', help = "type of analysis"),
				  make_option("--sample_name", default = NA, type = 'numeric', help = "sample name"))
				  
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

if (as.numeric(opt$type)==1) {

	log2_ = read.csv(file=paste0("cnvaccess/log2/", opt$sample_name, ".txt"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
	log2_ = subset(log2_, log2_$chr<=22 & !is.na(log2_$log2))
	log2_ = winsorize(data = log2_[,c("chr", "start", "log2"),drop=FALSE], method = "mad", tau = 2.05, k = 15, verbose = FALSE)
	pdf(file=paste0("cnvaccess/report/log2/", opt$sample_name, ".pdf"), width=14, height=5)
	par(mar=c(6.1, 6.5, 4.1, 1.1))
	plot(log2_$log2,
		 col = "grey35",
		 pch = 16,
		 cex = .8,
		 ylim = c(-3.5,3.5),
		 las = 1,
		 xlab = "",
		 ylab = "",
		 frame.plot = FALSE,
		 axes = FALSE)
	points(c(0-.05*length(log2_$log2),length(log2_$log2)+.01*length(log2_$log2)), c(0,0), type="l", col="red")
	axis(2, at=NULL, cex.axis=1.5, las=1, lwd=1.5)
	mtext(side=1, text="Chromosome", line=4, cex=1.5)
	mtext(side=2, text=expression(Log[2]~"Ratio"), line=4, cex=1.5)
	title(main=opt$sample_name, cex.main=1.5)
	start = end = NULL
	for (ii in sort(as.numeric(unique(log2_$chr)))) {
		inx = which(log2_$chr==ii)
		start = c(start, min(inx))
		end = c(end, max(inx))
	}
	axis(1, at=c(1, end), labels=rep("", length(end)+1), cex.axis=1.5, lwd=1.5, tcl=.5)
	axis(1, at=.5*(start+end), labels=sort(as.numeric(unique(log2_$chr))), tcl=-.5, lwd=0, lwd.ticks=1.25, tcl=-.5, cex.axis=1.25)
	dev.off()

}