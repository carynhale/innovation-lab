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

'absolute_' <- function(rho, psi, gamma=1, x)
{
	rho = ifelse(is.na(rho), 1, rho)
	psi = ifelse(is.na(psi), 2, psi)
	return(invisible(((((2^(x/gamma))*(rho*psi+(1-rho)*2)) - ((1-rho)*2))/rho)))
}

'prune_' <- function(x, n=10)
{
	cnm = matrix(NA, nrow=nrow(x), ncol=nrow(x))
	for (j in 1:nrow(x)) {
		cnm[,j] = abs(2^x[j,"log2"] - 2^x[,"log2"])
	}
	cnt = hclust(as.dist(cnm), "average")
	cnc = cutree(tree=cnt, k=n)
	for (j in unique(cnc)) {
		indx = which(cnc==j)
		if (length(indx)>2) {
			mcl = mean(x[indx,"log2"])
			scl = sd(x[indx,"log2"])
			ind = which(x[indx,"log2"]<(mcl+1.96*scl) & x[indx,"log2"]>(mcl-1.96*scl))
			x[indx[ind],"log2"] = mean(x[indx[ind],"log2"])
		} else {
			x[indx,"log2"] = mean(x[indx,"log2"])
		}
	}
	return(x)
}

if (as.numeric(opt$type)==1) {

	log2_ = read.csv(file=paste0("cnvaccess/log2/", opt$sample_name, "-ontarget.txt"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
	log2_ = subset(log2_, log2_$chr<=22 & !is.na(log2_$log2))
	log2_ = winsorize(data = log2_[,c("chr", "start", "log2"),drop=FALSE], method = "mad", tau = 2.05, k = 15, verbose = FALSE)
	pdf(file=paste0("cnvaccess/report/log2/", opt$sample_name, "-ontarget.pdf"), width=14, height=5)
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

} else if (as.numeric(opt$type)==2) {

	log2_ = read.csv(file=paste0("cnvaccess/log2/", opt$sample_name, "-ontarget.txt"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
	log2_ = subset(log2_, log2_$chr<=22 & !is.na(log2_$log2))
	log2_ = winsorize(data = log2_[,c("chr", "start", "log2"),drop=FALSE], method = "mad", tau = 2.15, k = 10, verbose = FALSE)
	segmented_ = pcf(data=log2_, kmin = 2, gamma = 5,
					 normalize = FALSE, fast = FALSE, assembly = "hg19", digits = 4,
           			 return.est = FALSE, save.res = FALSE, file.names = NULL, verbose = FALSE)[,2:7,drop=FALSE]
	colnames(log2_) = c("chr", "pos", "log2")
    colnames(segmented_) = c("chr", "arm", "start", "end", "n", "log2")
	save(log2_, segmented_, file=paste0("cnvaccess/report/segmented/", opt$sample_name, "-ontarget.RData"))
	
	segmented_ = prune_(x=segmented_, n=5)

	pdf(file=paste0("cnvaccess/report/segmented/", opt$sample_name, "-ontarget.pdf"), width=14, height=5)
	par(mar=c(6.1, 6.5, 4.1, 1.1))
	plot(log2_$log2,
		 col = "grey80",
		 pch = 16,
		 cex = .8,
		 ylim = c(-4,4),
		 las = 1,
		 xlab = "",
		 ylab = "",
		 frame.plot = FALSE,
		 axes = FALSE)
	segmented_[,"n"] = cumsum(segmented_[,"n"])
	for (j in 2:nrow(segmented_)) {
 		lines(x=c(segmented_[j-1,"n"], segmented_[j,"n"]), y=rep(segmented_[j,"log2"],2), lty=1, lwd=3.0, col="red")
 	}
 	lines(x=c(1, segmented_[1,"n"]), y=rep(segmented_[1,"log2"],2), lty=1, lwd=3.0, col="red")
	points(c(0-.05*length(log2_$log2),length(log2_$log2)+.01*length(log2_$log2)), c(0,0), type="l", col="red")
	axis(2, at=seq(-3, 3, by=1), labels=seq(-3,3, by=1), cex.axis=1.5, las=1, lwd=1.5)
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

} else if (as.numeric(opt$type)==3) {

	load(file=paste0("cnvaccess/report/segmented/", opt$sample_name, "-ontarget.RData"))
	segmented_ = prune_(x=segmented_, n=5)
	
	rho = seq(from=.05, to=1, length=250)
	psi = seq(from=1.5, to=5, length=250)
	SSE = matrix(NA, nrow=250, ncol=250)
	for (i in 1:length(rho)) {
		for (j in 1:length(psi)) {
			x = absolute_(rho[i], psi[j], gamma=.55, x=segmented_$log2)
			SSE[j,i] = (segmented_$end - segmented_$start) %*% (round(x)-x)^2
		}
	}
	
	pdf(file=paste0("cnvaccess/report/ASCAT/", opt$sample_name, "-ontarget.pdf"), width=7, height=7)
	par(mar=c(6.1, 6.5, 4.1, 1.1))
	image(log2(SSE), col = (colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(25)), useRaster = TRUE, axes = FALSE)
	axis(side=1, at=seq(from=0, to=1, length=5), labels=signif(seq(from=1.5, to=5, length=5), 2))
	axis(side=2, at=seq(from=0, to=1, length=5), labels=signif(seq(from=.05, to=1, length=5), 2), las=1)
	mtext(side=1, text="Ploidy", line=4, cex=1.5)
    mtext(side=2, text="Purity", line=4, cex=1.5)
	box(lwd=1.5)
	dev.off()
	
} else if (as.numeric(opt$type)==4) {

	log2_ = read.csv(file=paste0("cnvaccess/log2/", opt$sample_name, "-offtarget.txt"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
	log2_ = subset(log2_, log2_$chromosome %in% c(1:22) & !is.na(log2_$log2))
	log2_ = winsorize(data = log2_[,c("chromosome", "start", "log2"),drop=FALSE], method = "mad", tau = 1.5, k = 25, verbose = FALSE)
	pdf(file=paste0("cnvaccess/report/log2/", opt$sample_name, "-offtarget.pdf"), width=14, height=5)
	par(mar=c(6.1, 6.5, 4.1, 1.1))
	plot(log2_$log2,
		 col = "grey35",
		 pch = 16,
		 cex = .5,
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

} else if (as.numeric(opt$type)==5) {

	log2_ = read.csv(file=paste0("cnvaccess/log2/", opt$sample_name, "-offtarget.txt"), header=TRUE, sep="\t", stringsAsFactors=FALSE)[,c("chromosome", "start", "end", "log2"),drop=FALSE]
	log2_ = subset(log2_, log2_$chromosome %in% c(1:22) & !is.na(log2_$log2))
	log2_ = winsorize(data = log2_[,c("chromosome", "start", "log2"),drop=FALSE], method = "mad", tau = 1.5, k = 25, verbose = FALSE)
	segmented_ = pcf(data=log2_, kmin = 100, gamma = 60,
					 normalize = FALSE, fast = FALSE, assembly = "hg19", digits = 4,
           			 return.est = FALSE, save.res = FALSE, file.names = NULL, verbose = FALSE)[,2:7,drop=FALSE]
	colnames(log2_) = c("chr", "pos", "log2")
    colnames(segmented_) = c("chr", "arm", "start", "end", "n", "log2")
	save(log2_, segmented_, file=paste0("cnvaccess/report/segmented/", opt$sample_name, "-offtarget.RData"))
	
	segmented_ = prune_(x=segmented_, n=5)

	pdf(file=paste0("cnvaccess/report/segmented/", opt$sample_name, "-offtarget.pdf"), width=14, height=5)
	par(mar=c(6.1, 6.5, 4.1, 1.1))
	plot(log2_$log2,
		 col = "grey80",
		 pch = 16,
		 cex = .5,
		 ylim = c(-4,4),
		 las = 1,
		 xlab = "",
		 ylab = "",
		 frame.plot = FALSE,
		 axes = FALSE)
	segmented_[,"n"] = cumsum(segmented_[,"n"])
	for (j in 2:nrow(segmented_)) {
 		lines(x=c(segmented_[j-1,"n"], segmented_[j,"n"]), y=rep(segmented_[j,"log2"],2), lty=1, lwd=3.0, col="red")
 	}
 	lines(x=c(1, segmented_[1,"n"]), y=rep(segmented_[1,"log2"],2), lty=1, lwd=3.0, col="red")
	points(c(0-.05*length(log2_$log2),length(log2_$log2)+.01*length(log2_$log2)), c(0,0), type="l", col="red")
	axis(2, at=seq(-3, 3, by=1), labels=seq(-3,3, by=1), cex.axis=1.5, las=1, lwd=1.5)
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
