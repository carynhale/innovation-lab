#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("copynumber"))
suppressPackageStartupMessages(library("colorspace"))
suppressPackageStartupMessages(library("ASCAT"))
suppressPackageStartupMessages(library("GAP"))

'plot_log2_' <- function(x, title = "")
{
   	par(mar=c(5, 5, 4, 2)+.1)
   	data("CytoBand")
   	end = NULL
   	for (i in 1:23) {
   		end = c(end, max(CytoBand[CytoBand[,1]==i,"End"]))
   	}
   	end = cumsum(end)
   	start = c(1, end[1:22]+1)
   	CytoBand = cbind(start, end)
   	index = NULL
   	for (i in 1:23) {
   		index = c(index, seq(from = CytoBand[i, "start"], to=CytoBand[i, "end"], length=sum(x$chromosome==i)))
   	}
	plot(index, x$log2, type="p", pch=".", cex=1.95, col=c("grey80", "lightblue")[x$chromosome%%2 + 1], axes=FALSE, frame=FALSE, xlab="", ylab="", main="", ylim=c(-4.5,4.5))
  	axis(1, at = c(start, end[length(end)]), labels=rep("", length(start)+1), tcl=.5)
  	axis(1, at = .5*(start+end), labels=c(1:22, "X"), tcl=-.5, lwd=0, lwd.ticks=1, tcl=-.25)
  	axis(2, at = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4), cex.axis = 1, las = 1)
	mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
	points(c(0-.05*max(index),max(index)+.01*max(index)), c(0,0), type="l", col="black", lwd=1)
	title(main = title, cex.main=.75, font.main=1)
}

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--in_file", default = NA, type = 'character', help = "input file name"))
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

outfile_on_target = gsub("cnr", "log2", gsub(".cnr", ".ontarget.pdf", opt$in_file, fixed=TRUE), fixed=TRUE)
outfile_off_target = gsub("cnr", "log2", gsub(".cnr", ".offtarget.pdf", opt$in_file, fixed=TRUE), fixed=TRUE)

data = read.table(file=opt$in_file, header=TRUE, sep="\t", comment.char="#", stringsAsFactors=FALSE)
data = subset(data, data[,"depth"]!=0)

if (nrow(data)==0) {
	system(paste0("touch ", outfile_on_target))
	system(paste0("touch ", outfile_off_target))
} else {
	data[,"chromosome"] = gsub(pattern="chr", replacement="", x=data[,"chromosome"], fixed=TRUE)
	data[data[,"chromosome"]=="X", "chromosome"] = 23
	data[data[,"chromosome"]=="Y", "chromosome"] = 24
	data[,"chromosome"] = as.numeric(data[,"chromosome"])
	data = subset(data, data[,"chromosome"]<=23)
	
	if (sum(data$gene=="-")>0) {
		flag = 1
	} else if (sum(data$gene=="Antitarget")>0) {
		flag = 2
	}
	
	if (flag==1) {
		ontarget = subset(data, data$gene=="-")
	} else if (flag==2) {
		ontarget = subset(data, data$gene!="Antitarget")
	}
	
	pdf(file=outfile_on_target, width=10, height=4.25)
	plot_log2_(x=ontarget, title=gsub("cnvkit/cnr/", "", gsub(".cnr", "", opt$in_file, fixed=TRUE), fixed=TRUE))
	dev.off()
	
	if (flag==1) {
		offtarget = subset(data, data$gene!="-")
	} else if (flag==2) {
		offtarget = subset(data, data$gene=="Antitarget")
	}
	
	tmp = offtarget[,c("chromosome", "start", "log2"),drop=FALSE]
	tmp = winsorize(data=tmp, tau=3.5, k=25, verbose=FALSE, return.outliers=TRUE)
	offtarget[tmp$wins.outliers[,3]!=0,"log2"] = NA
	pdf(file=outfile_off_target, width=10, height=4.25)
	plot_log2_(x=offtarget, title=gsub("cnvkit/cnr/", "", gsub(".cnr", "", opt$in_file, fixed=TRUE), fixed=TRUE))
	dev.off()
}
