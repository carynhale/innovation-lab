#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(
				 make_option("--type", default = NA, type = 'character', help = "type of analysis"),
				 make_option("--sample_name", default = NA, type = 'character', help = "sample name")
				)
				  
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

if (as.numeric(opt$type)==1) {

	suppressPackageStartupMessages(library("copynumber"))
	suppressPackageStartupMessages(library("colorspace"))
	suppressPackageStartupMessages(library("ASCAT"))
	suppressPackageStartupMessages(library("GAP"))

	data = read.csv(file=paste0("cnvkit/cnr/", opt$sample_name, ".cnr"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
	data = subset(data, data[,"depth"]!=0)
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
		data = subset(data, data$gene=="-")
		data = data[,c("chromosome", "start", "log2"),drop=FALSE]
		colnames(data) = c("chrom", "pos", "log2")
	} else if (flag==2) {
		data = subset(data, data$gene!="Antitarget")
		data = data[,c("chromosome", "start", "log2"),drop=FALSE]
		colnames(data) = c("chrom", "pos", "log2")
	}
	
	'plot_log2_' <- function(x, title = " ")
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
   			index = c(index, seq(from = CytoBand[i, "start"], to=CytoBand[i, "end"], length=sum(x$chrom==i)))
   		}
		plot(index, x$log2, type="p", pch=".", cex=1.95, col="grey80", axes=FALSE, frame=FALSE, xlab="", ylab="", main="", ylim=c(-4.5,4.5))
  		axis(1, at = c(start, end[length(end)]), labels=rep("", length(start)+1), tcl=.5)
  		axis(1, at = .5*(start+end), labels=c(1:22, "X"), lwd=0, lwd.ticks=1, tcl=-.35)
  		axis(2, at = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4), cex.axis = 1, las = 1, lwd = 1.15, lwd.ticks = 0.95, line = -.5)
		mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
		points(c(0-.05*max(index),max(index)+.01*max(index)), c(0,0), type="l", col="black", lwd=1)
		title(main = title, cex.main=.75, font.main=1)
	}
	
	pdf(file=paste0("cnvkit/log2/", opt$sample_name, ".ontarget.pdf"), width=10, height=4.25)
	plot_log2_(x=data, title=opt$sample_name)
	dev.off()
	
} else if (as.numeric(opt$type)==2) {

	suppressPackageStartupMessages(library("copynumber"))
	suppressPackageStartupMessages(library("colorspace"))
	suppressPackageStartupMessages(library("ASCAT"))
	suppressPackageStartupMessages(library("GAP"))

	data = read.csv(file=paste0("cnvkit/cnr/", opt$sample_name, ".cnr"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
	data = subset(data, data[,"depth"]!=0)
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
		data = subset(data, data$gene!="-")
		data = data[,c("chromosome", "start", "log2"),drop=FALSE]
		colnames(data) = c("chrom", "pos", "log2")	
	} else if (flag==2) {
		data = subset(data, data$gene=="Antitarget")
		data = data[,c("chromosome", "start", "log2"),drop=FALSE]
		colnames(data) = c("chrom", "pos", "log2")
	}
	
	'plot_log2_' <- function(x, title = " ")
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
   			index = c(index, seq(from = CytoBand[i, "start"], to=CytoBand[i, "end"], length=sum(x$chrom==i)))
   		}
		plot(index, x$log2, type="p", pch=".", cex=1.95, col="grey80", axes=FALSE, frame=FALSE, xlab="", ylab="", main="", ylim=c(-4.5,4.5))
  		axis(1, at = c(start, end[length(end)]), labels=rep("", length(start)+1), tcl=.5)
  		axis(1, at = .5*(start+end), labels=c(1:22, "X"), lwd=0, lwd.ticks=1, tcl=-.35)
  		axis(2, at = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4), cex.axis = 1, las = 1, lwd = 1.15, lwd.ticks = 0.95, line = -.5)
		mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
		points(c(0-.05*max(index),max(index)+.01*max(index)), c(0,0), type="l", col="black", lwd=1)
		title(main = title, cex.main=.75, font.main=1)
	}
	
	tmp = winsorize(data=data, tau=3.5, k=25, verbose=FALSE, return.outliers=TRUE)
	data[tmp$wins.outliers[,3]!=0,"log2"] = NA
	pdf(file=paste0("cnvkit/log2/", opt$sample_name, ".offtarget.pdf"), width=10, height=4.25)
	plot_log2_(x=data, title=opt$sample_name)
	dev.off()

} else if (as.numeric(opt$type)==3) {

	suppressPackageStartupMessages(library("copynumber"))
	suppressPackageStartupMessages(library("colorspace"))
	suppressPackageStartupMessages(library("ASCAT"))
	suppressPackageStartupMessages(library("GAP"))
	
	data = read.csv(file=paste0("cnvkit/cnr/", opt$sample_name, ".cnr"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
	CN = data[,c("chromosome", "start", "log2"),drop=FALSE]
	colnames(CN) = c("Chromosome", "Position", "Log2Ratio")
	CN[,"Chromosome"] = gsub(pattern="chr", replacement="", x=CN[,"Chromosome"], fixed=TRUE)
	CN[CN[,"Chromosome"]=="X","Chromosome"] = 23
	CN[CN[,"Chromosome"]=="Y","Chromosome"] = 24
	CN[,"Chromosome"] = as.numeric(CN[,"Chromosome"])
	CN[CN[,"Log2Ratio"]<(-4) | CN[,"Log2Ratio"]>(4),"Log2Ratio"] = 0
	CN = subset(CN, CN[,"Chromosome"]<=23)
	tmp = pcf(data=winsorize(data=CN, method="mad", tau=2.5, k=10, verbose=FALSE), kmin = 10, gamma=40, fast=FALSE, verbose=FALSE)[,2:7,drop=FALSE]
	colnames(tmp) = c("Chromosome", "Arm", "Start", "End", "N", "Log2Ratio")
	save(CN, tmp, file=paste0("cnvkit/totalcopy/", opt$sample_name, ".RData"))
	
	'plot_log2_' <- function(x, y, title = "", alpha=NA, psi=NA)
	{
   		par(mar=c(5, 5, 4, 2)+.1)
   		data("CytoBand")
   		end = NULL
		for (j in 1:23) {
			end = c(end, max(CytoBand$End[CytoBand$Chromosome==j]))
		}
		end = cumsum(end)
		start = rep(0, 23)
		start[2:23] = end[1:22]+1
		for (j in 1:23) {
			y[y[,"Chromosome"]==j,"Start"] = y[y[,"Chromosome"]==j,"Start"] + start[j]
			y[y[,"Chromosome"]==j,"End"] = y[y[,"Chromosome"]==j,"End"] + start[j]
			x[x[,"chrom"]==j,"pos"] = x[x[,"chrom"]==j,"pos"] + start[j]
		}
		plot(x[,"pos"], x[,"Log2Ratio"], type="p", pch=".", cex=1, col="grey80", axes=FALSE, frame=FALSE, xlab="", ylab="", main="", ylim=c(-4.5,4.5))
		for (j in 1:nrow(y)) {
 			lines(x=c(y[j,"Start"], y[j,"End"]), y=rep(y[j,"Log2Ratio"],2), lty=1, lwd=1.75, col="red")
 		}
 		axis(1, at = c(start, end[length(end)]), labels=rep("", length(start)+1), tcl=.5)
  		axis(1, at = .5*(start+end), labels=c(1:22, "X"), lwd=0, lwd.ticks=1, tcl=-.35)
  		axis(2, at = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4), cex.axis = 1, las = 1, lwd = 1.15, lwd.ticks = 0.95, line = -.5)
		mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
		points(c(0-.05*max(x[,"pos"]),max(x[,"pos"])+.01*max(x[,"pos"])), c(0,0), type="l", col="black", lwd=1)
		title(main = paste0(title, " | alpha = ", signif(alpha, 3), " | psi = ", signif(psi, 3)), cex.main=.75, font.main=1)
	}
	
	'prune_' <- function(x, n=10)
	{
		cnm = matrix(NA, nrow=nrow(x), ncol=nrow(x))
		for (j in 1:nrow(x)) {
			cnm[,j] = abs(2^x[j,"Log2Ratio"] - 2^x[,"Log2Ratio"])
		}
		cnt = hclust(as.dist(cnm), "average")
		cnc = cutree(tree=cnt, k=n)
		for (j in unique(cnc)) {
			indx = which(cnc==j)
			if (length(indx)>2) {
				mcl = mean(x[indx,"Log2Ratio"])
				scl = sd(x[indx,"Log2Ratio"])
				ind = which(x[indx,"Log2Ratio"]<(mcl+1.96*scl) & x[indx,"Log2Ratio"]>(mcl-1.96*scl))
				x[indx[ind],"Log2Ratio"] = mean(x[indx[ind],"Log2Ratio"])
			} else {
				x[indx,"Log2Ratio"] = mean(x[indx,"Log2Ratio"])
			}
		}
		return(x)
	}

	CN = winsorize(data=CN[,c("Chromosome","Position","Log2Ratio")], tau=2.5, k=15, verbose=FALSE)
	tmp = prune_(x=tmp, n=10)
	
	pdf(file=paste0("cnvkit/segmented/", opt$sample_name, ".pdf"), width=10, height=4.25)
	file_names = dir(path="facets/cncf", pattern=opt$sample_name, full.names=TRUE)
	file_names = file_names[grep(".Rdata", file_names, fixed=TRUE)]
	if (length(file_names)==1) {
		load(file_names)
		alpha = fit$purity
		psi = fit$ploidy
	} else {
		alpha = NA
		psi = NA
	}
	plot_log2_(x=CN, y=tmp, title = opt$sample_name, alpha=alpha, psi=psi)
	dev.off()

} else if (as.numeric(opt$type)==4) {

	suppressPackageStartupMessages(library("copynumber"))
	suppressPackageStartupMessages(library("colorspace"))
	suppressPackageStartupMessages(library("ASCAT"))
	suppressPackageStartupMessages(library("GAP"))
	
	load(paste0("cnvkit/totalcopy/", opt$sample_name, ".RData"))
	file_names = dir(path="facets/cncf", pattern=opt$sample_name, full.names=TRUE)
	file_names = file_names[grep(".Rdata", file_names, fixed=TRUE)]
	if (length(file_names)==1) {
		load(file_names)
		alpha = ifelse(is.na(fit$purity), 1, fit$purity)
		psi = ifelse(is.na(fit$ploidy), 2, fit$ploid)
	} else {
		alpha = 1
		psi = 2
	}
	
	'prune_' <- function(x, n=10)
	{
		cnm = matrix(NA, nrow=nrow(x), ncol=nrow(x))
		for (j in 1:nrow(x)) {
			cnm[,j] = abs(2^x[j,"Log2Ratio"] - 2^x[,"Log2Ratio"])
		}
		cnt = hclust(as.dist(cnm), "average")
		cnc = cutree(tree=cnt, k=n)
		for (j in unique(cnc)) {
			indx = which(cnc==j)
			if (length(indx)>2) {
				mcl = mean(x[indx,"Log2Ratio"])
				scl = sd(x[indx,"Log2Ratio"])
				ind = which(x[indx,"Log2Ratio"]<(mcl+1.96*scl) & x[indx,"Log2Ratio"]>(mcl-1.96*scl))
				x[indx[ind],"Log2Ratio"] = mean(x[indx[ind],"Log2Ratio"])
			} else {
				x[indx,"Log2Ratio"] = mean(x[indx,"Log2Ratio"])
			}
		}
		return(x)
	}
	
	tmp = prune_(x=tmp, n=10)
	qt = round((((2^(tmp[,"Log2Ratio"])) * (alpha*psi + 2*(1-alpha))) - 2*(1-alpha))/alpha)
	qt[is.na(qt)] = 2
	qt[is.infinite(qt)] = 2
	cat5 = rep(0, length(qt))
	if (round(psi)==1 | round(psi)==2) {
		cat5t = c(0, 1, 3, 7)
	} else if (round(psi)==3) {
		cat5t = c(0, 1, 4, 9)
	} else if (round(psi)==4) {
		cat5t = c(0, 1, 5, 10)
	} else if (round(psi)==5) {
		cat5t = c(0, 2, 6, 12)
	} else if (round(psi)>=6) {
		cat5t = c(0, 2, 7, 15)
	} else {
		cat5t = c(0, 1, 3, 7)
	}
	cat5[qt <= cat5t[2]] = -1
	cat5[qt <= cat5t[1]] = -2
	cat5[qt >= cat5t[3]] = 1
	cat5[qt >= cat5t[4]] = 2
	tmp = cbind(tmp, "Cat5"=cat5)
	save(CN, tmp, file=paste0("cnvkit/called/", opt$sample_name, ".RData"))
	
} else if (as.numeric(opt$type)==5) {

	suppressPackageStartupMessages(library("RColorBrewer"))
	suppressPackageStartupMessages(library("GenomicRanges"))
	suppressPackageStartupMessages(library("plyr"))
	suppressPackageStartupMessages(library("dplyr"))
	suppressPackageStartupMessages(library("stringr"))
	suppressPackageStartupMessages(library("tidyr"))
	suppressPackageStartupMessages(library("magrittr"))
	suppressPackageStartupMessages(library("foreach"))
	suppressPackageStartupMessages(library("rtracklayer"))
	suppressPackageStartupMessages(library("grid"))
	suppressPackageStartupMessages(library("rlist"))

	sample_names = unlist(strsplit(opt$sample_name, split=" ", fixed=TRUE))
	genes = read.csv(file="~/share/reference/annotation_gene_lists/annotation_impact_468.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) %>%
			filter(Chromosome %in% as.character(c(1:22, "X", "Y"))) %>%
			filter(!duplicated(Gene_Symbol)) %>%
			arrange(as.integer(Chromosome), Start, End)

	genes_granges = genes %$%
					GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), Gene_Symbol = Gene_Symbol)
	mm = lapply(1:length(sample_names), function(i, sample_names, genes, genes_granges) {
	    load(paste0("cnvkit/called/", sample_names[i], ".RData"))
		tmp[tmp[,"Chromosome"]==23,"Chromosome"] = "X"
		tmp[tmp[,"Chromosome"]==24,"Chromosome"] = "Y"
		tmp_granges = tmp %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End))
		mcols(tmp_granges) = tmp %>% select(Cat5)
		fo = findOverlaps(tmp_granges, genes_granges)
		x = mcols(genes_granges)[subjectHits(fo),]
		y = mcols(tmp_granges)[queryHits(fo),]
		df = data.frame("Gene_Symbol"=x, "Cat5"=y)
		df = df %>%
			 group_by(Gene_Symbol) %>%
			 top_n(1, abs(Cat5))
		z = as.numeric(df$Cat5)
		names(z) = as.character(df$Gene_Symbol)
		z = z[names(z) %in% genes[,1]]
		res = rep(NA, nrow(genes))
		names(res) = genes[,1]
		res[names(z)] = z
		return(res)
	}, sample_names, genes, genes_granges)

	bygene = do.call(cbind, mm)
	colnames(bygene) = sample_names
	bygene = cbind(genes, bygene) %>%
		 	 arrange(as.integer(Chromosome), Start, End)

	save(bygene, file="cnvkit/summary/bygene.RData")
	write.table(bygene, file="cnvkit/summary/bygene.txt", sep="\t", col.names=TRUE, row.names=FALSE, na="", quote=FALSE)

}
