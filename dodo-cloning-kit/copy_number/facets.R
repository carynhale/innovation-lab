#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
				make_option("--option", default = NA, type = 'integer', help = "type of analysis"),
				make_option("--sample_name", default = NA, type = 'character', help = "sample name"),
				make_option("--snp_nbhd", default = 250, type = 'integer', help = "window size"),
				make_option("--pre_cval", default = 50, type = 'integer', help = "pre-processing critical value"),
				make_option("--cval1", default = 150, type = 'integer', help = "critical value for estimating diploid log Ratio"),
				make_option("--cval2", default = 50, type = 'integer', help = "starting critical value for segmentation (increases by 10 until success)"),
				make_option("--max_cval", default = 5000, type = 'integer', help = "maximum critical value for segmentation (increases by 10 until success)"),
				make_option("--min_nhet", default = 25, type = 'integer', help = "minimum number of heterozygote snps in a segment used for bivariate t-statistic during clustering of segment"),
				make_option("--het_threshold", default = 0.25, type = 'double', help = "AF threshold for heterozygous SNPs"),
				make_option("--diplogr", default = NULL, type = 'double', help = "override diploid log-ratio"),
				make_option("--ndepth_max", default = 1000, type = 'integer', help = "normal depth max"),
				make_option("--use_emcncf2", default = F, action = 'store_true', help = "use emcncf version 2"),
				make_option("--gene_loc_file", default = '~/share/reference/IMPACT410_genes_for_copynumber.txt', type = 'character', help = "file containing gene locations"),
				make_option("--genome", default = 'b37', type = 'character', help = "genome of counts file"),
				make_option("--out_prefix", default = NULL, help = "output prefix"),
				make_option("--centromereFile", default = NULL, type = "character", action = "store", help ="centromere file"),
				make_option("--pqLine", default = F, action = "store_true", help = "draw pq centromere line"),
				make_option("--outFile", default = NULL, help = "output file"),
				make_option("--mysqlHost", default = '10.0.200.48', help = "MySQL server hostname"),
				make_option("--mysqlPort", default = 38493, help = "MySQL server port"),
				make_option("--mysqlUser", default = 'embl', help = "MySQL server username"),
				make_option("--mysqlPassword", default = NULL, help = "MySQL server password"),
				make_option("--mysqlDb", default = 'homo_sapiens_core_75_37', help = "MySQL server database"),
				make_option("--genesFile", default = NULL, help = "list of genes to include (hgnc symbols)"),
				make_option("--annotFile", default = "~/share/reference/annotation_gene_lists/geneCN.txt", help = "file with annotations to replace MySQL sever query"),
				make_option("--includeChrY", action="store_true", default=F, help="Include Chromosome Y (drop by default)"),
				make_option("--sampleColumnPostFix", default="_AWN_LRR_threshold$", help="Postfix of columns that represent samples"),
				make_option("--summaryConfig", default=NULL, help="Summary config yaml for ordering samples and replacing IDs")
			   )

parser <- OptionParser(usage = "%prog [options] [tumor-normal base counts file]", option_list = optList)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options


if (as.numeric(opt$option)==1) {

	suppressPackageStartupMessages(library("RColorBrewer"))
	suppressPackageStartupMessages(library("plyr"))
	suppressPackageStartupMessages(library("dplyr"))
	suppressPackageStartupMessages(library("tidyr"))
	suppressPackageStartupMessages(library("stringr"))
	suppressPackageStartupMessages(library("magrittr"))
	suppressPackageStartupMessages(library("facets"))
	suppressPackageStartupMessages(library("pctGCdata"))
	suppressPackageStartupMessages(library("foreach"))
	suppressPackageStartupMessages(library("GAP"))

	snpPileupFile = arguments$args[1]
	tumorName = snpPileupFile %>% sub('.*/', '', .) %>% sub('_.*', '', .)
	normalName = snpPileupFile %>% sub('.*/', '', .) %>% sub('^.*_', '', .) %>% sub('\\..*', '', .)
	switch(opt$genome,
       b37={
            facetsGenome = 'hg19'
       },
       GRCh37={
            facetsGenome = 'hg19'
       },
       hg19={
            facetsGenome = 'hg19'
       },
       mm9={
            facetsGenome = 'mm9'
       },
       mm10={
            facetsGenome = 'mm10'
       },
       GRCm38={
            facetsGenome = 'mm10'
       },
       {
           stop(paste("Invalid Genome",opt$genome))
       })
	buildData = installed.packages()["facets",]
	cat("#Module Info\n")
	for (fi in c("Package","LibPath","Version","Built")) {
    	cat("#",paste(fi,":",sep=""),buildData[fi],"\n")
	}
	version = buildData["Version"]
	cat("\n")

	set.seed(1)
	snpmat = readSnpMatrix(snpPileupFile)
	preOut = snpmat %>%
			 preProcSample(snp.nbhd = opt$snp_nbhd,
			 			   het.thresh = opt$het_threshold,
			 			   cval = opt$pre_cval,
			 			   gbuild = facetsGenome,
			 			   ndepthmax = opt$ndepth_max)
	if (!is.null(opt$diplogr)) {
	    cval = opt$cval2
	    out2 = preOut %>%
	    	   procSample(cval = cval,
	    	   			  min.nhet = opt$min_nhet,
	    	   			  dipLogR = opt$diplogr)
	    if (opt$use_emcncf2) {
	        fit = out2 %>% emcncf2
	    } else {
	        fit = out2 %>% emcncf
	    }
	} else {
	    out1 = preOut %>%
	    	   procSample(cval = opt$cval1,
	    	   			  min.nhet = opt$min_nhet)
		cval = opt$cval2
	    success = FALSE
	    while (!success && cval < opt$max_cval) {
 	       out2 = preOut %>%
 	       		  procSample(cval = cval,
 	       		  			 min.nhet = opt$min_nhet,
 	       		  			 dipLogR = out1$dipLogR)
 	       print(str_c("attempting to run emncf() with cval2 = ", cval))
 	       fit = tryCatch({
 	           if (opt$use_emcncf2) {
 	               out2 %>% emcncf2
 	           } else {
 	               out2 %>% emcncf
 	           }
 	       }, error = function(e) {
 	           print(paste("Error:", e))
 	           return(NULL)
 	       })
 	       if (!is.null(fit)) {
 	           success = TRUE
 	       } else {
 	           cval = cval + 100
 	       }
 	   }
 	   if (!success) {
 	       stop("Failed to segment data\n")
 	   }
	}

	'formatSegmentOutput' <- function(out,sampID) {
		seg = list()
		seg$ID = rep(sampID,nrow(out$out))
		seg$chrom = out$out$chr
		seg$loc.start = rep(NA,length(seg$ID))
		seg$loc.end = seg$loc.start
		seg$num.mark = out$out$num.mark
		seg$seg.mean = out$out$cnlr.median
		for (i in 1:nrow(out$out)) {
			lims = range(out$jointseg$maploc[(out$jointseg$chrom==out$out$chr[i] & out$jointseg$seg==out$out$seg[i])],na.rm=T)
			seg$loc.start[i] = lims[1]
			seg$loc.end[i] = lims[2]
		}
		return(invisible(as.data.frame(seg)))
	}
	id = paste(tumorName, normalName, sep = '_')
	out2$IGV = formatSegmentOutput(out2, id)
	save(out2, fit, file = str_c(opt$out_prefix, ".Rdata"), compress=TRUE)

	ff = str_c(opt$out_prefix, ".out")
	cat("# Version =", version, "\n", file = ff, append = TRUE)
	cat("# Input =", basename(snpPileupFile), "\n", file = ff, append = TRUE)
	cat("# tumor =", tumorName, "\n", file = ff, append = TRUE)
	cat("# normal =", normalName, "\n", file = ff, append = TRUE)
	cat("# snp.nbhd =", opt$snp_nbhd, "\n", file = ff, append = TRUE)
	cat("# cval1 =", opt$cval1, "\n", file = ff, append = TRUE)
	cat("# cval2 =", cval, "\n", file = ff, append = TRUE)
	cat("# min.nhet =", opt$min_nhet, "\n", file = ff, append = TRUE)
	cat("# genome =", opt$genome, "\n", file = ff, append = TRUE)
	cat("# Purity =", fit$purity, "\n", file = ff, append = TRUE)
	cat("# Ploidy =", fit$ploidy, "\n", file = ff, append = TRUE)
	cat("# dipLogR =", fit$dipLogR, "\n", file = ff, append = TRUE)
	cat("# dipt =", fit$dipt, "\n", file = ff, append = TRUE)
	cat("# loglik =", fit$loglik, "\n", file = ff, append = TRUE)
	tab = cbind(select(out2$IGV, ID:num.mark), select(fit$cncf, -start, -end, -chrom, -num.mark))
	write.table(tab, str_c(opt$out_prefix, ".txt"), row.names = F, quote = F, sep = '\t')
	write.table(out2$IGV, str_c(opt$out_prefix, '.seg'), row.names = F, quote = F, sep = '\t')

} else if (as.numeric(opt$option)==2) {

	suppressPackageStartupMessages(library("RColorBrewer"))
	suppressPackageStartupMessages(library("plyr"))
	suppressPackageStartupMessages(library("dplyr"))
	suppressPackageStartupMessages(library("tidyr"))
	suppressPackageStartupMessages(library("stringr"))
	suppressPackageStartupMessages(library("magrittr"))
	suppressPackageStartupMessages(library("facets"))
	suppressPackageStartupMessages(library("pctGCdata"))
	suppressPackageStartupMessages(library("foreach"))
	suppressPackageStartupMessages(library("GAP"))

	'plot_log2_' <- function(x, y, n=10, purity=NA, ploidy=NA, title = "")
	{
		cna = x$jointseg %>%
			  select(chrom, pos = maploc, log2 = cnlr)
		seg = y$cncf %>%
			  select(chrom, start = start, end = end, log2 = cnlr.median, n=num.mark)
		seg = prune_(x=seg, n) %>%
			  mutate(n = cumsum(n))
		purity = ifelse(is.na(purity), 1, purity)
		ploidy = ifelse(is.na(ploidy), 2, ploidy)
		data(CytoBand)
   		par(mar=c(5, 5, 4, 2)+.1)
		end = NULL
   		for (i in 1:23) {
   			end = c(end, max(CytoBand[CytoBand[,1]==i,"End"]))
   		}
   		end = cumsum(end)
   		start = c(1, end[1:22]+1)
   		CytoBand = cbind(start, end)
   		index = NULL
   		for (i in 1:23) {
   			index = c(index, seq(from = CytoBand[i, "start"], to=CytoBand[i, "end"], length=sum(cna$chrom==i)))
   		}
		plot(index, cna$log2, type="p", pch=".", cex=1.95, col="grey80", axes=FALSE, frame=FALSE, xlab="", ylab="", main="", ylim=c(-4.5,4.5))
 		for (j in 1:nrow(seg)) {
 			if (j == 1) {
 				lines(x=c(1, index[seg[j,"n"]]), y=rep(seg[j,"log2"],2), lty=1, lwd=2.75, col="red")
 			} else {
 				lines(x=c(index[seg[j-1,"n"]], index[seg[j,"n"]]), y=rep(seg[j,"log2"],2), lty=1, lwd=2.75, col="red")
 			}
  		}
  		axis(side=1, at=c(CytoBand[,"start"],CytoBand[nrow(CytoBand),"end"]), labels=rep("", nrow(CytoBand)+1), tcl=.5)
		axis(side=1, at=apply(CytoBand[,c("start", "end"),drop=FALSE], 1, mean), labels=c(1:22, "X"), lwd=0, lwd.ticks=1, tcl=-.35)
  		axis(2, at = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4), cex.axis = 1, las = 1, lwd = 1.15, lwd.ticks = 0.95, line = -.5)
		mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
		points(c(0-.05*max(index),max(index)+.01*max(index)), c(0,0), type="l", col="black", lwd=1)
		title(main = paste0(title, " | alpha = ", signif(purity, 3), " | psi = ", signif(ploidy, 3)), cex.main=.75, font.main=1)
	}

	'plot_cncf_' <- function(x, emfit = NULL, clustered = FALSE, plot.type = c("em","naive","both","none"), sname = NULL)
	{
	    def.par = par(no.readonly = TRUE)
	    plot.type = match.arg(plot.type)
	    if (plot.type=="none") {
	    	layout(matrix(1:2, ncol=1))
	    }
	    if (plot.type=="em") {
	    	layout(matrix(rep(1:4, c(9,9,6,1)), ncol=1))
	    }
	    if (plot.type=="naive") {
	    	layout(matrix(rep(1:4, c(9,9,6,1)), ncol=1))
	    }
	    if (plot.type=="both") {
	    	layout(matrix(rep(1:6, c(9,9,6,1,6,1)), ncol=1))
	    }
	    par(mar=c(0.25+1,5,0.25+1,1), mgp=c(1.75, 0.6, 0), oma=c(3,0,1.25,0))
	    jseg = x$jointseg
	    chrbdry = which(diff(jseg$chrom) != 0)
	    if (missing(emfit)) {
    	    out = x$out
    	    if (plot.type=="em" | plot.type=="both") {
    	        warning("emfit is missing; plot.type set to naive")
    	        plot.type = "naive"
    	    }
    	} else {
    	    out = emfit$cncf
    	    out$tcn = x$out$tcn
    	    out$lcn = x$out$lcn
    	    out$cf = x$out$cf
    	}
    	if (clustered) {
    	    cnlr.median = out$cnlr.median.clust
    	    mafR = out$mafR.clust
    	    mafR[is.na(mafR)] = out$mafR[is.na(mafR)]
    	} else {
    	    cnlr.median = out$cnlr.median
    	    mafR = out$mafR
    	}
    	mafR = abs(mafR)
    	chrcol = 1+rep(out$chr-2*floor(out$chr/2), out$num.mark)
    	nn = cumsum(table(jseg$chrom[is.finite(jseg$cnlr)]))
    	segbdry = cumsum(c(0,out$num.mark))
    	segstart = segbdry[-length(segbdry)]
    	segend = segbdry[-1]
    	
    	plot(jseg$cnlr[is.finite(jseg$cnlr)], type="p", pch=".", cex=4, col = c("grey","lightblue","azure4","slateblue")[chrcol], axes=FALSE, frame=FALSE, xlab="", ylab="", main="", ylim=c(-4.5,4.5))
  		axis(2, at = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4), cex.axis = 1.25, cex.lab = 1.25, las = 1, lwd = 1.15, lwd.ticks = 0.95)
    	mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 2.75, cex = .85)
    	abline(v=c(1, chrbdry), lwd=1, lty=3, col="brown")
    	abline(v=length(jseg$cnlr[is.finite(jseg$cnlr)]), lwd=1, lty=3, col="brown")
    	abline(h=median(jseg$cnlr, na.rm=TRUE), col = "green2")
    	abline(h=x$dipLogR, col = "magenta4")
    	segments(segstart, cnlr.median, segend, cnlr.median, lwd=2, col = 2)
		chromlevels = x$chromlevels
    	if (is.null(chromlevels)) {
    		chromlevels = 1:length(nn)
    	}
    	axis(side=1, at=(nn+c(0,nn[-length(nn)]))/2, labels=rep("", length(chromlevels)))
		box(lwd=2)
    	
  	    plot(jseg$valor[is.finite(jseg$cnlr)], type="p", pch=".", cex=4, col = c("grey","lightblue","azure4","slateblue")[chrcol], axes=FALSE, frame=FALSE, xlab="", ylab="", main="", ylim=c(-4.5,4.5))
  	    axis(2, at = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4), cex.axis = 1.25, cex.lab = 1.25, las = 1, lwd = 1.15, lwd.ticks = 0.95)
    	mtext(side = 2, text = expression(Log~"OR"), line = 2.75, cex = .85)
    	abline(v=c(1, chrbdry), lwd=1, lty=3, col="brown")
    	abline(v=length(jseg$cnlr[is.finite(jseg$cnlr)]), lwd=1, lty=3, col="brown")
    	abline(h=0, lwd=1, lty=1, col="black")
    	segments(segstart, sqrt(mafR), segend, sqrt(mafR), lwd=2.15, col = 2)
    	segments(segstart, -sqrt(mafR), segend, -sqrt(mafR), lwd=2, col = 2)
    	chromlevels = x$chromlevels
    	if (is.null(chromlevels)) {
    		chromlevels = 1:length(nn)
    	}
    	axis(side=1, at=(nn+c(0,nn[-length(nn)]))/2, labels=rep("", length(chromlevels)))
		box(lwd=2)
  	    
	    cfpalette = c(colorRampPalette(c("white", "steelblue"))(10),"bisque2")
	    if (plot.type=="naive" | plot.type=="both") {
		    out$tcn[out$tcn > 10] = 9 + log10(out$tcn[out$tcn > 10])
    	    ii = which(out$lcn > 5)
    	    if (length(ii)>0) {
    	    	out$lcn[ii] = 5 + log10(out$lcn[ii])
    	    }
    	    plot(c(0,length(jseg$cnlr)), c(0,max(out$tcn)), type="n", axes=FALSE, frame=FALSE, xlab="", ylab="", main="", ylim=c(0,8))
        	axis(2, at = c(0, 2, 4, 6, 8), labels = c(0, 2, 4, 6, 8), cex.axis = 1.25, cex.lab = 1.25, las = 1, lwd = 1.15, lwd.ticks = 0.95)
        	axis(2, at = c(1, 3, 5, 7), labels=rep("", 4), lwd=0, lwd.ticks=1, tcl=-.35)
    		mtext(side = 2, text = expression("CN"), line = 2.75, cex = .85, srt=90)
    		abline(v=c(1, chrbdry), lwd=1, lty=3, col="brown")
    		abline(v=length(jseg$cnlr[is.finite(jseg$cnlr)]), lwd=1, lty=3, col="brown")
        	segments(segstart, out$lcn, segend, out$lcn, lwd=1.75, col=2)
        	segments(segstart, out$tcn, segend, out$tcn, lwd=1.75, col=1)
        	chromlevels = x$chromlevels
	    	if (is.null(chromlevels)) {
	    		chromlevels = 1:length(nn)
	    	}
	    	axis(side=1, at=(nn+c(0,nn[-length(nn)]))/2, labels=rep("", length(chromlevels)))
			box(lwd=2)
    	}
    	if (plot.type=="em" | plot.type=="both") {
        	out$tcn.em[out$tcn.em > 10] = 9 + log10(out$tcn.em[out$tcn.em > 10])
        	ii = which(out$lcn.em > 5)
        	if (length(ii)>0) {
        		out$lcn.em[ii] = 5 + log10(out$lcn.em[ii])
        	}
    	    plot(c(0,length(jseg$cnlr)), c(0,max(out$tcn)), type="n", axes=FALSE, frame=FALSE, xlab="", ylab="", main="", ylim=c(0,8))
        	axis(2, at = c(0, 2, 4, 6, 8), labels = c(0, 2, 4, 6, 8), cex.axis = 1.25, cex.lab = 1.25, las = 1, lwd = 1.15, lwd.ticks = 0.95)
    		mtext(side = 2, text = expression("CN"), line = 2.75, cex = .85, srt=90)
    		abline(v=c(1, chrbdry), lwd=1, lty=3, col="brown")
    		abline(v=length(jseg$cnlr[is.finite(jseg$cnlr)]), lwd=1, lty=3, col="brown")
        	segments(segstart, out$lcn.em, segend, out$lcn.em, lwd=1.75, col=2)
        	segments(segstart, out$tcn.em, segend, out$tcn.em, lwd=1.75, col=1)
        	chromlevels = x$chromlevels
    		if (is.null(chromlevels)) {
    			chromlevels = 1:length(nn)
    		}
    		axis(side=1, at=(nn+c(0,nn[-length(nn)]))/2, labels=rep(" ", length(chromlevels)), cex.axis=1.25, cex.lab=1.25, line=0)
    		axis(side=1, at=(nn+c(0,nn[-length(nn)]))/2, labels=chromlevels, cex.axis=1.25, cex.lab=1.25, line=.65, lwd=-1)
			box(lwd=2)
    	}
    	par(def.par)
	}

	'psi' <- function (x, z)
	{
	    xwin = x
	    xwin[x < -z] = -z
	    xwin[x > z] = z
	    return(xwin)
	}

	'numericChrom' <- function (chrom) 
	{
	    if (!is.numeric(chrom)) {
	        if (is.factor(chrom)) {
	            chrom = as.character(chrom)
	        }
	        chrx = c(which(chrom == "x"), which(chrom == "X"))
	        chrom[chrx] = 23
	        chry = c(which(chrom == "y"), which(chrom == "Y"))
	        chrom[chry] = 24
	        chrom = as.numeric(chrom)
	    }
	    return(chrom)
	}

	'numericArms' <- function (chrom, char.arms)
	{
	    p.arm = which(char.arms == "p")
	    q.arm = which(char.arms == "q")
	    arms = rep(NA, length(char.arms))
	    arms[p.arm] = chrom[p.arm] * 2 - 1
	    arms[q.arm] = chrom[q.arm] * 2
	    return(arms)
	}

	'madWins' <- function (x, tau, k, digits) 
	{
	    xhat = medianFilter(x, k)
	    d = x - xhat
	    SD = mad(d)
	    z = tau * SD
	    xwin = xhat + psi(d, z)
	    outliers = rep(0, length(x))
	    outliers[round(x, digits) > round(xwin, digits)] = 1
	    outliers[round(x, digits) < round(xwin, digits)] = -1
	    return(list(ywin = xwin, sdev = SD, outliers = outliers))
	}

	'medianFilter' <- function (x, k) 
	{
	    n = length(x)
	    filtWidth = 2 * k + 1
	    if (filtWidth > n) {
	        if (n == 0) {
	            filtWidth = 1
	        } else if (n%%2 == 0) {
	            filtWidth = n - 1
	        } else {
            	filtWidth = n
        	}
    	}
    	runMedian = runmed(x, k = filtWidth, endrule = "median")
    	return(runMedian)
	}


	'winsorize' <- function (data, pos.unit = "bp", arms = NULL, method = "mad",  tau = 2.5, k = 25, gamma = 40, iter = 1, assembly = "hg19", 
    						 digits = 4, return.outliers = FALSE, save.res = FALSE, file.names = NULL, verbose = TRUE) 
	{
    	stopifnot(pos.unit %in% c("bp", "kbp", "mbp"))
    	stopifnot(method %in% c("mad", "pcf"))
    	if (!assembly %in% c("hg19", "hg18", "hg17", "hg16", "mm7", "mm8", "mm9")) {
    	    stop("assembly must be one of hg19, hg18, hg17 or hg16", call. = FALSE)
    	}
    	stopifnot(class(data) %in% c("matrix", "data.frame", "character"))
    	isfile <- class(data) == "character"
    	if (!isfile) {
    	    stopifnot(ncol(data) >= 3)
    	    chrom = data[, 1]
    	    pos = data[, 2]
    	    nSample = ncol(data) - 2
    	    sample.names = colnames(data)[-c(1:2)]
    	} else {
	        f = file(data, "r")
	        head = scan(f, nlines = 1, what = "character", quiet = TRUE, sep = "\t")
	        if (length(head) < 3) {
	            stop("Data in file must have at least 3 columns", call. = FALSE)
        	}
        	sample.names = head[-c(1:2)]
        	nSample = length(sample.names)
        	chrom.pos = read.table(file = data, sep = "\t", header = TRUE, colClasses = c(rep(NA, 2), rep("NULL", nSample)), as.is = TRUE)
	        chrom = chrom.pos[, 1]
    	    pos = chrom.pos[, 2]
    	}
    	if (is.factor(chrom)) {
    	    chrom = as.character(chrom)
    	}
    	num.chrom = numericChrom(chrom)
    	nProbe = length(num.chrom)
    	if (!is.numeric(pos)) {
    	    stop("input in data column 2 (posistions) must be numeric", call. = FALSE)
    	}
    	if (is.null(arms)) {
    	    arms = getArms(num.chrom, pos, pos.unit, get(assembly))
    	} else {
        	if (length(arms) != nProbe) {
            	stop("'arms' must be the same length as number of rows in data", call. = FALSE)
        	}
    	}
    	num.arms = numericArms(num.chrom, arms)
    	arm.list = unique(num.arms)
    	nArm <- length(arm.list)
    	if (!save.res) {
        	wins.data = matrix(nrow = 0, ncol = nSample)
        	if (return.outliers) {
            	wins.outliers = matrix(nrow = 0, ncol = nSample)
        	}
    	} else {
        if (is.null(file.names)) {
            dir.res = "Wins_res"
            if (!dir.res %in% dir()) {
                dir.create(dir.res)
            }
            file.names = c(paste(dir.res, "/", "wins.data.txt", sep = ""), paste(dir.res, "/", "wins.outliers.txt", sep = ""))
        } else {
            if (length(file.names) < 2) {
                stop("'file.names' must be of length 2", call. = FALSE)
            }
        }
    }
    for (c in 1:nArm) {
        probe.c = which(num.arms == arm.list[c])
        wins.data.c = matrix(nrow = length(probe.c), ncol = 0)
        if (return.outliers || save.res) {
            wins.outliers.c = matrix(nrow = length(probe.c), ncol = 0)
        }
        if (!isfile) {
            arm.data = data[probe.c, -c(1:2), drop = FALSE]
        } else {
            arm.data = read.table(f, nrows = length(probe.c), sep = "\t", colClasses = c(rep("NULL", 2), rep("numeric", nSample)))
        }
        if (any(!sapply(arm.data, is.numeric))) {
            stop("input in data columns 3 and onwards (copy numbers) must be numeric", call. = FALSE)
        }
        for (i in 1:nSample) {
            y = arm.data[, i]
            na = is.na(y)
            use.y = y[!na]
            ywins = rep(NA, length(y))
            outliers = rep(NA, length(y))
            wins = switch(method, mad = madWins(use.y, tau = tau, k = k, digits = digits), pcf = pcfWins(use.y, tau = tau, k = k, gamma = gamma, iter = iter, digits = digits))
            ywins[!na] = wins$ywin
            outliers[!na] = wins$outliers
            ywins = round(ywins, digits = digits)
            wins.data.c = cbind(wins.data.c, ywins)
            if (return.outliers || save.res) {
                wins.outliers.c = cbind(wins.outliers.c, outliers)
            }
        }
        if (!save.res) {
            wins.data = rbind(wins.data, wins.data.c)
            if (return.outliers) {
                wins.outliers = rbind(wins.outliers, wins.outliers.c)
            }
        } else {
            if (c == 1) {
                wd = file(file.names[1], "w")
                wo = file(file.names[2], "w")
            }
            write.table(data.frame(chrom[probe.c], pos[probe.c], 
                wins.data.c, stringsAsFactors = FALSE), file = wd, 
                col.names = if (c == 1) 
                  c("chrom", "pos", sample.names)
                else FALSE, row.names = FALSE, quote = FALSE, 
                sep = "\t")
            write.table(data.frame(chrom[probe.c], pos[probe.c], 
                wins.outliers.c, stringsAsFactors = FALSE), file = wo, 
                col.names = if (c == 1) 
                  c("chrom", "pos", sample.names)
                else FALSE, row.names = FALSE, quote = FALSE, 
                sep = "\t")
        }
        if (verbose) {
            chr = unique(chrom[probe.c])
            a = unique(arms[probe.c])
            cat(paste("winsorize finished for chromosome arm ", chr, a, sep = ""), "\n")
        }
    }
    if (isfile) {
        close(f)
    }
    if (!save.res) {
        wins.data = data.frame(chrom, pos, wins.data, stringsAsFactors = FALSE)
        colnames(wins.data) = c("chrom", "pos", sample.names)
        if (return.outliers) {
            wins.outliers = data.frame(chrom, pos, wins.outliers, stringsAsFactors = FALSE)
            colnames(wins.outliers) = c("chrom", "pos", sample.names)
            return(list(wins.data = wins.data, wins.outliers = wins.outliers))
        } else {
            return(wins.data)
        }
    } else {
        	close(wd)
        	close(wo)
        	cat(paste("winsorized data were saved in file", file.names[1]), sep = "\n")
        	cat(paste("outliers were saved in file", file.names[2]), sep = "\n")
    	    return(invisible(NULL))
    	}
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
	
	load(paste0("facets/cncf/", opt$sample_name, ".Rdata"))

	if (!dir.exists("facets/plots")) {
		dir.create("facets/plots")
	}
	if (!dir.exists("facets/plots/log2")) {
		dir.create("facets/plots/log2")
	}
	if (!dir.exists("facets/plots/cncf")) {
		dir.create("facets/plots/cncf")
	}
	if (!dir.exists("facets/plots/bychr")) {
		dir.create("facets/plots/bychr")
	}
	tumorName <- strsplit(opt$sample_name, "_")[[1]][1]
	normalName <- strsplit(opt$sample_name, "_")[[1]][2]

	pdf(file = paste0("facets/plots/log2/", opt$sample_name, ".pdf"), width=10, height=4.25)
	plot_log2_(x=out2, y=fit, purity=fit$purity, ploidy=fit$ploidy, title = opt$sample_name)
	dev.off()

	pdf(file = paste0("facets/plots/cncf/", opt$sample_name, ".pdf"), width=10, height=9)
	plot_cncf_(out2, fit)
	dev.off()

	df <- left_join(out2$jointseg, out2$out) %>%
		  mutate(chrom = as.character(chrom)) %>%
		  mutate(chrom = ifelse(chrom=="23", "X", chrom)) %>%
		  mutate(chrom = ifelse(chrom=="24", "Y", chrom))

	for (chr in unique(df$chrom)) {
	     chdf <- filter(df, chrom == chr)
	     if (nrow(chdf) > 0) {
	     	
	     	if (!dir.exists(paste0("facets/plots/bychr/", opt$sample_name))) {
	     		dir.create(paste0("facets/plots/bychr/", opt$sample_name))
	     	}
	        pdf(paste0("facets/plots/bychr/", opt$sample_name, "/chromosome_", chr, ".pdf"), height = 5, width = 5)
	        par(mar=c(5, 5, 4, 2)+.1)
	        plot(chdf$cnlr, type="p", pch=20, col="grey80", xlab="", ylab="", main = "", ylim=c(-4,5), axes=FALSE, frame=TRUE)
	        points(chdf$cnlr.median.clust, pch=20, col="red")
	        abline(h=0, col="goldenrod3", lty=1, lwd=1)
	        axis(2, at = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4), cex.axis = 1, las = 1)
			mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
			axis(1, at = NULL, labels = NULL, cex.axis = 0.85, las = 1)
	    	rect(xleft=1-1e10, xright=nrow(chdf)+10e3, ybottom=4, ytop=6, col="lightgrey", border="black", lwd=1.5)
			title(main = opt$sample_name, line=-1.5, cex.main=.75, font.main=1)
	    	box(lwd=1.5)

	        if (!is.null(opt$centromereFile)) {
	            cen <- read.table(opt$centromereFile, sep = '\t')
	            for (j in unique(cen[,1])) {
	                pos <- cen[which(cen[,1]==j)[1],3]
	                index <- which(chdf$chrom == j & chdf$maploc > pos)[1]
	                if (opt$pqLine && !is.na(index)) {
	                    abline(v = index, col = "darkgrey", lty = 3)
	                }
	            }
	        }
	        dev.off()
	    }
	}

} else if (as.numeric(opt$option)==3) {
	
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
	suppressPackageStartupMessages(library("RMySQL"))
	
	facetsFiles <- arguments$args
	
	if (is.null(opt$annotFile)) {
		connect <- function() dbConnect(MySQL(), host = opt$mysqlHost, port = opt$mysqlPort, user = opt$mysqlUser, password = opt$mysqlPassword, dbname = opt$mysqlDb)
		cat('Connecting to ensembl ... ')
		mydb <- connect()
		on.exit(dbDisconnect(mydb))
		query <- "select r.name as chrom,
		g.seq_region_start as start,
		g.seq_region_end as end,
		x.display_label as hgnc,
		k.band as band
		from gene as g
		join seq_region as r on g.seq_region_id = r.seq_region_id
		join xref as x on g.display_xref_id = x.xref_id
		left join karyotype k on g.seq_region_id = k.seq_region_id
		and ((g.seq_region_start >= k.seq_region_start and g.seq_region_start <= k.seq_region_end)
		or (g.seq_region_end >= k.seq_region_start and g.seq_region_end <= k.seq_region_end))
		where x.external_db_id = 1100;"
		repeat {
			rs <- try(dbSendQuery(mydb, query), silent = T)
			if (is(rs, "try-error")) {
				cat("Lost connection to mysql db ... ")
				mydb <- connect()
				cat("reconnected\n")
			} else {
				break
			}
		}
		genes <- dbFetch(rs, -1)
	} else {
		genes = read.csv(file=opt$annotFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
	}
	cat(paste("Found", nrow(genes), "records\n"))

	genes %<>% filter(chrom %in% as.character(c(1:22, "X", "Y"))) %>%
			   filter(!duplicated(hgnc)) %>%
			   arrange(as.integer(chrom), start, end)

	if (!is.null(opt$genesFile)) {
		g <- scan(opt$genesFile, what = 'character')
		genes %<>% filter(hgnc %in% g)
		absentGenes <- g[!g %in% genes$hgnc]
		if (length(absentGenes) > 0) {
			print("Unable to find", length(absentGenes), "in database\n");
			cat(absentGenes, sep = '\n');
		}
	}

	cat(paste("Filtering to", nrow(genes), "records\n"))

	genesGR <- genes %$% GRanges(seqnames = chrom, ranges = IRanges(start, end), band = band, hgnc = hgnc)
			
	mm <- lapply(facetsFiles, function(f) {
	    load(f)
	    tab <- fit$cncf
		tab$chrom[which(tab$chrom==23)] <- "X"
		tab$chrom[which(tab$chrom==24)] <- "Y"
	
		tabGR <- tab %$% GRanges(seqnames = chrom, ranges = IRanges(start, end))
		mcols(tabGR) <- tab %>% select(cnlr.median:lcn.em)

		fo <- findOverlaps(tabGR, genesGR)

		df <- as.data.frame(cbind(mcols(genesGR)[subjectHits(fo),], mcols(tabGR)[queryHits(fo),]))
		df %<>% group_by(hgnc) %>% top_n(1, abs(cnlr.median))

		ploidy <- table(df$tcn.em)
		ploidy <- as.numeric(names(ploidy)[which.max(ploidy)])

		df$GL <- 0
		df$GL[df$tcn.em < ploidy] <- -1
		df$GL[df$tcn.em == 0] <- -2
		df$GL[df$tcn.em > ploidy] <- 1
		df$GL[df$tcn.em >= ploidy + 4] <- 2

		load(f)
		noise <- median(abs(out2$jointseg$cnlr-  unlist(apply(out2$out[,c("cnlr.median", "num.mark")], 1, function(x) {rep(x[1], each=x[2])}))))

		lrr <- sort(out2$jointseg$cnlr)
		if (noise <= 0.2) {
			lrr <- lrr[round(0.25*length(lrr)):round(0.75*length(lrr))]
		} else if ( noise <= 0.3 ) {
			lrr <- lrr[round(0.275*length(lrr)):round(0.725*length(lrr))]
		} else {
			lrr <- lrr[round(0.3*length(lrr)):round(0.7*length(lrr))]
		}

		df$GL2 <- 0
		df$GL2[df$cnlr.median < median(lrr)-(2.5*sd(lrr))] <- -1
		df$GL2[df$cnlr.median < median(lrr)-(7*sd(lrr))] <- -2
		df$GL2[df$cnlr.median > median(lrr)+(2*sd(lrr))] <- 1
		df$GL2[df$cnlr.median > median(lrr)+(6*sd(lrr))] <- 2
		df %>% select(hgnc, GL, GL2) %>% ungroup
	})
	names(mm) <- facetsFiles
	for (f in facetsFiles) {
		n <- sub('\\..*', '', sub('.*/', '', f))
		colnames(mm[[f]])[2:3] <- paste(n, c("EM", "LRR_threshold"), sep="_")
	}

	mm <- left_join(genes, join_all(mm, type = 'full', by="hgnc")) %>% arrange(as.integer(chrom), start, end)
	write.table(mm, file=opt$outFile, sep="\t", row.names=F, na="", quote=F)
	
} else if (as.numeric(opt$option)==4) {

	geneCN <- arguments$args[1]
    outFile <- arguments$args[2]

	if ("showtext" %in% rownames(installed.packages())) {
		suppressPackageStartupMessages(library("showtext"))
		font.add("DejaVuSansMono", "DejaVuSansMono.ttf")
		showtext.auto()
		fontfamily = "DejaVuSansMono"
	} else {
		fontfamily = "serif"
	}

	'plot_heatmap_' <- function(facets_tab,
								plot_file,
								sample_column_postfix,
								fontfamily,
								sample_names=NULL,
								summaryConfig=NULL,
								col=c("blue", "lightblue", "white", "darksalmon", "red"),
								zlim=c(-2,2)) {
		mm = facets_tab
		if (is.null(sample_names)) {
			sample_names = list(sort(colnames(mm)[sapply(colnames(mm), function(x) {grepl(paste(sample_column_postfix,"$",sep=""), x)})]))
		}
		if (!is.null(summaryConfig) && ("sample_order" %in% names(summaryConfig))) {
			sample_names[[1]] = sample_names[[1]][charmatch(summaryConfig$sample_order, sample_names[[1]])]
		}
		chrsep = cumsum(rle(mm$chrom)$lengths)
		chrmid = c(0,chrsep[-length(chrsep)]) + (rle(mm$chrom)$lengths/2)
		pdf(plot_file, width=20, height=5 + .8*length(sample_names[[1]]))
		par(mfrow=c(length(sample_names),1), mar=c(8,1*(max(sapply(sample_names,nchar))-nchar(sample_column_postfix)),1,2))
		lapply(sample_names, function(x, mm) {
			mm2 = mm[,rev(x)]
        	image(as.matrix(mm2), col=col, xaxt='n', yaxt='n', zlim=zlim)
        	box()
        	for (i in (chrsep*2)-1) {
        		abline(v=i/((max(chrsep)-1)*2), col="grey")
        	}
        	for (i in seq(-1, max(((2*(ncol(mm2)-1))+1),1), 2)) {
        		abline(h=i/(2*(ncol(mm2)-1)), col="white", lwd=2)
        	}
			axis(1,at=chrmid/(max(chrsep)-1), label=rle(mm$chrom)$values, cex.axis=0.8, tick=F, family=fontfamily)
			axis(2,at=seq(0,1,1/max((ncol(mm2)-1),1)), label=sub(sample_column_postfix, "", colnames(mm2)), las=2, cex.axis=1, tick=F, family=fontfamily)
		}, mm)
    	par(xpd=FALSE)
		legend("bottom", inset=c(0,-.1), legend=c("Homozygous deletion", "Loss", "Gain", "Amplification"), fill=c("blue", "lightblue", "darksalmon", "red"), xpd=TRUE, ncol=2)
		par(xpd=TRUE)
		dev.off()
	}

	geneCN_tab = read.table(geneCN, sep="\t", header=T, stringsAsFactors=F)
	if (!opt$includeChrY) {
    	geneCN_tab = geneCN_tab[geneCN_tab$chrom != "Y",]
	}
	if (!is.null(opt$summaryConfig)) {
    	library(yaml)
    	summaryConfig = yaml.load_file("summary_config.yaml")
		if ("sample_rename" %in% names(summaryConfig)) {
        	for (c in names(summaryConfig$sample_rename)) {
            	colnames(geneCN_tab) = sub(c, summaryConfig$sample_rename[c], colnames(geneCN_tab))
        	}
    	}
	} else {
    	summaryConfig = NULL
	}
	plot_heatmap_(geneCN_tab, outFile, opt$sampleColumnPostFix, fontfamily, summaryConfig=summaryConfig)

} else if (as.numeric(opt$option)==5) {

	suppressPackageStartupMessages(library("dplyr"))
	
	facetsFiles = arguments$args
	Df = data.frame()
	for (facetsFile in facetsFiles) {
    	load(facetsFile)
    	tumorName = facetsFile %>% sub('.*/', '', .) %>% sub('_.*', '', .) %>% sub('\\..*', '', .)
    	normalName = facetsFile %>% sub('.*/', '', .) %>% sub('^.*_', '', .) %>% sub('\\..*', '', .)
    	n = paste(tumorName, normalName, sep = '_')
    	Df[n, 'tumorName'] = tumorName
    	if (tumorName != normalName) {
    	    Df[n, 'normalName'] = normalName
    	}
    	Df[n, 'purity'] = fit$purity
    	Df[n, 'ploidy'] = fit$ploidy
    	Df[n, 'dipLogR'] = fit$dipLogR
	}
	Df = mutate(Df, bad = purity <= 0.3 | is.na(purity))
	write.table(Df, file = opt$outFile, sep = '\t', quote = F, row.names = F)

}
