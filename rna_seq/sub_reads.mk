include innovation-lab/Makefile.inc

LOGDIR = log/sub_reads.$(NOW)

sub_reads : $(foreach sample,$(SAMPLES),sumreads/$(sample).sumreads_bygene.txt) \
	    




..DUMMY := $(shell mkdir -p version; \
	     /home/brownd7/share/usr/env/r-rsubread-2.4.3/bin/R --version &> version/sub_reads.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: sub_reads
