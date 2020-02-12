include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/align_fastq.$(NOW)
PHONY += marianas




.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
