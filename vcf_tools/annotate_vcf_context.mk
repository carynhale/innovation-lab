include innovation-lab/Makefile.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/annotate_vcf_context.$(NOW)

annotate_vcf_context : 

vcf/.vcf : vcf/.vcf
	$(call RUN,-s 16G -m 20G,"set -o pipefail && \
							  ")
							  
	
..DUMMY := $(shell mkdir -p version; \
			 R --version >> version/annotate_vcf_context.txt)
.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: annotate_vcf_context