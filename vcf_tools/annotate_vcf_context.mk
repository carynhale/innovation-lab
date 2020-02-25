include innovation-lab/Makefile.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/annotate_vcf_context.$(NOW)

annotate_vcf_context : vcf/external.txt

ANNOTATION_ENV = $(HOME)/share/usr/env/r-bsgenome.hsapiens.ucsc.hg19-1.4.0

vcf/external.txt : vcf/external.vcf
	$(call RUN,-s 24G -m 48G -v $(ANNOTATION_ENV),"set -o pipefail && \
							                       $(RSCRIPT) $(SCRIPTS_DIR)/vcf_tools/annotate_vcf_context.R \
                                                   --file_in vcf/external.vcf \
                                                   --file_out vcf/external.txt \
                                                   --genome_build 'hg19' \
                                                   --ensembl_gene 'ensgene.RData'")
							  
	
..DUMMY := $(shell mkdir -p version; \
			 R --version >> version/annotate_vcf_context.txt)
.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: annotate_vcf_context
