include innovation-lab/Makefile.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/annotate_vcf_context.$(NOW)

annotate_vcf_context : vcf/target_vcf.txt

vcf/target_vcf.txt : vcf/target_vcf.vcf
	$(call RUN,-s 24G -m 48G -v $(ANNOTATION_ENV),"set -o pipefail && \
							                       $(RSCRIPT) $(SCRIPTS_DIR)/vcf_tools/annotate_vcf_context.R \
                                                   --file_in $(<) \
                                                   --file_out $(@) \
                                                   --ensembl_gene $(ENSEMBL)")
							  
	
..DUMMY := $(shell mkdir -p version; \
			 R --version >> version/annotate_vcf_context.txt)
.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: annotate_vcf_context
