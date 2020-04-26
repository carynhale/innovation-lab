include innovation-lab/Makefile.inc

LOGDIR ?= log/absolute.$(NOW)

abs : $(foreach set,$(SAMPLE_SETS),absolute/$(set).vcf)
#		   $(foreach sample,$(SAMPLES),absolute/$(sample).txt)
		   

define run-absolute
absolute/$1.vcf : summary/mutation_summary.txt
	$$(call RUN,-c -s 8G -m 12G -v $(ABSOLUTE_ENV),"set -o pipefail && \
												   	mkdir -p absolute && \
												   	$(RSCRIPT) $(SCRIPTS_DIR)/clonality/absolute.R \
												   	--option 1 \
												   	--file_name $$(<) \
												   	--sample_set $1")
												  
endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call run-absolute,$(set))))
		
#define run-sufam
#absolute/$1.txt : $(wildcard $(foreach set,$(SAMPLE_SETS),absolute/$(set).vcf)) bam/$1.bam
#	$$(call RUN,-c -s 6G -m 8G -v $(ABSOLUTE_ENV) -w 12:00:00,"set -o pipefail && \
#															   $(RSCRIPT) $(SCRIPTS_DIR)/clonality/absolute.R \
#															   --option 2 \
#															   --sample_name $1 \
#															   --sample_set $(SAMPLE_SETS) \
#															   --ref_fasta $(REF_FASTA)")
#												  
#endef
#$(foreach sample,$(SAMPLES),\
#		$(eval $(call run-sufam,$(sample))))

..DUMMY := $(shell mkdir -p version; \
			 ~/share/usr/env/cntu-0.0.1/bin/R --version > version/absolute.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: abs
