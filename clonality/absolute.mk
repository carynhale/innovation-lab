include innovation-lab/Makefile.inc

LOGDIR ?= log/absolute.$(NOW)

absolute : $(foreach set,$(SAMPLE_SETS),absolute/$(set)/$(set).vcf)

define run-absolute
absolute/%/%.vcf : summary/mutation_summary.txt
	$$(call RUN,-c -s 8G -m 12G -v $(ABSOLUTE_ENV),"set -o pipefail && \
												   	mkdir -p absolute && \
												   	mkdir -p absolute/$$(*) && \
												   	$(RSCRIPT) $(SCRIPTS_DIR)/clonality/absolute.R \
												   	--option 1 \
												   	--file_name $$(@) \
												   	--sample_set $$(*)")
												  
endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call run-absolute,$(set))))

..DUMMY := $(shell mkdir -p version; \
			 ~/share/usr/env/cntu-0.0.1/bin/R --version > version/absolute.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: absolute
