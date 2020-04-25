include innovation-lab/Makefile.inc

LOGDIR ?= log/absolute.$(NOW)

absolute : $(foreach set,$(SAMPLE_SETS),absolute/$(set)/$(set).vcf) \
		   $(foreach set,$(SAMPLE_SETS),absolute/$(set)/$(set).taskcomplete)
		   
CMD = IFS='_' read -ra SAMPLE_SET <<< $1; for i in $SAMPLE_SET; do echo $i; done

define run-absolute
absolute/$1/$1.vcf : summary/mutation_summary.txt
	$$(call RUN,-c -s 8G -m 12G -v $(ABSOLUTE_ENV),"set -o pipefail && \
												   	mkdir -p absolute && \
												   	mkdir -p absolute/$1 && \
												   	$(RSCRIPT) $(SCRIPTS_DIR)/clonality/absolute.R \
												   	--option 1 \
												   	--file_name $$(<) \
												   	--sample_set $1")
												  
endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call run-absolute,$(set))))
		
define run-sufam
absolute/$1/$1.taskcomplete : absolute/$1/$1.vcf
	$$(call RUN,-c -s 6G -m 8G -v $(ABSOLUTE_ENV),"set -o pipefail && \
												   $(CMD) && \
												   touch absolute/$1/$1.taskcomplete")
												  
endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call run-sufam,$(set))))


..DUMMY := $(shell mkdir -p version; \
			 ~/share/usr/env/cntu-0.0.1/bin/R --version > version/absolute.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: absolute
