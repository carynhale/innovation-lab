include innovation-lab/Makefile.inc

LOGDIR ?= log/onco_fuse.$(NOW)

FUSION_CALLERS = arriba starfusion defuse

oncofuse :  $(foreach caller,$(FUSION_CALLERS),oncofuse/$(caller).txt)
	   
define oncofuse-annotate
oncofuse/$1.txt : $1/summary.txt
	$$(call RUN,-c -n 1 -s 4G -m 8G,"set -o pipefail && \
					 $(RSCRIPT) $(SCRIPTS_DIR)/rna_seq/summarize_oncofuse.R --fusion_caller $1")
					 
endef
$(foreach caller,$(FUSION_CALLERS),\
		$(eval $(call oncofuse-annotate,$(caller))))

..DUMMY := $(shell mkdir -p version; \
	     R --version &> version/onco_fuse.txt;)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: oncofuse
