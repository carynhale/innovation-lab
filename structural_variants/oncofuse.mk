include innovation-lab/Makefile.inc

LOGDIR ?= log/onco_fuse.$(NOW)

FUSION_CALLERS = arriba starfusion defuse

oncofuse :  $(foreach caller,$(FUSION_CALLERS),oncofuse/$(caller).txt) \
	    $(foreach caller,$(FUSION_CALLERS),oncofuse/$(caller)_annotated_oncofuse.txt)
	   
define oncofuse-annotate
oncofuse/$1.txt : $1/summary.txt
	$$(call RUN,-c -n 1 -s 4G -m 8G,"set -o pipefail && \
					 $(RSCRIPT) $(SCRIPTS_DIR)/rna_seq/summarize_oncofuse.R --fusion_caller $1")
					 
oncofuse/$1_annotated_oncofuse.txt : oncofuse/$1.txt
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(ONCOFUSE_ENV),"set -o pipefail && \
							    $$(ONCOFUSE) oncofuse/$1.txt coord - $$(@)")
					 
endef
$(foreach caller,$(FUSION_CALLERS),\
		$(eval $(call oncofuse-annotate,$(caller))))

..DUMMY := $(shell mkdir -p version; \
	     R --version &> version/onco_fuse.txt;)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: oncofuse
