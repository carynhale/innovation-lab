include innovation-lab/Makefile.inc

LOGDIR ?= log/onco_fuse.$(NOW)

FUSION_CALLERS = arriba starfusion defuse fusioncatcher

oncofuse :  $(foreach caller,$(FUSION_CALLERS),oncofuse/$(caller)/summary.txt) \
	    $(foreach caller,$(FUSION_CALLERS),oncofuse/$(caller)/oncofuse.txt)
	   
define oncofuse-annotate
oncofuse/$1/summary.txt : $1/summary.txt
	$$(call RUN,-c -n 1 -s 4G -m 8G,"set -o pipefail && \
					 $(RSCRIPT) $(SCRIPTS_DIR)/rna_seq/summarize_oncofuse.R --fusion_caller $1")
					 
oncofuse/$1/oncofuse.txt : oncofuse/$1/summary.txt
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(ONCOFUSE_ENV),"set -o pipefail && \
							    $$(ONCOFUSE) $$(<) coord - $$(@)")
					 
endef
$(foreach caller,$(FUSION_CALLERS),\
		$(eval $(call oncofuse-annotate,$(caller))))

..DUMMY := $(shell mkdir -p version; \
	     R --version &> version/onco_fuse.txt;)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: oncofuse
