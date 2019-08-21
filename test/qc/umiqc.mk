include modules/Makefile.inc

LOGDIR ?= log/umi_qc.$(NOW)
PHONY += marianas

umi_qc : $(foreach sample,$(SAMPLES),marianas/$(sample)/umi-info.RData) \
		 $(foreach sample,$(SAMPLES),marianas/$(sample)/umi-composite.pdf) \

define umi-frequencies
marianas/$1/umi-info.RData : marianas/$1/timestamp
	$$(call RUN,-c -n 1 -s 4G -m 6G ,"set -o pipefail && \
									 $(RSCRIPT) modules/test/qc/umiqc.R --type 0 --sample_name $1")

marianas/$1/umi-composite.pdf : marianas/$1/umi-info.RData	
	$$(call RUN,-c -n 1 -s 4G -m 6G ,"set -o pipefail && \
									 $(RSCRIPT) modules/test/qc/umiqc.R --type 1 --sample_name $1")
																           
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call umi-frequencies,$(sample))))


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
