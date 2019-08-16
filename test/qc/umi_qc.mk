include modules/Makefile.inc

LOGDIR ?= log/umi_qc.$(NOW)
PHONY += marianas

umi_qc : $(foreach sample,$(SAMPLES),marianas/$(sample)/umi-info.RData)

define umi-frequencies
marianas/$1/umi-info.RData : marianas/$1/umi-frequencies.txt
	$$(call RUN,-c -n 1 -s 4G -m 6G ,"set -o pipefail && \
									 $(RSCRIPT) modules/test/qc/umi_qc.R --type 0 --sample_name $1")
																		           
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call umi-frequencies,$(sample))))


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
