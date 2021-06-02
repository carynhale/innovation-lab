include innovation-lab/Makefile.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/msi_sensor.$(NOW)

MSISENSOR_OPTS ?= -d $(REF_MSI) $(if $(TARGETS_FILE),-e $(TARGETS_FILE))

msi_sensor : $(foreach pair,$(SAMPLE_PAIRS),msisensor/$(pair).msi) \
	     summary/msi_metrics.txt
	     
define msisensor-tumor-normal
msisensor/$1_$2.msi : bam/$1.bam bam/$2.bam
	$$(call RUN,-c -n 8 -s 1G -m 2G -v $(MSISENSOR_ENV),"msisensor-pro msi $$(MSISENSOR_OPTS) -t $$(<) -n $$(<<) -b 8 -o $$(@)")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call msisensor-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
	
summary/msi_metrics.txt : $(foreach pair,$(SAMPLE_PAIRS),msisensor/$(pair).msi)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/summary/msi_sensor.R --option 1 --sample_pairs '$(SAMPLE_PAIRS)'")


..DUMMY := $(shell mkdir -p version; \
	     R --version > version/msi_sensor.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: msi_sensor
