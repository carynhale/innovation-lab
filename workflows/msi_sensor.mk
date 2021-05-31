include innovation-lab/Makefile.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/msi_sensor.$(NOW)

MSISENSOR_OPTS ?= -d $(REF_MSI) $(if $(TARGETS_FILE),-e $(TARGETS_FILE))

msi_sensor : $(foreach pair,$(SAMPLE_PAIRS),msisensor/$(pair).msi)
	     
define msisensor-tumor-normal
msisensor/$1_$2.msi : bam/$1.bam bam/$2.bam
	$$(call RUN,-c -n 8 -s 1G -m 2G -v $(MSISENSOR_ENV),"msisensor-pro msi $$(MSISENSOR_OPTS) -t $$(<) -n $$(<<) -b 8 -o $$(@)")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call msisensor-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))


..DUMMY := $(shell mkdir -p version; \
	     R --version > version/msi_sensor.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: msi_sensor
