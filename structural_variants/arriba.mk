include innovation-lab/Makefile.inc

LOGDIR ?= log/arriba.$(NOW)

arriba : $(foreach sample,$(SAMPLES),arriba/$(sample)/$(sample).fusion.tsv) \
		 arriba/summary.txt

define run-arriba
arriba/$1/$1.fusion.tsv : bam/$1.bam
	$$(call RUN,-c -n 8 -s 1G -m 2G -v $(STARFUSION_ENV) -w 36:00:00,"set -o pipefail && \
																	  ")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call run-arriba,$(sample))))
		
arriba/summary.txt : $(foreach sample,$(SAMPLES),arriba/$(sample)/$(sample).fusion.tsv)
	echo "" > arriba/summary.txt; \
	for i in $(SAMPLES); do \
		sed -e "1d" arriba/$$i/$$i.fusion.tsv | sed "s/$$/\t$$i/" >> arriba/summary.txt; \
	done

..DUMMY := $(shell mkdir -p version; \
			 share/usr/arriba_v2.0.0/arriba -h &> version/arriba.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: arriba
