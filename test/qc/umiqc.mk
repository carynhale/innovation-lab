include modules/Makefile.inc

LOGDIR ?= log/umi_qc.$(NOW)
PHONY += marinanas metrics metrics/summary

umi_qc : $(foreach sample,$(SAMPLES),marianas/$(sample)/family-sizes.txt) \
		 metrics/summary/umi_frequencies.tsv \
		 metrics/summary/umi_composite.tsv \
		 metrics/summary/umi_families.tsv

define family-size-metric
marianas/$1/family-sizes.txt : marianas/$1/second-pass-alt-alleles.txt
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
									   cd marianas/$1 && \
									   source ../../modules/test/qc/umiqc.sh $(UMI_QC_BED_FILE_A) $(UMI_QC_BED_FILE_B) $1")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call family-size-metric,$(sample))))
		 
metrics/summary/umi_frequencies.tsv : $(wildcard marianas/$(SAMPLES)/umi-frequencies.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"$(RSCRIPT) modules/test/qc/umiqc.R --type 1 --sample_names '$(SAMPLES)'")
	
metrics/summary/umi_composite.tsv : $(wildcard marianas/$(SAMPLES)/composite-umi-frequencies.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"$(RSCRIPT) modules/test/qc/umiqc.R --type 2 --sample_names '$(SAMPLES)'")

metrics/summary/umi_families.tsv : $(wildcard marianas/$(SAMPLES)/family-sizes.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"$(RSCRIPT) modules/test/qc/umiqc.R --type 3 --sample_names '$(SAMPLES)'")

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
