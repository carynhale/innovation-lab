include modules/Makefile.inc

LOGDIR ?= log/umi_qc.$(NOW)
PHONY += metrics metrics/summary

umi_qc : metrics/summary/umi-frequencies.tsv
		 
metrics/summary/umi-frequencies.tsv : $(wildcard marianas/$(SAMPLES)/composite-umi-frequencies.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"$(RSCRIPT) modules/test/qc/umiqc.R --type 1 --sample_names $(SAMPLES)")

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
