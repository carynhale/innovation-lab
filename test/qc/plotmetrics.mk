include modules/Makefile.inc

LOGDIR ?= log/plot_metrics.$(NOW)
PHONY += metrics metrics/report

plot_metrics : metrics/report/umi_frequencies.pdf

metrics/report/umi_frequencies.pdf : metrics/summary/umi_frequencies.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G,"$(RSCRIPT) modules/test/qc/plotmetrics.R --type 1")
	
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
