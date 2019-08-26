include modules/Makefile.inc

LOGDIR ?= log/plot_metrics.$(NOW)
PHONY += metrics metrics/report

plot_metrics : metrics/report/umi_frequencies.pdf \
			   metrics/report/umi_family_types_probe-A.pdf \
			   metrics/report/umi_family_types_probe-B.pdf \
			   metrics/report/umi_family_types_probe-AB.pdf

metrics/report/umi_frequencies.pdf : metrics/summary/umi_frequencies.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G,"$(RSCRIPT) modules/test/qc/plotmetrics.R --type 1")
	
metrics/report/umi_family_types_probe-A.pdf : metrics/summary/umi_family_types.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G,"$(RSCRIPT) modules/test/qc/plotmetrics.R --type 2")

metrics/report/umi_family_types_probe-B.pdf : metrics/summary/umi_family_types.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G,"$(RSCRIPT) modules/test/qc/plotmetrics.R --type 3")

metrics/report/umi_family_types_probe-AB.pdf : metrics/summary/umi_family_types.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G,"$(RSCRIPT) modules/test/qc/plotmetrics.R --type 4")

	
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
