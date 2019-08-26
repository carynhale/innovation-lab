include modules/Makefile.inc

LOGDIR ?= log/plot_metrics.$(NOW)
PHONY += metrics metrics/report

plot_metrics : metrics/report/umi_frequencies.pdf \
			   metrics/report/umi_family_types_probe-A.pdf \
			   metrics/report/umi_family_types_probe-B.pdf \
			   metrics/report/umi_family_types_probe-AB.pdf \
			   metrics/report/umi_family_sizes_all.pdf \
			   metrics/report/umi_family_sizes_duplex.pdf \
			   metrics/report/umi_family_sizes_simplex.pdf

metrics/report/umi_frequencies.pdf : metrics/summary/umi_frequencies.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G,"$(RSCRIPT) modules/test/qc/plotmetrics.R --type 1")
	
metrics/report/umi_family_types_probe-A.pdf : metrics/summary/umi_family_types.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G,"$(RSCRIPT) modules/test/qc/plotmetrics.R --type 2")

metrics/report/umi_family_types_probe-B.pdf : metrics/summary/umi_family_types.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G,"$(RSCRIPT) modules/test/qc/plotmetrics.R --type 3")

metrics/report/umi_family_types_probe-AB.pdf : metrics/summary/umi_family_types.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G,"$(RSCRIPT) modules/test/qc/plotmetrics.R --type 4")

metrics/report/umi_family_sizes_all.pdf : metrics/summary/umi_families.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G,"$(RSCRIPT) modules/test/qc/plotmetrics.R --type 5")

metrics/report/umi_family_sizes_duplex.pdf : metrics/summary/umi_families.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G,"$(RSCRIPT) modules/test/qc/plotmetrics.R --type 6")

metrics/report/umi_family_sizes_simplex.pdf : metrics/summary/umi_families.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G,"$(RSCRIPT) modules/test/qc/plotmetrics.R --type 7")
	
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
