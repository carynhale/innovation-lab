include modules/Makefile.inc

LOGDIR ?= log/plot_metrics.$(NOW)
PHONY += metrics metrics/report

plot_metrics : metrics/report/umi_frequencies.pdf \
			   metrics/report/umi_family_types_probe-A.pdf \
			   metrics/report/umi_family_types_probe-B.pdf \
			   metrics/report/umi_family_types_probe-AB.pdf \
			   metrics/report/umi_family_sizes_all.pdf \
			   metrics/report/umi_family_sizes_duplex.pdf \
			   metrics/report/umi_family_sizes_simplex.pdf \
			   metrics/report/mean_standard_target_coverage-dedup.pdf \
			   metrics/report/mean_standard_target_coverage-nodedup.pdf \
			   metrics/report/mean_unfiltered_target_coverage.pdf \
			   metrics/report/mean_duplex_target_coverage.pdf \
			   metrics/report/mean_simplex_target_coverage.pdf \
			   metrics/report/aligment_summary.pdf \
			   metrics/report/insert_size_summary.pdf \
			   metrics/report/insert_size_distribution.pdf \
			   metrics/report/read_alignment_summary.pdf \
			   metrics/report/non_reference_calls.pdf

metrics/report/umi_frequencies.pdf : metrics/summary/umi_frequencies.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G -v $(SUPERHEAT_ENV),"$(RSCRIPT) modules/test/qc/plotmetrics.R --type 1 && \
														  gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dFirstPage=2 -dLastPage=2 -sOutputFile=metrics/report/umi_frequencies-2.pdf metrics/report/umi_frequencies.pdf && \
														  rm metrics/report/umi_frequencies.pdf && \
														  mv metrics/report/umi_frequencies-2.pdf metrics/report/umi_frequencies.pdf")
	
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
	
metrics/report/mean_standard_target_coverage-dedup.pdf : metrics/summary/metrics_hs.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G,"$(RSCRIPT) modules/test/qc/plotmetrics.R --type 8")

metrics/report/mean_standard_target_coverage-nodedup.pdf : metrics/summary/metrics_hs.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G,"$(RSCRIPT) modules/test/qc/plotmetrics.R --type 9")

metrics/report/mean_unfiltered_target_coverage.pdf : metrics/summary/metrics_hs.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G,"$(RSCRIPT) modules/test/qc/plotmetrics.R --type 10")

metrics/report/mean_duplex_target_coverage.pdf : metrics/summary/metrics_hs.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G,"$(RSCRIPT) modules/test/qc/plotmetrics.R --type 11")

metrics/report/mean_simplex_target_coverage.pdf : metrics/summary/metrics_hs.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G,"$(RSCRIPT) modules/test/qc/plotmetrics.R --type 12")
	
metrics/report/aligment_summary.pdf : metrics/summary/metrics_idx.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G,"$(RSCRIPT) modules/test/qc/plotmetrics.R --type 13")
	
metrics/report/insert_size_summary.pdf : metrics/summary/metrics_insert.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G,"$(RSCRIPT) modules/test/qc/plotmetrics.R --type 14")
	
metrics/report/insert_size_distribution.pdf : metrics/summary/metrics_insert_distribution.tsv
	$(call RUN, -c -n 1 -s 12G -m 16G,"$(RSCRIPT) modules/test/qc/plotmetrics.R --type 15")
	
metrics/report/read_alignment_summary.pdf : metrics/summary/metrics_ts.tsv
	$(call RUN, -c -n 1 -s 12G -m 16G,"$(RSCRIPT) modules/test/qc/plotmetrics.R --type 20")
	
metrics/report/non_reference_calls.pdf : $(wildcard metrics/standard/$(SAMPLES)-pileup.txt) $(wildcard metrics/simplex/$(SAMPLES)-pileup.txt) $(wildcard metrics/duplex/$(SAMPLES)-pileup.txt)
	$(call RUN, -c -n 1 -s 48G -m 72G -v $(SUPERHEAT_ENV),"$(RSCRIPT) modules/test/qc/plotmetrics.R --type 21 --sample_names '$(SAMPLES)'")


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
