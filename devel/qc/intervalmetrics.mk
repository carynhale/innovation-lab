include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/interval_metrics.$(NOW)
PHONY += metrics metrics/standard metrics/unfiltered metrics/duplex metrics/simplex metrics/summary

POOL_A_TARGET_FILE ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-A.sorted.list
POOL_B_TARGET_FILE ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-B.sorted.list
ONTARGET_FILE_A ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-A.sorted.bed
ONTARGET_FILE_B ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-B.sorted.bed
OFFTARGET_FILE ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-AB.offtarget.bed

interval_metrics : $(foreach sample,$(SAMPLES),metrics/standard/$(sample).idx_stats.txt) \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).aln_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).insert_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).oxog_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).probe-A.hs_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).probe-B.hs_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).probe-A.hs_metrics-nodedup.txt) \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).probe-B.hs_metrics-nodedup.txt) \
				   metrics/standard/metrics_idx.tsv \
				   metrics/standard/metrics_aln.tsv \
				   metrics/standard/metrics_insert.tsv \
				   metrics/standard/metrics_insert_distribution.tsv \
				   metrics/standard/metrics_oxog.tsv \
				   metrics/standard/metrics_hs.tsv \
				   $(foreach sample,$(SAMPLES),metrics/unfiltered/$(sample).idx_stats.txt) \
				   $(foreach sample,$(SAMPLES),metrics/unfiltered/$(sample).aln_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/unfiltered/$(sample).insert_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/unfiltered/$(sample).probe-A.hs_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/unfiltered/$(sample).probe-B.hs_metrics.txt) \
				   metrics/unfiltered/metrics_idx.tsv \
				   metrics/unfiltered/metrics_aln.tsv \
				   metrics/unfiltered/metrics_insert.tsv \
				   metrics/unfiltered/metrics_insert_distribution.tsv \
				   metrics/unfiltered/metrics_hs.tsv \
				   $(foreach sample,$(SAMPLES),metrics/duplex/$(sample).idx_stats.txt) \
				   $(foreach sample,$(SAMPLES),metrics/duplex/$(sample).aln_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/duplex/$(sample).insert_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/duplex/$(sample).probe-A.hs_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/duplex/$(sample).probe-B.hs_metrics.txt) \
				   metrics/duplex/metrics_idx.tsv \
				   metrics/duplex/metrics_aln.tsv \
				   metrics/duplex/metrics_insert.tsv \
				   metrics/duplex/metrics_insert_distribution.tsv \
				   metrics/duplex/metrics_hs.tsv \
				   $(foreach sample,$(SAMPLES),metrics/simplex/$(sample).idx_stats.txt) \
				   $(foreach sample,$(SAMPLES),metrics/simplex/$(sample).aln_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/simplex/$(sample).insert_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/simplex/$(sample).probe-A.hs_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/simplex/$(sample).probe-B.hs_metrics.txt) \
				   metrics/simplex/metrics_idx.tsv \
				   metrics/simplex/metrics_aln.tsv \
				   metrics/simplex/metrics_insert.tsv \
				   metrics/simplex/metrics_insert_distribution.tsv \
				   metrics/simplex/metrics_hs.tsv \
				   metrics/summary/metrics_idx.tsv \
				   metrics/summary/metrics_aln.tsv \
				   metrics/summary/metrics_insert.tsv \
				   metrics/summary/metrics_insert_distribution.tsv \
				   metrics/summary/metrics_hs.tsv \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).A.ontarget.txt) \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).B.ontarget.txt) \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).AB.offtarget.txt) \
				   metrics/summary/metrics_ts.tsv \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample)-pileup.txt) \
				   $(foreach sample,$(SAMPLES),metrics/simplex/$(sample)-pileup.txt) \
				   $(foreach sample,$(SAMPLES),metrics/duplex/$(sample)-pileup.txt)



		
metrics/standard/metrics_idx.tsv : $(wildcard metrics/standard/$(SAMPLES).idx_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 1 --sample_names '$(SAMPLES)'")
		
metrics/standard/metrics_aln.tsv : $(wildcard metrics/standard/$(SAMPLES).aln_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 2 --sample_names '$(SAMPLES)'")
	
metrics/standard/metrics_insert.tsv : $(wildcard metrics/standard/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 3 --sample_names '$(SAMPLES)'")
	
metrics/standard/metrics_insert_distribution.tsv : $(wildcard metrics/standard/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 16G -m 24G,"set -o pipefail && \
									   $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 4 --sample_names '$(SAMPLES)'")
	
metrics/standard/metrics_oxog.tsv : $(wildcard metrics/standard/$(SAMPLES).oxog_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 5 --sample_names '$(SAMPLES)'")

metrics/standard/metrics_hs.tsv : $(wildcard metrics/standard/$(SAMPLES).probe-A.hs_metrics.txt) $(wildcard metrics/standard/$(SAMPLES).probe-B.hs_metrics.txt) $(wildcard metrics/standard/$(SAMPLES).probe-A.hs_metrics-nodedup.txt) $(wildcard metrics/standard/$(SAMPLES).probe-B.hs_metrics-nodedup.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 6 --sample_names '$(SAMPLES)'")
	
metrics/unfiltered/metrics_idx.tsv : $(wildcard metrics/unfiltered/$(SAMPLES).idx_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 7 --sample_names '$(SAMPLES)'")
		
metrics/unfiltered/metrics_aln.tsv : $(wildcard metrics/unfiltered/$(SAMPLES).aln_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 8 --sample_names '$(SAMPLES)'")
	
metrics/unfiltered/metrics_insert.tsv : $(wildcard metrics/unfiltered/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 9 --sample_names '$(SAMPLES)'")
	
metrics/unfiltered/metrics_insert_distribution.tsv : $(wildcard metrics/unfiltered/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 16G -m 24G,"set -o pipefail && \
									   $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 10 --sample_names '$(SAMPLES)'")
	
metrics/unfiltered/metrics_hs.tsv : $(wildcard metrics/unfiltered/$(SAMPLES).probe-A.hs_metrics.txt) $(wildcard metrics/unfiltered/$(SAMPLES).probe-B.hs_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 11 --sample_names '$(SAMPLES)'")
	
metrics/duplex/metrics_idx.tsv : $(wildcard metrics/duplex/$(SAMPLES).idx_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 12 --sample_names '$(SAMPLES)'")
		
metrics/duplex/metrics_aln.tsv : $(wildcard metrics/duplex/$(SAMPLES).aln_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 13 --sample_names '$(SAMPLES)'")
	
metrics/duplex/metrics_insert.tsv : $(wildcard metrics/duplex/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 14 --sample_names '$(SAMPLES)'")
	
metrics/duplex/metrics_insert_distribution.tsv : $(wildcard metrics/duplex/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 16G -m 24G,"set -o pipefail && \
									   $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 15 --sample_names '$(SAMPLES)'")
	
metrics/duplex/metrics_hs.tsv : $(wildcard metrics/duplex/$(SAMPLES).probe-A.hs_metrics.txt) $(wildcard metrics/duplex/$(SAMPLES).probe-B.hs_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 16 --sample_names '$(SAMPLES)'")
	
metrics/simplex/metrics_idx.tsv : $(wildcard metrics/simplex/$(SAMPLES).idx_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 17 --sample_names '$(SAMPLES)'")
		
metrics/simplex/metrics_aln.tsv : $(wildcard metrics/simplex/$(SAMPLES).aln_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 18 --sample_names '$(SAMPLES)'")
	
metrics/simplex/metrics_insert.tsv : $(wildcard metrics/simplex/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 19 --sample_names '$(SAMPLES)'")
	
metrics/simplex/metrics_insert_distribution.tsv : $(wildcard metrics/simplex/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 16G -m 24G,"set -o pipefail && \
									   $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 20 --sample_names '$(SAMPLES)'")
	
metrics/simplex/metrics_hs.tsv : $(wildcard metrics/simplex/$(SAMPLES).probe-A.hs_metrics.txt) $(wildcard metrics/simplex/$(SAMPLES).probe-B.hs_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 21 --sample_names '$(SAMPLES)'")
	
metrics/summary/metrics_idx.tsv : metrics/standard/metrics_idx.tsv metrics/unfiltered/metrics_idx.tsv metrics/duplex/metrics_idx.tsv metrics/simplex/metrics_idx.tsv
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 22")
		
metrics/summary/metrics_aln.tsv : metrics/standard/metrics_aln.tsv metrics/unfiltered/metrics_aln.tsv metrics/duplex/metrics_aln.tsv metrics/simplex/metrics_aln.tsv
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 23")
	
metrics/summary/metrics_insert.tsv : metrics/standard/metrics_insert.tsv metrics/unfiltered/metrics_insert.tsv metrics/duplex/metrics_insert.tsv metrics/simplex/metrics_insert.tsv
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 24")
	
metrics/summary/metrics_insert_distribution.tsv : metrics/standard/metrics_insert_distribution.tsv metrics/unfiltered/metrics_insert_distribution.tsv metrics/duplex/metrics_insert_distribution.tsv metrics/simplex/metrics_insert_distribution.tsv
	$(call RUN, -c -n 1 -s 16G -m 24G,"set -o pipefail && \
									   $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 25")
	
metrics/summary/metrics_hs.tsv : metrics/standard/metrics_hs.tsv metrics/unfiltered/metrics_hs.tsv metrics/duplex/metrics_hs.tsv metrics/simplex/metrics_hs.tsv
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 26")

metrics/summary/metrics_ts.tsv : $(wildcard metrics/standard/$(SAMPLES).A.ontarget.txt) $(wildcard metrics/standard/$(SAMPLES).B.ontarget.txt) $(wildcard metrics/standard/$(SAMPLES).AB.offtarget.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 27 --sample_names '$(SAMPLES)'")
	
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
