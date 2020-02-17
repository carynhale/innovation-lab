include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/interval_metrics.$(NOW)
PHONY += metrics metrics/standard metrics/unfiltered metrics/duplex metrics/simplex metrics/summary

POOL_A_TARGET_FILE ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-A.sorted.list
POOL_B_TARGET_FILE ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-B.sorted.list
ONTARGET_FILE_A ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-A.sorted.bed
ONTARGET_FILE_B ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-B.sorted.bed
OFFTARGET_FILE ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-AB.offtarget.bed





# 			 metrics/standard/metrics_oxog.tsv \
# 			 metrics/standard/metrics_hs.tsv \



# 			 metrics/unfiltered/metrics_hs.tsv \



# 			 metrics/duplex/metrics_hs.tsv \



# 			 metrics/simplex/metrics_hs.tsv \

# 			 metrics/summary/metrics_idx.tsv \
# 			 metrics/summary/metrics_aln.tsv \
# 			 metrics/summary/metrics_insert.tsv \
# 			 metrics/summary/metrics_insert_distribution.tsv \
# 			 metrics/summary/metrics_hs.tsv \
# 			 metrics/summary/metrics_ts.tsv



		

	



									   	
metrics/standard/metrics_oxog.tsv : $(wildcard metrics/standard/$(SAMPLES).oxog_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 5 --sample_names '$(SAMPLES)'")

metrics/standard/metrics_hs.tsv : $(wildcard metrics/standard/$(SAMPLES).probe-A.hs_metrics.txt) $(wildcard metrics/standard/$(SAMPLES).probe-B.hs_metrics.txt) $(wildcard metrics/standard/$(SAMPLES).probe-A.hs_metrics-nodedup.txt) $(wildcard metrics/standard/$(SAMPLES).probe-B.hs_metrics-nodedup.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 6 --sample_names '$(SAMPLES)'")
	
	
	

	
metrics/unfiltered/metrics_hs.tsv : $(wildcard metrics/unfiltered/$(SAMPLES).probe-A.hs_metrics.txt) $(wildcard metrics/unfiltered/$(SAMPLES).probe-B.hs_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 11 --sample_names '$(SAMPLES)'")
	
		

	
	

	
metrics/duplex/metrics_hs.tsv : $(wildcard metrics/duplex/$(SAMPLES).probe-A.hs_metrics.txt) $(wildcard metrics/duplex/$(SAMPLES).probe-B.hs_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 16 --sample_names '$(SAMPLES)'")
	
		

	
	

	
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
