include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/interval_metrics.$(NOW)
PHONY += metrics metrics/standard metrics/unfiltered metrics/simplex metrics/duplex

POOL_A_TARGET_FILE ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-A.sorted.list
POOL_B_TARGET_FILE ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-B.sorted.list

interval_metrics : $(foreach sample,$(SAMPLES),metrics/standard/$(sample).idx_stats.txt) \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).aln_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).insert_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).oxog_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).probe-A.hs_metrics.txt) \
				   $(foreach sample,$(SAMPLES),metrics/standard/$(sample).probe-B.hs_metrics.txt) \
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
				   metrics/unfiltered/metrics_hs.tsv


define picard-metrics-standard
metrics/standard/$1.idx_stats.txt : bam/$1-standard.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) BamIndexStats \
									   I=$$(<) \
									   TMP_DIR=$(TMPDIR) \
									   > $$(@)")
									   
metrics/standard/$1.aln_metrics.txt : bam/$1-standard.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) CollectAlignmentSummaryMetrics \
									   R=$(REF_FASTA) \
									   I=$$(<) \
									   O=$$(@) \
									   TMP_DIR=$(TMPDIR)")

metrics/standard/$1.insert_metrics.txt : bam/$1-standard.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) CollectInsertSizeMetrics \
									   I=$$(<) \
									   O=$$(@) \
									   H=metrics/standard/$1.insert_metrics.pdf \
									   M=0.5 \
									   TMP_DIR=$(TMPDIR)")
												
metrics/standard/$1.oxog_metrics.txt : bam/$1-standard.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) CollectOxoGMetrics \
									   R=$(REF_FASTA) \
									   I=$$(<) \
									   O=$$(@) \
									   TMP_DIR=$(TMPDIR)")

metrics/standard/$1.probe-A.hs_metrics.txt : bam/$1-standard.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) CalculateHsMetrics \
									   R=$(REF_FASTA) \
									   I=$$(<) \
									   O=$$(@) \
									   BAIT_INTERVALS=$(POOL_A_TARGET_FILE) \
									   TARGET_INTERVALS=$(POOL_A_TARGET_FILE) \
									   TMP_DIR=$(TMPDIR)")
												
metrics/standard/$1.probe-B.hs_metrics.txt : bam/$1-standard.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx12G -jar $$(PICARD_JAR) CalculateHsMetrics \
									   R=$(REF_FASTA) \
									   I=$$(<) \
									   O=$$(@) \
									   BAIT_INTERVALS=$(POOL_B_TARGET_FILE) \
									   TARGET_INTERVALS=$(POOL_B_TARGET_FILE) \
									   TMP_DIR=$(TMPDIR)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call picard-metrics-standard,$(sample))))
		
define picard-metrics-unfiltered
metrics/unfiltered/$1.idx_stats.txt : bam/$1-unfiltered.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) BamIndexStats \
									   I=$$(<) \
									   TMP_DIR=$(TMPDIR) \
									   > $$(@)")
									   
metrics/unfiltered/$1.aln_metrics.txt : bam/$1-unfiltered.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) CollectAlignmentSummaryMetrics \
									   R=$(REF_FASTA) \
									   I=$$(<) \
									   O=$$(@) \
									   TMP_DIR=$(TMPDIR)")

metrics/unfiltered/$1.insert_metrics.txt : bam/$1-unfiltered.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) CollectInsertSizeMetrics \
									   I=$$(<) \
									   O=$$(@) \
									   H=metrics/unfiltered/$1.insert_metrics.pdf \
									   M=0.5 \
									   TMP_DIR=$(TMPDIR)")
												
metrics/unfiltered/$1.probe-A.hs_metrics.txt : bam/$1-unfiltered.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx8G -jar $$(PICARD_JAR) CalculateHsMetrics \
									   R=$(REF_FASTA) \
									   I=$$(<) \
									   O=$$(@) \
									   BAIT_INTERVALS=$(POOL_A_TARGET_FILE) \
									   TARGET_INTERVALS=$(POOL_A_TARGET_FILE) \
									   TMP_DIR=$(TMPDIR)")
												
metrics/unfiltered/$1.probe-B.hs_metrics.txt : bam/$1-unfiltered.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx12G -jar $$(PICARD_JAR) CalculateHsMetrics \
									   R=$(REF_FASTA) \
									   I=$$(<) \
									   O=$$(@) \
									   BAIT_INTERVALS=$(POOL_B_TARGET_FILE) \
									   TARGET_INTERVALS=$(POOL_B_TARGET_FILE) \
									   TMP_DIR=$(TMPDIR)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call picard-metrics-unfiltered,$(sample))))		


metrics/standard/metrics_idx.tsv : $(wildcard metrics/standard/$(SAMPLES).idx_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"$(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 1 --sample_names '$(SAMPLES)'")
		
metrics/standard/metrics_aln.tsv : $(wildcard metrics/standard/$(SAMPLES).aln_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"$(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 2 --sample_names '$(SAMPLES)'")
	
metrics/standard/metrics_insert.tsv : $(wildcard metrics/standard/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"$(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 3 --sample_names '$(SAMPLES)'")
	
metrics/standard/metrics_insert_distribution.tsv : $(wildcard metrics/standard/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 16G -m 24G,"$(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 4 --sample_names '$(SAMPLES)'")
	
metrics/standard/metrics_oxog.tsv : $(wildcard metrics/standard/$(SAMPLES).oxog_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"$(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 5 --sample_names '$(SAMPLES)'")

metrics/standard/metrics_hs.tsv : $(wildcard metrics/standard/$(SAMPLES).probe-A.hs_metrics.txt) $(wildcard metrics/standard/$(SAMPLES).probe-B.hs_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"$(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 6 --sample_names '$(SAMPLES)'")
	
metrics/unfiltered/metrics_idx.tsv : $(wildcard metrics/unfiltered/$(SAMPLES).idx_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"$(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 7 --sample_names '$(SAMPLES)'")
		
metrics/unfiltered/metrics_aln.tsv : $(wildcard metrics/unfiltered/$(SAMPLES).aln_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"$(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 8 --sample_names '$(SAMPLES)'")
	
metrics/unfiltered/metrics_insert.tsv : $(wildcard metrics/unfiltered/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"$(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 9 --sample_names '$(SAMPLES)'")
	
metrics/unfiltered/metrics_insert_distribution.tsv : $(wildcard metrics/unfiltered/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 16G -m 24G,"$(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 10 --sample_names '$(SAMPLES)'")
	
metrics/unfiltered/metrics_hs.tsv : $(wildcard metrics/unfiltered/$(SAMPLES).probe-A.hs_metrics.txt) $(wildcard metrics/unfiltered/$(SAMPLES).probe-B.hs_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G,"$(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 11 --sample_names '$(SAMPLES)'")
	

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
