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
				   metrics/standard/idx_metrics.tsv \
				   metrics/standard/aln_metrics.tsv \
				   metrics/standard/insert_metrics.tsv

define picard-metrics
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
		$(eval $(call picard-metrics,$(sample))))

metrics/standard/idx_metrics.tsv : $(wildcard metrics/standard/$(SAMPLES).idx_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"$(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 1 --sample_names '$(SAMPLES)'")
		
metrics/standard/aln_metrics.tsv : $(wildcard metrics/standard/$(SAMPLES).aln_stats.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"$(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 2 --sample_names '$(SAMPLES)'")
	
metrics/standard/insert_metrics.tsv : $(wildcard metrics/standard/$(SAMPLES).insert_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 16G,"$(RSCRIPT) modules/test/qc/intervalmetrics.R --metric_type 3 --sample_names '$(SAMPLES)'")


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
