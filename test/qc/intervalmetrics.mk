include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/interval_metrics.$(NOW)
PHONY += metrics metrics/standard metrics/unfiltered metrics/simplex metrics/duplex

POOL_A_TARGET_FILE ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-A.sorted.list
POOL_B_TARGET_FILE ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-B..sorted.list

interval_metrics : $(foreach sample,$(SAMPLES),metrics/standard/$(sample).idx_stats.txt)

define picard-metrics
metrics/standard/$1.idx_stats.txt : bam/%-standard.bam
	$$(call RUN, -c -n 1 -s 12G -m 18G -w 1440,"java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx16G -jar $$(PICARD_JAR) BamIndexStats \
												I=$$(<) \
												O=$$(@) \
												TMP_DIR=$(TMPDIR)")

#metrics/standard/$1.aln_metrics.txt : bam/%-standard.bam
#	$$(call RUN, -c -n 1 -s 12G -m 18G -w 1440,"java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx16G -jar $$(PICARD_JAR) CollectAlignmentSummaryMetrics \
#												R=$(REF_FASTA) \
#												I=$$(<) \
#												O=$$(@) \
#												TMP_DIR=$(TMPDIR)")
#
#metrics/standard/$1.insert_metrics.txt : bam/%-standard.bam
#	$$(call RUN, -c -n 1 -s 12G -m 18G -w 1440,"java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx16G -jar $$(PICARD_JAR) CollectInsertSizeMetrics \
#												I=$$(<) \
#												O=$$(@) \
#												H=metrics/standard/$1.insert_metrics.pdf \
#												M=0.5 \
#												TMP_DIR=$(TMPDIR)")
#												
#metrics/standard/$1.insert_metrics.txt : bam/%-standard.bam
#	$$(call RUN, -c -n 1 -s 12G -m 18G -w 1440,"java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx16G -jar $$(PICARD_JAR) CollectOxoGMetrics \
#												R=$(REF_FASTA) \
#												I=$$(<) \
#												O=$$(@) \
#												TMP_DIR=$(TMPDIR)")
#
#metrics/standard/$1.pool-A.hs_metrics.txt : bam/%-standard.bam
#	$$(call RUN, -c -n 1 -s 12G -m 18G -w 1440,"java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx16G -jar $$(PICARD_JAR) CollectHsMetrics \
#												R=$(REF_FASTA) \
#												I=$$(<) \
#												O=$$(@) \
#												BAIT_INTERVALS=$(POOL_A_TARGET_FILE) \
#												TARGET_INTERVALS=$(POOL_A_TARGET_FILE) \
#												TMP_DIR=$(TMPDIR)")
#												
#metrics/standard/$1.pool-B.hs_metrics.txt : bam/%-standard.bam
#	$$(call RUN, -c -n 1 -s 12G -m 18G -w 1440,"java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx16G -jar $$(PICARD_JAR) CollectHsMetrics \
#												R=$(REF_FASTA) \
#												I=$$(<) \
#												O=$$(@) \
#												BAIT_INTERVALS=$(POOL_B_TARGET_FILE) \
#												TARGET_INTERVALS=$(POOL_B_TARGET_FILE) \
#												TMP_DIR=$(TMPDIR)")
#
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call picard-metrics,$(sample))))


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
