include innovation-lab/Makefile.inc

LOGDIR ?= log/rnaseq_metrics.$(NOW)

REF_FLAT ?= $(HOME)/share/lib/resource_files/refFlat_hg19.txt

rnaseq_metrics : $(foreach sample,$(SAMPLES),metrics/$(sample).taskcomplete)

define rnaseq-metrics
metrics/$1.taskcomplete : bam/$1.bam
	$$(call RUN,-c -s 6G -m 12G,"set -o pipefail && \
                                 for i in $(REF_FLAT); do \
                                 echo $i; \
                                 done;")
									 
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call rnaseq-metrics,$(sample))))

..DUMMY := $(shell mkdir -p version; \
			 echo "picard" >> version/rnaseq_metrics.txt; \
			 $(PICARD) CollectRnaSeqMetrics --version &>> version/rnaseq_metrics.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY:


#$$(COLLECT_RNASEQ_METRICS) \
#                                 INPUT=$$(<) \
#                                 OUTPUT=$$(@) \
#                                 REF_FLAT=$$(REF_FLAT) \
#                                 CHART_OUTPUT=metrics/$1.pdf \
#                                 STRAND_SPECIFICITY=NONE"