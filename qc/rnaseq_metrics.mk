include innovation-lab/Makefile.inc

LOGDIR ?= log/rnaseq_metrics.$(NOW)

REF_FLAT ?= $(HOME)/share/lib/resource_files/refFlat_hg19.txt
RIBOSOMAL_INTERVALS ?= /ifs/work/bergerm1/RNAseq/ref/Homo_sapiens.GRCh37.75.rRNA.interval_list

rnaseq_metrics : $(foreach sample,$(SAMPLES),metrics/$(sample).txt)

define rnaseq-metrics
metrics/$1.txt : bam/$1.bam
	$$(call RUN,-c -s 6G -m 12G,"set -o pipefail && \
                                 $$(COLLECT_RNASEQ_METRICS) \
                                 INPUT=$$(<) \
                                 OUTPUT=$$(@) \
                                 REF_FLAT=$$(REF_FLAT) \
                                 RIBOSOMAL_INTERVALS=$$(RIBOSOMAL_INTERVALS) \
                                 CHART_OUTPUT=metrics/$1.pdf \
                                 STRAND_SPECIFICITY=NONE")
									 
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call rnaseq-metrics,$(sample))))

..DUMMY := $(shell mkdir -p version; \
			 echo "picard" >> version/rnaseq_metrics.txt; \
			 $(PICARD) CollectRnaSeqMetrics --version &>> version/rnaseq_metrics.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY:
