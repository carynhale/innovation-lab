include innovation-lab/Makefile.inc

LOGDIR ?= log/rnaseq_metrics.$(NOW)

REF_FLAT ?= $(HOME)/share/lib/resource_files/refFlat_ensembl.v75.txt
RIBOSOMAL_INTERVALS ?= $(HOME)/share/lib/resource_files/Homo_sapiens.GRCh37.75.rRNA.interval_list
STRAND_SPECIFICITY ?= NONE
TARGETS_LIST ?= $(HOME)/share/lib/resource_files/IDT_exome_research_targets.list

rnaseq_metrics : $(foreach sample,$(SAMPLES),metrics/$(sample)_rnaseq_metrics.txt) \
		 $(foreach sample,$(SAMPLES),metrics/$(sample)_alignment_metrics.txt) \
		 $(foreach sample,$(SAMPLES),metrics/$(sample)_insert_metrics.txt) \
		 $(foreach sample,$(SAMPLES),metrics/$(sample)_hs_metrics.txt) \
		 summary/rnaseq_metrics.txt \
		 summary/alignment_metrics.txt \
		 summary/insert_metrics.txt \
		 summary/hs_metrics.txt

define rnaseq-metrics
metrics/$1_rnaseq_metrics.txt : bam/$1.bam
	$$(call RUN,-c -s 6G -m 12G,"set -o pipefail && \
				     $$(COLLECT_RNASEQ_METRICS) \
				     INPUT=$$(<) \
				     OUTPUT=$$(@) \
				     REF_FLAT=$$(REF_FLAT) \
				     RIBOSOMAL_INTERVALS=$$(RIBOSOMAL_INTERVALS) \
				     CHART_OUTPUT=metrics/$1_rnaseq_metrics.pdf \
				     STRAND_SPECIFICITY=$$(STRAND_SPECIFICITY)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call rnaseq-metrics,$(sample))))
		
define alignment-metrics
metrics/$1_alignment_metrics.txt : bam/$1.bam
	$$(call RUN, -c -n 1 -s 6G -m 12G,"set -o pipefail && \
					   $$(COLLECT_ALIGNMENT_METRICS) \
					   REFERENCE_SEQUENCE=$$(REF_FASTA) \
					   INPUT=$$(<) \
					   OUTPUT=$$(@)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call alignment-metrics,$(sample))))
		
define insert-metrics
metrics/$1_insert_metrics.txt : bam/$1.bam
	$$(call RUN,-c -s 6G -m 12G,"set -o pipefail && \
				     $$(COLLECT_INSERT_METRICS) \
				     INPUT=$$(<) \
				     OUTPUT=$$(@) \
				     HISTOGRAM_FILE=metrics/$1_insert_metrics.pdf")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call insert-metrics,$(sample))))
		
define hs-metrics
metrics/$1_hs_metrics.txt : bam/$1.bam
	$$(call RUN,-c -n 1 -s 6G -m 12G,"set -o pipefail && \
					  $$(CALC_HS_METRICS) \
					  R=$(REF_FASTA) \
					  I=$$(<) \
					  O=$$(@) \
					  BAIT_INTERVALS=$(TARGETS_LIST) \
					  TARGET_INTERVALS=$(TARGETS_LIST) \
					  TMP_DIR=$(TMPDIR)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call hs-metrics,$(sample))))


summary/rnaseq_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample)_rnaseq_metrics.txt)
	$(call RUN, -c -n 1 -s 4G -m 6G,"set -o pipefail && \
					 $(RSCRIPT) $(SCRIPTS_DIR)/qc/rnaseq_metrics.R --option 1 --sample_names '$(SAMPLES)'")

summary/alignment_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample)_alignment_metrics.txt)
	$(call RUN, -c -n 1 -s 4G -m 6G,"set -o pipefail && \
					 $(RSCRIPT) $(SCRIPTS_DIR)/qc/rnaseq_metrics.R --option 2 --sample_names '$(SAMPLES)'")
									 
summary/insert_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample)_insert_metrics.txt)
	$(call RUN, -c -n 1 -s 4G -m 6G,"set -o pipefail && \
					 $(RSCRIPT) $(SCRIPTS_DIR)/qc/rnaseq_metrics.R --option 3 --sample_names '$(SAMPLES)'")

summary/hs_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample)_hs_metrics.txt)
	$(call RUN, -c -n 1 -s 4G -m 6G,"set -o pipefail && \
					 $(RSCRIPT) $(SCRIPTS_DIR)/qc/rnaseq_metrics.R --option 4 --sample_names '$(SAMPLES)'")

..DUMMY := $(shell mkdir -p version; \
	     echo "picard" >> version/rnaseq_metrics.txt; \
	     $(PICARD) CollectRnaSeqMetrics --version &>> version/rnaseq_metrics.txt; \
             R --version >> version/rnaseq_metrics.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: rnaseq_metrics
