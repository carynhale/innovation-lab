include modules/Makefile.inc

LOGDIR ?= log/rnaseq_fusion.$(NOW)
PHONY += fastq defuse fusion_catcher

RNASEQ_FUSION_WORKFLOW += fastq
RNASEQ_FUSION_WORKFLOW += defuse
RNASEQ_FUSION_WORKFLOW += fusion_catcher

rnaseq_fusion_workflow : $(RNASEQ_FUSION_WORKFLOW)

include modules/fastq_tools/merge_split_fastq.mk
include modules/structural_variants/defuse.mk
include modules/structural_variants/fusioncatcher.mk

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
