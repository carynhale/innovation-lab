include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvaccess.$(NOW)
PHONY += cnvkit cnvkit/cnn cnvkit/cnn/tumor cnvkit/cnn/normal cnvkit/reference cnvkit/cnr cnvkit/log2 cnvkit/segmented cnvkit/called cnvkit/summary

CNV_KIT_WORKFLOW += cnvkit_coverage
CNV_KIT_WORKFLOW += cnvkit_reference

cnv_kit_workflow : $(CNV_KIT_WORKFLOW)

include modules/test/copy_number/cnvkitcoverage.mk
include modules/copy_number/cnvkitreference.mk
#include modules/copy_number/cnvkitfix.mk
#include modules/copy_number/cnvkitplot.mk
#include modules/copy_number/cnvkitsegment.mk
#include modules/copy_number/cnvkitsummary.mk

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
