include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnv_access.$(NOW)
PHONY += cnvaccess

CNV_ACCESS_WORKFLOW += cnvaccess_coverage
CNV_ACCESS_WORKFLOW += cnvaccess_fix
CNV_ACCESS_WORKFLOW += cnvaccess_plot

cnv_access_workflow : $(CNV_ACCESS_WORKFLOW)

include modules/test/copy_number/cnvaccesscoverage.mk
include modules/test/copy_number/cnvaccessfix.mk
include modules/test/copy_number/cnvaccessplot.mk

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
