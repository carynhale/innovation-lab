include modules/Makefile.inc

LOGDIR ?= log/fusion_summary.$(NOW)
PHONY += fusion_summary

fusion_summary : summary/tsv/fusion_summary.tsv \
				 summary/fusion_summary.xlsx

summary/tsv/fusion_summary.tsv : $(wildcard $(foreach sample,$(TUMOR_SAMPLES),defuse/$(sample).taskcomplete)) $(wildcard $(foreach sample,$(TUMOR_SAMPLES),fusion_catcher/$(sample)/out/taskcomplete)) $(wildcard $(foreach sample,$(TUMOR_SAMPLES),integrate_rnaseq/breakpoints/$(sample).breakpoints.tsv))
	$(call RUN,-n 1 -s 4G -m 4G,"set -o pipefail && \
								 do something here ...")
				     			 
summary/fusion_summary.xlsx : summary/tsv/fusion_summary.tsv
	$(call RUN,-n 1 -s 4G -m 4G,"set -o pipefail && \
								 do something here ...")

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
