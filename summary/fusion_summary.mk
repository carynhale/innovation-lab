include modules/Makefile.inc

LOGDIR ?= log/fusion_summary.$(NOW)
PHONY += fusion_summary

fusion_summary : summary/fusion_summary.tsv

summary/fusion_summary.tsv : integrate_rnaseq/summary.tsv fusion_catcher/summary.tsv defuse/summary.tsv
	$(call RUN,-n 1 -s 4G -m 4G,"set -o pipefail && \
								 mkdir -p summary && \
								 do something here ...")
				     			 
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
