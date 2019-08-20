include modules/Makefile.inc

LOGDIR ?= log/genome_summary.$(NOW)
PHONY += genome_stats summary summary/tsv

LST_SCORE ?= $(wildcard $(foreach set,$(SAMPLE_PAIRS),genome_stats/$(set).lst))
GENOME_ALTERED ?= $(wildcard $(foreach set,$(SAMPLE_PAIRS),genome_stats/$(set).fga))
NTAI_SCORE ?= $(wildcard $(foreach set,$(SAMPLE_PAIRS),genome_stats/$(set).ntai))
MYRIAD_SCORE ?= $(wildcard $(foreach set,$(SAMPLE_PAIRS),genome_stats/$(set).mrs))

genome_summary : genome_stats/lst_score.tsv genome_stats/genome_altered.tsv genome_stats/ntai_score.tsv genome_stats/myriad_score.tsv summary/tsv/genome_summary.tsv summary/genome_summary.xlsx

genome_stats/lst_score.tsv : $(wildcard $(foreach set,$(SAMPLE_PAIRS),genome_stats/$(set).lst))
	$(call RUN,-n 1 -s 4G -m 4G,"set -o pipefail && \
								 cat $(LST_SCORE) > genome_stats/lst_score.tsv")
				     			 
genome_stats/genome_altered.tsv : $(wildcard $(foreach set,$(SAMPLE_PAIRS),genome_stats/$(set).fga))
	$(call RUN,-n 1 -s 4G -m 4G,"set -o pipefail && \
								 cat $(GENOME_ALTERED) > genome_stats/genome_altered.tsv")
				     			 
genome_stats/ntai_score.tsv : $(wildcard $(foreach set,$(SAMPLE_PAIRS),genome_stats/$(set).ntai))
	$(call RUN,-n 1 -s 4G -m 4G,"set -o pipefail && \
				     			 cat $(NTAI_SCORE) > genome_stats/ntai_score.tsv")
				     			 
genome_stats/myriad_score.tsv : $(wildcard $(foreach set,$(SAMPLE_PAIRS),genome_stats/$(set).mrs))
	$(call RUN,-n 1 -s 4G -m 4G,"set -o pipefail && \
				     			 cat $(MYRIAD_SCORE) > genome_stats/myriad_score.tsv")

summary/tsv/genome_summary.tsv : genome_stats/lst_score.tsv genome_stats/genome_altered.tsv genome_stats/ntai_score.tsv genome_stats/myriad_score.tsv
	$(call RUN,-n 1 -s 4G -m 4G,"set -o pipefail && \
				     			 $(RSCRIPT) modules/summary/genomesummary.R")
				     			 
summary/genome_summary.xlsx : summary/tsv/genome_summary.tsv
	$(call RUN,-n 1 -s 4G -m 4G,"set -o pipefail && \
								 python modules/summary/genome_summary_excel.py")

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
