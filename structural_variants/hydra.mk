include modules/Makefile.inc
include modules/config/gatk.inc

LOGDIR = log/hydra.$(NOW)

HYDRA = $(HOME)/share/usr/bin/hydra
override HYDRA_OPTS ?= -mld 500 -mn 1500
BAM_TO_FASTQ = $(HOME)/share/usr/bin/bamToFastq
BAM_TO_BED = /opt/common/bedtools/bedtools-2.17.0/bin/bamToBed
DEDUP_DISCORDANTS = $(HOME)/share/usr/bin/dedupDiscordants.py
PAIR_DISCORDANTS = $(HOME)/share/usr/bin/pairDiscordants.py


.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all

all : $(foreach sample,$(SAMPLES),hydra/breaks/$(sample).breaks)

hydra/bam/%.disc.bam : bam/%.bam
	$(call RUN,,"$(SAMTOOLS) view -bF 2 $< > $@")

hydra/bed/%.disc.bedpe : hydra/bam/%.disc.bam
	$(call RUN,,"$(BAM_TO_BED) -i $< -tag NM | $(PAIR_DISCORDANTS) -i stdin -m hydra -z 800 > $@")

hydra/bed/%.disc.dedup.bedpe : hydra/bed/%.disc.bedpe
	$(call RUN,,"$(DEDUP_DISCORDANTS) -i $< -s 3 > $@")

hydra/breaks/%.breaks : hydra/bed/%.disc.dedup.bedpe
	$(call RUN,,"$(HYDRA) -in $< -out $@ $(HYDRA_OPTS)")
