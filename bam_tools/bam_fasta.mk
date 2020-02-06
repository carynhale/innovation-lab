include modules/Makefile.inc

LOGDIR ?= log/bam_fasta.$(NOW)

bam_fasta : $(foreach sample,$(SAMPLES),unmapped_reads/$(sample).fasta)

define bam-to-fasta
unmapped_reads/%.fasta : unmapped_reads/%.bam
	$$(call RUN,-c -n 4 -s 4G -m 6G,"$(SAMTOOLS2) fasta $$(<) > unmapped_reads/$$(*).fasta")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call bam-to-fasta,$(sample))))


.DUMMY := $(shell mkdir -p version; $(SAMTOOLS2) &> version/samtools.txt )
.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY: bam_fasta
