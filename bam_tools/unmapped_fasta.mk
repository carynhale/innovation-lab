include innovation-lab/Makefile.inc

LOGDIR ?= log/unmapped_fasta.$(NOW)

unmapped_fasta : $(foreach sample,$(SAMPLES),unmapped_reads/$(sample).fasta)

define extract-unmapped-reads
unmapped_reads/%.fasta : unmapped_reads/%.bam
	$$(call RUN,-c -n 4 -s 4G -m 6G,"$(SAMTOOLS2) fasta $$(<) > unmapped_reads/$$(*).fasta")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call extract-unmapped-reads,$(sample))))


..DUMMY := $(shell mkdir -p version; $(SAMTOOLS2) &> version/unmapped_fasta.txt )
.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY: unmapped_fasta
