include innovation-lab/Makefile.inc

LOGDIR ?= log/subsample_fastq.$(NOW)

subsample_fastq : $(foreach sample,$(SAMPLES),FASTQ_DOWNSAMPLE/fastq/$(sample).1.fastq.gz)

fastq/%.1.fastq.gz fastq/%.2.fastq.gz : %.bam
	$(call RUN,-n 4 -s 4G -m 9G,"$(SAMTOOLS) sort -T $(<D)/$* -O bam -n -@ 4 -m 6G $< | \
				     $(SAMTOOLS) fastq -f 1 -1 >(gzip -c > fastq/$*.1.fastq.gz) -2 >(gzip -c > fastq/$*.2.fastq.gz) -")

..DUMMY := $(shell mkdir -p version; \
	     $(SAMTOOLS) --version > version/extract_fastq.txt)
.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY: extract_fastq
