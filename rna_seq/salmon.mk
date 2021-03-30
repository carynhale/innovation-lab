include innovation-lab/Makefile.inc

LOGDIR = log/salmon.$(NOW)

salmon : $(foreach sample,$(SAMPLES),salmon/$(sample)/$(sample).1.fastq.gz)

define bam-to-fastq
salmon/$1/$1.1.fastq : bam/$1.bam
	$$(call RUN,-n 4 -s 4G -m 9G,"set -o pipefail && \
				      mkdir -p salmon/$1 && \
				      $$(SAMTOOLS) sort -T salmon/$1/$1 -O bam -n -@ 4 -m 6G $$(<) | \
				      bedtools bamtofastq -i - -fq salmon/$1/$1.1.fastq -fq2 salmon/$1/$1.2.fastq && \
				      gzip salmon/$1/$1.1.fastq && \
				      gzip salmon/$1/$1.2.fastq")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call bam-to-fastq,$(sample))))

..DUMMY := $(shell mkdir -p version; \
	     $(SAMTOOLS) --version > version/salmon.txt;)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: salmon
