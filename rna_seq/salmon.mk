include innovation-lab/Makefile.inc

LOGDIR = log/salmon.$(NOW)

salmon : $(foreach sample,$(SAMPLES),salmon/$(sample)/$(sample).1.fastq.gz) \
	 $(foreach sample,$(SAMPLES),kallisto/$(sample)/quant.sf)

define bam-to-fastq
salmon/$1/$1.1.fastq.gz : bam/$1.bam
	$$(call RUN,-n 4 -s 4G -m 9G,"set -o pipefail && \
				      mkdir -p salmon/$1 && \
				      $$(SAMTOOLS) sort -T salmon/$1/$1 -O bam -n -@ 4 -m 6G $$(<) | \
				      bedtools bamtofastq -i - -fq salmon/$1/$1.1.fastq -fq2 salmon/$1/$1.2.fastq && \
				      gzip salmon/$1/$1.1.fastq && \
				      gzip salmon/$1/$1.2.fastq")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call bam-to-fastq,$(sample))))

define fastq-to-salmon
salmon/$1/quant.sf : salmon/$1/$1.1.fastq.gz
	$$(call RUN,-c -n 12 -s 2G -m 3G -v $(SALMON_ENV),"set -o pipefail && \
							   salmon quant \
							   -i $$(SALMON_INDEX) \
							   -l A \
							   -1 salmon/$1/$1.1.fastq.gz \
							   -2 salmon/$1/$1.2.fastq.gz \
							   -p 12 \
							   --validateMappings \
							   -o salmon/$1/quant.sf")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call fastq-to-salmon,$(sample))))
	

..DUMMY := $(shell mkdir -p version; \
	     $(SAMTOOLS) --version > version/salmon.txt; \
	     ~/share/usr/env/salmon-1.4.0/bin/salmon version >> version/salmon.txt; \)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: salmon
