include innovation-lab/Makefile.inc

LOGDIR = log/kallisto.$(NOW)

kallisto : $(foreach sample,$(SAMPLES),kallisto/$(sample)/$(sample).1.fastq.gz) \
		   $(foreach sample,$(SAMPLES),kallisto/$(sample)/taskcomplete)

define bam-to-fastq
kallisto/$1/$1.1.fastq.gz : bam/$1.bam
	$$(call RUN,-n 4 -s 4G -m 9G,"set -o pipefail && \
								  mkdir -p kallisto/$1 && \
								  $$(SAMTOOLS) sort -T kallisto/$1/$1 -O bam -n -@ 4 -m 6G $$(<) | \
								  bedtools bamtofastq -i - -fq >(gzip -c > kallisto/$1/$1.1.fastq.gz) -fq2 >(gzip -c > kallisto/$1/$1.2.fastq.gz)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call bam-to-fastq,$(sample))))

define fastq-to-kallisto
kallisto/$1/taskcomplete : kallisto/$1/$1.1.fastq.gz
	$$(call RUN,-c -n 8 -s 2G -m 3G -v $(KALLISTO_ENV),"set -o pipefail && \
														kallisto -quant \
														-i $$(KALLISTO_INDEX) \
														-o kallisto/$1 \
														kallisto/$1/$1.1.fastq.gz kallisto/$1/$1.2.fastq.gz \
														--bias \
														-b 100 \
														-t 8 \
														echo $1 > kallisto/$1/taskcomplete")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call fastq-to-kallisto,$(sample))))

..DUMMY := $(shell mkdir -p version; \
			 $(SAMTOOLS) --version > version/kallisto.txt; \
			 ./share/usr/env/kallisto-0.46.2/bin/kallisto version >> version/kallisto.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: kallisto
