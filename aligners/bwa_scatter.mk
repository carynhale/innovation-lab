include innovation-lab/Makefile.inc
include innovation-lab/config/gatk.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/bwa_scatter.$(NOW)

bwa_scatter : $(foreach sample,$(SAMPLES),bwa_scatter/$(sample)/$(sample)_R1.fastq.gz) \
	      $(foreach sample,$(SAMPLES),bwa_scatter/$(sample)/$(sample)_R2.fastq.gz) \
	      $(foreach sample,$(SAMPLES),bwa_scatter/$(sample)/$(sample)_R1.taskcomplete)
	       
CHUNKS = 10

define merge-fastq
bwa_scatter/$1/$1_R1.fastq.gz : $$(foreach split,$2,$$(word 1, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 4G -m 6G,"zcat $$(^) | gzip -c > $$(@)")
	
bwa_scatter/$1/$1_R2.fastq.gz : $$(foreach split,$2,$$(word 2, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 4G -m 6G,"zcat $$(^) | gzip -c > $$(@)")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call merge-fastq,$(sample),$(split.$(sample)))))

define split-fastq
bwa_scatter/$1/$1_R1.taskcomplete : bwa_scatter/$1/$1_R1.fastq.gz
	$$(call RUN,-c -n $(CHUNKS) -s 1G -m 2G -v $(FASTQ_SPLITTER_ENV),"set -o pipefail && \
									  cd bwa_scatter/$1 && \
									  FASTQ_CHUNKS=($(shell seq 1 $(CHUNKS))) && \
									  ALL_CHUNKS=() && \
									  $(foreach chunk,$(FASTQ_CHUNKS[@]), ALL_CHUNKS+=" -o "$1.1.$(chunk).fastq.gz) && \
									  fastqsplitter -i $1_R1.fastq.gz $(ALL_CHUNKS[@]) && \
									  touch $1_R1.taskcomplete && \
									  cd ../..")
endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call split-fastq,$(sample))))


..DUMMY := $(shell mkdir -p version; \
	     echo "picard" >> version/bwa_scatter.txt; \
	     $(PICARD) SortSam --version &>> version/bwa_scatter.txt; \
	     $(PICARD) MarkIlluminaAdapters --version &>> version/bwa_scatter.txt; \
	     $(PICARD) SamToFastq --version &>> version/bwa_scatter.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: bwa_scatter

