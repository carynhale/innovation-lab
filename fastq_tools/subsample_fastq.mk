include innovation-lab/Makefile.inc

LOGDIR ?= log/subsample_fastq.$(NOW)

SEED = 1
N = 1 2 3 4 5
    
subsample_fastq : $(foreach sample,$(SAMPLES),FASTQ_DOWNSAMPLE/fastq/$(sample)/$(sample)_R1.fastq.gz) \
		  $(foreach sample,$(SAMPLES),FASTQ_DOWNSAMPLE/fastq/$(sample)/$(sample)_R2.fastq.gz) \
		  $(foreach sample,$(SAMPLES), \
		  	$(foreach n,$(N),FASTQ_DOWNSAMPLE/fastq/$(sample)/$(sample)_R1--$(n).fastq.gz)) \
		  $(foreach sample,$(SAMPLES), \
		  	$(foreach n,$(N),FASTQ_DOWNSAMPLE/fastq/$(sample)/$(sample)_R2--$(n).fastq.gz))

define merge-fastq
FASTQ_DOWNSAMPLE/fastq/$1/$1_R1.fastq.gz : $$(foreach split,$2,$$(word 1, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 4G -m 6G -w 72:00:00,"zcat $$(^) | gzip -c > $$(@)")
	
FASTQ_DOWNSAMPLE/fastq/$1/$1_R2.fastq.gz : $$(foreach split,$2,$$(word 2, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 4G -m 6G -w 72:00:00,"zcat $$(^) | gzip -c > $$(@)")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call merge-fastq,$(sample),$(split.$(sample)))))

define sample-fastq
FASTQ_DOWNSAMPLE/fastq/$1/$1_R1--$2.fastq.gz : FASTQ_DOWNSAMPLE/fastq/$1/$1_R1.fastq.gz
	$$(call RUN, -c -s 12G -m 24G -v $(SEQTK_ENV),"set -o pipefail && \
						       $$(SEQTK) sample -s $(SEED) $$(<) '$${READS.$2}' > FASTQ_DOWNSAMPLE/fastq/$1/$1_R1--$2.fastq && \
						       gzip FASTQ_DOWNSAMPLE/fastq/$1/$1_R1--$2.fastq")

FASTQ_DOWNSAMPLE/fastq/$1/$1_R2--$2.fastq.gz : FASTQ_DOWNSAMPLE/fastq/$1/$1_R2.fastq.gz
	$$(call RUN, -c -s 12G -m 24G -v $(SEQTK_ENV),"set -o pipefail && \
						       $$(SEQTK) sample -s $(SEED) $$(<) '$${READS.$2}' > FASTQ_DOWNSAMPLE/fastq/$1/$1_R2--$2.fastq && \
						       gzip FASTQ_DOWNSAMPLE/fastq/$1/$1_R2--$2.fastq")

endef
$(foreach sample,$(SAMPLES), \
	$(foreach n,$(N), \
		$(eval $(call sample-fastq,$(sample),$(n)))))


..DUMMY := $(shell mkdir -p version; \
	     $(HOME)/share/usr/env/seqtk-1.3/bin/seqtk &> version/subsample_fastq.txt)
.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY: subsample_fastq
