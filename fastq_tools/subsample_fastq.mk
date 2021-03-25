include innovation-lab/Makefile.inc

LOGDIR ?= log/subsample_fastq.$(NOW)

SEED = 1
TARGET_READS = 1000000 \
	       2000000 \
	       3000000 \
	       4000000 \
	       5000000 \
	       7000000 \
	       10000000 \
	       15000000 \
	       20000000 \
	       25000000 \
	       30000000 \
	       40000000 \
	       50000000

subsample_fastq : $(foreach sample,$(SAMPLES),FASTQ_DOWNSAMPLE/$(sample)/$(sample)_R1.fastq.gz) \
		  $(foreach sample,$(SAMPLES), \
		  	$(foreach targetreads,$(TARGET_READS),FASTQ_DOWNSAMPLE/$(sample)/$(sample)_R1--$(targetreads).fastq.gz))

define copy-fastq
FASTQ_DOWNSAMPLE/$1/$1_R1.fastq.gz : $3
	$$(call RUN,-c -n 1 -s 2G -m 4G,"set -o pipefail && \
					 mkdir -p FASTQ_DOWNSAMPLE/$1 && \
					 $(RSCRIPT) $(SCRIPTS_DIR)/fastq_tools/copy_fastq.R \
					 --sample_name $1 \
					 --directory_name FASTQ_DOWNSAMPLE \
					 --fastq_files '$$^'")

endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)),$(eval $(call copy-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))

define sample-fastq
FASTQ_DOWNSAMPLE/$1/$1_R1--$2.fastq.gz : FASTQ_DOWNSAMPLE/$1/$1_R1.fastq.gz
	$$(call RUN, -c -s 12G -m 24G -v $(SEQTK_ENV),"set -o pipefail && \
						       $$(SEQTK) sample -s $(SEED) FASTQ_DOWNSAMPLE/$1/$1_R1.fastq.gz $2 > FASTQ_DOWNSAMPLE/$1/$1_R1--$2.fastq.gz")

endef
$(foreach sample,$(SAMPLES), \
	$(foreach targetreads,$(TARGET_READS), \
		$(eval $(call sample-fastq,$(sample),$(targetreads)))))


..DUMMY := $(shell mkdir -p version; \
	     $(HOME)/share/usr/env/seqtk-1.3/bin/seqtk &> version/subsample_fastq.txt)
.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY: subsample_fastq
