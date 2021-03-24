include innovation-lab/Makefile.inc

LOGDIR ?= log/subsample_fastq.$(NOW)

subsample_fastq : $(foreach sample,$(SAMPLES),FASTQ_DOWNSAMPLE/fastq/$(sample).taskcomplete)

SEED = 1
TARGETREADS = 1000000 \
	      2000000 \
	      5000000 \
	      10000000 \
	      15000000 \
	      20000000 \
	      25000000 \
	      30000000 \
	      40000000 \
	      50000000

define copy-fastq
FASTQ_DOWNSAMPLE/fastq/$1_R1--0.fastq.gz : $3
	$$(call RUN,-c -n 1 -s 2G -m 4G,"set -o pipefail && \
					 mkdir -p FASTQ_DOWNSAMPLE/fastq && \
					 $(RSCRIPT) $(SCRIPTS_DIR)/fastq_tools/copy_fastq.R \
					 --sample_name $1 \
					 --directory_name FASTQ_DOWNSAMPLE/fastq \
					 --fastq_files '$$^' && \
					 mv FASTQ_DOWNSAMPLE/fastq/$1_R1.fastq.gz FASTQ_DOWNSAMPLE/fastq/$1_R1--0.fastq.gz && \
					 mv FASTQ_DOWNSAMPLE/fastq/$1_R2.fastq.gz FASTQ_DOWNSAMPLE/fastq/$1_R2--0.fastq.gz")

endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)),$(eval $(call copy-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))

define sample-fastq
FASTQ_DOWNSAMPLE/fastq/$1.taskcomplete : FASTQ_DOWNSAMPLE/fastq/$1_R1--0.fastq.gz
	$$(call RUN, -c -n 1 -s 6G -m 12G -v $(SEQTK_ENV),"set -o pipefail && \
							   for x in $$(TARGETREADS); do $$(SEQTK) sample -s $(SEED) FASTQ_DOWNSAMPLE/fastq/$1_R1--0.fastq.gz $(x) > FASTQ_DOWNSAMPLE/fastq/$1_R1--$(x).fastq.gz; done && \
							   echo $1 >  FASTQ_DOWNSAMPLE/fastq/$1.taskcomplete")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call sample-fastq,$(sample))))


..DUMMY := $(shell mkdir -p version; \
	     $(HOME)/share/usr/env/seqtk-1.3/bin/seqtk &> version/subsample_fastq.txt)
.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY: subsample_fastq
