include innovation-lab/Makefile.inc

LOGDIR ?= log/msk_access.$(NOW)

fgbio_access : $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_R1.fastq.gz)

define copy-fastq
fgbio/$1/$1_R1.fastq.gz : $3
	$$(call RUN,-c -n 1 -s 2G -m 4G,"set -o pipefail && \
									 mkdir -p fgbio/$1 && \
								     $(RSCRIPT) $(SCRIPTS_DIR)/fastq_tools/copy_fastq.R --sample_name $1 --directory_name 'fgbio' --fastq_files '$$^'")

endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)),$(eval $(call copy-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))



..DUMMY := $(shell mkdir -p version)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: fgbio_access
