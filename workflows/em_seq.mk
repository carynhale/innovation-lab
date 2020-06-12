include innovation-lab/Makefile.inc

LOGDIR ?= log/em_seq.$(NOW)

em_seq : $(foreach sample,$(SAMPLES),emseq/$(sample)/$(sample)_R1.fastq.gz) \

define copy-fastq
emseq/$1/$1_R1.fastq.gz : $3
	$$(call RUN,-c -n 1 -s 2G -m 4G,"set -o pipefail && \
									 mkdir -p emseq/$1 && \
								     $(RSCRIPT) $(SCRIPTS_DIR)/fastq_tools/copy_fastq.R \
								     --sample_name $1 \
								     --directory_name emseq \
								     --fastq_files '$$^'")

endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)),$(eval $(call copy-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))
 
..DUMMY := $(shell mkdir -p version; \
			 R --version >> version.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: em_seq
