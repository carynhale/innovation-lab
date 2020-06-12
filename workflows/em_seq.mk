include innovation-lab/Makefile.inc

LOGDIR ?= log/em_seq.$(NOW)

em_seq : $(foreach sample,$(SAMPLES),emseq/$(sample)/$(sample)_R1.fastq.gz) \
		 $(foreach sample,$(SAMPLES),emseq/$(sample)/$(sample).fastq) \
		 $(foreach sample,$(SAMPLES),emseq/$(sample)/$(sample).txt) \
		 $(foreach sample,$(SAMPLES),emseq/$(sample)/qual.txt) \
		 $(foreach sample,$(SAMPLES),emseq/$(sample)/seq.txt)

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
	
define extract-fastq
emseq/$1/$1.fastq : emseq/$1/$1_R1.fastq.gz
	$$(call RUN,-c -n 1 -s 4G -m 6G,"set -o pipefail && \
									 zcat $$(<) > $$(@)")
									 
emseq/$1/$1.txt : emseq/$1/$1.fastq
	$$(call RUN,-c -n 1 -s 4G -m 6G,"set -o pipefail && \
									 $(SCRIPTS_DIR)/fastq_tools/extract_fastq.sh 1 $$(<) $$(@)")
									 
emseq/$1/qual.txt : emseq/$1/$1.txt
	$$(call RUN,-c -n 1 -s 4G -m 6G,"set -o pipefail && \
									 $(SCRIPTS_DIR)/fastq_tools/extract_fastq.sh 1 $$(<) $$(@)")

emseq/$1/seq.txt : emseq/$1/$1.txt
	$$(call RUN,-c -n 1 -s 4G -m 6G,"set -o pipefail && \
									 $(SCRIPTS_DIR)/fastq_tools/extract_fastq.sh 2 $$(<) $$(@)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call extract-fastq,$(sample))))
		
..DUMMY := $(shell mkdir -p version; \
			 R --version >> version.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: em_seq
