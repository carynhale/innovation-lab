include innovation-lab/Makefile.inc

LOGDIR ?= log/fastqc.$(NOW)

fast_qc : $(foreach sample,$(SAMPLES),fastqc/$(sample)/taskcomplete.txt)

define fast-qc
fastqc/$1/taskcomplete.txt : $3
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(FASTQC_ENV),"set -o pipefail && \
							  mkdir -p fastqc/$1 && \
							  $(FASTQC) \
							  -o fastqc/$1 \
							  '$$^' && \
							  touch fastqc/$1/taskcomplete.txt")

endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)),$(eval $(call copy-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))

..DUMMY := $(shell mkdir -p version; \
	     $(FASTQC_ENV)/bin/$(FASTQC) --version &> version/fastqc.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: fast_qc
