include innovation-lab/Makefile.inc

LOGDIR ?= log/fusion_catcher.$(NOW)

CACHE = $(HOME)/share/usr/env/fusioncatcher-1.2.0/share/fusioncatcher-1.20/db/current

fusion_catcher : $(foreach sample,$(SAMPLES),fusioncatcher/$(sample)/$(sample).1.fastq.gz) \
		 		 $(foreach sample,$(SAMPLES),fusioncatcher/$(sample)/$(sample).2.fastq.gz) \
		 		 $(foreach sample,$(SAMPLES),fusioncatcher/$(sample)/out/taskcomplete) \
		 		 fusioncatcher/summary.txt

define merged-fastq
fusioncatcher/$1/$1.1.fastq.gz : $$(foreach split,$2,$$(word 1, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 2G -m 4G,"set -o pipefail && \
									 mkdir -p fusioncatcher/$1 && \
									 cp $$(^) $$(@)")
fusioncatcher/$1/$1.2.fastq.gz : $$(foreach split,$2,$$(word 2, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 2G -m 4G,"set -o pipefail && \
									 mkdir -p fusioncatcher/$1 && \
									 cp $$(^) $$(@)")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call merged-fastq,$(sample),$(split.$(sample)))))

define fusion-catcher
fusioncatcher/$1/out/taskcomplete : fusioncatcher/$1/$1.1.fastq.gz fusioncatcher/$1/$1.2.fastq.gz
	$$(call RUN,-c -n 8 -s 2G -m 3G -v $(FUSIONCATCHER_ENV) -w 72:00:00,"set -o pipefail && \
																		 mkdir -p fusioncatcher/$1/out && \
																		 fusioncatcher.py \
																		 -i fusioncatcher/$1 \
																		 -o fusioncatcher/$1/out \
																		 -d $$(CACHE) \
																		 -p 8 && \
																		 echo $1 > fusioncatcher/$1/out/taskcomplete")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call fusion-catcher,$(sample))))
		
fusioncatcher/summary.txt : $(wildcard $(foreach sample,$(SAMPLES),fusioncatcher/$(sample)/out/taskcomplete))
	for i in $(SAMPLES); do \
		sed \"s/$$/\t$$i/\" fusioncatcher/$$i/out/final-list_candidate-fusion-genes.hg19.txt > fusioncatcher/$$i/out/final-list_candidate-fusion-genes.hg19.tmp; \
		cat fusioncatcher/$$i/out/final-list_candidate-fusion-genes.hg19.tmp >> fusioncatcher/summary.txt; \
	done

..DUMMY := $(shell mkdir -p version; \
			 $(PYTHON) --version &> version/fusioncatcher.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: fusion_catcher
