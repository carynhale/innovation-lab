include modules/Makefile.inc

LOGDIR ?= log/merge_split_fastq.$(NOW)
.PHONY : fastq

fastq: $(foreach sample,$(SAMPLES),fastq/$(sample).1.fastq.gz fastq/$(sample).2.fastq.gz)

define merged-fastq
fastq/$1.%.fastq.gz : $$(foreach split,$2,fastq/$$(split).%.fastq.gz)
	$$(call RUN,-c -n 1 -s 2G -m 4G,"zcat $$(^) | gzip -c > $$(@)")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call merged-fastq,$(sample),$(split.$(sample)))))

define merged-fastq2
fastq/$1.1.fastq.gz : $$(foreach split,$2,$$(word 1, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 2G -m 4G,"zcat $$(^) | gzip -c > $$(@)")
fastq/$1.2.fastq.gz : $$(foreach split,$2,$$(word 2, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 2G -m 4G,"zcat $$(^) | gzip -c > $$(@)")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call merged-fastq2,$(sample),$(split.$(sample)))))


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
