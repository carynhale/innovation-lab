include innovation-lab/Makefile.inc

LOGDIR ?= log/merge_split_fastq.$(NOW)

fastq: $(foreach sample,$(SAMPLES),fastq/$(sample).1.fastq.gz fastq/$(sample).2.fastq.gz)

define merged-fastq
fastq/$1.1.fastq.gz : $$(foreach split,$2,$$(word 1, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 4G -m 6G -w 72:00:00,"zcat $$(^) | gzip -c > $$(@)")
fastq/$1.2.fastq.gz : $$(foreach split,$2,$$(word 2, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 4G -m 6G -w 72:00:00,"zcat $$(^) | gzip -c > $$(@)")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call merged-fastq,$(sample),$(split.$(sample)))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : fastq
