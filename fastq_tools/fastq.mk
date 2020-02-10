include innovation-lab/Makefile.inc

LOGDIR ?= log/fastq.$(NOW)

VPATH ?= unprocessed_bam
EXTRACT_TOOL ?= $(PICARD)
FASTQ_TRIMMER = $(PERL) $(SCRIPTS_DIR)/fastq_tools/trim_fastq.pl

ifeq ($(TRIM_READS),true)
   FASTQ_FILTER := trim
   TRIM_LENGTH ?= 150
   TRIM_OPTS ?= -l $(TRIM_LENGTH)
endif

ifeq ($(MERGE_SPLIT_FASTQ),true)
fastq: $(foreach sample,$(SAMPLES),fastq/$(sample).1.fastq.gz)
else
fastq: $(foreach split,$(UNSPLIT_SAMPLES),fastq/$(split).1.fastq.gz)
endif

ifdef FASTQ_FILTER
fastq/%.1.fastq.gz fastq/%.2.fastq.gz : unprocessed_fastq/%.1.$(FASTQ_FILTER).fastq.gz unprocessed_fastq/%.2.$(FASTQ_FILTER).fastq.gz
	$(INIT) ln $< fastq/$*.1.fastq.gz && ln $(word 2,$(^)) fastq/$*.2.fastq.gz && cp  $< fastq/$*.1.fastq.gz && cp $(word 2,$^) fastq/$*.2.fastq.gz
else
fastq/%.1.fastq.gz fastq/%.2.fastq.gz : unprocessed_fastq/%.1.fastq.gz unprocessed_fastq/%.2.fastq.gz
	$(INIT) ln $< fastq/$*.1.fastq.gz && ln $(word 2,$(^)) fastq/$*.2.fastq.gz && cp  $< fastq/$*.1.fastq.gz && cp $(word 2,$^) fastq/$*.2.fastq.gz
endif

fastq/%.1.fastq.gz fastq/%.2.fastq.gz : unprocessed_bam/%.nsorted.bam
	$(call RUN,-s 10G -m 15G,"TEMP=`mktemp`; mkfifo \$${TEMP}_1; mkfifo \$${TEMP}_2; \
							  gzip < \$${TEMP}_1 > fastq/$*.1.fastq.gz & \
							  gzip < \$${TEMP}_2 > fastq/$*.2.fastq.gz & \
							  $(call PICARD,SamToFastq,9G) QUIET=true I=$< FASTQ=\$${TEMP}_1 SECOND_END_FASTQ=\$${TEMP}_2")

%.nsorted.bam : %.bam
	$(call RUN,-s 20G -m 25G,"$(call PICARD,SortSam,19G) I=$< O=$@ SO=queryname")

unprocessed_fastq/%.trim.fastq.gz : unprocessed_fastq/%.fastq.gz
	$(call RUN,-s 2G -m 3G,"zcat $< | $(FASTQ_TRIMMER) $(TRIM_OPTS) | gzip -c > $@ ")

unprocessed_fastq/%.readtrim.1.fastq.gz unprocessed_fastq/%.readtrim.2.fastq.gz : %.bam %.read_len
	$(call RUN,-s 22G -m 32G,"NUM_READS=`awk '{ sum += $$1 } END { print sum }' $(word 2,$^)`; \
							  MAX_LENGTH=`sort -k 2 $(word 2,$^) | awk -v nreads="$$NUM_READS" '$$1 / nreads > 0.4 { print $$2 }' | head -1`; \
							  if [ "$$MAX_LENGTH" = "" ]; then MAX_LENGTH=`cut -d' ' -f 2 $(word 2,$^) | head -1`; fi; \
							  TEMP=`mktemp`; mkfifo \$${TEMP}_1; mkfifo \$${TEMP}_2; \
							  gzip < \$${TEMP}_1 > fastq/$*.readtrim.1.fastq.gz & \
							  gzip < \$${TEMP}_2 > fastq/$*.readtrim.2.fastq.gz & \
							  $(call SAM_TO_FASTQ_MEM,21G) I=$< FASTQ=\$${TEMP}_1 SECOND_END_FASTQ=\$${TEMP}_2 READ1_MAX_BASES_TO_WRITE=\$$MAX_LENGTH READ2_MAX_BASES_TO_WRITE=\$$MAX_LENGTH")

%.read_len : %.bam
	$(call RUN,-s 4G -m 6G,"$(SAMTOOLS) view $< | awk '{ print length($$10) }' | sort -n | uniq -c | sort -rn | sed 's/^ \+//' | awk ' > $@")

define merged-fastq
unprocessed_fastq/$1.%.fastq.gz : $$(foreach split,$2,unprocessed_fastq/$$(split).%.fastq.gz)
	$$(INIT) zcat $$(^) | gzip > $$@ 2> $$(LOG)
unprocessed_fastq/$1.%.fastq.gz : $$(foreach split,$2,unprocessed_fastq/$$(split).%.fastq)
	$$(INIT) cat $$(^) | gzip > $$@ 2> $$(LOG)
endef
$(foreach sample,$(SPLIT_SAMPLES),$(eval $(call merged-fastq,$(sample),$(split.$(sample)))))

.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY : fastq
