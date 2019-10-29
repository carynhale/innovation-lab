include modules/Makefile.inc
include modules/fastq_tools/merge_split_fastq.mk

LOGDIR ?= log/fusion_catcher.$(NOW)
PHONY += FUSION_CATCHER_WORKFLOW

fusion_catcher : $(foreach sample,$(SAMPLES),fusioncatcher/$(sample)/$(sample).taskcomplete)

FUSION_CATCHER_WORKFLOW += fastq
FUSION_CATCHER_WORKFLOW += fusion_catcher

FUSIONCATCHER = $(HOME)/share/usr/fusioncatcher/bin/fusioncatcher
FUSIONCATCHER_OPTS = -d $(HOME)/share/usr/fusioncatcher/data/current --extract-buffer-size=35000000000

define fusion-catcher
fusion_catcher/%/%.taskcomplete : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$(call RUN,-n 8 -s 1G -m 4G,"$(FUSIONCATCHER) $(FUSIONCATCHER_OPTS) -p 8 -o $(@D)/$* -i $<$(,)$(<<) && touch $@")
	
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call fusion-catcher,$(sample))))
					

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)


define defuse-single-sample
fastq/%.1.fastq fastq/%.2.fastq : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$$(call RUN,-c -s 4G -m 9G,"gzip -d $$(<) $$(<<)")

defuse/tables/%.results.filtered.tsv : fastq/%.1.fastq fastq/%.2.fastq
	$$(call RUN,-c -n 10 -s 3G -m 4G -w 540,"$$(PERL) $$(DEFUSE_SCRIPTS)/defuse_run.pl -c $$(CONFIG) -d $$(DEFUSE75) -o defuse/$$* --res defuse/tables/$$*.results.tsv -resfil defuse/tables/$$*.results.filtered.tsv -1 fastq/$$*.1.fastq -2 fastq/$$*.2.fastq -p 10 -s direct")
	
defuse/%.taskcomplete : defuse/tables/%.results.filtered.tsv
	$$(call RUN,-c -s 1G -m 3G,"echo $$* > defuse/$$*.taskcomplete")
endef