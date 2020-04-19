include innovation-lab/Makefile.inc

LOGDIR ?= log/defuse.$(NOW)

defuse : $(foreach sample,$(SAMPLES),defuse/$(sample).1.fastq) \
		 $(foreach sample,$(SAMPLES),defuse/$(sample).2.fastq)
#		 $(foreach sample,$(SAMPLES),defuse/$(sample).results.filtered.tsv) \
#		 $(foreach sample,$(SAMPLES),defuse/$(sample).taskcomplete) \
#		 defuse/summary.tsv

define merged-fastq
defuse/$1.1.fastq : $$(foreach split,$2,$$(word 1, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 2G -m 4G,"zcat $$(^) > $$(@)")
defuse/$1.2.fastq : $$(foreach split,$2,$$(word 2, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 2G -m 4G,"zcat $$(^) > $$(@)")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call merged-fastq,$(sample),$(split.$(sample)))))

#defuse/%.results.filtered.tsv : defuse/%.1.fastq defuse/%.2.fastq
#	$$(call RUN,-c -n 10 -s 3G -m 4G -w 540,"$$(PERL) $$(DEFUSE_SCRIPTS)/defuse_run.pl \
#											 -c $$(CONFIG) -d $$(DEFUSE75) -o defuse/$$(*) \
#											 --res defuse/$$(*).results.tsv \
#											 -resfil defuse/$$(*).results.filtered.tsv \
#											 -1 defuse/$$(*).1.fastq \
#											 -2 defuse/$$(*).2.fastq -p 10 -s direct")
#	
#defuse/%.taskcomplete : defuse/%.results.filtered.tsv
#	$$(call RUN,-c -s 1G -m 2G,"echo $$(*) > defuse/$$(*).taskcomplete")
#
#endef
#
#$(foreach sample,$(SAMPLES),\
#		$(eval $(call defuse-single-sample,$(sample))))
		
# defuse/summary.tsv : $(wildcard $(foreach sample,$(TUMOR_SAMPLES),defuse/$(sample).taskcomplete))
#	$(call RUN,-c -n 1 -s 6G -m 8G,"set -o pipefail && \
#				     			    $(RSCRIPT) $(SCRIPTS_DIR)/structural_variants/defuse.R --samples '$(TUMOR_SAMPLES)'")		
		
..DUMMY := $(shell mkdir -p version; \
			 $(PERL) --version > version/defuse.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: defuse

