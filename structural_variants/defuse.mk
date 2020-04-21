include innovation-lab/Makefile.inc

LOGDIR ?= log/defuse.$(NOW)

defuse : $(foreach sample,$(SAMPLES),defuse/$(sample).1.fastq) \
		 $(foreach sample,$(SAMPLES),defuse/$(sample).2.fastq) \
		 $(foreach sample,$(SAMPLES),defuse/$(sample).results.filtered.tsv) \
		 $(foreach sample,$(SAMPLES),defuse/$(sample).taskcomplete)
		 
DEFUSE_CONFIG = innovation-lab/config/defuse.inc
DEFUSE_E75 = /home/brownd7/share/lib/resource_files/defuse/homo_sapiens/Ensembl/Grch37.p13/Sequence/defuse_e75/

define merged-fastq
defuse/$1.1.fastq : $$(foreach split,$2,$$(word 1, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 2G -m 4G,"zcat $$(^) > $$(@)")
defuse/$1.2.fastq : $$(foreach split,$2,$$(word 2, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 2G -m 4G,"zcat $$(^) > $$(@)")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call merged-fastq,$(sample),$(split.$(sample)))))

define run-defuse
defuse/%.results.filtered.tsv : defuse/%.1.fastq defuse/%.2.fastq
	$$(call RUN,-c -n 10 -s 2G -m 3G -w 72:00:00 -v $(DEFUSE_ENV),"set -o pipefail && \
																   mkdir -p defuse && \
																   $$(DEFUSE) \
																   --config $$(DEFUSE_CONFIG) \
																   --dataset $$(DEFUSE_E75) \
																   --output defuse/$$(*) \
																   --res defuse/$$(*).results.tsv \
																   --rescla defuse/$$(*).results.classify.tsv \
																   --resfil defuse/$$(*).results.filtered.tsv \
																   -1 defuse/$$(*).1.fastq \
																   -2 defuse/$$(*).2.fastq \
																   -s direct \
																   -p 10")
	
defuse/%.taskcomplete : defuse/%.results.filtered.tsv
	$$(call RUN,-c -s 1G -m 2G,"echo $$(*) > defuse/$$(*).taskcomplete")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call run-defuse,$(sample))))
		
# defuse/summary.tsv : $(wildcard $(foreach sample,$(TUMOR_SAMPLES),defuse/$(sample).taskcomplete))
#	$(call RUN,-c -n 1 -s 6G -m 8G,"set -o pipefail && \
#				     			    $(RSCRIPT) $(SCRIPTS_DIR)/structural_variants/defuse.R --samples '$(TUMOR_SAMPLES)'")		
		
..DUMMY := $(shell mkdir -p version; \
			 $(PERL) --version > version/defuse.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: defuse
