include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvkit.$(NOW)
PHONY += cnvkit cnvkit/cnn cnvkit/cnn/tumor cnvkit/cnn/normal cnvkit/reference cnvkit/cnr cnvkit/log2 cnvkit/segmented cnvkit/called cnvkit/summary


cnvkit : $(foreach sample,$(TUMOR_SAMPLES),cnvkit/cnn/tumor/$(sample).targetcoverage.cnn) \
		 $(foreach sample,$(TUMOR_SAMPLES),cnvkit/cnn/tumor/$(sample).antitargetcoverage.cnn) \
		 $(foreach sample,$(NORMAL_SAMPLES),cnvkit/cnn/normal/$(sample).targetcoverage.cnn) \
		 $(foreach sample,$(NORMAL_SAMPLES),cnvkit/cnn/normal/$(sample).antitargetcoverage.cnn) \
		 cnvkit/reference/combined_reference.cnr \
		 $(foreach sample,$(TUMOR_SAMPLES),cnvkit/cnr/$(sample).cnr) \
		 $(foreach sample,$(TUMOR_SAMPLES),cnvkit/log2/$(sample).ontarget.pdf)

define cnvkit-tumor-cnn
cnvkit/cnn/tumor/%.targetcoverage.cnn : bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(ONTARGET_FILE) -o cnvkit/cnn/tumor/$$(*).targetcoverage.cnn")

cnvkit/cnn/tumor/%.antitargetcoverage.cnn : bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(OFFTARGET_FILE) -o cnvkit/cnn/tumor/$$(*).antitargetcoverage.cnn")
endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvkit-tumor-cnn,$(sample))))
		
define cnvkit-normal-cnn
cnvkit/cnn/normal/%.targetcoverage.cnn : bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(ONTARGET_FILE) -o cnvkit/cnn/normal/$$(*).targetcoverage.cnn")

cnvkit/cnn/normal/%.antitargetcoverage.cnn : bam/%.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"cnvkit.py coverage -p 4 -q 0 $$(<) $$(OFFTARGET_FILE) -o cnvkit/cnn/normal/$$(*).antitargetcoverage.cnn")
endef
 $(foreach sample,$(NORMAL_SAMPLES),\
		$(eval $(call cnvkit-normal-cnn,$(sample))))

cnvkit/reference/combined_reference.cnr : $(wildcard cnvkit/cnn/normal/$(NORMAL_SAMPLES).targetcoverage.cnn) $(wildcard cnvkit/cnn/normal/$(NORMAL_SAMPLES).antitargetcoverage.cnn)
	$(call RUN,-n 1 -s 24G -m 32G,"cnvkit.py reference cnvkit/cnn/normal/*.cnn -f $(REF_FASTA) --no-edge -o cnvkit/reference/combined_reference.cnr")

define cnvkit-cnr
cnvkit/cnr/%.cnr : cnvkit/cnn/tumor/%.targetcoverage.cnn cnvkit/cnn/tumor/%.antitargetcoverage.cnn cnvkit/reference/combined_reference.cnr
	$$(call RUN,-c -s 6G -m 8G,"cnvkit.py fix $$(<) $$(<<) cnvkit/reference/combined_reference.cnr -o cnvkit/cnr/$$(*).cnr")
	
endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvkit-cnr,$(sample))))

define cnvkit-plot
cnvkit/log2/%.ontarget.pdf : cnvkit/cnr/%.cnr
	$$(call RUN,-c -v $(ASCAT_ENV) -s 4G -m 6G,"$(RSCRIPT) $(SCRIPTS_DIR)/copy_number/cnvkit.R --type 'plot-log2' --sample_name $$(*)")
endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvkit-plot,$(sample))))
		
include modules/copy_number/cnvkitsegment.mk
include modules/copy_number/cnvkitsummary.mk

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
