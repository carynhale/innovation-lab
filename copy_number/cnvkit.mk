include innovation-lab/Makefile.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/cnv_kit.$(NOW)

cnv_kit : $(foreach sample,$(TUMOR_SAMPLES),cnvkit/cnn/tumor/$(sample).targetcoverage.cnn) \
	  $(foreach sample,$(TUMOR_SAMPLES),cnvkit/cnn/tumor/$(sample).antitargetcoverage.cnn) \
	  $(foreach sample,$(NORMAL_SAMPLES),cnvkit/cnn/normal/$(sample).targetcoverage.cnn) \
	  $(foreach sample,$(NORMAL_SAMPLES),cnvkit/cnn/normal/$(sample).antitargetcoverage.cnn)
	  
ONTARGET_FILE = $(HOME)/share/lib/bed_files/MSK-IMPACT-v3_cnvkit_ontarget.bed
OFFTARGET_FILE = $(HOME)/share/lib/bed_files/MSK-IMPACT-v3_cnvkit_offtarget.bed

define cnvkit-tumor-cnn
cnvkit/cnn/tumor/$1.targetcoverage.cnn : bam/$1.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"set -o pipefail && \
					 cnvkit.py coverage -p 4 -q 0 $$(<) $$(ONTARGET_FILE) -o cnvkit/cnn/tumor/$1.targetcoverage.cnn")

cnvkit/cnn/tumor/$1.antitargetcoverage.cnn : bam/$1.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"set -o pipefail && \
					 cnvkit.py coverage -p 4 -q 0 $$(<) $$(OFFTARGET_FILE) -o cnvkit/cnn/tumor/$1.antitargetcoverage.cnn")
endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvkit-tumor-cnn,$(sample))))
		
define cnvkit-normal-cnn
cnvkit/cnn/normal/$1.targetcoverage.cnn : bam/$1.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"set -o pipefail && \
					 cnvkit.py coverage -p 4 -q 0 $$(<) $$(ONTARGET_FILE) -o cnvkit/cnn/normal/$1.targetcoverage.cnn")

cnvkit/cnn/normal/$1.antitargetcoverage.cnn : bam/$1.bam
	$$(call RUN,-c -n 4 -s 6G -m 8G,"set -o pipefail && \
					 cnvkit.py coverage -p 4 -q 0 $$(<) $$(OFFTARGET_FILE) -o cnvkit/cnn/normal/$1.antitargetcoverage.cnn")
endef
 $(foreach sample,$(NORMAL_SAMPLES),\
		$(eval $(call cnvkit-normal-cnn,$(sample))))

cnvkit/reference/combined_reference.cnr : $(wildcard cnvkit/cnn/normal/$(NORMAL_SAMPLES).targetcoverage.cnn) $(wildcard cnvkit/cnn/normal/$(NORMAL_SAMPLES).antitargetcoverage.cnn)
	$(call RUN,-n 1 -s 24G -m 32G,"set -o pipefail && \
								   sleep 30 && \
								   cnvkit.py reference cnvkit/cnn/normal/*.cnn -f $(REF_FASTA) --no-edge -o cnvkit/reference/combined_reference.cnr")

define cnvkit-cnr
cnvkit/cnr/%.cnr : cnvkit/cnn/tumor/%.targetcoverage.cnn cnvkit/cnn/tumor/%.antitargetcoverage.cnn cnvkit/reference/combined_reference.cnr
	$$(call RUN,-c -s 6G -m 8G,"set -o pipefail && \
				    cnvkit.py fix $$(<) $$(<<) cnvkit/reference/combined_reference.cnr -o cnvkit/cnr/$$(*).cnr")
	
endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvkit-cnr,$(sample))))

define cnvkit-plot
cnvkit/log2/%.ontarget.pdf : cnvkit/cnr/%.cnr
	$$(call RUN,-c -v $(ASCAT_ENV) -s 4G -m 6G,"set -o pipefail && \
						    $(RSCRIPT) $(SCRIPTS_DIR)/copy_number/cnvkit.R --type '1' --sample_name $$(*)")
	
cnvkit/log2/%.offtarget.pdf : cnvkit/cnr/%.cnr
	$$(call RUN,-c -v $(ASCAT_ENV) -s 4G -m 6G,"set -o pipefail && \
						    $(RSCRIPT) $(SCRIPTS_DIR)/copy_number/cnvkit.R --type '2' --sample_name $$(*)")
	
endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvkit-plot,$(sample))))
		
define cnvkit-totalcopy
cnvkit/segmented/%.pdf cnvkit/totalcopy/%.RData : cnvkit/cnr/%.cnr
	$$(call RUN,-c -v $(ASCAT_ENV) -s 6G -m 12G,"set -o pipefail && \
						     mkdir -p cnvkit/segmented && \
						     mkdir -p cnvkit/totalcopy && \
						     $(RSCRIPT) $(SCRIPTS_DIR)/copy_number/cnvkit.R --type '3' --sample_name $$(*)")
												 
cnvkit/called/%.RData : cnvkit/totalcopy/%.RData
	$$(call RUN,-c -v $(ASCAT_ENV) -s 6G -m 12G,"set -o pipefail && \
						     mkdir -p cnvkit/called && \
						     $(RSCRIPT) $(SCRIPTS_DIR)/copy_number/cnvkit.R --type '4' --sample_name $$(*)")

endef
 $(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call cnvkit-totalcopy,$(sample))))
		

cnvkit/summary/bygene.txt : $(foreach sample,$(TUMOR_SAMPLES),cnvkit/called/$(sample).RData)
	$(call RUN,-c -s 24G -m 48G,"set -o pipefail && \
				     mkdir -p cnvkit/summary && \
				     $(RSCRIPT) $(SCRIPTS_DIR)/copy_number/cnvkit.R --type '5' --sample_name '$(TUMOR_SAMPLES)'")
		
		
..DUMMY := $(shell mkdir -p version; \
	     python $(CNVKIT_ENV)/bin/cnvkit.py version &> version/cnvkit.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: cnv_kit
