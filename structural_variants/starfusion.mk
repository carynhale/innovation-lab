include innovation-lab/Makefile.inc

LOGDIR ?= log/star_fusion.$(NOW)

CTAT_LIB ?= $(HOME)/share/lib/ref_files/CTAT_GRCh37

star_fusion : $(foreach sample,$(SAMPLES),starfusion/$(sample)/$(sample).1.fastq) \
			  $(foreach sample,$(SAMPLES),starfusion/$(sample)/$(sample).2.fastq) \
			  $(foreach sample,$(SAMPLES),starfusion/$(sample)/taskcomplete)

define merged-fastq
starfusion/$1/$1.1.fastq : $$(foreach split,$2,$$(word 1, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 2G -m 4G,"set -o pipefail && \
									 mkdir -p starfusion && \
									 mkdir -p starfusion/$1 && \
									 zcat $$(^) > $$(@)")
starfusion/$1/$1.2.fastq : $$(foreach split,$2,$$(word 2, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 2G -m 4G,"set -o pipefail && \
									 mkdir -p starfusion && \
									 mkdir -p starfusion/$1 && \
									 zcat $$(^) > $$(@)")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call merged-fastq,$(sample),$(split.$(sample)))))

define star-fusion
starfusion/$1/taskcomplete : starfusion/$1/$1.1.fastq starfusion/$1/$1.2.fastq
	$(call RUN,-c -n 10 -s 2G -m 3G -v $(STAR_FUSION_ENV),"$(STAR_FUSION) \
														  --left_fq starfusion/$1/$1.1.fastq \
														  --right_fq starfusion/$1/$1.2.fastq \
														  --CPU 10 \
														  --output_dir starfusion/$1 \
														  --genome_lib_dir $$(CTAT_LIB) && \
														  touch echo $1 > starfusion/$1/taskcomplete")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call star-fusion,$(sample))))

..DUMMY := $(shell mkdir -p version; \
			 ~/share/usr/env/starfusion-1.6.0/bin/STAR-Fusion --help > version/starfusion.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: star_fusion
