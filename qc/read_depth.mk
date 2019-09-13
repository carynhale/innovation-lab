include modules/Makefile.inc
include modules/config/gatk.inc

LOGDIR = log/read_depth.$(NOW)

EXOME ?= false

ifeq ($(EXOME),true)
READ_DEPTH_ARGS += -L $(EXOME_BED)
endif

.PHONY: all

all: $(foreach sample,$(SAMPLES),gatk/read_depth/$(sample).read_depth)

gatk/read_depth/%.read_depth : %.bam
	$(call RUN,-s 8G -m 12G,"$(call GATK_MEM,7G) -T DepthOfCoverage -R $(REF_FASTA) $(READ_DEPTH_ARGS) -o $@ -I $<")

