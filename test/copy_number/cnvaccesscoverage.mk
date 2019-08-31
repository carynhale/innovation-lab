include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvaccess_coverage.$(NOW)
PHONY += cnvaccess cnvaccess/cnn


cnvaccess_coverage : $(foreach sample,$(SAMPLES),cnvaccess/cov/$(sample).probe-A.txt)

R_COVERAGE ?= modules/test/copy_number/cnvaccesscoverage.R
EXOME_DEPTH_ENV ?= $(HOME)/usr/anaconda-envs/exomedepth-1.1.12

define cnvaccess-coverage
cnvaccess/cov/$1.probe-A.txt : bam/$1-standard.bam
	$$(call RUN,-c -n 1 -s 12G -m 26G -v $(EXOME_DEPTH_ENV),"$(RSCRIPT) $(R_COVERAGE) --sample_name $1 --probe 'A'")
	
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call cnvaccess-coverage,$(sample))))

.PHONY: $(PHONY)
