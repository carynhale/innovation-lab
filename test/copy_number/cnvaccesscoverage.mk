include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvaccess_coverage.$(NOW)
PHONY += cnvaccess cnvaccess/cnn


cnvaccess_coverage : $(foreach sample,$(SAMPLES),cnvaccess/cov/$(sample).probe-A.txt) \
					 $(foreach sample,$(SAMPLES),cnvaccess/cov/$(sample).probe-B.txt) \
					 $(foreach sample,$(SAMPLES),cnvaccess/cov/$(sample).probe-AB.txt)

R_COVERAGE ?= modules/test/copy_number/cnvaccesscoverage.R
EXOME_DEPTH_ENV ?= $(HOME)/share/usr/anaconda-envs/exomedepth-1.1.12

define cnvaccess-coverage
cnvaccess/cov/$1.probe-A.txt : bam/$1-standard.bam
	$$(call RUN,-c -n 1 -s 6G -m 12G -v $(EXOME_DEPTH_ENV),"$(RSCRIPT) $(R_COVERAGE) --sample_name $1 --probe 'A'")
	
cnvaccess/cov/$1.probe-B.txt : bam/$1-standard.bam
	$$(call RUN,-c -n 1 -s 6G -m 12G -v $(EXOME_DEPTH_ENV),"$(RSCRIPT) $(R_COVERAGE) --sample_name $1 --probe 'B'")
	
cnvaccess/cov/$1.probe-AB.txt : bam/$1-standard.bam
	$$(call RUN,-c -n 1 -s 12G -m 24G -v $(EXOME_DEPTH_ENV),"$(RSCRIPT) $(R_COVERAGE) --sample_name $1 --probe 'NA'")
	
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call cnvaccess-coverage,$(sample))))

.PHONY: $(PHONY)
