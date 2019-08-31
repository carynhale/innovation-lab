include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvaccess_coverage.$(NOW)
PHONY += cnvaccess cnvaccess/cnn


cnvaccess_coverage : $(foreach sample,$(SAMPLES),cnvaccess/cov/$(sample).probe-A.txt)

R_COVERAGE ?= modules/test/copy_number/cnvaccesscoverage.R

define cnvaccess-coverage
cnvaccess/cov/%.probe-A.txt : bam/%-standard.bam
	$$(call RUN,-c -n 1 -s 12G -m 26G,"$(RSCRIPT) $(R_COVERAGE) --sample_name $1 --probe 'A'")
	
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call cnvaccess-coverage,$(sample))))

.PHONY: $(PHONY)
