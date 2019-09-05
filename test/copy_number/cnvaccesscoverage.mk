include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvaccess_coverage.$(NOW)
PHONY += cnvaccess cnvaccess/cov cnvaccess/cnn

cnvaccess_coverage : $(foreach sample,$(SAMPLES),cnvaccess/cov/$(sample).probe-A.txt) \
					 $(foreach sample,$(SAMPLES),cnvaccess/cov/$(sample).probe-B.txt) \
					 $(foreach sample,$(SAMPLES),cnvaccess/cnn/$(sample).antitargetcoverage.cnn) \
					 $(foreach sample,$(SAMPLES),cnvaccess/cnn/$(sample).targetcoverage.cnn)
					 

R_COVERAGE ?= modules/test/copy_number/cnvaccesscoverage.R
R_FIX ?= modules/test/copy_number/cnvaccessfix.R
R_PLOT ?= modules/test/copy_number/cnvaccessplot.R
EXOME_DEPTH_ENV ?= $(HOME)/share/usr/anaconda-envs/exomedepth-1.1.12
OFFTARGET_FILE ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-AB.offtarget.bed
REFERENCE_FILE ?= $(HOME)/share/reference/cnvkit_reference/MSK-ACCESS-v1_0-noprobe.cnr

define cnvaccess-coverage
cnvaccess/cov/$1.probe-A.txt : bam/$1-standard.bam
	$$(call RUN,-c -n 1 -s 6G -m 12G -v $(EXOME_DEPTH_ENV),"set -o pipefail && \
															$(RSCRIPT) $(R_COVERAGE) --sample_name $1 --probe 'A'")
	
cnvaccess/cov/$1.probe-B.txt : bam/$1-standard.bam
	$$(call RUN,-c -n 1 -s 6G -m 12G -v $(EXOME_DEPTH_ENV),"set -o pipefail && \
															$(RSCRIPT) $(R_COVERAGE) --sample_name $1 --probe 'B'")
	
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call cnvaccess-coverage,$(sample))))
		
define cnvaccess-cnvkit-cov
cnvaccess/cnn/$1.antitargetcoverage.cnn : bam/$1-standard.bam
		$$(call RUN,-c -n 4 -s 6G -m 8G,"set -o pipefail && \
										 cnvkit.py coverage -p 4 -q 0 $$(<) $$(OFFTARGET_FILE) -o cnvaccess/cnn/$1.antitargetcoverage.cnn")
	
cnvaccess/cnn/$1.targetcoverage.cnn : bam/$1-standard.bam
		$$(call RUN,-n 1 -s 1G -m 2G,"set -o pipefail && \
									  touch cnvaccess/cnn/$1.targetcoverage.cnn")
	
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call cnvaccess-cnvkit-cov,$(sample))))

.PHONY: $(PHONY)
