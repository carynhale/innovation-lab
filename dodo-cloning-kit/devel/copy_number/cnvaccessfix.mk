include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvaccess_fix.$(NOW)
PHONY += cnvaccess cnvaccess/log2

cnvaccess_fix : $(foreach sample,$(SAMPLES),cnvaccess/log2/$(sample)-ontarget.txt) \
				$(foreach sample,$(SAMPLES),cnvaccess/log2/$(sample)-offtarget.txt)

R_COVERAGE ?= modules/test/copy_number/cnvaccesscoverage.R
R_FIX ?= modules/test/copy_number/cnvaccessfix.R
R_PLOT ?= modules/test/copy_number/cnvaccessplot.R
EXOME_DEPTH_ENV ?= $(HOME)/share/usr/anaconda-envs/exomedepth-1.1.12
OFFTARGET_FILE ?= $(HOME)/share/reference/target_panels/MSK-ACCESS-v1_0-probe-AB.offtarget.bed
REFERENCE_FILE ?= $(HOME)/share/reference/cnvkit_reference/MSK-ACCESS-v1_0-noprobe.cnr

define cnvaccess-fix
cnvaccess/log2/$1-ontarget.txt : cnvaccess/cov/$1.probe-A.txt cnvaccess/cov/$1.probe-B.txt
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(RSCRIPT) $(R_FIX) --sample_name $1")
	
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call cnvaccess-fix,$(sample))))
		
define cnvaccess-cnvkit-fix
cnvaccess/log2/$1-offtarget.txt : cnvaccess/cnn/$1.targetcoverage.cnn cnvaccess/cnn/$1.antitargetcoverage.cnn
	$$(call RUN,-c -s 6G -m 8G,"set -o pipefail && \
								cnvkit.py fix $$(<) $$(<<) $(REFERENCE_FILE) -o cnvaccess/log2/$1-offtarget.txt")
	
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call cnvaccess-cnvkit-fix,$(sample))))

.PHONY: $(PHONY)
