include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvaccess_fix.$(NOW)
PHONY += cnvaccess cnvaccess/log2

cnvaccess_fix : $(foreach sample,$(SAMPLES),cnvaccess/log2/$(sample).txt)

R_COVERAGE ?= modules/test/copy_number/cnvaccesscoverage.R
R_FIX ?= modules/test/copy_number/cnvaccessfix.R
R_PLOT ?= modules/test/copy_number/cnvaccessplot.R
EXOME_DEPTH_ENV ?= $(HOME)/share/usr/anaconda-envs/exomedepth-1.1.12

define cnvaccess-fix
cnvaccess/log2/$1.txt : cnvaccess/cov/$1.probe-A.txt cnvaccess/cov/$1.probe-B.txt
	$$(call RUN,-c -n 1 -s 8G -m 16G,"$(RSCRIPT) $(R_FIX) --sample_name $1")
	
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call cnvaccess-fix,$(sample))))

.PHONY: $(PHONY)
