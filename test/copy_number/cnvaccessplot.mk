include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/cnvaccess_plot.$(NOW)
PHONY += cnvaccess cnvaccess/report cnvaccess/report/log2 cnvaccess/report/ASCAT

cnvaccess_plot : $(foreach sample,$(SAMPLES),cnvaccess/report/log2/$(sample).pdf) \
				 $(foreach sample,$(SAMPLES),cnvaccess/report/segmented/$(sample).RData) \
				 $(foreach sample,$(SAMPLES),cnvaccess/report/ASCAT/$(sample).pdf) \
				 $(foreach sample,$(SAMPLES),cnvaccess/report/log2/$(sample)-2.pdf)
				 

R_COVERAGE ?= modules/test/copy_number/cnvaccesscoverage.R
R_FIX ?= modules/test/copy_number/cnvaccessfix.R
R_PLOT ?= modules/test/copy_number/cnvaccessplot.R
EXOME_DEPTH_ENV ?= $(HOME)/share/usr/anaconda-envs/exomedepth-1.1.12

define cnvaccess-plot
cnvaccess/report/log2/$1.pdf : cnvaccess/log2/$1.txt
	$$(call RUN,-c -n 1 -s 4G -m 6G -v $(ASCAT_ENV),"$(RSCRIPT) $(R_PLOT) --type 1 --sample_name $1")
	
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call cnvaccess-plot,$(sample))))
		
define cnvaccess-segment
cnvaccess/report/segmented/$1.RData : cnvaccess/log2/$1.txt
	$$(call RUN,-c -n 1 -s 4G -m 6G -v $(ASCAT_ENV),"$(RSCRIPT) $(R_PLOT) --type 2 --sample_name $1")
	
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call cnvaccess-segment,$(sample))))
		
define cnvaccess-ascat
cnvaccess/report/ASCAT/$1.pdf : cnvaccess/report/segmented/$1.RData
	$$(call RUN,-c -n 1 -s 4G -m 6G -v $(ASCAT_ENV),"$(RSCRIPT) $(R_PLOT) --type 3 --sample_name $1")
	
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call cnvaccess-ascat,$(sample))))


define cnvaccess-plot-cnvkit
cnvaccess/report/log2/$1-2.pdf : cnvaccess/cnr/$1.txt
	$$(call RUN,-c -n 1 -s 4G -m 6G -v $(ASCAT_ENV),"$(RSCRIPT) $(R_PLOT) --type 4 --sample_name $1")
	
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call cnvaccess-plot-cnvkit,$(sample))))


.PHONY: $(PHONY)
