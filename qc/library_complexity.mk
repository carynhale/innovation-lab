include innovation-lab/Makefile.inc

LOGDIR ?= log/library_complexity.$(NOW)

picard_metrics : $(foreach sample,$(SAMPLES),metrics/$(sample)_library_complexity.txt) \
		 summary/complexity_summary.txt \
		 summary/complexity_metrics.txt

define library-complexity
metrics/$1_library_complexity.txt : bam/$1.bam
	$$(call RUN,-c -s 12G -m 24G,"set -o pipefail && \
				     $$(CALC_LIBRARY_COMPLEXITY) \
				     INPUT=$$(<) \
				     OUTPUT=$$(@) \
				     MIN_IDENTICAL_BASES=5 \
				     MAX_DIFF_RATE=0.03")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call library-complexity,$(sample))))
		
summary/complexity_summary.txt : $(foreach sample,$(SAMPLES),metrics/$(sample)_library_complexity.txt)
	$(call RUN, -c -n 1 -s 12G -m 24G,"set -o pipefail && \
					   $(RSCRIPT) $(SCRIPTS_DIR)/qc/complexity_metrics.R --option 1 --sample_names '$(SAMPLES)'")
					 
summary/complexity_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample)_library_complexity.txt)
	$(call RUN, -c -n 1 -s 12G -m 24G,"set -o pipefail && \
					   $(RSCRIPT) $(SCRIPTS_DIR)/qc/complexity_metrics.R --option 2 --sample_names '$(SAMPLES)'")


..DUMMY := $(shell mkdir -p version; \
	     echo "picard" >> version/library_complexity.txt; \
	     $(PICARD) EstimateLibraryComplexity --version &>> version/library_complexity.txt; \
             R --version >> version/library_complexity.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: picard_metrics
