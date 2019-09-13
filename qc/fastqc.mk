include modules/Makefile.inc

FASTQC_SUMMARY_PLOT = $(RSCRIPT) $(SCRIPTS_DIR)/qc/fastqc_summary.R

LOGDIR ?= log/fastqc.$(NOW)

.PHONY: fastqc
.SECONDARY: 

fastqc : $(foreach sample,$(SAMPLES),fastqc/$(sample)_fastqc/summary.txt) fastqc/all_summary.txt

fastqc/%_fastqc.zip : bam/%.bam
	$(call RUN,-N $*_fastqc -s 4G -m 12G,"$(FASTQC) -o fastqc $^")

fastqc/%_fastqc/summary.txt : fastqc/%_fastqc.zip
	$(INIT) $(UNZIP) -o -d fastqc $< &> $(LOG) && touch $@

fastqc/all_summary.txt : $(foreach sample,$(SAMPLES),fastqc/$(sample)_fastqc/summary.txt)
	$(INIT) $(FASTQC_SUMMARY_PLOT) --outPrefix fastqc/all_summary $^ &> $(LOG)
