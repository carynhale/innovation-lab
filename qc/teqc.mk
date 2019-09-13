include modules/Makefile.inc

LOGDIR ?= teqc/log

teqc : teqc_report/index.html

teqc/%.Rdata : bam/%.bam bam/%.bam.bai
	$(call INIT_MEM,12G,14G) $(RSCRIPT) $(SCRIPTS_DIR)/qc/teqc.R --ref=$(REF) --outFile $@ $< $(TARGETS_FILE) &> $(LOGDIR)/$(@F).log

teqc_report/index.html : $(foreach sample,$(SAMPLES),teqc/$(sample).Rdata)
	$(call INIT_MEM,12G,14G) $(MKDIR) teqc_report; $(RSCRIPT) $(SCRIPTS_DIR)/qc/teqc_report.R --outDir=$(@D) $^

.PHONY: teqc
.DELETE_ON_ERROR:
.SECONDARY:
