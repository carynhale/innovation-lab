include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR = log/deseq.$(NOW)

DESEQ_RNW = $(SCRIPTS_DIR)/rnaseq/deseq.Rnw
SWEAVE = $(RSCRIPT) $(SCRIPTS_DIR)/runtime/Sweave.R

DESEQ_CONDITION ?= condition
DESEQ_REF_CONDITION ?= ref
PHENO_FILE ?= pheno.txt

deseq_results.txt : sumreads/geneCounts.txt
	mkdir -p graphics; $(SWEAVE) $(DESEQ_RNW) --condition $(DESEQ_CONDITION) --refCondition $(DESEQ_REF_CONDITION) --outFile $@ $< $(PHENO_FILE)

.PHONY : all
.DELETE_ON_ERROR: 
.SECONDARY: 