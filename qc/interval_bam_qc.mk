include modules/Makefile.inc

INTERVAL_FILE ?= intervals.bed
VPATH ?= bam

TEQC = $(SCRIPTS_DIR)/qc/teqc.R
INTERVAL_BAM_QC = $(SCRIPTS_DIR)/qc/interval_bam_qc.R
VARIANT_EVAL_REPORT = $(SCRIPTS_DIR)/qc/variant_eval_gatk_report.R

LOGDIR ?= log/interval_qc.$(NOW)

all : interval_qc/coverage/index.html
rdata : $(foreach sample,$(SAMPLES),amplicon_qc/rdata/$(sample).Rdata)

interval_qc/rdata/%.Rdata : %.bam
	$(call RUN,-s 8G -m 18G,"$(RSCRIPT) $(TEQC) --outFile $@ --ref $(REF) $< $(INTERVAL_FILE)")

interval_qc/variantEval.grp : $(foreach sample,$(SAMPLES),gatk/vcf/$(sample).variants.vcf)
	$(call RUN,-n 4 -s 1G -m 1.5G,"$(call GATK_MEM,4G) -T VariantEval -nt 4 -o $@ -R $(REF_FASTA) --stratIntervals $(INTERVAL_FILE) --dbsnp $(DBSNP) $(foreach vcf,$^, --eval:$(call strip-suffix,$(notdir $(vcf))) $(vcf)) -ST IntervalStratification -ST Filter")

interval_qc/coverage/index.html : $(foreach sample,$(SAMPLES),interval_qc/rdata/$(sample).Rdata)
	$(call RUN,-s 2G -m 4G,"$(RSCRIPT) $(INTERVAL_BAM_QC) --outDir $(@D) $^")

interval_qc/variant_eval/index.html : interval_qc/variantEval.grp
	$(call RUN,-s 2G -m 4G,"$(RSCRIPT) $(VARIANT_EVAL_REPORT) --outDir $(@D) $< &> $(LOG)")

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: rdata all
