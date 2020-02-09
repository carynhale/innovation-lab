include modules/Makefile.inc
include modules/bam_tools/process_bam.mk

LOGDIR = log/integrate_rnaseq.$(NOW)
.PHONY: integrate_rnaseq

INTEGRATE_MINW ?= 2.0
INTEGRATE_LARGENUM ?= 4
INTEGRATE_OPTS = -minW $(INTEGRATE_MINW) -largeNum $(INTEGRATE_LARGENUM)
INTEGRATE_ONCOFUSE = $(RSCRIPT) $(SCRIPTS_DIR)/structural_variants/integrateOncofuse.R
INTEGRATE_ONCOFUSE_OPTS = --oncofuseJar $(ONCOFUSE_JAR) --oncofuseTissueType $(ONCOFUSE_TISSUE_TYPE) --java $(JAVA_BIN) \
						  --mysqlHost $(EMBL_MYSQLDB_HOST) --mysqlPort $(EMBL_MYSQLDB_PORT) --mysqlUser $(EMBL_MYSQLDB_USER) \
						  $(if $(EMBL_MYSQLDB_PW),--mysqlPassword $(EMBL_MYSQLDB_PW)) --mysqlDb $(EMBL_MYSQLDB_DB)
ONCOFUSE_TISSUE_TYPE ?= EPI
INTEGRATE_TO_USV = python $(SCRIPTS_DIR)/structural_variants/integrate2usv.py

integrate_rnaseq: integrate_rnaseq/summary.tsv

define init-integrate
integrate_rnaseq/reads/%.reads.txt integrate_rnaseq/sum/%.sum.tsv integrate_rnaseq/exons/%.exons.tsv integrate_rnaseq/breakpoints/%.breakpoints.tsv : bam/%.bam bam/%.bam.bai
	$$(call RUN,-s 8G -m 40G,"mkdir -p integrate_rnaseq/reads && \
							  mkdir -p integrate_rnaseq/sum && \
							  mkdir -p integrate_rnaseq/exons && \
							  mkdir -p integrate_rnaseq/breakpoints && \
							  $$(INTEGRATE) fusion $$(INTEGRATE_OPTS) \
							  -reads integrate_rnaseq/reads/$$(*).reads.txt \
							  -sum integrate_rnaseq/sum/$$(*).sum.tsv \
							  -ex integrate_rnaseq/exons/$$(*).exons.tsv \
							  -bk integrate_rnaseq/breakpoints/$$(*).breakpoints.tsv \
							  $$(REF_FASTA) $$(INTEGRATE_ANN) $$(INTEGRATE_BWTS) $$(<) $$(<)")

endef
$(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call init-integrate,$(sample))))

define init-oncofuse
integrate_rnaseq/oncofuse/%.oncofuse.txt : integrate_rnaseq/sum/%.sum.tsv integrate_rnaseq/exons/%.exons.tsv integrate_rnaseq/breakpoints/%.breakpoints.tsv
	$$(call RUN,-s 7G -m 10G,"$$(INTEGRATE_ONCOFUSE) $$(INTEGRATE_ONCOFUSE_OPTS) \
							  --ref $$(REF) \
					  		  --sumFile $$(<) \
							  --exonsFile $$(<<) \
							  --breakpointsFile $$(<<<) \
							  --outPrefix $$(@D)/$$(*)")
		
endef
$(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call init-oncofuse,$(sample))))

integrate_rnaseq/summary.tsv : $(wildcard $(foreach sample,$(TUMOR_SAMPLES),integrate_rnaseq/oncofuse/$(sample).oncofuse.txt))
	$(call RUN,-c -n 1 -s 6G -m 8G,"set -o pipefail && \
				     			    $(RSCRIPT) $(SCRIPTS_DIR)/structural_variants/integraternaseq.R --samples '$(TUMOR_SAMPLES)'")

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
