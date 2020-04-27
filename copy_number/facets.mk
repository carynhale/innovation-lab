include innovation-lab/Makefile.inc

LOGDIR ?= log/facets.$(NOW)

facets : $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).txt) \
		 $(foreach pair,$(SAMPLE_PAIRS),facets/plots/log2/$(pair).pdf) \
		 facets/summary/bygene.txt \
		 facets/summary/bygene.pdf \
		 facets/summary/summary.tsv

RUN_FACETS = $(RSCRIPT) $(SCRIPTS_DIR)/copy_number/facets.R
MERGE_TN = python $(SCRIPTS_DIR)/copy_number/facets_merge_tn.py
FACETS_PRE_CVAL ?= 50
FACETS_CVAL1 ?= 150
FACETS_CVAL2 ?= 50
FACETS_MIN_NHET ?= 25
FACETS_SNP_NBHD ?= 250
FACETS_NDEPTH_MAX ?= 1000
FACETS_HET_THRESHOLD ?= 0.25
FACETS_GATK_VARIANTS ?= false
FACETS_OPTS = --genome $(REF) \
			  $(if $(facets_cval1.$1),--cval1 $(facets_diplogr.$1),--cval1 $(FACETS_CVAL1)) \
			  $(if $(facets_cval2.$1),--cval2 $(facets_diplogr.$1),--cval2 $(FACETS_CVAL1)) \
			  --het_threshold $(FACETS_HET_THRESHOLD) \
			  --min_nhet $(FACETS_MIN_NHET) \
			  --snp_nbhd $(FACETS_SNP_NBHD) \
			  $(if $(facets_pre_cval.$1),--pre_cval $(facets_pre_cval.$1),--pre_cval $(FACETS_PRE_CVAL)) \
			  --ndepth_max $(FACETS_NDEPTH_MAX) \
			  --use_emcncf2 \
			  $(if $(facets_diplogr.$1),--diplogr $(facets_diplogr.$1)) \
			  $(if $(facets_purity.$1),--purity $(facets_purity.$1))
SNP_PILEUP = snp-pileup
SNP_PILEUP_OPTS = -A --min-map-quality=15 --min-base-quality=15 --gzip --max-depth=15000
FACETS_DBSNP = $(if $(TARGETS_FILE),facets/vcf/targets_dbsnp.vcf,$(DBSNP))
CONVERT_BASECOUNT ?= false
FACETS_UNION_GATK_DBSNP ?= false
ifeq ($(FACETS_UNION_GATK_DBSNP),true)
FACETS_SNP_VCF = facets/vcf/dbsnp_het_gatk.snps.vcf
else
FACETS_SNP_VCF = $(FACETS_DBSNP)
endif

FACETS_GENE_CN_OPTS = $(if $(GENES_FILE),--genesFile $(GENES_FILE)) \
					  --mysqlHost $(EMBL_MYSQLDB_HOST) --mysqlPort $(EMBL_MYSQLDB_PORT) \
					  --mysqlUser $(EMBL_MYSQLDB_USER) $(if $(EMBL_MYSQLDB_PW),--mysqlPassword $(EMBL_MYSQLDB_PW)) \
					  --mysqlDb $(EMBL_MYSQLDB_DB)
FACETS_PLOT_GENE_CN_OPTS = --sampleColumnPostFix '_LRR_threshold'



facets/vcf/dbsnp_het_gatk.snps.vcf : $(FACETS_DBSNP) $(foreach sample,$(SAMPLES),gatk/vcf/$(sample).variants.snps.het.pass.vcf)
	$(call RUN,-c -s 4G -m 6G,"set -o pipefail && \
							   $(call GATK_MEM,3G) $(if $(TARGETS_FILE),-L $(TARGETS_FILE)) -T CombineVariants --minimalVCF $(foreach i,$^, --variant $i) -R $(REF_FASTA) -o $@")

%.het.vcf : %.vcf
	$(call RUN,-c -s 9G -m 12G,"set -o pipefail && \
								$(call GATK_MEM,8G) -V $< -T VariantFiltration -R $(REF_FASTA) --genotypeFilterName 'hom' --genotypeFilterExpression 'isHet == 0' -o $@")

facets/vcf/targets_dbsnp.vcf : $(TARGETS_FILE)
	$(INIT) $(BEDTOOLS) intersect -header -u -a $(DBSNP) -b $< > $@

ifeq ($(CONVERT_BASECOUNT),true)
CONVERT_BC_TO_SNP_PILEUP = python $(SCRIPTS_DIR)/copy_number/convert_basecount_to_snp_pileup.py
facets/pileup/%.gz : facets/base_count/%.bc.gz
	$(call RUN,-s 12G -m 14G,"set -o pipefail && \
							  $(CONVERT_BC_TO_SNP_PILEUP) $< | gzip -c > $@")
else
define snp-pileup-tumor-normal
facets/pileup/$1_$2.gz : bam/$1.bam bam/$2.bam $$(FACETS_SNP_VCF)
	$$(call RUN,-c -s 8G -m 20G,"set -o pipefail && \
								 rm -f $$@ && $$(SNP_PILEUP) $$(SNP_PILEUP_OPTS) $$(<<<) $$@ $$(<<) $$(<)")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call snp-pileup-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
endif

facets/cncf/%.txt facets/cncf/%.Rdata : facets/pileup/%.gz
	$(call RUN,-c -v $(FACETS_ENV) -s 8G -m 60G,"set -o pipefail && \
												 $(RUN_FACETS) --option 1 $(call FACETS_OPTS,$*) --out_prefix $(@D)/$* $<")

facets/plots/log2/%.pdf : facets/cncf/%.Rdata
	$(call RUN,-v $(FACETS_ENV) -s 4G -m 6G,"set -o pipefail && \
											 $(RUN_FACETS) --option 2 --centromereFile $(CENTROMERE_TABLE) --sample_name $(*)")

facets/summary/bygene.txt : $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).Rdata)
	$(call RUN,-c -s 8G -m 30G,"set -o pipefail && \
								$(RUN_FACETS) --option 3 $(FACETS_GENE_CN_OPTS) --outFile $@ $^")

facets/summary/bygene.pdf : facets/summary/bygene.txt
	$(call RUN,-s 8G -m 10G,"set -o pipefail && \
							 $(RUN_FACETS) --option 4 $(FACETS_PLOT_GENE_CN_OPTS) $< $@")
							 
facets/summary/summary.tsv : $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).Rdata)
	$(call RUN,-c -s 8G -m 12G,"set -o pipefail && \
								$(RUN_FACETS) --option 5 --outFile $@ $^")

include modules/variant_callers/gatk.mk
include modules/bam_tools/process_bam.mk

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : $(PHONY)
