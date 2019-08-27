include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR = log/cluster_samples.$(NOW)
PHONE += marianas metrics/summary metrics/report

cluster_samples : $(foreach sample,$(SAMPLES),marianas/$(sample)/$(sample)-snps.vcf) \
				  $(foreach sample,$(SAMPLES),marianas/$(sample)/$(sample)-snps-filtered.vcf)

DBSNP_SUBSET ?= $(HOME)/share/reference/dbsnp_tseq_intersect.bed
CLUSTER_VCF ?= $(RSCRIPT) modules/test/qc/clustersample.R

define genotype-snps
marianas/$1/$1-snps.vcf : marianas/$1/$1.standard.bam
	$$(call RUN,-n 4 -s 2.5G -m 3G,"$(call GATK_MEM,8G) -T UnifiedGenotyper -nt 4 -R $(REF_FASTA) --dbsnp $(DBSNP) \
									-I $$(<) \
									-L $(DBSNP_SUBSET) \
									-o $$(@) \
									--output_mode EMIT_ALL_SITES")
									
marianas/$1/$1-snps-filtered.vcf : marianas/$1/$1-snps.vcf
	$$(call RUN,-n 1 -s 6G -m 12G,"grep '^#' $$(<) > $$(@) && \
								   grep -e '0/1' -e '1/1' $$(<) >> $$(@)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call genotype-snps,$(sample))))


include modules/vcf_tools/vcftools.mk

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all
