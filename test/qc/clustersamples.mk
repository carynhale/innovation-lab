include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR = log/cluster_samples.$(NOW)
PHONE += marianas metrics/summary metrics/report

cluster_samples : $(foreach sample,$(SAMPLES),marianas/$(sample)/$(sample)-snps.vcf)

DBSNP_SUBSET ?= $(HOME)/share/reference/dbsnp_tseq_intersect.bed
CLUSTER_VCF ?= $(RSCRIPT) modules/test/qc/clustersample.R

marianas/%/%-snps.vcf : bam/%-duplex.bam
	$(call RUN,-n 4 -s 2.5G -m 3G,"$(call GATK_MEM,8G) -T UnifiedGenotyper -nt 4 -R $(REF_FASTA) --dbsnp $(DBSNP) $(foreach bam,$(filter %.bam,$^),-I $(bam) ) -L $(DBSNP_SUBSET) -o $@ --output_mode EMIT_ALL_SITES")

include modules/vcf_tools/vcftools.mk

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all
