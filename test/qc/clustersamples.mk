include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR = log/cluster_samples.$(NOW)
PHONE += marianas metrics/summary metrics/report

cluster_samples : $(foreach sample,$(SAMPLES),marianas/$(sample)/$(sample)-snps.vcf) \
				  metrics/summary/snps.vcf \
				  metrics/summary/snps-filtered.vcf \
				  metrics/summary/snps-filtered.tsv
#				  metrics/report/snp_clustering.pdf

DBSNP_SUBSET = $(HOME)/share/reference/dbsnp_tseq_intersect.bed
CLUSTER_VCF = $(RSCRIPT) modules/test/qc/clustersamples.R

define genotype-snps
marianas/$1/$1-snps.vcf : marianas/$1/$1.standard.bam
	$$(call RUN,-n 4 -s 2.5G -m 3G,"$(call GATK_MEM,8G) -T UnifiedGenotyper -nt 4 -R $(REF_FASTA) --dbsnp $(DBSNP) \
									-I $$(<) \
									-L $(DBSNP_SUBSET) \
									-o $$(@) \
									--output_mode EMIT_ALL_SITES")
									
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call genotype-snps,$(sample))))
		
metrics/summary/snps.vcf : $(foreach sample,$(SAMPLES),marianas/$(sample)/$(sample)-snps.vcf)
	$(call RUN,-s 16G -m 20G,"$(call GATK_MEM,14G) -T CombineVariants $(foreach vcf,$^,--variant $(vcf) ) \
							  -o $@ --genotypemergeoption UNSORTED -R $(REF_FASTA)")
							  
metrics/summary/snps-filtered.vcf : metrics/summary/snps.vcf
	$(INIT) grep '^#' $< > $@ && grep -e '0/1' -e '1/1' $< >> $@
	
metrics/summary/snps-filtered.tsv : metrics/summary/snps-filtered.vcf
	$(INIT) $(CLUSTER_VCF)

include modules/vcf_tools/vcftools.mk

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all
