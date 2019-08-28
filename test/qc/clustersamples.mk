include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR = log/cluster_samples.$(NOW)
PHONE += metrics metrics/summary metrics/report

cluster_samples : $(foreach sample,$(SAMPLES),metrics/standard/$(sample)-snps.vcf) \
				  metrics/summary/snps_combined-standard.vcf \
				  metrics/summary/snps_filtered-standard.vcf \
				  metrics/summary/snps_filtered-standard.tsv
#				  metrics/report/snps_clustering-standard.pdf

DBSNP_SUBSET = $(HOME)/share/reference/dbsnp_tseq_intersect.bed
CLUSTER_VCF = $(RSCRIPT) modules/test/qc/clustersamples.R

define genotype-snps-standard
metrics/standard/$1-snps.vcf : bam/$1-standard.bam
	$$(call RUN,-n 4 -s 2.5G -m 3G,"$(call GATK_MEM,8G) -T UnifiedGenotyper -nt 4 -R $(REF_FASTA) --dbsnp $(DBSNP) \
									-I $$(<) \
									-L $(DBSNP_SUBSET) \
									-o $$(@) \
									--output_mode EMIT_ALL_SITES")
									
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call genotype-snps-standard,$(sample))))
		
define genotype-snps-unfiltered
metrics/unfiltered/$1-snps.vcf : bam/$1-unfiltered.bam
	$$(call RUN,-n 4 -s 2.5G -m 3G,"$(call GATK_MEM,8G) -T UnifiedGenotyper -nt 4 -R $(REF_FASTA) --dbsnp $(DBSNP) \
									-I $$(<) \
									-L $(DBSNP_SUBSET) \
									-o $$(@) \
									--output_mode EMIT_ALL_SITES \
									--min_base_quality_score 10")
									
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call genotype-snps-unfiltered,$(sample))))
		
define genotype-snps-simplex
metrics/simplex/$1-snps.vcf : bam/$1-simplex.bam
	$$(call RUN,-n 4 -s 2.5G -m 3G,"$(call GATK_MEM,8G) -T UnifiedGenotyper -nt 4 -R $(REF_FASTA) --dbsnp $(DBSNP) \
									-I $$(<) \
									-L $(DBSNP_SUBSET) \
									-o $$(@) \
									--output_mode EMIT_ALL_SITES \
									--min_base_quality_score 10")
									
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call genotype-snps-simplex,$(sample))))

define genotype-snps-duplex
metrics/duplex/$1-snps.vcf : bam/$1-duplex.bam
	$$(call RUN,-n 4 -s 2.5G -m 3G,"$(call GATK_MEM,8G) -T UnifiedGenotyper -nt 4 -R $(REF_FASTA) --dbsnp $(DBSNP) \
									-I $$(<) \
									-L $(DBSNP_SUBSET) \
									-o $$(@) \
									--output_mode EMIT_ALL_SITES \
									--min_base_quality_score 10")
									
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call genotype-snps-duplex,$(sample))))


metrics/summary/snps_combined-standard.vcf : $(foreach sample,$(SAMPLES),metrics/standard/$(sample)-snps.vcf)
	$(call RUN,-s 16G -m 20G,"$(call GATK_MEM,14G) -T CombineVariants $(foreach vcf,$^,--variant $(vcf) ) \
							  -o $@ --genotypemergeoption UNSORTED -R $(REF_FASTA)")
							  
metrics/summary/snps_filtered-standard.vcf : metrics/summary/snps_combined-standard.vcf
	$(INIT) grep '^#' $< > $@ && grep -e '0/1' -e '1/1' $< >> $@
	
metrics/summary/snps_filtered-standard.tsv : metrics/summary/snps_filtered-standard.vcf
	$(INIT) $(CLUSTER_VCF) --library 'STANDARD'
	
metrics/summary/snps_combined-unfiltered.vcf : $(foreach sample,$(SAMPLES),metrics/unfiltered/$(sample)-snps.vcf)
	$(call RUN,-s 16G -m 20G,"$(call GATK_MEM,14G) -T CombineVariants $(foreach vcf,$^,--variant $(vcf) ) \
							  -o $@ --genotypemergeoption UNSORTED -R $(REF_FASTA)")
							  
metrics/summary/snps_filtered-unfiltered.vcf : metrics/summary/snps_combined-unfiltered.vcf
	$(INIT) grep '^#' $< > $@ && grep -e '0/1' -e '1/1' $< >> $@
	
metrics/summary/snps_filtered-unfiltered.tsv : metrics/summary/snps_filtered-unfiltered.vcf
	$(INIT) $(CLUSTER_VCF) --library 'UNFILTERED'
	
metrics/summary/snps_combined-simplex.vcf : $(foreach sample,$(SAMPLES),metrics/simplex/$(sample)-snps.vcf)
	$(call RUN,-s 16G -m 20G,"$(call GATK_MEM,14G) -T CombineVariants $(foreach vcf,$^,--variant $(vcf) ) \
							  -o $@ --genotypemergeoption UNSORTED -R $(REF_FASTA)")
							  
metrics/summary/snps_filtered-simplex.vcf : metrics/summary/snps_combined-simplex.vcf
	$(INIT) grep '^#' $< > $@ && grep -e '0/1' -e '1/1' $< >> $@
	
metrics/summary/snps_filtered-simplex.tsv : metrics/summary/snps_filtered-simplex.vcf
	$(INIT) $(CLUSTER_VCF) --library 'SIMPLEX'

metrics/summary/snps_combined-duplex.vcf : $(foreach sample,$(SAMPLES),metrics/duplex/$(sample)-snps.vcf)
	$(call RUN,-s 16G -m 20G,"$(call GATK_MEM,14G) -T CombineVariants $(foreach vcf,$^,--variant $(vcf) ) \
							  -o $@ --genotypemergeoption UNSORTED -R $(REF_FASTA)")
							  
metrics/summary/snps_filtered-duplex.vcf : metrics/summary/snps_combined-duplex.vcf
	$(INIT) grep '^#' $< > $@ && grep -e '0/1' -e '1/1' $< >> $@
	
metrics/summary/snps_filtered-duplex.tsv : metrics/summary/snps_filtered-duplex.vcf
	$(INIT) $(CLUSTER_VCF) --library 'SIMPLEX'
	
metrics/report/snps_clustering-standard.pdf : metrics/summary/snps_filtered-standard.tsv
	$(call RUN, -c -n 1 -s 12G -m 16G -v $(SUPERHEAT_ENV),"$(RSCRIPT) modules/test/qc/plotmetrics.R --type 16 && \
									   					   gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dFirstPage=2 -dLastPage=2 -sOutputFile=metrics/report/snp_clustering-2.pdf metrics/report/snp_clustering.pdf && \
									   					   rm metrics/report/snp_clustering.pdf && \
									   					   mv metrics/report/snp_clustering-2.pdf metrics/report/snp_clustering.pdf")

include modules/vcf_tools/vcftools.mk

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: $(PHONY)
