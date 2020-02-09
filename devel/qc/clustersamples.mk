include modules/Makefile.inc
include modules/config/gatk.inc

LOGDIR = log/cluster_samples.$(NOW)
PHONE += metrics metrics/summary metrics/report

cluster_samples : $(foreach sample,$(SAMPLES),metrics/standard/$(sample)-snps.vcf) \
				  metrics/summary/snps_combined-standard.vcf \
				  metrics/summary/snps_filtered-standard.vcf \
				  metrics/summary/snps_filtered-standard.tsv \
				  $(foreach sample,$(SAMPLES),metrics/unfiltered/$(sample)-snps.vcf) \
				  metrics/summary/snps_combined-unfiltered.vcf \
				  metrics/summary/snps_filtered-unfiltered.vcf \
				  metrics/summary/snps_filtered-unfiltered.tsv \
				  $(foreach sample,$(SAMPLES),metrics/simplex/$(sample)-snps.vcf) \
				  metrics/summary/snps_combined-simplex.vcf \
				  metrics/summary/snps_filtered-simplex.vcf \
				  metrics/summary/snps_filtered-simplex.tsv \
				  $(foreach sample,$(SAMPLES),metrics/duplex/$(sample)-snps.vcf) \
				  metrics/summary/snps_combined-duplex.vcf \
				  metrics/summary/snps_filtered-duplex.vcf \
				  metrics/summary/snps_filtered-duplex.tsv \
				  metrics/report/snps_clustering-standard.pdf \
				  metrics/report/snps_clustering-unfiltered.pdf \
				  metrics/report/snps_clustering-simplex.pdf \
				  metrics/report/snps_clustering-duplex.pdf

DBSNP_SUBSET = $(HOME)/share/reference/dbsnp_tseq_intersect.bed
CLUSTER_VCF = $(RSCRIPT) modules/test/qc/clustersamples.R
POOL_AB_INTERVAL ?= /home/${USER}/share/reference/target_panels/MSK-ACCESS-v1_0-probe-AB.sorted.list

define genotype-snps-standard
metrics/standard/$1-snps.vcf : bam/$1-standard.bam
	$$(call RUN,-n 4 -s 2.5G -m 3G,"set -o pipefail && \
									$(call GATK_MEM,8G) -T UnifiedGenotyper -nt 4 -R $(REF_FASTA) --dbsnp $(DBSNP) \
									-I $$(<) \
									-L $(DBSNP_SUBSET) \
									-o $$(@) \
									--output_mode EMIT_ALL_SITES \
									--min_base_quality_score 30 \
									--standard_min_confidence_threshold_for_calling 20 \
									-L $(POOL_AB_INTERVAL)")
									
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call genotype-snps-standard,$(sample))))
		
metrics/summary/snps_combined-standard.vcf : $(foreach sample,$(SAMPLES),metrics/standard/$(sample)-snps.vcf)
	$(call RUN,-s 16G -m 20G,"set -o pipefail && \
							  $(call GATK_MEM,14G) -T CombineVariants $(foreach vcf,$^,--variant $(vcf) ) \
							  -o $@ --genotypemergeoption UNSORTED -R $(REF_FASTA) \
							  --disable_auto_index_creation_and_locking_when_reading_rods")
							  
metrics/summary/snps_filtered-standard.vcf : metrics/summary/snps_combined-standard.vcf
	$(INIT) grep '^#' $< > $@ && grep -e '0/1' -e '1/1' $< >> $@
	
metrics/summary/snps_filtered-standard.tsv : metrics/summary/snps_filtered-standard.vcf
	$(INIT) $(CLUSTER_VCF) --library 'STANDARD'
	
metrics/report/snps_clustering-standard.pdf : metrics/summary/snps_filtered-standard.tsv
	$(call RUN, -c -n 1 -s 12G -m 16G -v $(SUPERHEAT_ENV),"set -o pipefail && \
														   $(RSCRIPT) modules/test/qc/plotmetrics.R --type 15 && \
									   					   gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dFirstPage=2 -dLastPage=2 -sOutputFile=metrics/report/snps_clustering-standard-2.pdf metrics/report/snps_clustering-standard.pdf && \
									   					   rm metrics/report/snps_clustering-standard.pdf && \
									   					   mv metrics/report/snps_clustering-standard-2.pdf metrics/report/snps_clustering-standard.pdf")
		
define genotype-snps-unfiltered
metrics/unfiltered/$1-snps.vcf : bam/$1-unfiltered.bam
	$$(call RUN,-n 4 -s 2.5G -m 3G,"set -o pipefail && \
									$(call GATK_MEM,8G) -T UnifiedGenotyper -nt 4 -R $(REF_FASTA) --dbsnp $(DBSNP) \
									-I $$(<) \
									-L $(DBSNP_SUBSET) \
									-o $$(@) \
									--output_mode EMIT_ALL_SITES \
									--allow_potentially_misencoded_quality_scores \
									--min_base_quality_score 10 \
									--standard_min_confidence_threshold_for_calling 20 \
									-L $(POOL_AB_INTERVAL)")
									
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call genotype-snps-unfiltered,$(sample))))
		
metrics/summary/snps_combined-unfiltered.vcf : $(foreach sample,$(SAMPLES),metrics/unfiltered/$(sample)-snps.vcf)
	$(call RUN,-s 16G -m 20G,"set -o pipefail && \
							  $(call GATK_MEM,14G) -T CombineVariants $(foreach vcf,$^,--variant $(vcf) ) \
							  -o $@ --genotypemergeoption UNSORTED -R $(REF_FASTA)")
							  
metrics/summary/snps_filtered-unfiltered.vcf : metrics/summary/snps_combined-unfiltered.vcf
	$(INIT) grep '^#' $< > $@ && grep -e '0/1' -e '1/1' $< >> $@
	
metrics/summary/snps_filtered-unfiltered.tsv : metrics/summary/snps_filtered-unfiltered.vcf
	$(INIT) $(CLUSTER_VCF) --library 'UNFILTERED'
	
metrics/report/snps_clustering-unfiltered.pdf : metrics/summary/snps_filtered-unfiltered.tsv
	$(call RUN, -c -n 1 -s 12G -m 16G -v $(SUPERHEAT_ENV),"set -o pipefail && \
														   $(RSCRIPT) modules/test/qc/plotmetrics.R --type 16 && \
									   					   gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dFirstPage=2 -dLastPage=2 -sOutputFile=metrics/report/snps_clustering-unfiltered-2.pdf metrics/report/snps_clustering-unfiltered.pdf && \
									   					   rm metrics/report/snps_clustering-unfiltered.pdf && \
									   					   mv metrics/report/snps_clustering-unfiltered-2.pdf metrics/report/snps_clustering-unfiltered.pdf")
		
define genotype-snps-simplex
metrics/simplex/$1-snps.vcf : bam/$1-simplex.bam
	$$(call RUN,-n 4 -s 2.5G -m 3G,"set -o pipefail && \
									$(call GATK_MEM,8G) -T UnifiedGenotyper -nt 4 -R $(REF_FASTA) --dbsnp $(DBSNP) \
									-I $$(<) \
									-L $(DBSNP_SUBSET) \
									-o $$(@) \
									--output_mode EMIT_ALL_SITES \
									--allow_potentially_misencoded_quality_scores \
									--min_base_quality_score 10 \
									--standard_min_confidence_threshold_for_calling 20 \
									-L $(POOL_AB_INTERVAL)")
									
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call genotype-snps-simplex,$(sample))))

metrics/summary/snps_combined-simplex.vcf : $(foreach sample,$(SAMPLES),metrics/simplex/$(sample)-snps.vcf)
	$(call RUN,-s 16G -m 20G,"set -o pipefail && \
							  $(call GATK_MEM,14G) -T CombineVariants $(foreach vcf,$^,--variant $(vcf) ) \
							  -o $@ --genotypemergeoption UNSORTED -R $(REF_FASTA)")
							  
metrics/summary/snps_filtered-simplex.vcf : metrics/summary/snps_combined-simplex.vcf
	$(INIT) grep '^#' $< > $@ && grep -e '0/1' -e '1/1' $< >> $@
	
metrics/summary/snps_filtered-simplex.tsv : metrics/summary/snps_filtered-simplex.vcf
	$(INIT) $(CLUSTER_VCF) --library 'SIMPLEX'
	
metrics/report/snps_clustering-simplex.pdf : metrics/summary/snps_filtered-simplex.tsv
	$(call RUN, -c -n 1 -s 12G -m 16G -v $(SUPERHEAT_ENV),"set -o pipefail && \
														   $(RSCRIPT) modules/test/qc/plotmetrics.R --type 17 && \
									   					   gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dFirstPage=2 -dLastPage=2 -sOutputFile=metrics/report/snps_clustering-simplex-2.pdf metrics/report/snps_clustering-simplex.pdf && \
									   					   rm metrics/report/snps_clustering-simplex.pdf && \
									   					   mv metrics/report/snps_clustering-simplex-2.pdf metrics/report/snps_clustering-simplex.pdf")
	
define genotype-snps-duplex
metrics/duplex/$1-snps.vcf : bam/$1-duplex.bam
	$$(call RUN,-n 4 -s 2.5G -m 3G,"set -o pipefail && \
									$(call GATK_MEM,8G) -T UnifiedGenotyper -nt 4 -R $(REF_FASTA) --dbsnp $(DBSNP) \
									-I $$(<) \
									-L $(DBSNP_SUBSET) \
									-o $$(@) \
									--output_mode EMIT_ALL_SITES \
									--allow_potentially_misencoded_quality_scores \
									--min_base_quality_score 10 \
									--standard_min_confidence_threshold_for_calling 20 \
									-L $(POOL_AB_INTERVAL)")
									
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call genotype-snps-duplex,$(sample))))
		
metrics/summary/snps_combined-duplex.vcf : $(foreach sample,$(SAMPLES),metrics/duplex/$(sample)-snps.vcf)
	$(call RUN,-s 16G -m 20G,"set -o pipefail && \
							  $(call GATK_MEM,14G) -T CombineVariants $(foreach vcf,$^,--variant $(vcf) ) \
							  -o $@ --genotypemergeoption UNSORTED -R $(REF_FASTA)")
							  
metrics/summary/snps_filtered-duplex.vcf : metrics/summary/snps_combined-duplex.vcf
	$(INIT) grep '^#' $< > $@ && grep -e '0/1' -e '1/1' $< >> $@
	
metrics/summary/snps_filtered-duplex.tsv : metrics/summary/snps_filtered-duplex.vcf
	$(INIT) $(CLUSTER_VCF) --library 'DUPLEX'
	
metrics/report/snps_clustering-duplex.pdf : metrics/summary/snps_filtered-duplex.tsv
	$(call RUN, -c -n 1 -s 12G -m 16G -v $(SUPERHEAT_ENV),"set -o pipefail && \
														   $(RSCRIPT) modules/test/qc/plotmetrics.R --type 18 && \
									   					   gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dFirstPage=2 -dLastPage=2 -sOutputFile=metrics/report/snps_clustering-duplex-2.pdf metrics/report/snps_clustering-duplex.pdf && \
									   					   rm metrics/report/snps_clustering-duplex.pdf && \
									   					   mv metrics/report/snps_clustering-duplex-2.pdf metrics/report/snps_clustering-duplex.pdf")
	


include modules/vcf_tools/vcftools.mk

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: $(PHONY)
