include modules/Makefile.inc
include modules/config/gatk.inc

LOGDIR = log/cluster_samples.$(NOW)

cluster_samples : $(foreach sample,$(SAMPLES),metrics/snps/$(sample).vcf) \
				  metrics/summary/snps_combined.vcf \
				  metrics/summary/snps_filtered.vcf

#metrics/summary/snps_filtered.tsv
#metrics/report/snps_clustering.pdf

DBSNP_SUBSET = $(HOME)/share/lib/bed_files/dbsnp_137.b37_subset.bed
CLUSTER_VCF = $(RSCRIPT) $(SCRIPTS_DIR)/qc/cluster_access.R

define genotype-snps
metrics/snps/$1.vcf : bam/$1.bam
	$$(call RUN,-n 4 -s 2.5G -m 3G,"set -o pipefail && \
									$$(call GATK_CMD,8G) -T UnifiedGenotyper -nt 4 \
									-I $$(<) \
									-R $$(REF_FASTA) \
									-o $$(@) \
									--output_mode EMIT_ALL_SITES \
									--min_base_quality_score 30 \
									--standard_min_confidence_threshold_for_calling 20 \
									-L $$(DBSNP_SUBSET)")
									
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call genotype-snps,$(sample))))
		
metrics/summary/snps_combined.vcf : $(foreach sample,$(SAMPLES),metrics/snps/$(sample).vcf)
	$(call RUN,-s 16G -m 20G,"set -o pipefail && \
							  $(call GATK_MEM,14G) -T CombineVariants $(foreach vcf,$^,--variant $(vcf) ) \
							  -o $@ --genotypemergeoption UNSORTED -R $(REF_FASTA) \
							  --disable_auto_index_creation_and_locking_when_reading_rods")
							  
metrics/summary/snps_filtered.vcf : metrics/summary/snps_combined.vcf
	$(INIT) grep '^#' $< > $@ && grep -e '0/1' -e '1/1' $< >> $@
	
#metrics/summary/snps_filtered.tsv : metrics/summary/snps_filtered.vcf
#	$(INIT) $(CLUSTER_VCF) --switch 1
#	
#metrics/report/snps_clustering.pdf : metrics/summary/snps_filtered.tsv
#	$(call RUN, -c -n 1 -s 12G -m 16G -v $(SUPERHEAT_ENV),"set -o pipefail && \
#														   $(CLUSTER_VCF) --switch 2 && \
#									   					   gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dFirstPage=2 -dLastPage=2 -dAutoRotatePages=/None -sOutputFile=metrics/report/snps_clustering-2.pdf metrics/report/snps_clustering.pdf && \
#									   					   rm metrics/report/snps_clustering.pdf && \
#									   					   mv metrics/report/snps_clustering-2.pdf metrics/report/snps_clustering.pdf")
#		

..DUMMY := $(shell mkdir -p version; \
			 echo "gatk3" > version/cluster_samples.txt; \
			 $(GATK) --version >> version/cluster_samples.txt; \
			 R --version >> version/cluster_samples.txt)
.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: cluster_samples
