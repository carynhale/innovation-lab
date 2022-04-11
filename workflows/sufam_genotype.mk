include innovation-lab/Makefile.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/sufam_genotype.$(NOW)

sufam_genotype : $(foreach sample,$(SAMPLES),sufam_genotype/$(sample).txt)
#		 summary/summary_genotype.txt

MPILEUP_PARAMETERS = "--count-orphans \
		     --ignore-RG \
		     --min-MQ 1 \
		     --max-depth 250000 \
		     --max-idepth 250000"

define sufam-genotype
sufam_genotype/$1.txt : bam/$1.bam vcf/$1.vcf
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(SUFAM_ENV),"set -o pipefail && \
							 mkdir -p sufam_genotype && \
							 $$(SUFAM) \
							 --sample_name $1 \
							 --format sufam \
							 --mpileup-parameters $$(MPILEUP_PARAMETERS) \
							 $$(REF_FASTA) \
							 vcf/$1.vcf \
							 bam/$1.bam \
							 2> sufam_genotype/$1.log \
							 > sufam_genotype/$1.txt")
									 
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call sufam-genotype,$(sample))))
		
summary/summary_genotype.txt : $(foreach sample,$(SAMPLES),sufam_genotype/$(sample).txt)
	$(call RUN, -c -n 1 -s 12G -m 24G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/summary/sufam_summary.R --sample_names '$(SAMPLES)'")

		
..DUMMY := $(shell mkdir -p version; \
	     $(SUFAM_ENV)/bin/$(SAMTOOLS) --version > version/sufam_genotype.txt; \
	     $(SUFAM_ENV)/bin/$(SUFAM) --version &>> version/sufam_genotype.txt; \
	     R --version >> version/sufam_genotype.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: sufam_genotype
