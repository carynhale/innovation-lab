include innovation-lab/Makefile.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/sufam_genotype.$(NOW)

sufam_genotype : $(foreach sample,$(SAMPLES),sufam_genotype/$(sample).txt) \
		 summary/summary_genotype.txt

MPILEUP_PARAMETERS = --count-orphans \
		     --ignore-RG \
		     --min-MQ 1 \
		     --max-depth 250000 \
		     --max-idepth 250000

define sufam-genotype
sufam_genotype/$1.txt : bam/
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(SUFAM_ENV),"set -o pipefail && \
							 mkdir -p sufam_genotype && \
								      cd genotype_variants && \
								      $$(GENOTYPE_VARIANTS) small_variants all \
								      -i ../maf/$1.maf \
								      -r $$(REF_FASTA) \
								      -p $1 \
								      -b ../bam/$1_cl_aln_srt_MD_IR_FX_BR.bam \
								      -d ../bam/$1_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.bam \
								      -s ../bam/$1_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-simplex.bam \
								      -g $$(GBCMS_PATH) \
								      -fd $$(FILTER_DUPLICATES) \
								      -fc $$(FRAGMENT_COUNT) \
								      -mapq $$(MAPPING_QUALITY) \
								      -t $$(THREADS) && \
								      echo 'task completed!' > $1.taskcomplete")
									 
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call sufam-genotype,$(sample))))
		
summary/summary_genotype.txt : $(foreach sample,$(SAMPLES),genotype_variants/$(sample).taskcomplete)
	$(call RUN, -c -n 1 -s 12G -m 24G,"set -o pipefail && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/summary/genotype_access.R --sample_names '$(SAMPLES)'")

		
..DUMMY := $(shell mkdir -p version; \
	     $(SUFAM_ENV)/bin/$(SAMTOOLS) --version > version/sufam_genotype.txt; \
	     $(SUFAM_ENV)/bin/sufam --version &>> version/sufam_genotype.txt; \
	     R --version >> version/sufam_genotype.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: sufam_genotype
