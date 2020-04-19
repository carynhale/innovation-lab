include innovation-lab/Makefile.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/genotype_access.$(NOW)

genotype_access : $(foreach sample,$(SAMPLES),genotype_variants/$(sample).taskcomplete)

REF_FASTA = /ifs/depot/resources/dmp/data/pubdata/hg-fasta/VERSIONS/hg19/Homo_sapiens_assembly19.fasta
GBCMS_PATH ?= $(HOME)/share/usr/bin/GetBaseCountsMultiSample
FILTER_DUPLICATES ?= 0
FRAGMENT_COUNT ?= 1
MAPPING_QUALITY ?= 20
THREADS ?= 10
VERBOSITY ?= INFO

define genotype-access
genotype_variants/$1.taskcomplete : bam/$1_cl_aln_srt_MD_IR_FX_BR.bam bam/$1_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX.bam bam/$1_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-simplex.bam bam/$1_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.bam
	$$(call RUN,-c -n 10 -s 1G -m 2G -v $(GENOTYPE_VARIANTS_ENV),"set -o pipefail && \
									  							  mkdir -p genotype_variants && \
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
		$(eval $(call genotype-access,$(sample))))
		
..DUMMY := $(shell mkdir -p version; \
			 $(GBCMS_PATH) --help &> version/genotype_access.txt; \
			 echo 'genotype_variants' >> version/genotype_access.txt; \
			 $(GENOTYPE_VARIANTS_ENV)/bin/$(GENOTYPE_VARIANTS) --version >> version/genotype_access.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: genotype_access
