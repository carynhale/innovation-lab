include innovation-lab/Makefile.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/genotype_access.$(NOW)

genotype_access : $(foreach sample,$(SAMPLES),)

REF_FASTA = /ifs/depot/resources/dmp/data/pubdata/hg-fasta/VERSIONS/hg19/Homo_sapiens_assembly19.fasta
GBCMS_PATH ?= $(HOME)/share/usr/bin/GetBaseCountsMultiSample
FILTER_DUPLICATE ?= 0
FRAGMENT_COUNT ?= 1
MAPPING_QUALITY ?= 20
THREADS ?= 10
VERBOSITY ?= INFO


define genotype-access
waltz/$1-pileup.txt.gz : bam/$1.bam
	$$(call RUN,-c -n 10 -s 4G -m 6G,"set -o pipefail && \
									  mkdir -p genotype_variants && \
									 ")
									 
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call waltz-genotype,$(sample))))
		
..DUMMY := $(shell mkdir -p version; \
			 $(GBCMS_PATH) --help &> version/genotype_access.txt; \
			 $(GENOTYPE_VARIANTS) --version >> version/genotype_access.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: genotype_access
