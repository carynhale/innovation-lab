include innovation-lab/Makefile.inc

LOGDIR ?= log/facets.$(NOW)

facets : facets/targets/targets.vcf \
		 $(foreach pair,$(SAMPLE_PAIRS),facets/pileup/$(pair).txt.gz)

facets/targets/targets.vcf : $(TARGETS_FILE)
	$(INIT) $(BEDTOOLS) intersect -header -u -a $(DBSNP) -b $< > $@
	
define snp-pileup-tumor-normal
facets/pileup/$1_$2.txt.gz : bam/$1.bam bam/$2.bam facets/targets/targets.vcf
	$$(call RUN,-c -s 8G -m 16G -v $(ABSOLUTE_ENV),"set -o pipefail && \
													rm -f $$@ && \
													snp-pileup \
													-A \
													--min-map-quality=15 \
													--min-base-quality=15 \
													--gzip \
													--max-depth=15000 \
													$$(<<<) \
													$$@ $$(<<) \
													$$(<)")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call snp-pileup-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))


..DUMMY := $(shell mkdir -p version; \
			 $(BEDTOOLS) --version > version/facets.txt; \
			 ~/share/usr/env/cntu-0.0.1/bin/snp-pileup --help >> version/facets.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: facets
