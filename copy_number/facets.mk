include innovation-lab/Makefile.inc
include innovation-lab/config/facets.inc

LOGDIR ?= log/facets.$(NOW)

facets : facets/targets/targets.vcf \
	 $(foreach pair,$(SAMPLE_PAIRS),facets/pileup/$(pair).txt.gz) \
	 $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).RData)
#	 $(foreach pair,$(SAMPLE_PAIRS),facets/plots/log2/$(pair).pdf)

facets/targets/targets.vcf : $(TARGETS_FILE)
	$(call RUN,-c -s 8G -m 16G,"set -o pipefail && \
				    $(BEDTOOLS) intersect -header -u -a $(DBSNP) -b $(<) > $(@)")
	
define snp-pileup-tumor-normal
facets/pileup/$1_$2.txt.gz : bam/$1.bam bam/$2.bam facets/targets/targets.vcf
	$$(call RUN,-c -s 8G -m 16G -v $(FACETS_ENV),"set -o pipefail && \
						      rm -f $$@ && \
						      snp-pileup \
						      -A \
						      --min-map-quality=15 \
						      --min-base-quality=15 \
						      --gzip \
						      --max-depth=15000 \
						      $$(<<<) \
						      $$@ \
						      $$(<<) \
						      $$(<)")

endef
$(foreach pair,$(SAMPLE_PAIRS),\
	       $(eval $(call snp-pileup-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

define facets-tumor-normal
facets/cncf/$1_$2.RData : facets/pileup/$1_$2.txt.gz
	$$(call RUN,-c -s 8G -m 16G -v $(FACETS_ENV),"set -o pipefail && \
						      $(RSCRIPT) $(SCRIPTS_DIR)/copy_number/facets.R \
						      --option 1 \
						      --sample_name $1_$2 \
						      $(call FACETS_OPTS,$*)")
													
facets/plots/log2/$1_$2.pdf : facets/cncf/$1_$2.RData
	$$(call RUN,-c -s 8G -m 16G -v $(FACETS_ENV),"set -o pipefail && \
						      $(RSCRIPT) $(SCRIPTS_DIR)/copy_number/facets.R \
						      --option 2 \
						      --sample_name $1_$2")

endef
$(foreach pair,$(SAMPLE_PAIRS),\
	       $(eval $(call facets-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))


..DUMMY := $(shell mkdir -p version; \
	     $(BEDTOOLS) --version > version/facets.txt; \
	     $(FACETS_ENV)/bin/snp-pileup --help >> version/facets.txt; \
	     $(FACETS_ENV)/bin/R --version >> version/facets.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: facets
