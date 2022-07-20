include innovation-lab/Makefile.inc
include innovation-lab/config/gatk.inc

LOGDIR ?= log/mutect.$(NOW)

MUTECT_MAX_ALT_IN_NORMAL ?= 500
MUTECT_MAX_ALT_IN_NORMAL_FRACTION ?= 0.05
MUTECT_FILTERS = DuplicateRead \
		 FailsVendorQualityCheck \
		 NotPrimaryAlignment \
		 MappingQualityUnavailable \
		 UnmappedRead \
		 BadMate \
		 BadCigar
MUTECT_OPTS = --enable_extended_output \
	      --max_alt_alleles_in_normal_count $(MUTECT_MAX_ALT_IN_NORMAL) \
	      --max_alt_allele_in_normal_fraction $(MUTECT_MAX_ALT_IN_NORMAL_FRACTION) \
	      -R $(REF_FASTA) \
	      --dbsnp $(DBSNP) \
	      $(foreach ft,$(MUTECT_FILTERS),-rf $(ft))
	       
BED_SPLIT = 500
BED_CHUNKS = $(shell seq 1 $(BED_SPLIT))

mutect : mutect/bed/taskcomplete \
	 $(foreach pair,$(SAMPLE_PAIRS), \
		  	$(foreach n,$(BED_CHUNKS),mutect/$(pair)/$(pair)--$(n).vcf))
#	 $(foreach pair,$(SAMPLE_PAIRS), \
#		  	$(foreach n,$(BED_CHUNKS),mutect/$(pair)/$(pair)--$(n).maf)) \ 		
#	 $(foreach pair,$(SAMPLE_PAIRS), \
#		  	$(foreach n,$(BED_CHUNKS),mutect/$(pair)/$(pair)--$(n).tsv)) \
#	 $(foreach pair,$(SAMPLE_PAIRS),mutect/$(pair)/$(pair).txt)


mutect/bed/taskcomplete : $(TARGETS_FILE)
	$(call RUN, -c -n 1 -s 2G -m 4G,"set -o pipefail && \
					 $(RSCRIPT) $(SCRIPTS_DIR)/bed_tools/split_bed.R --n_splits $(BED_SPLIT) --bed_file $(TARGETS_FILE)")

define mutect-tumor-normal
mutect/$1_$2/$1_$2--$3.vcf : bam/$1.bam bam/$2.bam mutect/bed/taskcomplete
	$$(call RUN,-c -s 6G -m 9G,"set -o pipefail && \
				    $$(MUTECT) \
				    --tumor_sample_name $1\
				    --normal_sample_name $2 \
				    --intervals mutect/bed/$3.bed \
				    -I:tumor $$(<) \
				    -I:normal $$(<<) \
				    --out mutect/$1_$2/$1_$2--$3.txt \
				    -vcf $$(@) \
				    --coverage_file mutect/$1_$2/$1_$2--$3.wig")
endef
$(foreach chunk,$(BED_CHUNKS), \
	$(foreach pair,$(SAMPLE_PAIRS), \
			$(eval $(call mutect-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)),$(chunk)))))

..DUMMY := $(shell mkdir -p version; \
	     echo "$(MUTECT) &> version/mutect.txt")
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: mutect
