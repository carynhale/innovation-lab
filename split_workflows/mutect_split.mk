include innovation-lab/Makefile.inc
include innovation-lab/config/gatk.inc

LOGDIR ?= log/mutect_split.$(NOW)

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
	       
BED_SPLIT ?= $(MUTECT_SPLITS)
BED_CHUNKS = $(shell seq 1 $(BED_SPLIT))

mutect : mutect/bed_chunks/taskcomplete \
	 $(foreach pair,$(SAMPLE_PAIRS), \
		  	$(foreach n,$(BED_CHUNKS),mutect/$(pair)/$(pair)--$(n).vcf)) \
	 $(foreach pair,$(SAMPLE_PAIRS), \
		  	$(foreach n,$(BED_CHUNKS),mutect/$(pair)/$(pair)--$(n).ft.vcf)) \
	 $(foreach pair,$(SAMPLE_PAIRS), \
		  	$(foreach n,$(BED_CHUNKS),mutect/$(pair)/$(pair)--$(n).ft.maf)) \
	 $(foreach pair,$(SAMPLE_PAIRS),mutect/$(pair)/$(pair).ft.vcf) \
	 $(foreach pair,$(SAMPLE_PAIRS),mutect/$(pair)/$(pair).ft.maf)


mutect/bed_chunks/$(BED_SPLIT).bed : $(TARGETS_FILE)
	$(call RUN, -c -n 1 -s 2G -m 4G,"set -o pipefail && \
					 $(RSCRIPT) $(SCRIPTS_DIR)/bed_tools/split_bed.R --n_splits $(BED_SPLIT) --bed_file $(TARGETS_FILE)")

define mutect-tumor-normal
mutect/$1_$2/$1_$2--$3.vcf : bam/$1.bam bam/$2.bam mutect/bed/taskcomplete
	$$(call RUN,-c -s 6G -m 9G,"set -o pipefail && \
				    $$(MUTECT) \
				    $$(MUTECT_OPTS) \
				    --tumor_sample_name $1\
				    --normal_sample_name $2 \
				    --intervals mutect/bed/$3.bed \
				    -I:tumor $$(<) \
				    -I:normal $$(<<) \
				    --out mutect/$1_$2/$1_$2--$3.txt \
				    -vcf $$(@) \
				    --coverage_file mutect/$1_$2/$1_$2--$3.wig")
				    
mutect/$1_$2/$1_$2--$3.ft.vcf : mutect/$1_$2/$1_$2--$3.vcf
	$$(call RUN,-c -n 1 -s 1G -m 2G,"set -o pipefail && \
					 grep '^##' $$(<) > $$(@) && \
					 grep 'PASS\|#' $$(<) | grep -v '^##' >> $$(@)")
					 
mutect/$1_$2/$1_$2--$3.ft.maf : mutect/$1_$2/$1_$2--$3.ft.vcf
	$$(call RUN,-c -n 4 -s 1G -m 2G -v $(VCF2MAF_ENV) -w 72:00:00,"set -o pipefail && \
									$$(VCF2MAF) \
									--input-vcf $$(<) \
									--output-maf $$(@) \
									--tmp-dir $$(TMPDIR) \
									--tumor-id $1 \
									--normal-id $2 \
									--vep-path $$(VCF2MAF_ENV)/bin \
									--vep-data $$(HOME)/share/lib/resource_files/VEP/GRCh37/ \
									--vep-forks 4 \
									--ref-fasta $$(HOME)/share/lib/resource_files/VEP/GRCh37/homo_sapiens/99_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
									--filter-vcf $$(HOME)/share/lib/resource_files/VEP/GRCh37/homo_sapiens/99_GRCh37/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz \
									--species homo_sapiens \
									--ncbi-build GRCh37 \
									--maf-center MSKCC && \
									$$(RM) $$(TMPDIR)/$1_$2--$3.ft.vep.vcf")

endef
$(foreach chunk,$(BED_CHUNKS), \
	$(foreach pair,$(SAMPLE_PAIRS), \
			$(eval $(call mutect-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)),$(chunk)))))
			
			
define merge-splits
mutect/$1_$2/$1_$2.ft.vcf : $(foreach n,$(BED_CHUNKS),mutect/$1_$2/$1_$2--$(n).ft.vcf)
	$$(call RUN,-c -s 12G -m 18G,"set -o pipefail && \
				      grep '^#' mutect/$1_$2/$1_$2--1.ft.vcf > $$(@) && \
				      $(RSCRIPT) $(SCRIPTS_DIR)/vcf_tools/concat_vcf.R \
				      --vcf_in '$$(^)' \
				      --vcf_out $$(@)")
				      
mutect/$1_$2/$1_$2.ft.maf : $(foreach n,$(BED_CHUNKS),mutect/$1_$2/$1_$2--$(n).ft.maf)
	$$(call RUN,-c -s 12G -m 18G,"set -o pipefail && \
				      grep '^#' mutect/$1_$2/$1_$2--1.ft.maf > $$(@) && \
				      $(RSCRIPT) $(SCRIPTS_DIR)/maf_tools/concat_maf.R \
				      --maf_in '$$(^)' \
				      --maf_out $$(@)")

endef
$(foreach pair,$(SAMPLE_PAIRS),\
	       $(eval $(call merge-splits,$(tumor.$(pair)),$(normal.$(pair)))))
	       
..DUMMY := $(shell mkdir -p version; \
	     echo "$(MUTECT) &> version/mutect_split.txt")
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: mutect
