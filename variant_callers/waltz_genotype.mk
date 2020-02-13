include innovation-lab/Makefile.inc
include innovation-lab/config/waltz.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/waltz_genotype.$(NOW)

waltz_genotype : $(foreach sample,$(SAMPLES),waltz/$(sample)-STANDARD-pileup.txt.gz)
#				 $(foreach sample,$(SAMPLES),waltz/$(sample)-COLLAPSED-pileup.txt.gz) \
#				 $(foreach sample,$(SAMPLES),waltz/$(sample)-SIMPLEX-pileup.txt.gz) \
#				 $(foreach sample,$(SAMPLES),waltz/$(sample)-DUPLEX-pileup.txt.gz)

WALTZ_MIN_MAPQ ?= 20

define waltz-genotype
waltz/$1-STANDARD-pileup.txt.gz : bam/$1-STANDARD.bam
	$$(call RUN,-c -n 4 -s 4G -m 6G,"set -o pipefail && \
									 mkdir -p waltz && \
									 cd waltz && \
									 ln -sf ../bam/$1-STANDARD.bam $1-STANDARD.bam && \
									 ln -sf ../bam/$1-STANDARD.bai $1-STANDARD.bai && \
									 $$(call WALTZ_CMD,2G,8G) org.mskcc.juber.waltz.Waltz PileupMetrics $(WALTZ_MIN_MAPQ) $1-STANDARD.bam $(DMP_FASTA) $(TARGETS_FILE) && \
									 gzip $1-STANDARD-pileup.txt && \
									 gzip $1-STANDARD-pileup-without-duplicates.txt && \
									 gzip $1-STANDARD-intervals.txt && \
									 gzip $1-STANDARD-intervals-without-duplicates.txt && \
									 cd ..")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call waltz-genotype,$(sample))))


..DUMMY := $(shell mkdir -p version)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: waltz_genotype
