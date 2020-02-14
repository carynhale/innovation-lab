include innovation-lab/Makefile.inc
include innovation-lab/config/waltz.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/waltz_genotype.$(NOW)

waltz_genotype : $(foreach sample,$(SAMPLES),waltz/$(sample)-STANDARD-pileup.txt.gz) \
				 $(foreach sample,$(SAMPLES),waltz/$(sample)-COLLAPSED-pileup.txt.gz) \
				 $(foreach sample,$(SAMPLES),waltz/$(sample)-SIMPLEX-pileup.txt.gz) \
				 $(foreach sample,$(SAMPLES),waltz/$(sample)-DUPLEX-pileup.txt.gz) \
				 $(foreach sample,$(SAMPLES),waltz/noise_metrics_with_duplicates.txt) \
				 $(foreach sample,$(SAMPLES),waltz/noise_metrics_without_duplicates.txt)

WALTZ_MIN_MAPQ ?= 15
TARGETS_FILE_NOMSI ?= $(HOME)/share/lib/resource_files/MSK-ACCESS-v1_0-A-good-positions-noMSI.txt

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
									 
waltz/$1-COLLAPSED-pileup.txt.gz : bam/$1-COLLAPSED.bam
	$$(call RUN,-c -n 4 -s 4G -m 6G,"set -o pipefail && \
									 mkdir -p waltz && \
									 cd waltz && \
									 ln -sf ../bam/$1-COLLAPSED.bam $1-COLLAPSED.bam && \
									 ln -sf ../bam/$1-COLLAPSED.bai $1-COLLAPSED.bai && \
									 $$(call WALTZ_CMD,2G,8G) org.mskcc.juber.waltz.Waltz PileupMetrics $(WALTZ_MIN_MAPQ) $1-COLLAPSED.bam $(DMP_FASTA) $(TARGETS_FILE) && \
									 gzip $1-COLLAPSED-pileup.txt && \
									 gzip $1-COLLAPSED-pileup-without-duplicates.txt && \
									 gzip $1-COLLAPSED-intervals.txt && \
									 gzip $1-COLLAPSED-intervals-without-duplicates.txt && \
									 cd ..")
									 
waltz/$1-SIMPLEX-pileup.txt.gz : bam/$1-SIMPLEX.bam
	$$(call RUN,-c -n 4 -s 4G -m 6G,"set -o pipefail && \
									 mkdir -p waltz && \
									 cd waltz && \
									 ln -sf ../bam/$1-SIMPLEX.bam $1-SIMPLEX.bam && \
									 ln -sf ../bam/$1-SIMPLEX.bai $1-SIMPLEX.bai && \
									 $$(call WALTZ_CMD,2G,8G) org.mskcc.juber.waltz.Waltz PileupMetrics $(WALTZ_MIN_MAPQ) $1-SIMPLEX.bam $(DMP_FASTA) $(TARGETS_FILE) && \
									 gzip $1-SIMPLEX-pileup.txt && \
									 gzip $1-SIMPLEX-pileup-without-duplicates.txt && \
									 gzip $1-SIMPLEX-intervals.txt && \
									 gzip $1-SIMPLEX-intervals-without-duplicates.txt && \
									 cd ..")
									 
waltz/$1-DUPLEX-pileup.txt.gz : bam/$1-DUPLEX.bam
	$$(call RUN,-c -n 4 -s 4G -m 6G,"set -o pipefail && \
									 mkdir -p waltz && \
									 cd waltz && \
									 ln -sf ../bam/$1-DUPLEX.bam $1-DUPLEX.bam && \
									 ln -sf ../bam/$1-DUPLEX.bai $1-DUPLEX.bai && \
									 $$(call WALTZ_CMD,2G,8G) org.mskcc.juber.waltz.Waltz PileupMetrics $(WALTZ_MIN_MAPQ) $1-DUPLEX.bam $(DMP_FASTA) $(TARGETS_FILE) && \
									 gzip $1-DUPLEX-pileup.txt && \
									 gzip $1-DUPLEX-pileup-without-duplicates.txt && \
									 gzip $1-DUPLEX-intervals.txt && \
									 gzip $1-DUPLEX-intervals-without-duplicates.txt && \
									 cd ..")


endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call waltz-genotype,$(sample))))
		
waltz/noise_metrics_with_duplicates.txt : $(wildcard waltz/$(SAMPLES)-STANDARD-pileup.txt.gz) $(wildcard waltz/$(SAMPLES)-COLLAPSED-pileup.txt.gz) $(wildcard waltz/$(SAMPLES)-SIMPLEX-pileup.txt.gz) $(wildcard waltz/$(SAMPLES)-DUPLEX-pileup.txt.gz)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/waltz_metrics.R --type 1 --target_file $(TARGETS_FILE_NOMSI) --sample_names '$(SAMPLES)'")

waltz/noise_metrics_without_duplicates.txt : $(wildcard waltz/$(SAMPLES)-STANDARD-pileup.txt.gz) $(wildcard waltz/$(SAMPLES)-COLLAPSED-pileup.txt.gz) $(wildcard waltz/$(SAMPLES)-SIMPLEX-pileup.txt.gz) $(wildcard waltz/$(SAMPLES)-DUPLEX-pileup.txt.gz)
	$(call RUN, -c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  $(RSCRIPT) $(SCRIPTS_DIR)/qc/waltz_metrics.R --type 1 --target_file $(TARGETS_FILE_NOMSI) --sample_names '$(SAMPLES)'")


..DUMMY := $(shell mkdir -p version; \
			 $(JAVA8) -version &> version/waltz_genotype.txt; \
			 R --version >> version/waltz_genotype.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: waltz_genotype
