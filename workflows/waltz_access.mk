include innovation-lab/Makefile.inc
include innovation-lab/config/waltz.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/waltz_access.$(NOW)

waltz_access : $(foreach sample,$(SAMPLES),waltz/$(sample)-pileup.txt.gz) \
			   $(foreach sample,$(SAMPLES),waltz/$(sample)__aln_srt_IR_FX-pileup.txt.gz) \
			   $(foreach sample,$(SAMPLES),waltz/$(sample)__aln_srt_IR_FX-simplex-pileup.txt.gz) \
			   $(foreach sample,$(SAMPLES),waltz/$(sample)__aln_srt_IR_FX-duplex-pileup.txt.gz) \
			   waltz/noise_metrics_with_duplicates.txt

               #waltz/noise_metrics_without_duplicates.txt \
			   #waltz/noise_by_position_standard_with_duplicates.txt \
			   #waltz/noise_by_position_standard_without_duplicates.txt \
			   #waltz/noise_by_position_simplex_without_duplicates.txt \
			   #waltz/noise_by_position_duplex_without_duplicates.txt \
			   #waltz/noise_by_position.pdf \
			   #waltz/noise_by_position.png

WALTZ_MIN_MAPQ ?= 20
TARGETS_FILE_NOMSI ?= $(HOME)/share/lib/resource_files/MSK-ACCESS-v1_0-A-good-positions-noMSI.txt
SUPERHEAT_ENV = $(HOME)/share/usr/env/r-complexheatmap-2.2.0

define waltz-genotype
waltz/$1-pileup.txt.gz : bam/$1.bam
	$$(call RUN,-c -n 4 -s 4G -m 6G,"set -o pipefail && \
									 mkdir -p waltz && \
									 cd waltz && \
									 ln -sf ../bam/$1.bam $1.bam && \
									 ln -sf ../bam/$1.bai $1.bai && \
									 if [[ ! -f '.bed' ]]; then cut -f 4 $$(TARGETS_FILE) | paste -d '\t' $$(TARGETS_FILE) - > .bed; fi && \
									 $$(call WALTZ_CMD,2G,8G) org.mskcc.juber.waltz.Waltz PileupMetrics $$(WALTZ_MIN_MAPQ) $1.bam $$(DMP_FASTA) .bed && \
									 gzip $1-pileup.txt && \
									 gzip $1-pileup-without-duplicates.txt && \
									 gzip $1-intervals.txt && \
									 gzip $1-intervals-without-duplicates.txt && \
									 cd ..")
									 
waltz/$1__aln_srt_IR_FX-pileup.txt.gz : bam/$1__aln_srt_IR_FX.bam
	$$(call RUN,-c -n 4 -s 4G -m 6G,"set -o pipefail && \
									 mkdir -p waltz && \
									 cd waltz && \
									 ln -sf ../bam/$1__aln_srt_IR_FX.bam $1__aln_srt_IR_FX.bam && \
									 ln -sf ../bam/$1__aln_srt_IR_FX.bai $1__aln_srt_IR_FX.bai && \
									 if [[ ! -f '.bed' ]]; then cut -f 4 $$(TARGETS_FILE) | paste -d '\t' $$(TARGETS_FILE) - > .bed; fi && \
									 $$(call WALTZ_CMD,2G,8G) org.mskcc.juber.waltz.Waltz PileupMetrics $$(WALTZ_MIN_MAPQ) $1__aln_srt_IR_FX.bam $$(DMP_FASTA) .bed && \
									 gzip $1__aln_srt_IR_FX-pileup.txt && \
									 gzip $1__aln_srt_IR_FX-pileup-without-duplicates.txt && \
									 gzip $1__aln_srt_IR_FX-intervals.txt && \
									 gzip $1__aln_srt_IR_FX-intervals-without-duplicates.txt && \
									 cd ..")
									 
waltz/$1__aln_srt_IR_FX-simplex-pileup.txt.gz : bam/$1__aln_srt_IR_FX-simplex.bam
	$$(call RUN,-c -n 4 -s 4G -m 6G,"set -o pipefail && \
									 mkdir -p waltz && \
									 cd waltz && \
									 ln -sf ../bam/$1__aln_srt_IR_FX-simplex.bam $1__aln_srt_IR_FX-simplex.bam && \
									 ln -sf ../bam/$1__aln_srt_IR_FX-simplex.bai $1__aln_srt_IR_FX-simplex.bai && \
									 if [[ ! -f '.bed' ]]; then cut -f 4 $$(TARGETS_FILE) | paste -d '\t' $$(TARGETS_FILE) - > .bed; fi && \
									 $$(call WALTZ_CMD,2G,8G) org.mskcc.juber.waltz.Waltz PileupMetrics $$(WALTZ_MIN_MAPQ) $1__aln_srt_IR_FX-simplex.bam $$(DMP_FASTA) .bed && \
									 gzip $1__aln_srt_IR_FX-simplex-pileup.txt && \
									 gzip $1__aln_srt_IR_FX-simplex-pileup-without-duplicates.txt && \
									 gzip $1__aln_srt_IR_FX-simplex-intervals.txt && \
									 gzip $1__aln_srt_IR_FX-simplex-intervals-without-duplicates.txt && \
									 cd ..")
									 
waltz/$1__aln_srt_IR_FX-duplex-pileup.txt.gz : bam/$1__aln_srt_IR_FX-duplex.bam
	$$(call RUN,-c -n 4 -s 4G -m 6G,"set -o pipefail && \
									 mkdir -p waltz && \
									 cd waltz && \
									 ln -sf ../bam/$1__aln_srt_IR_FX-duplex.bam $1__aln_srt_IR_FX-duplex.bam && \
									 ln -sf ../bam/$1__aln_srt_IR_FX-duplex.bai $1__aln_srt_IR_FX-duplex.bai && \
									 if [[ ! -f '.bed' ]]; then cut -f 4 $$(TARGETS_FILE) | paste -d '\t' $$(TARGETS_FILE) - > .bed; fi && \
									 $$(call WALTZ_CMD,2G,8G) org.mskcc.juber.waltz.Waltz PileupMetrics $$(WALTZ_MIN_MAPQ) $1__aln_srt_IR_FX-duplex.bam $$(DMP_FASTA) .bed && \
									 gzip $1__aln_srt_IR_FX-duplex-pileup.txt && \
									 gzip $1__aln_srt_IR_FX-duplex-pileup-without-duplicates.txt && \
									 gzip $1__aln_srt_IR_FX-duplex-intervals.txt && \
									 gzip $1__aln_srt_IR_FX-duplex-intervals-without-duplicates.txt && \
									 cd ..")


endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call waltz-genotype,$(sample))))
		
waltz/noise_metrics_with_duplicates.txt : $(wildcard waltz/$(SAMPLES)-pileup.txt.gz) $(wildcard waltz/$(SAMPLES)__aln_srt_IR_FX-pileup.txt.gz) $(wildcard waltz/$(SAMPLES)__aln_srt_IR_FX-simplex-pileup.txt.gz) $(wildcard waltz/$(SAMPLES)__aln_srt_IR_FX-duplex-pileup.txt.gz)
	$(call RUN, -c -n 1 -s 12G -m 18G,"set -o pipefail && \
									   $(RSCRIPT) $(SCRIPTS_DIR)/qc/waltz_metrics.R --type 1 --target_file $(TARGETS_FILE_NOMSI) --sample_names '$(SAMPLES)'")

waltz/noise_metrics_without_duplicates.txt : $(wildcard waltz/$(SAMPLES)-pileup.txt.gz) $(wildcard waltz/$(SAMPLES)__aln_srt_IR_FX-pileup.txt.gz) $(wildcard waltz/$(SAMPLES)__aln_srt_IR_FX-simplex-pileup.txt.gz) $(wildcard waltz/$(SAMPLES)__aln_srt_IR_FX-duplex-pileup.txt.gz)
	$(call RUN, -c -n 1 -s 12G -m 18G,"set -o pipefail && \
									   $(RSCRIPT) $(SCRIPTS_DIR)/qc/waltz_metrics.R --type 2 --target_file $(TARGETS_FILE_NOMSI) --sample_names '$(SAMPLES)'")
									  
waltz/noise_by_position_standard_with_duplicates.txt : $(wildcard waltz/$(SAMPLES)-pileup.txt.gz)
	$(call RUN, -c -n 1 -s 12G -m 24G,"set -o pipefail && \
									   $(RSCRIPT) $(SCRIPTS_DIR)/qc/waltz_metrics.R --type 3 --target_file $(TARGETS_FILE_NOMSI) --sample_names '$(SAMPLES)'")
									  
waltz/noise_by_position_standard_without_duplicates.txt : $(wildcard waltz/$(SAMPLES)-pileup.txt.gz)
	$(call RUN, -c -n 1 -s 12G -m 24G,"set -o pipefail && \
									   $(RSCRIPT) $(SCRIPTS_DIR)/qc/waltz_metrics.R --type 4 --target_file $(TARGETS_FILE_NOMSI) --sample_names '$(SAMPLES)'")
									  
waltz/noise_by_position_simplex_without_duplicates.txt : $(wildcard waltz/$(SAMPLES)-aln_srt_IR_FX-simplex-pileup.txt.gz)
	$(call RUN, -c -n 1 -s 12G -m 24G,"set -o pipefail && \
									   $(RSCRIPT) $(SCRIPTS_DIR)/qc/waltz_metrics.R --type 5 --target_file $(TARGETS_FILE_NOMSI) --sample_names '$(SAMPLES)'")

waltz/noise_by_position_duplex_without_duplicates.txt : $(wildcard waltz/$(SAMPLES)-aln_srt_IR_FX-duplex-pileup.txt.gz)
	$(call RUN, -c -n 1 -s 12G -m 24G,"set -o pipefail && \
									   $(RSCRIPT) $(SCRIPTS_DIR)/qc/waltz_metrics.R --type 6 --target_file $(TARGETS_FILE_NOMSI) --sample_names '$(SAMPLES)'")

waltz/noise_by_position.pdf : waltz/noise_by_position_standard_with_duplicates.txt waltz/noise_by_position_standard_without_duplicates.txt waltz/noise_by_position_simplex_without_duplicates.txt waltz/noise_by_position_duplex_without_duplicates.txt
	$(call RUN, -c -n 1 -s 36G -m 48G -v $(SUPERHEAT_ENV),"set -o pipefail && \
									   					   $(RSCRIPT) $(SCRIPTS_DIR)/qc/waltz_metrics.R --type 7 --target_file $(TARGETS_FILE_NOMSI) --sample_names '$(SAMPLES)'")

waltz/noise_by_position.png : waltz/noise_by_position.pdf
	$(call RUN, -c -n 1 -s 4G -m 8G,"set -o pipefail && \
									 gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dFirstPage=2 -dLastPage=2 -sOutputFile=waltz/noise_by_position-2.pdf waltz/noise_by_position.pdf && \
									 mv waltz/noise_by_position-2.pdf waltz/noise_by_position.pdf && \
									 convert -antialias -background white -quality 100 -sharpen 0x1.0 -density 150 waltz/noise_by_position.pdf waltz/noise_by_position.png")


..DUMMY := $(shell mkdir -p version; \
			 $(JAVA8) -version &> version/waltz_access.txt; \
			 R --version >> version/waltz_access.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: waltz_access
