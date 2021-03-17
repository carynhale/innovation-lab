include innovation-lab/Makefile.inc
include innovation-lab/config/waltz.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/pileup_metrics.$(NOW)

pileup_metrics : $(foreach sample,$(SAMPLES),waltz/$(sample)_cl_aln_srt_MD_IR_FX2-pileup.txt.gz) \
		 $(foreach sample,$(SAMPLES),waltz/$(sample)_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX-pileup.txt.gz) \
		 $(foreach sample,$(SAMPLES),waltz/$(sample)_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_SIMPLEX-pileup.txt.gz) \
		 $(foreach sample,$(SAMPLES),waltz/$(sample)_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_DUPLEX-pileup.txt.gz)


WALTZ_MIN_MAPQ ?= 30
TARGETS_FILE_NOMSI ?= $(HOME)/share/lib/resource_files/MSK-ACCESS-v1_0-A-good-positions-noMSI.txt

define waltz-genotype
waltz/$1_cl_aln_srt_MD_IR_FX2-pileup.txt.gz : bam/$1_cl_aln_srt_MD_IR_FX2.bam
	$$(call RUN,-c -n 4 -s 4G -m 6G,"set -o pipefail && \
					 mkdir -p waltz && \
					 cd waltz && \
					 ln -sf ../bam/$1_cl_aln_srt_MD_IR_FX2.bam $1_cl_aln_srt_MD_IR_FX2.bam && \
					 ln -sf ../bam/$1_cl_aln_srt_MD_IR_FX2.bai $1_cl_aln_srt_MD_IR_FX2.bai && \
					 if [[ ! -f '.bed' ]]; then cut -f 4 $$(TARGETS_FILE_NOMSI) | paste -d '\t' $$(TARGETS_FILE_NOMSI) - > .bed; fi && \
					 $$(call WALTZ_CMD,2G,8G) org.mskcc.juber.waltz.Waltz PileupMetrics $$(WALTZ_MIN_MAPQ) $1_cl_aln_srt_MD_IR_FX2.bam $$(REF_FASTA) .bed && \
					 gzip $1_cl_aln_srt_MD_IR_FX2-pileup.txt && \
					 gzip $1_cl_aln_srt_MD_IR_FX2-pileup-without-duplicates.txt && \
					 gzip $1_cl_aln_srt_MD_IR_FX2-intervals.txt && \
					 gzip $1_cl_aln_srt_MD_IR_FX2-intervals-without-duplicates.txt && \
					 unlink $1_cl_aln_srt_MD_IR_FX2.bam && \
					 unlink $1_cl_aln_srt_MD_IR_FX2.bai && \
					 cd ..")
									 
waltz/$1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX-pileup.txt.gz : bam/$1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX.bam
	$$(call RUN,-c -n 4 -s 4G -m 6G,"set -o pipefail && \
					 mkdir -p waltz && \
					 cd waltz && \
					 ln -sf ../bam/$1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX.bam $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX.bam && \
					 ln -sf ../bam/$1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX.bai $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX.bai && \
					 if [[ ! -f '.bed' ]]; then cut -f 4 $$(TARGETS_FILE_NOMSI) | paste -d '\t' $$(TARGETS_FILE_NOMSI) - > .bed; fi && \
					 $$(call WALTZ_CMD,2G,8G) org.mskcc.juber.waltz.Waltz PileupMetrics $$(WALTZ_MIN_MAPQ) $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX.bam $$(REF_FASTA) .bed && \
					 gzip $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX-pileup.txt && \
					 gzip $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX-pileup-without-duplicates.txt && \
					 gzip $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX-intervals.txt && \
					 gzip $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX-intervals-without-duplicates.txt && \
					 unlink $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX.bam && \
					 unlink $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX.bai && \
					 cd ..")
									 
waltz/$1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_SIMPLEX-pileup.txt.gz : bam/$1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_SIMPLEX.bam
	$$(call RUN,-c -n 4 -s 4G -m 6G,"set -o pipefail && \
					 mkdir -p waltz && \
					 cd waltz && \
					 ln -sf ../bam/$1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_SIMPLEX.bam $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_SIMPLEX.bam && \
					 ln -sf ../bam/$1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_SIMPLEX.bai $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_SIMPLEX.bai && \
					 if [[ ! -f '.bed' ]]; then cut -f 4 $$(TARGETS_FILE_NOMSI) | paste -d '\t' $$(TARGETS_FILE_NOMSI) - > .bed; fi && \
					 $$(call WALTZ_CMD,2G,8G) org.mskcc.juber.waltz.Waltz PileupMetrics $$(WALTZ_MIN_MAPQ) $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_SIMPLEX.bam $$(REF_FASTA) .bed && \
					 gzip $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_SIMPLEX-pileup.txt && \
					 gzip $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_SIMPLEX-simplex-pileup-without-duplicates.txt && \
					 gzip $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_SIMPLEX-intervals.txt && \
					 gzip $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_SIMPLEX-intervals-without-duplicates.txt && \
					 unlink $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_SIMPLEX.bam && \
					 unlink $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_SIMPLEX.bai && \
					 cd ..")
									 
waltz/$1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_DUPLEX-pileup.txt.gz : bam/$1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_DUPLEX.bam
	$$(call RUN,-c -n 4 -s 4G -m 6G,"set -o pipefail && \
					 mkdir -p waltz && \
					 cd waltz && \
					 ln -sf ../bam/$1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_DUPLEX.bam $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_DUPLEX.bam && \
					 ln -sf ../bam/$1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_DUPLEX.bai $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_DUPLEX.bai && \
					 if [[ ! -f '.bed' ]]; then cut -f 4 $$(TARGETS_FILE_NOMSI) | paste -d '\t' $$(TARGETS_FILE_NOMSI) - > .bed; fi && \
					 $$(call WALTZ_CMD,2G,8G) org.mskcc.juber.waltz.Waltz PileupMetrics $$(WALTZ_MIN_MAPQ) $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_DUPLEX.bam $$(REF_FASTA) .bed && \
					 gzip $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_DUPLEX-pileup.txt && \
					 gzip $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_DUPLEX-pileup-without-duplicates.txt && \
					 gzip $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_DUPLEX-intervals.txt && \
					 gzip $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_DUPLEX-intervals-without-duplicates.txt && \
					 unlink $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_DUPLEX.bam && \
					 unlink $1_cl_aln_srt_MD_IR_FX2__grp_DC_MA_RG_IR_FX_DUPLEX.bai && \
					 cd ..")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call waltz-genotype,$(sample))))
		
..DUMMY := $(shell mkdir -p version; \
	     $(JAVA8) -version &> version/pileup_metrics.txt; \
	     R --version >> version/pileup_metrics.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: pileup_metrics
