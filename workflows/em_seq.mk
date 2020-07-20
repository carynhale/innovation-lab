include innovation-lab/Makefile.inc
include innovation-lab/genome_inc/b37.inc
include innovation-lab/config/waltz.inc

LOGDIR ?= log/em_seq.$(NOW)

em_seq : $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_R1.fastq.gz) \
		 $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_aln.bam) \
		 $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_aln_srt.bam) \
         $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_aln_srt_fx.bam) \
		 $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_aln_srt_fx__F1R2.bam) \
		 $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_aln_srt_fx__F2R1.bam) \
		 $(foreach sample,$(SAMPLES),bam/$(sample)_aln_srt_fx.bam) \
		 $(foreach sample,$(SAMPLES),bam/$(sample)_aln_srt_fx__F1R2.bam) \
		 $(foreach sample,$(SAMPLES),bam/$(sample)_aln_srt_fx__F2R1.bam) \
		 $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_aln_srt_fx.rrbs_summary_metrics) \
		 $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_aln_srt_fx__F1R2.rrbs_summary_metrics) \
		 $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_aln_srt_fx__F2R1.rrbs_summary_metrics) \
		 $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_aln_srt_fx.txt) \
		 $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_aln_srt_fx__F1R2.txt) \
		 $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_aln_srt_fx__F2R1.txt) \
		 $(foreach sample,$(SAMPLES),waltz/$(sample)_aln_srt_fx-pileup.txt.gz) \
		 $(foreach sample,$(SAMPLES),waltz/$(sample)_aln_srt_fx__F1R2-pileup.txt.gz) \
		 $(foreach sample,$(SAMPLES),waltz/$(sample)_aln_srt_fx__F2R1-pileup.txt.gz) \
		 summary/rrbs_metrics.txt \
		 summary/alignment_metrics.txt \
		 summary/noise_by_position.txt \
		 summary/noise_by_position__F1R2.txt \
		 summary/noise_by_position__F2R1.txt

REF_FASTA = $(REF_DIR)/IDT_oligo/idt_oligo.fasta
GENOME_FOLDER = $(REF_DIR)/IDT_oligo/

SAMTOOLS_THREADS = 8
SAMTOOLS_MEM_THREAD = 2G

WALTZ_MIN_MAPQ = 10

define copy-fastq
bismark/$1/$1_R1.fastq.gz : $3
	$$(call RUN,-c -n 1 -s 2G -m 4G,"set -o pipefail && \
								     $(RSCRIPT) $(SCRIPTS_DIR)/fastq_tools/copy_fastq.R \
								     --sample_name $1 \
								     --directory_name bismark \
								     --fastq_files '$$^'")

endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)),$(eval $(call copy-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))
	
define fastq-to-bam
bismark/$1/$1_aln.bam : bismark/$1/$1_R1.fastq.gz
	$$(call RUN,-c -s 8G -m 16G -v $(BISMARK_ENV),"set -o pipefail && \
												   cd bismark/$1 && \
												   bismark --fastq \
												   --genome_folder $$(GENOME_FOLDER) \
												   -1 $1_R1.fastq.gz \
												   -2 $1_R2.fastq.gz \
												   --output_dir . && \
												   mv $1_R1_bismark_bt2_pe.bam $1_aln.bam && \
												   mv $1_R1_bismark_bt2_PE_report.txt $1_aln.txt && \
												   cd ../..")

bismark/$1/$1_aln_srt.bam : bismark/$1/$1_aln.bam
	$$(call RUN,-c -n $(SAMTOOLS_THREADS) -s 1G -m $(SAMTOOLS_MEM_THREAD),"set -o pipefail && \
									  									   $$(SAMTOOLS) sort -@ $$(SAMTOOLS_THREADS) -m $$(SAMTOOLS_MEM_THREAD) $$(^) -o $$(@) -T $$(TMPDIR) && \
									  									   $$(SAMTOOLS) index $$(@) && \
									  									   cp bismark/$1/$1_aln_srt.bam.bai bismark/$1/$1_aln_srt.bai")

bismark/$1/$1_aln_srt_fx.bam : bismark/$1/$1_aln_srt.bam
    $$(call RUN,-c -n 1 -s 12G -m 16G,"set -o pipefail && \
									   $$(FIX_MATE) \
									   INPUT=$$(<) \
									   OUTPUT=$$(@) \
									   SORT_ORDER=coordinate \
									   COMPRESSION_LEVEL=0 \
									   CREATE_INDEX=true")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call fastq-to-bam,$(sample))))

define filter-bam
bismark/$1/$1_aln_srt_fx__F1R2.bam : bismark/$1/$1_aln_srt_fx.bam
	$$(call RUN,-c -n $(SAMTOOLS_THREADS) -s 1G -m $(SAMTOOLS_MEM_THREAD),"set -o pipefail && \
																		   $$(SAMTOOLS) view -b -f 144 $$(<) > bismark/$1/$1_aln_srt__F1R2.bam && \
																		   $$(SAMTOOLS) index bismark/$1/$1_aln_srt__F1R2.bam && \
																		   $$(SAMTOOLS) view -b -f 64 -F 16 $$(<) > bismark/$1/$1_aln_srt__F1F2.bam && \
																		   $$(SAMTOOLS) index bismark/$1/$1_aln_srt__F1F2.bam && \
																		   $$(SAMTOOLS) merge -f $$(@) bismark/$1/$1_aln_srt__F1R2.bam bismark/$1/$1_aln_srt__F1F2.bam && \
																		   $$(SAMTOOLS) index $$(@)")

bismark/$1/$1_aln_srt_fx__F2R1.bam : bismark/$1/$1_aln_srt_fx.bam
	$$(call RUN,-c -n $(SAMTOOLS_THREADS) -s 1G -m $(SAMTOOLS_MEM_THREAD),"set -o pipefail && \
																		   $$(SAMTOOLS) view -b -f 128 -F 16 $$(<) > bismark/$1/$1_aln_srt__F2R1.bam && \
																		   $$(SAMTOOLS) index bismark/$1/$1_aln_srt__F2R1.bam && \
																		   $$(SAMTOOLS) view -b -f 80 $$(<) > bismark/$1/$1_aln_srt__F2F1.bam && \
																		   $$(SAMTOOLS) index bismark/$1/$1_aln_srt__F2F1.bam && \
																		   $$(SAMTOOLS) merge -f $$(@) bismark/$1/$1_aln_srt__F2R1.bam bismark/$1/$1_aln_srt__F2F1.bam && \
																		   $$(SAMTOOLS) index $$(@)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call filter-bam,$(sample))))
		
define copy-bam
bam/$1_aln_srt_fx.bam : bismark/$1/$1_aln_srt_fx.bam
	$$(call RUN,-c -s 1G -m 2G,"set -o pipefail && \
								cp $$(<) $$(@) && \
								$$(SAMTOOLS) index $$(@) && \
								cp bam/$1_aln_srt_fx.bam.bai bam/$1_aln_srt_fx.bai")

bam/$1_aln_srt_fx__F1R2.bam : bismark/$1/$1_aln_srt_fx__F1R2.bam
	$$(call RUN,-c -s 1G -m 2G,"set -o pipefail && \
								cp $$(<) $$(@) && \
								$$(SAMTOOLS) index $$(@) && \
								cp bam/$1_aln_srt_fx__F1R2.bam.bai bam/$1_aln_srt_fx__F1R2.bai")

bam/$1_aln_srt_fx__F2R1.bam : bismark/$1/$1_aln_srt_fx__F2R1.bam
	$$(call RUN,-c -s 1G -m 2G,"set -o pipefail && \
								cp $$(<) $$(@) && \
								$$(SAMTOOLS) index $$(@) && \
								cp bam/$1_aln_srt_fx__F2R1.bam.bai bam/$1_aln_srt_fx__F2R1.bai")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call copy-bam,$(sample))))
		
define picard-metrics
bismark/$1/$1_aln_srt_fx.rrbs_summary_metrics : bismark/$1/$1_aln_srt_fx.bam
	$$(call RUN,-c -s 12G -m 16G,"set -o pipefail && \
								  $$(COLLECT_RRBS_METRICS) \
								  R=$$(REF_FASTA) \
								  I=$$(<) \
								  M=bismark/$1/$1_aln_srt_fx")
								  
bismark/$1/$1_aln_srt_fx.txt : bismark/$1/$1_aln_srt_fx.bam
	$$(call RUN,-c -s 12G -m 16G,"set -o pipefail && \
								  $$(COLLECT_ALIGNMENT_METRICS) \
								  R=$$(REF_FASTA) \
								  I=$$(<) \
								  O=$$(@)")

bismark/$1/$1_aln_srt_fx__F1R2.rrbs_summary_metrics : bismark/$1/$1_aln_srt_fx__F1R2.bam
	$$(call RUN,-c -s 12G -m 16G,"set -o pipefail && \
								  $$(COLLECT_RRBS_METRICS) \
								  R=$$(REF_FASTA) \
								  I=$$(<) \
								  M=bismark/$1/$1_aln_srt_fx__F1R2")

bismark/$1/$1_aln_srt_fx__F1R2.txt : bismark/$1/$1_aln_srt_fx__F1R2.bam
	$$(call RUN,-c -s 12G -m 16G,"set -o pipefail && \
								  $$(COLLECT_ALIGNMENT_METRICS) \
								  R=$$(REF_FASTA) \
								  I=$$(<) \
								  O=$$(@)")
								  
bismark/$1/$1_aln_srt_fx__F2R1.rrbs_summary_metrics : bismark/$1/$1_aln_srt_fx__F2R1.bam
	$$(call RUN,-c -s 12G -m 16G,"set -o pipefail && \
								  $$(COLLECT_RRBS_METRICS) \
								  R=$$(REF_FASTA) \
								  I=$$(<) \
								  M=bismark/$1/$1_aln_srt_fx__F2R1")

bismark/$1/$1_aln_srt_fx__F2R1.txt : bismark/$1/$1_aln_srt_fx__F2R1.bam
	$$(call RUN,-c -s 12G -m 16G,"set -o pipefail && \
								  $$(COLLECT_ALIGNMENT_METRICS) \
								  R=$$(REF_FASTA) \
								  I=$$(<) \
								  O=$$(@)")
								  
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call picard-metrics,$(sample))))
		
define waltz-genotype
waltz/$1_aln_srt_fx-pileup.txt.gz : bismark/$1/$1_aln_srt_fx.bam
	$$(call RUN,-c -n 4 -s 4G -m 6G,"set -o pipefail && \
									 mkdir -p waltz && \
									 cd waltz && \
									 ln -sf ../bismark/$1/$1_aln_srt_fx.bam $1_aln_srt_fx.bam && \
									 ln -sf ../bismark/$1/$1_aln_srt_fx.bai $1_aln_srt_fx.bai && \
									 if [[ ! -f '.bed' ]]; then cut -f 4 $$(TARGETS_FILE) | paste -d '\t' $$(TARGETS_FILE) - > .bed; fi && \
									 $$(call WALTZ_CMD,2G,8G) org.mskcc.juber.waltz.Waltz PileupMetrics $$(WALTZ_MIN_MAPQ) $1_aln_srt_fx.bam $$(REF_FASTA) .bed && \
									 gzip $1_aln_srt_fx-pileup.txt && \
									 gzip $1_aln_srt_fx-pileup-without-duplicates.txt && \
									 gzip $1_aln_srt_fx-intervals.txt && \
									 gzip $1_aln_srt_fx-intervals-without-duplicates.txt && \
									 cd ..")

waltz/$1_aln_srt_fx__F1R2-pileup.txt.gz : bismark/$1/$1_aln_srt_fx__F1R2.bam
	$$(call RUN,-c -n 4 -s 4G -m 6G,"set -o pipefail && \
									 mkdir -p waltz && \
									 cd waltz && \
									 ln -sf ../bismark/$1/$1_aln_srt_fx__F1R2.bam $1_aln_srt_fx__F1R2.bam && \
									 ln -sf ../bismark/$1/$1_aln_srt_fx__F1R2.bam.bai $1_aln_srt_fx__F1R2.bai && \
									 if [[ ! -f '.bed' ]]; then cut -f 4 $$(TARGETS_FILE) | paste -d '\t' $$(TARGETS_FILE) - > .bed; fi && \
									 $$(call WALTZ_CMD,2G,8G) org.mskcc.juber.waltz.Waltz PileupMetrics $$(WALTZ_MIN_MAPQ) $1_aln_srt_fx__F1R2.bam $$(REF_FASTA) .bed && \
									 gzip $1_aln_srt_fx__F1R2-pileup.txt && \
									 gzip $1_aln_srt_fx__F1R2-pileup-without-duplicates.txt && \
									 gzip $1_aln_srt_fx__F1R2-intervals.txt && \
									 gzip $1_aln_srt_fx__F1R2-intervals-without-duplicates.txt && \
									 cd ..")

waltz/$1_aln_srt_fx__F2R1-pileup.txt.gz : bismark/$1/$1_aln_srt_fx__F2R1.bam
	$$(call RUN,-c -n 4 -s 4G -m 6G,"set -o pipefail && \
									 mkdir -p waltz && \
									 cd waltz && \
									 ln -sf ../bismark/$1/$1_aln_srt_fx__F2R1.bam $1_aln_srt_fx__F2R1.bam && \
									 ln -sf ../bismark/$1/$1_aln_srt_fx__F2R1.bam.bai $1_aln_srt_fx__F2R1.bai && \
									 if [[ ! -f '.bed' ]]; then cut -f 4 $$(TARGETS_FILE) | paste -d '\t' $$(TARGETS_FILE) - > .bed; fi && \
									 $$(call WALTZ_CMD,2G,8G) org.mskcc.juber.waltz.Waltz PileupMetrics $$(WALTZ_MIN_MAPQ) $1_aln_srt_fx__F2R1.bam $$(REF_FASTA) .bed && \
									 gzip $1_aln_srt_fx__F2R1-pileup.txt && \
									 gzip $1_aln_srt_fx__F2R1-pileup-without-duplicates.txt && \
									 gzip $1_aln_srt_fx__F2R1-intervals.txt && \
									 gzip $1_aln_srt_fx__F2R1-intervals-without-duplicates.txt && \
									 cd ..")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call waltz-genotype,$(sample))))
		
summary/rrbs_metrics.txt : $(foreach sample,$(SAMPLES),waltz/$(sample)_aln_srt_fx-pileup.txt.gz) $(foreach sample,$(SAMPLES),waltz/$(sample)_aln_srt_fx__F1R2-pileup.txt.gz) $(foreach sample,$(SAMPLES),waltz/$(sample)_aln_srt_fx__F2R1-pileup.txt.gz)
	$(call RUN, -c -n 1 -s 12G -m 16G,"set -o pipefail && \
									   mkdir -p summary && \
									   $(RSCRIPT) $(SCRIPTS_DIR)/qc/emseq_metrics.R --option 1 --sample_names '$(SAMPLES)'")

summary/alignment_metrics.txt : $(foreach sample,$(SAMPLES),waltz/$(sample)_aln_srt_fx-pileup.txt.gz) $(foreach sample,$(SAMPLES),waltz/$(sample)_aln_srt_fx__F1R2-pileup.txt.gz) $(foreach sample,$(SAMPLES),waltz/$(sample)_aln_srt_fx__F2R1-pileup.txt.gz)
	$(call RUN, -c -n 1 -s 12G -m 16G,"set -o pipefail && \
									   mkdir -p summary && \
									   $(RSCRIPT) $(SCRIPTS_DIR)/qc/emseq_metrics.R --option 2 --sample_names '$(SAMPLES)'")
									   
summary/noise_by_position.txt : $(foreach sample,$(SAMPLES),waltz/$(sample)_aln_srt_fx-pileup.txt.gz)
	$(call RUN, -c -n 1 -s 12G -m 16G,"set -o pipefail && \
									   mkdir -p summary && \
									   $(RSCRIPT) $(SCRIPTS_DIR)/qc/emseq_metrics.R --option 3 --sample_names '$(SAMPLES)'")

summary/noise_by_position__F1R2.txt : $(foreach sample,$(SAMPLES),waltz/$(sample)_aln_srt_fx__F1R2-pileup.txt.gz)
	$(call RUN, -c -n 1 -s 12G -m 16G,"set -o pipefail && \
									   mkdir -p summary && \
									   $(RSCRIPT) $(SCRIPTS_DIR)/qc/emseq_metrics.R --option 4 --sample_names '$(SAMPLES)'")

summary/noise_by_position__F2R1.txt : $(foreach sample,$(SAMPLES),waltz/$(sample)_aln_srt_fx__F2R1-pileup.txt.gz)
	$(call RUN, -c -n 1 -s 12G -m 16G,"set -o pipefail && \
									   mkdir -p summary && \
									   $(RSCRIPT) $(SCRIPTS_DIR)/qc/emseq_metrics.R --option 5 --sample_names '$(SAMPLES)'")

..DUMMY := $(shell mkdir -p version; \
			 $(HOME)/share/usr/env/bismark-0.22.1/bin/bismark --version > version/em_seq.txt; \
			 $(SAMTOOLS) --version >> version/em_seq.txt; \
			 R --version >> version/em_seq.txt; \
			 $(JAVA8) -version &> version/em_seq.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: em_seq
