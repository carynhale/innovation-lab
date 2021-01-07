include innovation-lab/Makefile.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/bismark_bt2.$(NOW)

bismark : $(foreach sample,$(SAMPLES),bam/$(sample)_aln_srt_MD_FX.bam) \
	  $(foreach sample,$(SAMPLES),bam/$(sample)_aln_srt_MD_FX__F1R2.bam) \
	  $(foreach sample,$(SAMPLES),bam/$(sample)_aln_srt_MD_FX__F2R1.bam) \
	  $(foreach sample,$(SAMPLES),metrics/$(sample)_aln_srt_MD_FX.rrbs_summary_metrics) \
	  $(foreach sample,$(SAMPLES),metrics/$(sample)_aln_srt_MD_FX__F1R2.rrbs_summary_metrics) \
	  $(foreach sample,$(SAMPLES),metrics/$(sample)_aln_srt_MD_FX__F2R1.rrbs_summary_metrics) \
	  $(foreach sample,$(SAMPLES),metrics/$(sample)_aln_srt_MD_FX.txt) \
	  $(foreach sample,$(SAMPLES),metrics/$(sample)_aln_srt_MD_FX__F1R2.txt) \
	  $(foreach sample,$(SAMPLES),metrics/$(sample)_aln_srt_MD_FX__F2R1.txt) \
	  summary/rrbs_metrics.txt \
	  summary/alignment_metrics.txt

SAMTOOLS_THREADS = 8
SAMTOOLS_MEM_THREAD = 5G

BISMARK_PARALLEL = 8
BISMARK_THREADS = 40
BISMARK_MEM_THREAD = 4G

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
	$$(call RUN,-c -n $(BISMARK_THREADS) -s 1G -m $(BISMARK_MEM_THREAD) -v $(BISMARK_ENV) -w 72:00:00,"set -o pipefail && \
								    					   cd bismark/$1 && \
								    					   bismark \
								    					   --fastq \
								    					   --parallel $$(BISMARK_PARALLEL) \
								    					   --genome_folder $$(BISMARK_GENOME) \
								    					   -1 $1_R1.fastq.gz \
								    					   -2 $1_R2.fastq.gz \
								    					   --output_dir . && \
								    					   mv $1_R1_bismark_bt2_pe.bam $1_aln.bam && \
								    					   mv $1_R1_bismark_bt2_PE_report.txt $1_aln.txt && \
								    					   cd ../.. && \
													   $$(RM) bismark/$1/$1_R1.fastq.gz && \
													   $$(RM) bismark/$1/$1_R2.fastq.gz")

bismark/$1/$1_aln_srt.bam : bismark/$1/$1_aln.bam
	$$(call RUN,-c -n $(SAMTOOLS_THREADS) -s 1G -m $(SAMTOOLS_MEM_THREAD) -w 12:00:00,"set -o pipefail && \
											   $$(SAMTOOLS) sort -@ $$(SAMTOOLS_THREADS) -m $$(SAMTOOLS_MEM_THREAD) $$(^) -o $$(@) -T $$(TMPDIR) && \
											   $$(SAMTOOLS) index $$(@) && \
											   cp bismark/$1/$1_aln_srt.bam.bai bismark/$1/$1_aln_srt.bai && \
											   $$(RM) $$(^)")

bismark/$1/$1_aln_srt_MD.bam : bismark/$1/$1_aln_srt.bam
	$$(call RUN, -c -n 1 -s 36G -m 48G -w 12:00:00,"set -o pipefail && \
							$$(MARK_DUP) \
						       	INPUT=$$(<) \
						       	OUTPUT=$$(@) \
						       	METRICS_FILE=bismark/$1/$1_cl_aln_srt.txt \
						       	REMOVE_DUPLICATES=false \
						       	ASSUME_SORTED=true && \
						       	$$(SAMTOOLS) index $$(@) && \
						       	cp bismark/$1/$1_aln_srt_MD.bam.bai bismark/$1/$1_aln_srt_MD.bai && \
						       	$$(RM)  $$(^)")

bismark/$1/$1_aln_srt_MD_FX.bam : bismark/$1/$1_aln_srt_MD.bam
	$$(call RUN,-c -n 1 -s 36G -m 48G -w 24:00:00,"set -o pipefail && \
						      	$$(FIX_MATE) \
						      	INPUT=$$(<) \
						      	OUTPUT=$$(@) \
						      	SORT_ORDER=coordinate \
						      	COMPRESSION_LEVEL=0 \
						      	CREATE_INDEX=true && \
						      	$$(RM) $$(^)")
						      
endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call fastq-to-bam,$(sample))))

define filter-bam
bismark/$1/$1_aln_srt_MD_FX__F1.bam : bismark/$1/$1_aln_srt_MD_FX.bam
	$$(call RUN,-c -s 8G -m 12G -w 12:00:00,"set -o pipefail && \
						 $$(SAMTOOLS) view -b -f 144 $$(<) > $$(@) && \
						 $$(SAMTOOLS) index $$(@) &&")

bismark/$1/$1_aln_srt_MD_FX__R2.bam : bismark/$1/$1_aln_srt_MD_FX.bam
	$$(call RUN,-c -s 8G -m 12G -w 12:00:00,"set -o pipefail && \
						 $$(SAMTOOLS) view -b -f 64 -F 16 $$(<) > $$(@) && \
						 $$(SAMTOOLS) index $$(@)")

bismark/$1/$1_aln_srt_MD_FX__F1R2.bam : bismark/$1/$1_aln_srt_MD_FX__F1.bam bismark/$1/$1_aln_srt_MD_FX__R2.bam
	$$(call RUN,-c -s 8G -m 12G -w 12:00:00,"set -o pipefail && \
						 $$(SAMTOOLS) merge -f $$(@) $$(<) $$(<<) && \
						 $$(SAMTOOLS) index $$(@) && \
						 $$(RM) $$(<) && \
						 $$(RM) $$(<<)")
											   
bismark/$1/$1_aln_srt_MD_FX__F2.bam : bismark/$1/$1_aln_srt_MD_FX.bam
	$$(call RUN,-c -s 8G -m 12G -w 12:00:00,"set -o pipefail && \
						 $$(SAMTOOLS) view -b -f 128 -F 16 $$(<) > $$(@) && \
						 $$(SAMTOOLS) index $$(@)")

bismark/$1/$1_aln_srt_MD_FX__R1.bam : bismark/$1/$1_aln_srt_MD_FX.bam
	$$(call RUN,-c -s 8G -m 12G -w 12:00:00,"set -o pipefail && \
						 $$(SAMTOOLS) view -b -f 80 $$(<) > $$(@) && \
						 $$(SAMTOOLS) index $$(@)")

bismark/$1/$1_aln_srt_MD_FX__F2R1.bam : bismark/$1/$1_aln_srt_MD_FX__F2.bam bismark/$1/$1_aln_srt_MD_FX__R1.bam
	$$(call RUN,-c -s 8G -m 12G -w 12:00:00,"set -o pipefail && \
						 $$(SAMTOOLS) merge -f $$(@) $$(<) $$(<<) && \
						 $$(SAMTOOLS) index $$(@) && \
						 $$(RM) $$(<) && \
						 $$(RM) $$(<<)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call filter-bam,$(sample))))

define copy-bam
bam/$1_aln_srt_MD_FX.bam : bismark/$1/$1_aln_srt_MD_FX.bam
	$$(call RUN,-c -s 2G -m 4G,"set -o pipefail && \
				    cp $$(<) $$(@) && \
				    $$(SAMTOOLS) index $$(@) && \
				    cp bam/$1_aln_srt_MD_FX.bam.bai bam/$1_aln_srt_MD_FX.bai && \
				    $$(RM) $$(^)")

bam/$1_aln_srt_MD_FX__F1R2.bam : bismark/$1/$1_aln_srt_MD_FX__F1R2.bam
	$$(call RUN,-c -s 2G -m 4G,"set -o pipefail && \
				    cp $$(<) $$(@) && \
				    $$(SAMTOOLS) index $$(@) && \
				    cp bam/$1_aln_srt_MD_FX__F1R2.bam.bai bam/$1_aln_srt_MD_FX__F1R2.bai && \
				    $$(RM) $$(^)")

bam/$1_aln_srt_MD_FX__F2R1.bam : bismark/$1/$1_aln_srt_MD_FX__F2R1.bam
	$$(call RUN,-c -s 2G -m 4G,"set -o pipefail && \
				    cp $$(<) $$(@) && \
				    $$(SAMTOOLS) index $$(@) && \
				    cp bam/$1_aln_srt_MD_FX__F2R1.bam.bai bam/$1_aln_srt_MD_FX__F2R1.bai && \
				    $$(RM) $$(^)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call copy-bam,$(sample))))
		
define picard-metrics
metrics/$1_aln_srt_MD_FX.rrbs_summary_metrics : bam/$1_aln_srt_MD_FX.bam
	$$(call RUN,-c -s 12G -m 24G,"set -o pipefail && \
				      $$(COLLECT_RRBS_METRICS) \
				      R=$$(REF_FASTA) \
				      I=$$(<) \
				      M=metrics/$1_aln_srt_MD_FX")
								  
metrics/$1_aln_srt_MD_FX.txt : bam/$1_aln_srt_MD_FX.bam
	$$(call RUN,-c -s 12G -m 24G,"set -o pipefail && \
				      $$(COLLECT_ALIGNMENT_METRICS) \
				      R=$$(REF_FASTA) \
				      I=$$(<) \
				      O=$$(@)")

metrics/$1_aln_srt_MD_FX__F1R2.rrbs_summary_metrics : bam/$1_aln_srt_MD_FX__F1R2.bam
	$$(call RUN,-c -s 12G -m 24G,"set -o pipefail && \
				      $$(COLLECT_RRBS_METRICS) \
				      R=$$(REF_FASTA) \
				      I=$$(<) \
				      M=metrics/$1_aln_srt_MD_FX__F1R2")

metrics/$1_aln_srt_MD_FX__F1R2.txt : bam/$1_aln_srt_MD_FX__F1R2.bam
	$$(call RUN,-c -s 12G -m 24G,"set -o pipefail && \
				      $$(COLLECT_ALIGNMENT_METRICS) \
				      R=$$(REF_FASTA) \
				      I=$$(<) \
				      O=$$(@)")
								  
metrics/$1_aln_srt_MD_FX__F2R1.rrbs_summary_metrics : bam/$1_aln_srt_MD_FX__F2R1.bam
	$$(call RUN,-c -s 12G -m 24G,"set -o pipefail && \
				      $$(COLLECT_RRBS_METRICS) \
				      R=$$(REF_FASTA) \
				      I=$$(<) \
				      M=metrics/$1_aln_srt_MD_FX__F2R1")

metrics/$1_aln_srt_MD_FX__F2R1.txt : bam/$1_aln_srt_MD_FX__F2R1.bam
	$$(call RUN,-c -s 12G -m 24G,"set -o pipefail && \
				      $$(COLLECT_ALIGNMENT_METRICS) \
				      R=$$(REF_FASTA) \
				      I=$$(<) \
				      O=$$(@)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call picard-metrics,$(sample))))
		
summary/rrbs_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample)_aln_srt_MD_FX.rrbs_summary_metrics) $(foreach sample,$(SAMPLES),metrics/$(sample)_aln_srt_MD_FX__F1R2.rrbs_summary_metrics) $(foreach sample,$(SAMPLES),metrics/$(sample)_aln_srt_MD_FX__F2R1.rrbs_summary_metrics)
	$(call RUN, -c -n 1 -s 12G -m 16G,"set -o pipefail && \
					   mkdir -p summary && \
					   $(RSCRIPT) $(SCRIPTS_DIR)/qc/bismark_metrics.R --option 1 --sample_names '$(SAMPLES)'")

summary/alignment_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample)_aln_srt_MD_FX.txt) $(foreach sample,$(SAMPLES),metrics/$(sample)_aln_srt_MD_FX__F1R2.txt) $(foreach sample,$(SAMPLES),metrics/$(sample)_aln_srt_MD_FX__F2R1.txt)
	$(call RUN, -c -n 1 -s 12G -m 16G,"set -o pipefail && \
					   mkdir -p summary && \
					   $(RSCRIPT) $(SCRIPTS_DIR)/qc/bismark_metrics.R --option 2 --sample_names '$(SAMPLES)'")
									   
..DUMMY := $(shell mkdir -p version; \
	     $(HOME)/share/usr/env/bismark-0.22.1/bin/bismark --version > version/bismark_bt2.txt; \
	     $(SAMTOOLS) --version >> version/bismark_bt2.txt; \
	     R --version >> version/bismark_bt2.txt; \
	     $(JAVA8) -version &> version/bismark_bt2.txt)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: bismark
