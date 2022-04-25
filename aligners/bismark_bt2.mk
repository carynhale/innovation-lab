include innovation-lab/Makefile.inc
include innovation-lab/config/gatk.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/bismark_bt2.$(NOW)

bismark : $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_R1.fastq.gz) \
	  $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_R2.fastq.gz) \
	  $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_aln.bam) \
	  $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_aln_srt.bam) \
	  $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_aln_srt_MD.bam) \
	  $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_aln_srt_MD_RG.bam) \
	  $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_aln_srt_MD_RG.intervals) \
	  $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_aln_srt_MD_RG_IR.bam) \
	  $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_aln_srt_MD_RG_IR_FX.bam) \
#	  $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_aln_srt_RG_IR_FX__F1R2.bam) \
#	  $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_aln_srt_RG_IR_FX__F2R1.bam) \
#	  $(foreach sample,$(SAMPLES),bam/$(sample)_aln_srt_RG_IR_FX.bam) \
#	  $(foreach sample,$(SAMPLES),bam/$(sample)_aln_srt_RG_IR_FX__F1R2.bam) \
#	  $(foreach sample,$(SAMPLES),bam/$(sample)_aln_srt_RG_IR_FX__F2R1.bam) \
#	  $(foreach sample,$(SAMPLES),metrics/$(sample)_aln_srt_RG_IR_FX.rrbs_summary_metrics) \
#	  $(foreach sample,$(SAMPLES),metrics/$(sample)_aln_srt_RG_IR_FX__F1R2.rrbs_summary_metrics) \
#	  $(foreach sample,$(SAMPLES),metrics/$(sample)_aln_srt_RG_IR_FX__F2R1.rrbs_summary_metrics) \
#	  $(foreach sample,$(SAMPLES),metrics/$(sample)_aln_srt_RG_IR_FX.aln_metrics) \
#	  $(foreach sample,$(SAMPLES),metrics/$(sample)_aln_srt_RG_IR_FX__F1R2.aln_metrics) \
#	  $(foreach sample,$(SAMPLES),metrics/$(sample)_aln_srt_RG_IR_FX__F2R1.aln_metrics) \
#	  summary/rrbs_metrics.txt \
#	  summary/alignment_metrics.txt

SAMTOOLS_THREADS = 4
SAMTOOLS_MEM_THREAD = 4G

BISMARK_PARALLEL = 4
BISMARK_THREADS = 4
BISMARK_MEM_THREAD = 4G

GATK_THREADS = 8
GATK_MEM_THREAD = 2G

REF_FASTA = $(REF_DIR)/IDT_oligo/idt_oligo.fasta
BISMARK_GENOME = $(REF_DIR)/IDT_oligo/bismark_bt2

define merge-fastq
bismark/$1/$1_R1.fastq.gz : $$(foreach split,$2,$$(word 1, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 4G -m 6G,"set -o pipefail && \
					 zcat $$(^) | gzip -c > $$(@)")

bismark/$1/$1_R2.fastq.gz : $$(foreach split,$2,$$(word 2, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 4G -m 6G,"set -o pipefail && \
					 zcat $$(^) | gzip -c > $$(@)")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call merge-fastq,$(sample),$(split.$(sample)))))

define fastq-to-bam
bismark/$1/$1_aln.bam : bismark/$1/$1_R1.fastq.gz bismark/$1/$1_R2.fastq.gz
	$$(call RUN,-c -n $(BISMARK_THREADS) -s 2G -m $(BISMARK_MEM_THREAD) -v $(BISMARK_ENV),"set -o pipefail && \
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
											       cd ../..")

bismark/$1/$1_aln_srt.bam : bismark/$1/$1_aln.bam
	$$(call RUN,-c -n $(SAMTOOLS_THREADS) -s 2G -m $(SAMTOOLS_MEM_THREAD),"set -o pipefail && \
									       $$(SAMTOOLS) \
									       sort \
									       -@ $$(SAMTOOLS_THREADS) \
									       -m $$(SAMTOOLS_MEM_THREAD) \
									       $$(^) \
									       -o $$(@) \
									       -T $$(TMPDIR) && \
									       $$(SAMTOOLS) index $$(@) && \
									       cp bismark/$1/$1_aln_srt.bam.bai bismark/$1/$1_aln_srt.bai")
									       
bismark/$1/$1_aln_srt_MD.bam : bismark/$1/$1_aln_srt.bam
	$$(call RUN, -c -n 1 -s 8G -m 16G,"set -o pipefail && \
					   $$(MARK_DUP) \
					   INPUT=$$(<) \
					   OUTPUT=$$(@) \
					   METRICS_FILE=bismark/$1/$1_aln_srt.txt \
					   REMOVE_DUPLICATES=false \
					   ASSUME_SORTED=true && \
					   $$(SAMTOOLS) index $$(@) && \
					   cp bismark/$1/$1_aln_srt_MD.bam.bai bismark/$1/$1_aln_srt_MD.bai")
					   
bismark/$1/$1_aln_srt_MD_RG.bam : bismark/$1/$1_aln_srt_MD.bam
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
					  $$(ADD_RG) \
					  INPUT=$$(<) \
					  OUTPUT=$$(@) \
					  RGID=$1 \
					  RGLB=$1 \
					  RGPL=illumina \
					  RGPU=NA \
					  RGSM=$1 \
					  SORT_ORDER=coordinate \
					  COMPRESSION_LEVEL=0 && \
					  $$(SAMTOOLS) index $$(@) && \
					  cp bismark/$1/$1_aln_srt_MD_RG.bam.bai bismark/$1/$1_aln_srt_MD_RG.bai")

bismark/$1/$1_aln_srt_MD_RG.intervals : bismark/$1/$1_aln_srt_MD_RG.bam
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD) -v $(GATK_ENV),"set -o pipefail && \
										      $$(call GATK_CMD,16G) \
										      -T RealignerTargetCreator \
										      -I $$(^) \
										      -nt $$(GATK_THREADS) \
										      -R $$(REF_FASTA) \
										      -o $$(@)")

bismark/$1/$1_aln_srt_MD_RG_IR.bam : bismark/$1/$1_aln_srt_MD_RG.bam bismark/$1/$1_aln_srt_MD_RG.intervals
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD) -v $(GATK_ENV),"set -o pipefail && \
										      $$(call GATK_CMD,16G) \
										      -T IndelRealigner \
										      -I $$(<) \
										      -R $$(REF_FASTA) \
										      -targetIntervals $$(<<) \
										      -o $$(@)")

bismark/$1/$1_aln_srt_MD_RG_IR_FX.bam : bismark/$1/$1_aln_srt_MD_RG_IR.bam
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
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
bismark/$1/$1_aln_srt_RG_IR_FX__F1R2.bam : bismark/$1/$1_aln_srt_RG_IR_FX.bam
	$$(call RUN,-c -n $(SAMTOOLS_THREADS) -s 1G -m $(SAMTOOLS_MEM_THREAD),"set -o pipefail && \
									       $$(SAMTOOLS) view -b -f 144 $$(<) > bismark/$1/$1_aln_srt__F1R2.bam && \
									       $$(SAMTOOLS) index bismark/$1/$1_aln_srt__F1R2.bam && \
									       $$(SAMTOOLS) view -b -f 64 -F 16 $$(<) > bismark/$1/$1_aln_srt__F1F2.bam && \
									       $$(SAMTOOLS) index bismark/$1/$1_aln_srt__F1F2.bam && \
									       $$(SAMTOOLS) merge -f $$(@) bismark/$1/$1_aln_srt__F1R2.bam bismark/$1/$1_aln_srt__F1F2.bam && \
									       $$(SAMTOOLS) index $$(@)")

bismark/$1/$1_aln_srt_RG_IR_FX__F2R1.bam : bismark/$1/$1_aln_srt_RG_IR_FX.bam
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
bam/$1_aln_srt_RG_IR_FX.bam : bismark/$1/$1_aln_srt_RG_IR_FX.bam
	$$(call RUN,-c -s 1G -m 2G,"set -o pipefail && \
				    cp $$(<) $$(@) && \
				    $$(SAMTOOLS) index $$(@) && \
				    cp bam/$1_aln_srt_RG_IR_FX.bam.bai bam/$1_aln_srt_RG_IR_FX.bai")

bam/$1_aln_srt_RG_IR_FX__F1R2.bam : bismark/$1/$1_aln_srt_RG_IR_FX__F1R2.bam
	$$(call RUN,-c -s 1G -m 2G,"set -o pipefail && \
				    cp $$(<) $$(@) && \
				    $$(SAMTOOLS) index $$(@) && \
				    cp bam/$1_aln_srt_RG_IR_FX__F1R2.bam.bai bam/$1_aln_srt_RG_IR_FX__F1R2.bai")

bam/$1_aln_srt_RG_IR_FX__F2R1.bam : bismark/$1/$1_aln_srt_RG_IR_FX__F2R1.bam
	$$(call RUN,-c -s 1G -m 2G,"set -o pipefail && \
				    cp $$(<) $$(@) && \
				    $$(SAMTOOLS) index $$(@) && \
				    cp bam/$1_aln_srt_RG_IR_FX__F2R1.bam.bai bam/$1_aln_srt_RG_IR_FX__F2R1.bai")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call copy-bam,$(sample))))
		
define picard-metrics
metrics/$1_aln_srt_RG_IR_FX.rrbs_summary_metrics : bam/$1_aln_srt_RG_IR_FX.bam
	$$(call RUN,-c -s 12G -m 16G,"set -o pipefail && \
				      $$(COLLECT_RRBS_METRICS) \
				      R=$$(REF_FASTA) \
				      I=$$(<) \
				      M=metrics/$1_aln_srt_RG_IR_FX")
								  
metrics/$1_aln_srt_RG_IR_FX.aln_metrics : bam/$1_aln_srt_RG_IR_FX.bam
	$$(call RUN,-c -s 12G -m 16G,"set -o pipefail && \
				      $$(COLLECT_ALIGNMENT_METRICS) \
				      R=$$(REF_FASTA) \
				      I=$$(<) \
				      O=$$(@)")

metrics/$1_aln_srt_RG_IR_FX__F1R2.rrbs_summary_metrics : bam/$1_aln_srt_RG_IR_FX__F1R2.bam
	$$(call RUN,-c -s 12G -m 16G,"set -o pipefail && \
				      $$(COLLECT_RRBS_METRICS) \
				      R=$$(REF_FASTA) \
				      I=$$(<) \
				      M=metrics/$1_aln_srt_RG_IR_FX__F1R2")

metrics/$1_aln_srt_RG_IR_FX__F1R2.aln_metrics : bam/$1_aln_srt_RG_IR_FX__F1R2.bam
	$$(call RUN,-c -s 12G -m 16G,"set -o pipefail && \
				      $$(COLLECT_ALIGNMENT_METRICS) \
				      R=$$(REF_FASTA) \
				      I=$$(<) \
				      O=$$(@)")
								  
metrics/$1_aln_srt_RG_IR_FX__F2R1.rrbs_summary_metrics : bam/$1_aln_srt_RG_IR_FX__F2R1.bam
	$$(call RUN,-c -s 12G -m 16G,"set -o pipefail && \
				      $$(COLLECT_RRBS_METRICS) \
				      R=$$(REF_FASTA) \
				      I=$$(<) \
				      M=metrics/$1_aln_srt_RG_IR_FX__F2R1")

metrics/$1_aln_srt_RG_IR_FX__F2R1.aln_metrics : bam/$1_aln_srt_RG_IR_FX__F2R1.bam
	$$(call RUN,-c -s 12G -m 16G,"set -o pipefail && \
				      $$(COLLECT_ALIGNMENT_METRICS) \
				      R=$$(REF_FASTA) \
				      I=$$(<) \
				      O=$$(@)")
			  
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call picard-metrics,$(sample))))
		
summary/rrbs_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample)_aln_srt_RG_IR_FX.rrbs_summary_metrics) $(foreach sample,$(SAMPLES),metrics/$(sample)_aln_srt_RG_IR_FX__F1R2.rrbs_summary_metrics) $(foreach sample,$(SAMPLES),metrics/$(sample)_aln_srt_RG_IR_FX__F2R1.rrbs_summary_metrics)
	$(call RUN,-c -n 1 -s 12G -m 16G,"set -o pipefail && \
					  mkdir -p summary && \
					  $(RSCRIPT) $(SCRIPTS_DIR)/qc/bismark_metrics.R --option 1 --sample_names '$(SAMPLES)'")

summary/alignment_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample)_aln_srt_RG_IR_FX.aln_metrics) $(foreach sample,$(SAMPLES),metrics/$(sample)_aln_srt_RG_IR_FX__F1R2.aln_metrics) $(foreach sample,$(SAMPLES),metrics/$(sample)_aln_srt_RG_IR_FX__F2R1.aln_metrics)
	$(call RUN,-c -n 1 -s 12G -m 16G,"set -o pipefail && \
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
