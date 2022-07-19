include innovation-lab/Makefile.inc
include innovation-lab/config/gatk.inc
include innovation-lab/config/align.inc

LOGDIR ?= log/bwa_split.$(NOW)

FASTQ_SPLIT = 500
FASTQ_SEQ = $(shell seq 1 $(FASTQ_SPLIT))

bwa_split : $(foreach sample,$(SAMPLES),bwamem/$(sample)/$(sample)_R1.fastq.gz) \
	    $(foreach sample,$(SAMPLES),bwamem/$(sample)/$(sample)_R2.fastq.gz) \
	    $(foreach sample,$(SAMPLES), \
		  	$(foreach n,$(FASTQ_SEQ),bwamem/$(sample)/$(sample)--$(n)_R1.fastq.gz)) \
	    $(foreach sample,$(SAMPLES), \
		  	$(foreach n,$(FASTQ_SEQ),bwamem/$(sample)/$(sample)--$(n)_R2.fastq.gz))
#	    $(foreach sample,$(SAMPLES),bwamem/$(sample)/$(sample)--$(FASTQ_SPLIT)_R1.fastq.gz) \
#	    $(foreach sample,$(SAMPLES),bwamem/$(sample)/$(sample)--$(FASTQ_SPLIT)_R2.fastq.gz) \
#	    $(foreach sample,$(SAMPLES), \
#		  	$(foreach n,$(FASTQ_SEQ),bwamem/$(sample)/$(sample)--$(n)_aln.bam)) \
#	    $(foreach sample,$(SAMPLES), \
#		  	$(foreach n,$(FASTQ_SEQ),bwamem/$(sample)/$(sample)--$(n)_cl.fastq.gz)) \
#	    $(foreach sample,$(SAMPLES), \
#		  	$(foreach n,$(FASTQ_SEQ),bwamem/$(sample)/$(sample)--$(n)_cl_aln.bam)) \
#	    $(foreach sample,$(SAMPLES), \
#		  	$(foreach n,$(FASTQ_SEQ),bwamem/$(sample)/$(sample)--$(n)_cl_aln_srt.bam)) \
#	    $(foreach sample,$(SAMPLES), \
#		  	$(foreach n,$(FASTQ_SEQ),bwamem/$(sample)/$(sample)--$(n)_cl_aln_srt.intervals)) \
#	    $(foreach sample,$(SAMPLES), \
#		  	$(foreach n,$(FASTQ_SEQ),bwamem/$(sample)/$(sample)--$(n)_cl_aln_srt_IR.bam)) \
#	    $(foreach sample,$(SAMPLES), \
#		  	$(foreach n,$(FASTQ_SEQ),bwamem/$(sample)/$(sample)--$(n)_cl_aln_srt_IR_FX.bam)) \
#	    $(foreach sample,$(SAMPLES), \
#		  	$(foreach n,$(FASTQ_SEQ),bwamem/$(sample)/$(sample)--$(n)_cl_aln_srt_IR_FX.grp)) \
#	    $(foreach sample,$(SAMPLES), \
#		  	$(foreach n,$(FASTQ_SEQ),bwamem/$(sample)/$(sample)--$(n)_cl_aln_srt_IR_FX_BR.bam)) \
#	    $(foreach sample,$(SAMPLES),bwamem/$(sample)/$(sample)_cl_aln_srt_IR_FX_BR.bam) \
#	    $(foreach sample,$(SAMPLES),bwamem/$(sample)/$(sample)_cl_aln_srt_IR_FX_BR_MD.bam) \
#	    $(foreach sample,$(SAMPLES),bam/$(sample).bam)

SPLIT_THREADS = 8
SPLIT_MEM_THREAD = 1G

BWAMEM_THREADS = 4
BWAMEM_MEM_PER_THREAD = 1G

SAMTOOLS_THREADS = 4
SAMTOOLS_MEM_THREAD = 1G

GATK_THREADS = 4
GATK_MEM_THREAD = 2G

define merge-fastq
bwamem/$1/$1_R1.fastq.gz : $$(foreach split,$2,$$(word 1, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 4G -m 6G -w 24:00:00,"zcat $$(^) | gzip -c > $$(@)")
	
bwamem/$1/$1_R2.fastq.gz : $$(foreach split,$2,$$(word 2, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 4G -m 6G -w 24:00:00,"zcat $$(^) | gzip -c > $$(@)")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call merge-fastq,$(sample),$(split.$(sample)))))
		
define split-fastq
$(foreach n,$(FASTQ_SEQ),bwamem/$1/$1--$(n)_R1.fastq.gz) &: bwamem/$1/$1_R1.fastq.gz
	$$(call RUN,-c -n $(SPLIT_THREADS) -s 1G -m $(SPLIT_MEM_THREAD) -v $(FASTQ_SPLITTER_ENV),"set -o pipefail && \
												  $(SCRIPTS_DIR)/fastq_tools/split_fastq.sh \
												  $$(FASTQ_SPLIT) \
								   				  $$(<) \
												  bwamem/$1/$1 \
												  R1 \
												  -t $$(SPLIT_THREADS)")

$(foreach n,$(FASTQ_SEQ),bwamem/$1/$1--$(n)_R2.fastq.gz) &: bwamem/$1/$1_R2.fastq.gz
	$$(call RUN,-c -n $(SPLIT_THREADS) -s 1G -m $(SPLIT_MEM_THREAD) -v $(FASTQ_SPLITTER_ENV),"set -o pipefail && \
								   				  $(SCRIPTS_DIR)/fastq_tools/split_fastq.sh \
												  $$(FASTQ_SPLIT) \
												  $$(<) \
												  bwamem/$1/$1 \
												  R2 \
												  -t $$(SPLIT_THREADS)")
								   
endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call split-fastq,$(sample))))

define fastq-2-bam
bwamem/$1/$1--$2_aln.bam : bwamem/$1/$1--$(FASTQ_SPLIT)_R1.fastq.gz bwamem/$1/$1--$(FASTQ_SPLIT)_R2.fastq.gz
	$$(call RUN,-c -n 1 -s 2G -m 4G,"set -o pipefail && \
					 $$(FASTQ_TO_SAM) \
					 FASTQ=bwamem/$1/$1--$2_R1.fastq.gz \
					 FASTQ2=bwamem/$1/$1--$2_R2.fastq.gz \
					 OUTPUT=$$(@) \
					 SM=$1 \
					 LB=$1 \
					 PU=NA \
					 PL=illumina")
									       
bwamem/$1/$1--$2_cl.fastq.gz : bwamem/$1/$1--$2_aln.bam
	$$(call RUN,-c -n 1 -s 1G -m 2G,"set -o pipefail && \
					 $$(MARK_ADAPTERS) \
					 INPUT=$$(<) \
					 OUTPUT=/dev/stdout \
					 METRICS=bwamem/$1/$1--$2_adapter-metrics.txt | \
					 $$(SAM_TO_FASTQ) \
					 INPUT=/dev/stdin \
					 FASTQ=$$(@) \
					 INTERLEAVE=true \
					 CLIPPING_ATTRIBUTE=XT \
					 CLIPPING_ACTION=X \
					 CLIPPING_MIN_LENGTH=25")
									       
bwamem/$1/$1--$2_cl_aln.bam : bwamem/$1/$1--$2_cl.fastq.gz
	$$(call RUN,-c -n $(BWAMEM_THREADS) -s 1G -m $(BWAMEM_MEM_PER_THREAD),"set -o pipefail && \
									       $$(BWA) mem -p -M \
									       -R \"@RG\tID:$1\tLB:$1\tPL:illumina\tSM:$1\" \
									       -t $$(BWAMEM_THREADS) $$(REF_FASTA) $$(<) | \
									       $$(SAMTOOLS) view -bhS - > $$(@)")

bwamem/$1/$1--$2_cl_aln_srt.bam : bwamem/$1/$1--$2_cl_aln.bam
	$$(call RUN,-c -n $(SAMTOOLS_THREADS) -s 1G -m $(SAMTOOLS_MEM_THREAD),"set -o pipefail && \
									       $$(SAMTOOLS) sort $$(<) -o $$(@) && \
									       $$(SAMTOOLS) index $$(@) && \
									       cp bwamem/$1/$1--$2_cl_aln_srt.bam.bai bwamem/$1/$1--$2_cl_aln_srt.bai")

bwamem/$1/$1--$2_cl_aln_srt.intervals : bwamem/$1/$1--$2_cl_aln_srt.bam
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD) -v $(GATK_ENV),"set -o pipefail && \
										      $$(call GATK_CMD,8G) \
										      -T RealignerTargetCreator \
										      -I $$(^) \
										      -nt $$(GATK_THREADS) \
										      -R $$(REF_FASTA) \
										      -o $$(@) \
										      -known $$(KNOWN_INDELS)")
										      
bwamem/$1/$1--$2_cl_aln_srt_IR.bam : bwamem/$1/$1--$2_cl_aln_srt.bam bwamem/$1/$1--$2_cl_aln_srt.intervals
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD) -v $(GATK_ENV),"set -o pipefail && \
										      $$(call GATK_CMD,8G) \
										      -T IndelRealigner \
										      -I $$(<) \
										      -R $$(REF_FASTA) \
										      -targetIntervals $$(<<) \
										      -o $$(@) \
										      -known $$(KNOWN_INDELS)")
										      
bwamem/$1/$1--$2_cl_aln_srt_IR_FX.bam : bwamem/$1/$1--$2_cl_aln_srt_IR.bam
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
					  $$(FIX_MATE) \
					  INPUT=$$(<) \
					  OUTPUT=$$(@) \
					  SORT_ORDER=coordinate \
					  COMPRESSION_LEVEL=0 \
					  CREATE_INDEX=true")
										      
bwamem/$1/$1--$2_cl_aln_srt_IR_FX.grp : bwamem/$1/$1--$2_cl_aln_srt_IR_FX.bam
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD) -v $(GATK_ENV),"set -o pipefail && \
										      $$(call GATK_CMD,8G) \
										      -T BaseRecalibrator \
										      -R $$(REF_FASTA) \
										      -knownSites $$(DBSNP) \
										      -I $$(<) \
										      -o $$(@)")

bwamem/$1/$1--$2_cl_aln_srt_IR_FX_BR.bam : bwamem/$1/$1--$2_cl_aln_srt_IR_FX.bam bwamem/$1/$1--$2_cl_aln_srt_IR_FX.grp
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD) -v $(GATK_ENV),"set -o pipefail && \
										      $$(call GATK_CMD,8G) \
										      -T PrintReads \
										      -R $$(REF_FASTA) \
										      -I $$(<) \
										      -BQSR $$(<<) \
										      -o $$(@)")

endef
$(foreach sample,$(SAMPLES), \
	$(foreach n,$(FASTQ_SEQ), \
		$(eval $(call fastq-2-bam,$(sample),$(n)))))
		
define collect-dedup
bwamem/$1/$1_cl_aln_srt_IR_FX_BR.bam : $(foreach n,$(FASTQ_SEQ),bwamem/$1/$1--$(n)_cl_aln_srt_IR_FX_BR.bam)
	$$(call RUN,-c -n 4 -s 2G -m 4G -w 72:00:00,"set -o pipefail && \
						     $$(SAMTOOLS) \
						     merge \
						     -c -p \
						     --threads 4 \
						     $$(@) \
						     $$(^) && \
						     $$(SAMTOOLS) index $$(@) && \
						     cp bwamem/$1/$1_cl_aln_srt_IR_FX_BR.bam.bai bwamem/$1/$1_cl_aln_srt_IR_FX_BR.bai")

bwamem/$1/$1_cl_aln_srt_IR_FX_BR_MD.bam : bwamem/$1/$1_cl_aln_srt_IR_FX_BR.bam
	$$(call RUN, -c -n 1 -s 24G -m 36G -w 72:00:00,"set -o pipefail && \
							$$(MARK_DUP) \
							INPUT=$$(<) \
							OUTPUT=$$(@) \
							METRICS_FILE=bwamem/$1/$1_cl_aln_srt_IR_FX_BR_MD.txt \
							REMOVE_DUPLICATES=false \
							ASSUME_SORTED=true && \
							$$(SAMTOOLS) index $$(@) && \
							cp bwamem/$1/$1_cl_aln_srt_IR_FX_BR_MD.bam.bai bwamem/$1/$1_cl_aln_srt_IR_FX_BR_MD.bai")
							
bam/$1.bam : bwamem/$1/$1_cl_aln_srt_IR_FX_BR_MD.bam
	$$(call RUN, -c -n 1 -s 1G -m 2G,"set -o pipefail && \
					  cp $$(<) $$(@) && \
					  cp $$(<).bai $$(@).bai && \
					  cp $$(<).bai bam/$1.bai")

endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call collect-dedup,$(sample))))

		
..DUMMY := $(shell mkdir -p version; \
	     $(BWA) &> version/tmp.txt; \
	     head -3 version/tmp.txt | tail -2 > version/bwa_split.txt; \
	     rm version/tmp.txt; \
	     $(SAMTOOLS) --version >> version/bwa_split.txt; \
	     echo "gatk3" >> version/bwa_split.txt; \
	     $(GATK) --version >> version/bwa_split.txt; \
	     echo "picard" >> version/bwa_split.txt; \
	     $(PICARD) MarkDuplicates --version &>> version/bwa_split.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: bwa_split
