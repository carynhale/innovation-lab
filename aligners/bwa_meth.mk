include innovation-lab/Makefile.inc
include innovation-lab/config/gatk.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/bwa_meth.$(NOW)

bwa_meth : $(foreach sample,$(SAMPLES),bwameth/$(sample)/$(sample)_R1.fastq.gz) \
	   $(foreach sample,$(SAMPLES),bwameth/$(sample)/$(sample)_R2.fastq.gz) \
	   $(foreach sample,$(SAMPLES),bwameth/$(sample)/$(sample)_aln.bam) \
	   $(foreach sample,$(SAMPLES),bwameth/$(sample)/$(sample)_aln_srt.bam) \
	   $(foreach sample,$(SAMPLES),bwameth/$(sample)/$(sample)_aln_srt_MD.bam) \
	   $(foreach sample,$(SAMPLES),bwameth/$(sample)/$(sample)_aln_srt_MD_RG.bam) \
	   $(foreach sample,$(SAMPLES),bwameth/$(sample)/$(sample)_aln_srt_MD_RG.intervals) \
	   $(foreach sample,$(SAMPLES),bwameth/$(sample)/$(sample)_aln_srt_MD_RG_IR.bam) \
	   $(foreach sample,$(SAMPLES),bwameth/$(sample)/$(sample)_aln_srt_MD_RG_IR_FX.bam)

BWAMETH_GENOME = $(REF_DIR)/IDT_oligo/bwa_meth/idt_oligo.fasta

BWAMETH_THREADS = 12
BWAMETH_MEM_PER_THREAD = 2G

SAMTOOLS_THREADS = 8
SAMTOOLS_MEM_THREAD = 2G

GATK_THREADS = 8
GATK_MEM_THREAD = 2G

define merge-fastq
bwameth/$1/$1_R1.fastq.gz : $$(foreach split,$2,$$(word 1, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 4G -m 6G,"zcat $$(^) | gzip -c > $$(@)")
	
bwameth/$1/$1_R2.fastq.gz : $$(foreach split,$2,$$(word 2, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 4G -m 6G,"zcat $$(^) | gzip -c > $$(@)")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call merge-fastq,$(sample),$(split.$(sample)))))

define fastq-2-bam
bwameth/$1/$1_aln.bam : bwameth/$1/$1_R1.fastq.gz bwameth/$1/$1_R2.fastq.gz
	$$(call RUN,-c -n $(BWAMETH_THREADS) -s 1G -m $(BWAMETH_MEM_PER_THREAD) -v $(BWAMETH_ENV),"set -o pipefail && \
												   $$(BWAMETH) \
												   --threads $$(BWAMETH_THREADS) \
												   --reference $$(BWAMETH_GENOME) \
												   $$(<) $$(<<) | \
												   $$(SAMTOOLS) view -bhS - > $$@")

bwameth/$1/$1_aln_srt.bam : bwameth/$1/$1_aln.bam
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
					  $$(SORT_SAM) \
					  INPUT=$$(<) \
					  OUTPUT=$$(@) \
					  SORT_ORDER=coordinate \
					  $$(SAMTOOLS) \
					  index \
					  $$(@) && \
					  cp bwameth/$1/$1_aln_srt.bam.bai bwameth/$1/$1_aln_srt.bai")
					 
bwameth/$1/$1_aln_srt_MD.bam : bwameth/$1/$1_aln_srt.bam
	$$(call RUN, -c -n 1 -s 12G -m 24G,"set -o pipefail && \
					    $$(MARK_DUP) \
					    INPUT=$$(<) \
					    OUTPUT=$$(@) \
					    METRICS_FILE=bwameth/$1/$1_aln_srt.txt \
					    REMOVE_DUPLICATES=false \
					    ASSUME_SORTED=true && \
					    $$(SAMTOOLS) index $$(@) && \
					    cp bwameth/$1/$1_aln_srt_MD.bam.bai bwameth/$1/$1_aln_srt_MD.bai")

bwameth/$1/$1_aln_srt_MD_RG.bam : bwameth/$1/$1_aln_srt_MD.bam
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
					  cp bwameth/$1/$1_aln_srt_MD_RG.bam.bai bwameth/$1/$1_aln_srt_MD_RG.bai")

bwameth/$1/$1_aln_srt_MD_RG.intervals : bwameth/$1/$1_aln_srt_MD_RG.bam
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD) -v $(GATK_ENV),"set -o pipefail && \
										      $$(call GATK_CMD,16G) \
										      -T RealignerTargetCreator \
										      -I $$(^) \
										      -nt $$(GATK_THREADS) \
										      -R $$(BWAMETH_GENOME) \
										      -o $$(@)")

bwameth/$1/$1_aln_srt_MD_RG_IR.bam : bwameth/$1/$1_aln_srt_MD_RG.bam bwameth/$1/$1_aln_srt_MD_RG.intervals
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD) -v $(GATK_ENV),"set -o pipefail && \
										      $$(call GATK_CMD,16G) \
										      -T IndelRealigner \
										      -I $$(<) \
										      -R $$(BWAMETH_GENOME) \
										      -targetIntervals $$(<<) \
										      -o $$(@)")

bwameth/$1/$1_aln_srt_MD_RG_IR_FX.bam : bwameth/$1/$1_aln_srt_MD_RG_IR.bam
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
					  $$(FIX_MATE) \
					  INPUT=$$(<) \
					  OUTPUT=$$(@) \
					  SORT_ORDER=coordinate \
					  COMPRESSION_LEVEL=0 \
					  CREATE_INDEX=true")

endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call fastq-2-bam,$(sample))))


..DUMMY := $(shell mkdir -p version; \
	     $(HOME)/share/usr/env/bwameth-0.2.4/bin/bwameth.py --version > version/bwa_meth.txt; \
	     $(SAMTOOLS) --version >> version/bwa_meth.txt; \
	     R --version >> version/bwa_meth.txt; \
	     $(JAVA8) -version &> version/bwa_meth.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: bwa_meth
