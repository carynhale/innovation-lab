include innovation-lab/Makefile.inc
include innovation-lab/config/gatk.inc
include innovation-lab/config/align.inc

LOGDIR ?= log/bwa_split.$(NOW)

FASTQ_SPLIT = 500
FASTQ_SEQ = $(shell seq 1 $(FASTQ_SPLIT))

bwa_split : $(foreach sample,$(SAMPLES),bwamem/$(sample)/$(sample)_R1.fastq.gz) \
	    $(foreach sample,$(SAMPLES),bwamem/$(sample)/$(sample)_R2.fastq.gz) \
	    $(foreach sample,$(SAMPLES),bwamem/$(sample)/taskcomplete_r1.txt) \
	    $(foreach sample,$(SAMPLES),bwamem/$(sample)/taskcomplete_r2.txt) \
	    $(foreach sample,$(SAMPLES), \
		  	$(foreach n,$(FASTQ_SEQ),bwamem/$(sample)/$(sample)_aln--$(n).bam))
	    
BWAMEM_THREADS = 12
BWAMEM_MEM_PER_THREAD = 2G
BWA_ALN_OPTS ?= -M

SAMTOOLS_THREADS = 8
SAMTOOLS_MEM_THREAD = 2G

GATK_THREADS = 8
GATK_MEM_THREAD = 2G

define merge-fastq
bwamem/$1/$1_R1.fastq.gz : $$(foreach split,$2,$$(word 1, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 4G -m 6G -w 72:00:00,"zcat $$(^) | gzip -c > $$(@)")
	
bwamem/$1/$1_R2.fastq.gz : $$(foreach split,$2,$$(word 2, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 4G -m 6G -w 72:00:00,"zcat $$(^) | gzip -c > $$(@)")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call merge-fastq,$(sample),$(split.$(sample)))))
		

define split-fastq
bwamem/$1/taskcomplete_r1.txt : bwamem/$1/$1_R1.fastq.gz
	$$(call RUN,-c -n 12 -s 1G -m 2G -v $(FASTQ_SPLITTER_ENV),"set -o pipefail && \
								   $(SCRIPTS_DIR)/fastq_tools/split_fastq.sh \
								   $$(FASTQ_SPLIT) \
								   $$(<) \
								   bwamem/$1/ \
								   R1 \
								   -t 12 && \
								   echo 'taskcomplete!' > $$(@)")

bwamem/$1/taskcomplete_r2.txt : bwamem/$1/$1_R2.fastq.gz
	$$(call RUN,-c -n 12 -s 1G -m 2G -v $(FASTQ_SPLITTER_ENV),"set -o pipefail && \
								   $(SCRIPTS_DIR)/fastq_tools/split_fastq.sh \
								   $$(FASTQ_SPLIT) \
								   $$(<) \
								   bwamem/$1/ \
								   R2 \
								   -t 12 && \
								   echo 'taskcomplete!' > $$(@)")
								   
endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call split-fastq,$(sample))))

define fastq-2bam
bwamem/$1/$1_aln--$2.bam : bwamem/$1/taskcomplete_r1.txt bwamem/$1/taskcomplete_r2.txt
	$$(call RUN,-c -n $(BWAMEM_THREADS) -s 1G -m $(BWAMEM_MEM_PER_THREAD),"set -o pipefail && \
									       $$(BWA) mem -p -t $$(BWAMEM_THREADS) $$(BWA_ALN_OPTS) $$(REF_FASTA) \
									       bwamem/$1/$2_R1.fastq.gz bwamem/$1/$2_R2.fastq.gz | \
									       $$(SAMTOOLS) view -bhS - > $$(@)")
endef
$(foreach sample,$(SAMPLES), \
	$(foreach n,$(FASTQ_SEQ), \
		$(eval $(call fastq-2bam,$(sample),$(n)))))

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
