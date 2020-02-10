include innovation-lab/Makefile.inc
include innovation-lab/config/gatk.inc
include innovation-lab/config/align.inc

LOGDIR ?= log/bwa_mem.$(NOW)

ALIGNER := bwamem
SAMTOOLS_SORT_MEM = 2000000000
SEQ_PLATFORM = illumina
FASTQ_CHUNKS := 10
FASTQ_CHUNK_SEQ := $(shell seq 1 $(FASTQ_CHUNKS))
BWA_ALN_OPTS ?= -M
BWAMEM_REF_FASTA ?= $(REF_FASTA)
BWAMEM_THREADS = 8
BWAMEM_MEM_PER_THREAD = $(if $(findstring true,$(PDX)),4G,2G)
BWA_BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)

bwa_mem : $(BWA_BAMS) $(addsuffix .bai,$(BWA_BAMS))

bam/%.bam : bwamem/bam/%.bwamem.$(BAM_SUFFIX)
	$(call RUN,,"ln -f $(<) $(@)")

define align-split-fastq
bwamem/bam/$2.bwamem.bam : $3
	$$(call RUN,-n $$(BWAMEM_THREADS) -s 1G -m $$(BWAMEM_MEM_PER_THREAD),"$$(BWA) mem -t $$(BWAMEM_THREADS) $$(BWA_ALN_OPTS) -R \"@RG\tID:$2\tLB:$1\tPL:$${SEQ_PLATFORM}\tSM:$1\" $$(BWAMEM_REF_FASTA) $$^ | $$(SAMTOOLS) view -bhS - > $$@")
endef
$(foreach ss,$(SPLIT_SAMPLES),$(if $(fq.$(ss)),$(eval $(call align-split-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))

bwamem/bam/%.bwamem.bam : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	LBID=`echo "$*" | sed 's/_[A-Za-z0-9]\+//'`; \
	$(call RUN,-n $(BWAMEM_THREADS) -s 1G -m $(BWAMEM_MEM_PER_THREAD),"$(BWA) mem -t $(BWAMEM_THREADS) $(BWA_ALN_OPTS) -R \"@RG\tID:$*\tLB:$${LBID}\tPL:${SEQ_PLATFORM}\tSM:$${LBID}\" $(BWAMEM_REF_FASTA) $^ | $(SAMTOOLS) view -bhS - > $@")

bwamem/bam/%.bwamem.bam : fastq/%.fastq.gz
	LBID=`echo "$*" | sed 's/_[A-Za-z0-9]\+//'`; \
	$(call RUN,-n $(BWAMEM_THREADS) -s 1G -m $(BWAMEM_MEM_PER_THREAD),"$(BWA) mem -t $(BWAMEM_THREADS) $(BWA_ALN_OPTS) -R \"@RG\tID:$*\tLB:$${LBID}\tPL:${SEQ_PLATFORM}\tSM:$${LBID}\" $(BWAMEM_REF_FASTA) $^ | $(SAMTOOLS) view -bhS - > $@")

fastq/%.fastq.gz : fastq/%.fastq
	$(call RUN,,"gzip -c $< > $(@) && $(RM) $<")


..DUMMY := $(shell mkdir -p version && \
			 $(BWA) &> version/tmp.txt && \
			 head -3 version/tmp.txt | tail -2 > version/bwa_mem.txt && \
			 rm version/tmp.txt && \
			 $(SAMTOOLS) --version >> version/bwa_mem.txt && \
			 echo "gatk3" >> version/bwa_mem.txt && \
			 $(GATK) --version >> version/bwa_mem.txt && \
			 echo "picard" >> version/bwa_mem.txt && \
			 $(PICARD) MarkDuplicates --version >> version/bwa_mem.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: bwa_mem

include innovation-lab/bam_tools/process_bam.mk
include innovation-lab/fastq_tools/fastq.mk
include innovation-lab/aligners/align.mk
