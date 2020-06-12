include innovation-lab/Makefile.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/em_seq.$(NOW)

em_seq : $(foreach sample,$(SAMPLES),emseq/$(sample)/$(sample)_R1.fastq.gz) \
		 $(foreach sample,$(SAMPLES),emseq/$(sample)/$(sample)_aln.bam) \
		 $(foreach sample,$(SAMPLES),emseq/$(sample)/$(sample)_aln_srt.bam) \
		 $(foreach sample,$(SAMPLES),bam/$(sample).bam)

REF_FASTA = $(REF_DIR)/IDT_oligo/idt_oligo.fasta

BWA_ALN_OPTS ?= -M
BWAMEM_THREADS = 12
BWAMEM_MEM_PER_THREAD = 2G

SAMTOOLS_THREADS = 8
SAMTOOLS_MEM_THREAD = 2G
		 
define copy-fastq
emseq/$1/$1_R1.fastq.gz : $3
	$$(call RUN,-c -n 1 -s 2G -m 4G,"set -o pipefail && \
									 mkdir -p emseq/$1 && \
								     $(RSCRIPT) $(SCRIPTS_DIR)/fastq_tools/copy_fastq.R \
								     --sample_name $1 \
								     --directory_name emseq \
								     --fastq_files '$$^'")

endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)),$(eval $(call copy-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))
	
define fastq-to-bam
emseq/$1/$1_aln.bam : emseq/$1/$1_R1.fastq.gz
	$$(call RUN,-c -n $(BWAMEM_THREADS) -s 1G -m $(BWAMEM_MEM_PER_THREAD),"set -o pipefail && \
																		   $$(BWA) mem -t $$(BWAMEM_THREADS) $$(BWA_ALN_OPTS) \
																		   -R \"@RG\tID:$1\tLB:$1\tPL:$$(SEQ_PLATFORM)\tSM:$1\" $$(REF_FASTA) emseq/$1/$1_R1.fastq.gz | $$(SAMTOOLS) view -bhS - > $$(@)")

emseq/$1/$1_aln_srt.bam : emseq/$1/$1_aln.bam
	$$(call RUN,-c -n $(SAMTOOLS_THREADS) -s 1G -m $(SAMTOOLS_MEM_THREAD),"set -o pipefail && \
									  									   $$(SAMTOOLS) sort -@ $$(SAMTOOLS_THREADS) -m $$(SAMTOOLS_MEM_THREAD) $$(^) -o $$(@) -T $$(TMPDIR) && \
									  									   $$(SAMTOOLS) index $$(@) && \
									  									   cp emseq/$1/$1_aln_srt.bam.bai emseq/$1/$1_aln_srt.bai")
																		   
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call fastq-to-bam,$(sample))))
		
define copy-to-bam
bam/$1.bam : emseq/$1/$1_aln_srt.bam
	$$(call RUN, -c -s 2G -m 4G ,"set -o pipefail && \
								  cp $$(<) $$(@) && \
								  $$(SAMTOOLS) index $$(@) && \
								  cp bam/$1.bam.bai bam/$1.bai")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call copy-to-bam,$(sample))))
		
..DUMMY := $(shell mkdir -p version; \
			 $(BWA) &> version/tmp.txt; \
			 head -3 version/tmp.txt | tail -2 > version/em_seq.txt; \
			 rm version/tmp.txt; \
			 $(SAMTOOLS) --version >> version/em_seq.txt; \
			 R --version >> version/em_seq.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: em_seq
