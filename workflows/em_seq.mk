include innovation-lab/Makefile.inc

LOGDIR ?= log/em_seq.$(NOW)

em_seq : $(foreach sample,$(SAMPLES),emseq/$(sample)/$(sample)_R1.fastq.gz)
		 
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
emseq/$1/$1.bam : emseq/$1/$1_R1.fastq.gz
	$$(call RUN,-c -n $(BWAMEM_THREADS) -s 1G -m $(BWAMEM_MEM_PER_THREAD),"set -o pipefail && \
																		   $$(BWA) mem -t $$(BWAMEM_THREADS) $$(BWA_ALN_OPTS) \
																		   -R \"@RG\tID:$1\tLB:$1\tPL:$$(SEQ_PLATFORM)\tSM:$1\" $$(REF_FASTA) marianas/$1/$1_R1_umi-clipped.fastq.gz marianas/$1/$1_R2_umi-clipped.fastq.gz | $$(SAMTOOLS) view -bhS - > $$(@)")
																		   
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call fastq-to-bam,$(sample))))
		
..DUMMY := $(shell mkdir -p version; \
			 $(BWA) &> version/tmp.txt; \
			 head -3 version/tmp.txt | tail -2 > version/em_seq.txt; \
			 rm version/tmp.txt; \
			 $(SAMTOOLS) --version >> version/em_seq.txt; \
			 R --version >> version/em_seq.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: em_seq
