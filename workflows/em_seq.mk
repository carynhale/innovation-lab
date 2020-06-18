include innovation-lab/Makefile.inc
include innovation-lab/genome_inc/b37.inc
include innovation-lab/config/waltz.inc

LOGDIR ?= log/em_seq.$(NOW)

em_seq : $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_R1.fastq.gz) \
		 $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_aln.bam) \
		 $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_aln_srt.bam) \
		 $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_aln_srt__F1R1R2.bam) \
		 $(foreach sample,$(SAMPLES),bismark/$(sample)/$(sample)_aln_srt__F2R1R2.bam)

REF_FASTA = $(REF_DIR)/IDT_oligo/idt_oligo.fasta
GENOME_FOLDER = $(REF_DIR)/IDT_oligo/

SAMTOOLS_THREADS = 8
SAMTOOLS_MEM_THREAD = 2G

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
																		   
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call fastq-to-bam,$(sample))))
		
define filter-bam
bismark/$1/$1_aln_srt__F1R1R2.bam : bismark/$1/$1_aln_srt.bam
	$$(call RUN,-c -n $(SAMTOOLS_THREADS) -s 1G -m $(SAMTOOLS_MEM_THREAD),"set -o pipefail && \
																		   $$(SAMTOOLS) view -@ $$(SAMTOOLS_THREADS) -m $$(SAMTOOLS_MEM_THREAD) -b -f 144 $$(<) > bismak/$1/$1_aln_srt__F1R1.bam && \
																		   $$(SAMTOOLS) index bismark/$1/$1_aln_srt__F1R1.bam && \
																		   $$(SAMTOOLS) view -@ $$(SAMTOOLS_THREADS) -m $$(SAMTOOLS_MEM_THREAD) -b -f 64 -F 16 $$(<) > bismark/$1/$1_aln_srt__F1R2.bam && \
																		   $$(SAMTOOLS) index bismark/$1/$1_aln_srt__F1R2.bam && \
																		   $$(SAMTOOLS) merge -@ $$(SAMTOOLS_THREADS) -m $$(SAMTOOLS_MEM_THREAD) -f $$(@) bismark/$1/$1_aln_srt__F1R1.bam bismark/$1/$1_aln_srt__F1R2.bam && \
																		   $$(SAMTOOLS) index $$(@)")

bismark/$1/$1_aln_srt__F2R1R2.bam : bismark/$1/$1_aln_srt.bam
	$$(call RUN,-c -n $(SAMTOOLS_THREADS) -s 1G -m $(SAMTOOLS_MEM_THREAD),"set -o pipefail && \
																		   $$(SAMTOOLS) view -@ $$(SAMTOOLS_THREADS) -m $$(SAMTOOLS_MEM_THREAD) -b -f 128 -F 16 $$(<) > bismark/$1/$1_aln_srt__F2R1.bam && \
																		   $$(SAMTOOLS) index bismark/$1/$1_aln_srt__F2R1.bam && \
																		   $$(SAMTOOLS) view -@ $$(SAMTOOLS_THREADS) -m $$(SAMTOOLS_MEM_THREAD) -b -f 80 $$(<) > bismark/$1/$1_aln_srt__F2R2.bam && \
																		   $$(SAMTOOLS) index bismark/$1/$1_aln_srt__F2R2.bam && \
																		   $$(SAMTOOLS) merge -@ $$(SAMTOOLS_THREADS) -m $$(SAMTOOLS_MEM_THREAD) -f $$(@) bismark/$1/$1_aln_srt__F2R1.bam bismark/$1/$1_aln_srt__F2R2.bam && \
																		   $$(SAMTOOLS) index $$(@)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call filter-bam,$(sample))))

		
..DUMMY := $(shell mkdir -p version; \
			 $(HOME)/share/usr/env/bismark-0.22.1/bin/bismark --version > version/em_seq.txt; \
			 $(SAMTOOLS) --version >> version/em_seq.txt; \
			 R --version >> version/em_seq.txt; \
			 $(JAVA8) -version &> version/em_seq.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: em_seq
