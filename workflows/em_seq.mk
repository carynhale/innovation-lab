include innovation-lab/Makefile.inc
include innovation-lab/genome_inc/b37.inc
include innovation-lab/config/waltz.inc

LOGDIR ?= log/em_seq.$(NOW)

em_seq : $(foreach sample,$(SAMPLES),bwamem/$(sample)/$(sample)_R1.fastq.gz) \
		 $(foreach sample,$(SAMPLES),bwamem/$(sample)/$(sample)_aln.bam) \
		 $(foreach sample,$(SAMPLES),bwamem/$(sample)/$(sample)_aln_srt.bam) \
		 $(foreach sample,$(SAMPLES),bwamem/$(sample)/$(sample)_aln_srt_fx.bam)

REF_FASTA = $(REF_DIR)/IDT_oligo/idt_oligo.fasta

BWA_ALN_OPTS ?= -M -L 100,100
BWAMEM_THREADS = 12
BWAMEM_MEM_PER_THREAD = 2G

SAMTOOLS_THREADS = 8
SAMTOOLS_MEM_THREAD = 2G

WALTZ_MIN_MAPQ ?= 20
		 
define copy-fastq
bwamem/$1/$1_R1.fastq.gz : $3
	$$(call RUN,-c -n 1 -s 2G -m 4G,"set -o pipefail && \
								     $(RSCRIPT) $(SCRIPTS_DIR)/fastq_tools/copy_fastq.R \
								     --sample_name $1 \
								     --directory_name bwamem \
								     --fastq_files '$$^'")

endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)),$(eval $(call copy-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))
	
define fastq-to-bam
bwamem/$1/$1_aln.bam : bwamem/$1/$1_R1.fastq.gz
	$$(call RUN,-c -n $(BWAMEM_THREADS) -s 1G -m $(BWAMEM_MEM_PER_THREAD),"set -o pipefail && \
																		   $$(BWA) mem -t $$(BWAMEM_THREADS) $$(BWA_ALN_OPTS) \
																		   -R \"@RG\tID:$1\tLB:$1\tPL:$$(SEQ_PLATFORM)\tSM:$1\" $$(REF_FASTA) bwamem/$1/$1_R1.fastq.gz bwamem/$1/$1_R2.fastq.gz | $$(SAMTOOLS) view -bhS - > $$(@)")
									  									   
bwamem/$1/$1_aln_srt.bam : bwamem/$1/$1_aln.bam
	$$(call RUN,-c -n $(SAMTOOLS_THREADS) -s 1G -m $(SAMTOOLS_MEM_THREAD),"set -o pipefail && \
									  									   $$(SAMTOOLS) sort -@ $$(SAMTOOLS_THREADS) -m $$(SAMTOOLS_MEM_THREAD) $$(^) -o $$(@) -T $$(TMPDIR) && \
									  									   $$(SAMTOOLS) index $$(@) && \
									  									   cp bwamem/$1/$1_aln_srt.bam.bai bwamem/$1/$1_aln_srt.bai")
																		   
bwamem/$1/$1_aln_srt_fx.bam : bwamem/$1/$1_aln_srt.bam
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
		
define waltz-genotype
waltz/$1-pileup.txt.gz : bam/$1.bam
	$$(call RUN,-c -n 4 -s 4G -m 6G,"set -o pipefail && \
									 mkdir -p waltz && \
									 cd waltz && \
									 ln -sf ../bam/$1.bam $1.bam && \
									 ln -sf ../bam/$1.bai $1.bai && \
									 if [[ ! -f '.bed' ]]; then cut -f 4 $$(TARGETS_FILE) | paste -d '\t' $$(TARGETS_FILE) - > .bed; fi && \
									 $$(call WALTZ_CMD,2G,8G) org.mskcc.juber.waltz.Waltz PileupMetrics $$(WALTZ_MIN_MAPQ) $1.bam $$(REF_FASTA) .bed && \
									 gzip $1-pileup.txt && \
									 gzip $1-pileup-without-duplicates.txt && \
									 gzip $1-intervals.txt && \
									 gzip $1-intervals-without-duplicates.txt && \
									 cd ..")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call waltz-genotype,$(sample))))
		
..DUMMY := $(shell mkdir -p version; \
			 $(BWA) &> version/tmp.txt; \
			 head -3 version/tmp.txt | tail -2 > version/em_seq.txt; \
			 rm version/tmp.txt; \
			 $(SAMTOOLS) --version >> version/em_seq.txt; \
			 R --version >> version/em_seq.txt; \
			 $(JAVA8) -version &> version/em_seq.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: em_seq
