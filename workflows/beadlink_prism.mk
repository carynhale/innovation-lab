include innovation-lab/Makefile.inc
include innovation-lab/config/fgbio.inc
include innovation-lab/config/gatk.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/beadlink_prism.$(NOW)

beadlink_prism : $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_R1.fastq.gz) \
				 $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_fq.bam) \
				 $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_fq_srt.bam) \
				 $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl.fastq.gz) \
				 $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt.bam) \
				 $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD.bam) \
				 $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD.intervals) \
				 $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR.bam) \
				 $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX.bam) \
				 $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX__grp.bam) \
				 $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX__grp_DC.bam) \
				 $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX__grp_DC.duplex_umi_counts.txt) \
				 $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX__grp_DC_MA.bam) \
				 $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX__grp_DC_MA_RG.bam) \
				 $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX__grp_DC_MA_RG.intervals)
#				 $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX__grp_DC_MA_RG_IR.bam) \
#				 $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX__grp_DC_MA_RG_IR_FX.bam) \
#				 $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX__grp_DC_MA_RG_IR_FX-simplex.bam) \
#				 $(foreach sample,$(SAMPLES),fgbio/$(sample)/$(sample)_cl_aln_srt_MD_IR_FX__grp_DC_MA_RG_IR_FX-duplex.bam)

BWAMEM_THREADS = 12
BWAMEM_MEM_PER_THREAD = 2G

SAMTOOLS_THREADS = 8
SAMTOOLS_MEM_THREAD = 2G

GATK_THREADS = 8
GATK_MEM_THREAD = 2G

define copy-fastq
fgbio/$1/$1_R1.fastq.gz : $3
	$$(call RUN,-c -n 1 -s 2G -m 4G,"set -o pipefail && \
									 mkdir -p fgbio/$1 && \
								     $(RSCRIPT) $(SCRIPTS_DIR)/fastq_tools/copy_fastq.R \
								     --sample_name $1 \
								     --directory_name fgbio \
								     --fastq_files '$$^'")

endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)),$(eval $(call copy-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))

define fastq-2-bam
fgbio/$1/$1_fq.bam : fgbio/$1/$1_R1.fastq.gz
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $$(call FGBIO_CMD,2G,8G) \
									  FastqToBam \
									  --input fgbio/$1/$1_R1.fastq.gz fgbio/$1/$1_R2.fastq.gz \
									  --read-structures 8M+T 8M+T \
									  --output $$(@) \
									  --sample $1 \
									  --library $1 \
									  --platform illumina \
									  --platform-unit NA")

endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call fastq-2-bam,$(sample))))
	
define merge-bams
fgbio/$1/$1_fq_srt.bam : fgbio/$1/$1_fq.bam
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $$(SORT_SAM) \
									  INPUT=$$(<) \
									  OUTPUT=$$(@) \
									  SORT_ORDER=queryname")

fgbio/$1/$1_cl.fastq.gz : fgbio/$1/$1_fq_srt.bam
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $$(MARK_ADAPTERS) \
									  INPUT=$$(<) \
									  OUTPUT=/dev/stdout \
									  METRICS=fgbio/$1/$1_adapter-metrics.txt | \
									  $$(SAM_TO_FASTQ) \
									  INPUT=/dev/stdin \
									  FASTQ=$$(@) \
									  INTERLEAVE=true \
									  CLIPPING_ATTRIBUTE=XT \
									  CLIPPING_ACTION=X \
									  CLIPPING_MIN_LENGTH=25")

fgbio/$1/$1_cl_aln_srt.bam : fgbio/$1/$1_cl.fastq.gz fgbio/$1/$1_fq_srt.bam
	$$(call RUN,-c -n $(BWAMEM_THREADS) -s 1G -m $(BWAMEM_MEM_PER_THREAD),"set -o pipefail && \
																		   $$(BWA) mem -p -t $$(BWAMEM_THREADS) $$(REF_FASTA) $$(<) | \
																		   $$(MERGE_ALIGNMENTS) \
																		   UNMAPPED=$$(<<) \
																		   ALIGNED=/dev/stdin \
																		   OUTPUT=$$(@) \
																		   REFERENCE_SEQUENCE=$$(REF_FASTA) \
																		   SORT_ORDER=coordinate \
																		   MAX_GAPS=-1 \
																		   ORIENTATIONS=FR")

fgbio/$1/$1_cl_aln_srt_MD.bam : fgbio/$1/$1_cl_aln_srt.bam
	$$(call RUN, -c -n 1 -s 12G -m 18G,"set -o pipefail && \
										$$(MARK_DUP) \
										INPUT=$$(<) \
										OUTPUT=$$(@) \
										METRICS_FILE=fgbio/$1/$1_cl_aln_srt.txt \
										REMOVE_DUPLICATES=false \
										ASSUME_SORTED=true && \
										$$(SAMTOOLS) index $$(@) && \
									  	cp fgbio/$1/$1_cl_aln_srt_MD.bam.bai fgbio/$1/$1_cl_aln_srt_MD.bai")
									  	
fgbio/$1/$1_cl_aln_srt_MD.intervals : fgbio/$1/$1_cl_aln_srt_MD.bam
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD),"set -o pipefail && \
									   							   $$(call GATK_CMD,16G) \
									   							   -T RealignerTargetCreator \
									   							   -I $$(^) \
									   							   -nt $$(GATK_THREADS) \
									   							   -R $$(REF_FASTA) \
									   							   -o $$(@) \
									   							   -known $$(KNOWN_INDELS)")

fgbio/$1/$1_cl_aln_srt_MD_IR.bam : fgbio/$1/$1_cl_aln_srt_MD.bam fgbio/$1/$1_cl_aln_srt_MD.intervals
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD),"set -o pipefail && \
									   							   $$(call GATK_CMD,16G) \
							   							   		   -T IndelRealigner \
							   							   		   -I $$(<) \
							   							   		   -R $$(REF_FASTA) \
							   							   		   -targetIntervals $$(<<) \
							   							   		   -o $$(@) \
									   							   -known $$(KNOWN_INDELS)")
									   							   
fgbio/$1/$1_cl_aln_srt_MD_IR_FX.bam : fgbio/$1/$1_cl_aln_srt_MD_IR.bam
	$$(call RUN,-c -n 1 -s 12G -m 16G,"set -o pipefail && \
									   $$(FIX_MATE) \
									   INPUT=$$(<) \
									   OUTPUT=$$(@) \
									   SORT_ORDER=coordinate \
									   COMPRESSION_LEVEL=0 \
									   CREATE_INDEX=true")

endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call merge-bams,$(sample))))
	
define create-consensus
fgbio/$1/$1_cl_aln_srt_MD_IR_FX__grp.bam : fgbio/$1/$1_cl_aln_srt_MD_IR_FX.bam
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $$(call FGBIO_CMD,2G,8G) \
									  GroupReadsByUmi \
									  --strategy=paired \
									  --input=$$(<) \
									  --output=$$(@) \
									  --family-size-histogram=fgbio/$1/$1_cl_aln_srt_MD_IR_FX__grp.txt")
									  
fgbio/$1/$1_cl_aln_srt_MD_IR_FX__grp_DC.bam : fgbio/$1/$1_cl_aln_srt_MD_IR_FX__grp.bam
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD),"set -o pipefail && \
																   $$(call FGBIO_CMD,2G,16G) \
																   CallDuplexConsensusReads \
																   --min-reads 1 1 0 \
																   --input $$(<) \
																   --output $$(@) \
																   --threads $$(GATK_THREADS)")

fgbio/$1/$1_cl_aln_srt_MD_IR_FX__grp_DC.duplex_umi_counts.txt : fgbio/$1/$1_cl_aln_srt_MD_IR_FX__grp.bam
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $$(call FGBIO_CMD,2G,16G) \
									  CollectDuplexSeqMetrics \
									  --input $$(<) \
									  --output fgbio/$1/$1_cl_aln_srt_MD_IR_FX__grp_DC \
									  --duplex-umi-counts true")

endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call create-consensus,$(sample))))
	
define align-consensus
fgbio/$1/$1_cl_aln_srt_MD_IR_FX__grp_DC_MA.bam : fgbio/$1/$1_cl_aln_srt_MD_IR_FX__grp_DC.bam
	$$(call RUN,-c -n $(BWAMEM_THREADS) -s 1G -m $(BWAMEM_MEM_PER_THREAD),"set -o pipefail && \
																		   $$(SAM_TO_FASTQ) \
																		   INPUT=$$(<) \
																		   FASTQ=/dev/stdout \
																		   INTERLEAVE=true | \
																		   $$(BWA) mem -p -t $$(BWAMEM_THREADS) $$(REF_FASTA) /dev/stdin | \
																		   $$(MERGE_ALIGNMENTS) \
																		   UNMAPPED=$$(<) \
																		   ALIGNED=/dev/stdin \
																		   OUTPUT=$$(@) \
																		   REFERENCE_SEQUENCE=$$(REF_FASTA) \
																		   MAX_GAPS=-1 \
																		   ORIENTATIONS=FR")

fgbio/$1/$1_cl_aln_srt_MD_IR_FX__grp_DC_MA_RG.bam : fgbio/$1/$1_cl_aln_srt_MD_IR_FX__grp_DC_MA.bam
	$$(call RUN, -c -n 1 -s 12G -m 18G,"set -o pipefail && \
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
										cp fgbio/$1/$1_cl_aln_srt_MD_IR_FX__grp_DC_MA_RG.bam.bai fgbio/$1/$1_cl_aln_srt_MD_IR_FX__grp_DC_MA_RG.bai")

fgbio/$1/$1_cl_aln_srt_MD_IR_FX__grp_DC_MA_RG.intervals : fgbio/$1/$1_cl_aln_srt_MD_IR_FX__grp_DC_MA_RG.bam
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD),"set -o pipefail && \
									   							   $$(call GATK_CMD,16G) \
									   							   -T RealignerTargetCreator \
									   							   -I $$(^) \
									   							   -nt $$(GATK_THREADS) \
									   							   -R $$(REF_FASTA) \
									   							   -o $$(@) \
									   							   -known $$(KNOWN_INDELS) \
									   							   --allow_potentially_misencoded_quality_scores")

fgbio/$1/$1_cl_aln_srt_MD_IR_FX__grp_DC_MA_RG_IR.bam : fgbio/$1/$1_cl_aln_srt_MD_IR_FX__grp_DC_MA_RG.bam fgbio/$1/$1_cl_aln_srt_MD_IR_FX__grp_DC_MA_RG.intervals
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD),"set -o pipefail && \
									   							   $$(call GATK_CMD,16G) \
							   							   		   -T IndelRealigner \
							   							   		   -I $$(<) \
							   							   		   -R $$(REF_FASTA) \
							   							   		   -targetIntervals $$(<<) \
							   							   		   -o $$(@) \
									   							   -known $$(KNOWN_INDELS) \
									   							   --allow_potentially_misencoded_quality_scores")
									   							   
fgbio/$1/$1_cl_aln_srt_MD_IR_FX__grp_DC_MA_RG_IR_FX.bam : fgbio/$1/$1_cl_aln_srt_MD_IR_FX__grp_DC_MA_RG_IR.bam
	$$(call RUN,-c -n 1 -s 18G -m 24G,"set -o pipefail && \
									   $$(FIX_MATE) \
									   INPUT=$$(<) \
									   OUTPUT=$$(@) \
									   SORT_ORDER=coordinate \
									   COMPRESSION_LEVEL=0 \
									   CREATE_INDEX=true")
									   
endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call align-consensus,$(sample))))
	
define filter-consensus
fgbio/$1/$1_cl_aln_srt_MD_IR_FX__grp_DC_MA_RG_IR_FX-simplex.bam : fgbio/$1/$1_cl_aln_srt_MD_IR_FX__grp_DC_MA_RG_IR_FX.bam
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $$(call FGBIO_CMD,2G,16G) \
									  FilterConsensusReads \
									  --input $$(<) \
									  --output $$(@) \
									  --ref $$(REF_FASTA) \
									  --min-reads=3 3 0 \
									  --min-base-quality=30 \
									  --reverse-per-base-tags=true")

fgbio/$1/$1_cl_aln_srt_MD_IR_FX__grp_DC_MA_RG_IR_FX-duplex.bam : fgbio/$1/$1_cl_aln_srt_MD_IR_FX__grp_DC_MA_RG_IR_FX.bam
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $$(call FGBIO_CMD,2G,16G) \
									  FilterConsensusReads \
									  --input $$(<) \
									  --output $$(@) \
									  --ref $$(REF_FASTA) \
									  --min-reads=2 1 1 \
									  --min-base-quality=30 \
									  --reverse-per-base-tags=true")

endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call filter-consensus,$(sample))))

	
..DUMMY := $(shell mkdir -p version; \
			 $(JAVA8) -jar $(FGBIO) --help &> version/fgbio_access.txt; \
			 echo "picard" >> version/fgbio_access.txt; \
			 $(PICARD) SortSam --version &>> version/fgbio_access.txt; \
			 $(PICARD) MarkIlluminaAdapters --version &>> version/fgbio_access.txt; \
			 $(PICARD) SamToFastq --version &>> version/fgbio_access.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: fgbio_access
