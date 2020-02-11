include innovation-lab/Makefile.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/fix_bam.$(NOW)

fix_bam : $(foreach sample,$(SAMPLES),fixed_bam/$(sample).bam)

define fix-bam
unprocessed_bam/%.ubam : unprocessed_bam/%.bam
	$$(call RUN,-c -n 1 -s 12G -m 16G,"$$(REVERT_SAM) \
									   I=$$(<) \
									   O=unprocessed_bam/$$(*).ubam \
									   SANITIZE=true \
									   MAX_DISCARD_FRACTION=0.005 \
									   ATTRIBUTE_TO_CLEAR=XT \
									   ATTRIBUTE_TO_CLEAR=XN \
									   ATTRIBUTE_TO_CLEAR=AS \
									   ATTRIBUTE_TO_CLEAR=OC \
									   ATTRIBUTE_TO_CLEAR=OP \
									   SORT_ORDER=queryname \
									   RESTORE_ORIGINAL_QUALITIES=true \
									   REMOVE_DUPLICATE_INFORMATION=true \
									   REMOVE_ALIGNMENT_INFORMATION=true")
unprocessed_bam/%.fixed.bam : unprocessed_bam/%.ubam
	$$(call RUN, -c -n 1 -s 12G -m 16G,"$$(MERGE_ALIGNMENTS) \
										R=$$(DMP_FASTA) \
										UNMAPPED_BAM=$$(<) \
										ALIGNED_BAM=unprocessed_bam/$$(*).bam \
										O=unprocessed_bam/$$(*).fixed.bam \
										CREATE_INDEX=true \
										ADD_MATE_CIGAR=true \
										CLIP_ADAPTERS=true \
										CLIP_OVERLAPPING_READS=true \
										INCLUDE_SECONDARY_ALIGNMENTS=false \
										MAX_INSERTIONS_OR_DELETIONS=-1 && \
										rm -rf unprocessed_bam/$$(*).ubam")
unprocessed_bam/%.dedup.bam : unprocessed_bam/%.fixed.bam
	$$(call RUN, -c -n 1 -s 12G -m 16G,"$$(MARK_DUP) \
										I=$$(<) \
										O=unprocessed_bam/$$(*).dedup.bam \
										M=unprocessed_bam/$$(*).txt && \
										rm -rf unprocessed_bam/$$(*).fixed.bam")
fixed_bam/%.bam : unprocessed_bam/%.dedup.bam
	$$(call RUN, -c -n 1 -s 12G -m 16G,"$$(ADD_RG) \
										I=$$(<) \
										O=fixed_bam/$$(*).bam \
										RGID=$$(*) \
										RGLB=$$(*) \
										RGPL=illumina \
										RGPU=NA \
										RGSM=$$(*) && \
										samtools index fixed_bam/$$(*).bam && \
										cp fixed_bam/$$(*).bam.bai fixed_bam/$$(*).bai && \
										rm -rf unprocessed_bam/$$(*).dedup.bam && \
										rm -rf unprocessed_bam/$$(*).fixed.bai && \
										rm -rf unprocessed_bam/$$(*).dedup.bai && \
										rm -rf unprocessed_bam/$$(*).txt")
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call fix-bam,$(sample))))


..DUMMY := $(shell mkdir -p version && \
			 echo "picard" > version/fix_bam.txt && \
			 $(PICARD) MarkDuplicates --version >> version/fix_bam.txt && \
			 $(SAMTOOLS) --version >> version/fix_bam.txt)
.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY: fix_bam
