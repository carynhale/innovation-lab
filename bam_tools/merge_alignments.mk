include innovation-lab/Makefile.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/merge_alignments.$(NOW)

merge_alignments : $(foreach sample,$(SAMPLES),merged_alignments/$(sample).bam)

define merge-alignments
merged_alignments/%.ubam : bam/%.bam
	$$(call RUN,-c -n 1 -s 12G -m 16G,"set -o pipefail && \
					   $$(REVERT_SAM) \
					   INPUT=$$(<) \
					   OUTPUT=$$(@) \
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

merged_alignments/%.fixed.bam : merged_alignments/%.ubam
	$$(call RUN, -c -n 1 -s 12G -m 16G,"set -o pipefail && \
					    $$(MERGE_ALIGNMENTS) \
					    REFERENCE_SEQUENCE=$$(REF_FASTA) \
					    UNMAPPED_BAM=$$(<) \
					    ALIGNED_BAM=bam/$$(*).bam \
					    OUTPUT=$$(@) \
					    CREATE_INDEX=true \
					    ADD_MATE_CIGAR=true \
					    CLIP_ADAPTERS=true \
					    CLIP_OVERLAPPING_READS=true \
					    INCLUDE_SECONDARY_ALIGNMENTS=false \
					    MAX_INSERTIONS_OR_DELETIONS=-1 && \
					    $$(RM) $$(<)")
					    
merged_alignments/%.dedup.bam : merged_alignments/%.fixed.bam
	$$(call RUN, -c -n 1 -s 12G -m 16G,"set -o pipefail && \
					    $$(MARK_DUP) \
					    INPUT=$$(<) \
					    OUTPUT=$$(@) \
					    METRICS_FILE=merged_alignments/$$(*).txt && \
					    $(RM) $$(<)")
					    
merged_alignments/%.bam : merged_alignments/%.dedup.bam
	$$(call RUN, -c -n 1 -s 12G -m 16G,"set -o pipefail && \
					    $$(ADD_RG) \
					    INPUT=$$(<) \
					    OUTPUT=$$(@) \
					    RGID=$$(*) \
					    RGLB=$$(*) \
					    RGPL=illumina \
					    RGPU=NA \
					    RGSM=$$(*) && \
					    samtools index $$(@) && \
					    cp merged_alignments/$$(*).bam.bai merged_alignments/$$(*).bai && \
					    $$(RM) merged_alignments/$$(*).dedup.bam && \
					    $$(RM) merged_alignments/$$(*).fixed.bai && \
					    $$(RM) merged_alignments/$$(*).dedup.bai && \
					    $$(RM) merged_alignments/$$(*).txt")
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call merge-alignments,$(sample))))


..DUMMY := $(shell mkdir -p version; \
	     echo "picard" > version/merge_alignments.txt; \
	     $(PICARD) MarkDuplicates --version >> version/merge_alignments.txt; \
	     $(SAMTOOLS) --version >> version/merge_alignments.txt)
.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY: merge_alignments
