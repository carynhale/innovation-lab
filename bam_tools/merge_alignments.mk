include innovation-lab/Makefile.inc
include innovation-lab/genome_inc/b37.inc

LOGDIR ?= log/merge_alignments.$(NOW)

merge_alignments : $(foreach sample,$(SAMPLES),merged_alignments/$(sample).bam)

define fix-bam
merged_alignments/%.ubam : bam/%.bam
	$$(call RUN,-c -n 1 -s 12G -m 16G,"set -o pipefail && \
					   $$(REVERT_SAM) \
					   INPUT=$$(<) \
					   OUTPUT=merged_alignments/$$(*).ubam \
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
					    OUTPUT=merged_alignments/$$(*).fixed.bam \
					    CREATE_INDEX=true \
					    ADD_MATE_CIGAR=true \
					    CLIP_ADAPTERS=true \
					    CLIP_OVERLAPPING_READS=true \
					    INCLUDE_SECONDARY_ALIGNMENTS=false \
					    MAX_INSERTIONS_OR_DELETIONS=-1 && \
					    rm -rf merged_alignments/$$(*).ubam")
					    
merged_alignments/%.dedup.bam : merged_alignments/%.fixed.bam
	$$(call RUN, -c -n 1 -s 12G -m 16G,"set -o pipefail && \
					    $$(MARK_DUP) \
					    INPUT=$$(<) \
					    OUTPUT=merged_alignments/$$(*).dedup.bam \
					    METRICS_FILE=unprocessed_bam/$$(*).txt && \
					    rm -rf merged_alignments/$$(*).fixed.bam")
					    
fixed_bam/%.bam : unprocessed_bam/%.dedup.bam
	$$(call RUN, -c -n 1 -s 12G -m 16G,"set -o pipefail && \
					    $$(ADD_RG) \
					    INPUT=$$(<) \
					    OUTPUT=fixed_bam/$$(*).bam \
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


..DUMMY := $(shell mkdir -p version; \
	     echo "picard" > version/merge_alignments.txt; \
	     $(PICARD) MarkDuplicates --version >> version/merge_alignments.txt; \
	     $(SAMTOOLS) --version >> version/merge_alignments.txt)
.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY: merge_alignments
