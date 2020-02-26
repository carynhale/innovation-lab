include innovation-lab/Makefile.inc
include innovation-lab/config/align.inc

LOGDIR = log/star_align.$(NOW)

ALIGNER := star
BAM_NO_REALN = true
BAM_NO_RECAL = true
BAM_NO_FILTER = true
BAM_DUP_TYPE = none
SEQ_PLATFROM = illumina

STAR_OPTS = --genomeDir $(STAR_REF) \
            --outSAMtype BAM SortedByCoordinate \
			--twopassMode Basic \
            --outReadsUnmapped None \
            --chimSegmentMin 12 \
            --chimJunctionOverhangMin 12 \
            --alignSJDBoverhangMin 10 \
            --alignMatesGapMax 200000 \
            --alignIntronMax 200000 \
            --chimSegmentReadGapMax parameter 3 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --chimOutType WithinBAM \
			--quantMode GeneCounts

star : $(foreach sample,$(SAMPLES),star/$(sample).Aligned.sortedByCoord.out.bam) \
       $(foreach sample,$(SAMPLES),bam/$(sample).bam)

#$(foreach sample,$(SAMPLES),bam/$(sample).bam.bai) \
#$(foreach sample,$(SAMPLES),bam/$(sample).bai)

define align-split-fastq
star/$1.Aligned.sortedByCoord.out.bam : $3
	$$(call RUN,-n 4 -s 6G -m 10G,"set -o pipefail && \
                                   STAR $$(STAR_OPTS) \
                                   --outFileNamePrefix star/$1. \
                                   --runThreadN 4 \
                                   --outSAMattrRGline \"ID:$1\" \"LB:$1\" \"SM:$1\" \"PL:$${SEQ_PLATFORM}\" \
                                   --readFilesIn $$^ \
                                   --readFilesCommand zcat")

endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)),\
	$(eval $(call align-split-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))
    
    
define copy-bam
bam/$1.bam : star/$1.Aligned.sortedByCoord.out.bam
    $$(call RUN,-n 1 -s 2G -m 4G,"set -o pipefail && \
                                  cp $$(<) $$(@)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call copy-bam,$(sample))))

#bam/$1.bam.bai : bam/$1.bam
#    $$(call RUN,-n 1 -s 2G -m 4G,"set -o pipefail && \
#                                  $$(SAMTOOLS) index $$(<)")
#
#bam/$1.bai : bam/$1.bam.bai
#    $$(call RUN,-n 1 -s 2G -m 4G,"set -o pipefail && \
#                                  cp $$(<) $$(@)")


..DUMMY := $(shell mkdir -p version; \
             echo "STAR" > version/star_align.txt; \
			 STAR --version >> version/star_align.txt; \
             $(SAMTOOLS) --version >> version/star_align.txt)
.SECONDARY: 
.DELETE_ON_ERROR:
.PHONY: star

include innovation-lab/fastq_tools/fastq.mk
include innovation-lab/bam_tools/process_bam.mk
include innovation-lab/aligners/align.mk
