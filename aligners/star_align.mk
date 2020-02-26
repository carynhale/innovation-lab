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

star : $(foreach sample,$(SAMPLES),star/$(sample).taskcomplete)

#$(foreach sample,$(SAMPLES),bam/$(sample).bam) \
#$(foreach sample,$(SAMPLES),star/$(sample).Chimeric.out.junction)

define align-split-fastq
star/$1.taskcomplete : $3
	$$(call RUN,-n 4 -s 6G -m 10G,"STAR $$(STAR_OPTS) \
                                   --outFileNamePrefix star/$1. \
                                   --runThreadN 4 \
                                   --outSAMattrRGline \"ID:$1\" \"LB:$1\" \"SM:$1\" \"PL:$${SEQ_PLATFORM}\" \
                                   --readFilesIn $$^ \
                                   --readFilesCommand zcat && \
                                   touch $$@")
endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)),\
	$(eval $(call align-split-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))
    
#    star/$2.Aligned.sortedByCoord.out.bam : star/$2.taskcomplete
#    star/$2.Chimeric.out.junction : star/$2.taskcomplete
#
#
#    star/%.Aligned.sortedByCoord.out.bam star/%.Chimeric.out.junction : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
#        $(call RUN,-n 4 -s 6G -m 10G -v $(STAR_ENV),"STAR $(STAR_OPTS) \
#                                                     --outFileNamePrefix star/$*. --runThreadN 4 \
#                                                     --outSAMattrRGline \"ID:$*\" \"LB:$*\" \"SM:$*\" \"PL:${SEQ_PLATFORM}\" \
#                                                     --readFilesIn $^ --readFilesCommand zcat")
#
#    star/bam/%.star.sorted.bam : star/%.Aligned.sortedByCoord.out.bam
#        $(INIT) mv $< $@
#
#    bam/%.bam : star/bam/%.star.$(BAM_SUFFIX)
#        $(INIT) ln -f $(<) $(@) 
#
#    define merged-chimeric-junction
#    star/$1.Chimeric.out.junction : $$(foreach split,$$(split.$1),star/$$(split).Chimeric.out.junction)
#        $$(INIT) sort -V $$^ > $$@
#    endef
#    $(foreach sample,$(SAMPLES),\
#        $(eval $(call merged-chimeric-junction,$(sample))))

..DUMMY := $(shell mkdir -p version; \
             echo "STAR" > version/star_align.txt; \
			 STAR --version >> version/star_align.txt)
.SECONDARY: 
.DELETE_ON_ERROR:
.PHONY: star

include innovation-lab/fastq_tools/fastq.mk
include innovation-lab/bam_tools/process_bam.mk
include innovation-lab/aligners/align.mk
