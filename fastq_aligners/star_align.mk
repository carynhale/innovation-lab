include innovation-lab/Makefile.inc
include innovation-lab/config/align.inc

LOGDIR = log/star_align.$(NOW)

ALIGNER := star
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

star : $(foreach sample,$(SAMPLES),star/$(sample)/$(sample)_R1.fastq.gz) \
       $(foreach sample,$(SAMPLES),star/$(sample)/$(sample)_R2.fastq.gz) \
       $(foreach sample,$(SAMPLES),bam/$(sample).bam) \
       $(foreach sample,$(SAMPLES),bam/$(sample).bam.bai)
       
define merge-fastq
star/$1/$1_R1.fastq.gz : $$(foreach split,$2,$$(word 1, $$(fq.$$(split))))
	$$(call RUN,-c -n 12 -s 0.5G -m 1G -w 24:00:00 -v $(PIGZ_ENV),"set -o pipefail && \
								       $$(PIGZ) -cd -p 12 $$(^) | $$(PIGZ) -c -p 12 > $$(@)")
	
star/$1/$1_R2.fastq.gz : $$(foreach split,$2,$$(word 2, $$(fq.$$(split))))
	$$(call RUN,-c -n 12 -s 0.5G -m 1G -w 24:00:00 -v $(PIGZ_ENV),"set -o pipefail && \
								       $$(PIGZ) -cd -p 12 $$(^) | $$(PIGZ) -c -p 12 > $$(@)")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call merge-fastq,$(sample),$(split.$(sample)))))
       

define align-fastq
star/$1/$1.Aligned.sortedByCoord.out.bam : star/$1/$1_R1.fastq.gz star/$1/$1_R2.fastq.gz
	$$(call RUN,-n 12 -s 2G -m 4G,"set -o pipefail && \
				       STAR $$(STAR_OPTS) \
                                       --outFileNamePrefix star/$1/$1. \
                                       --runThreadN 12 \
                                       --outSAMattrRGline \"ID:$1\" \"LB:$1\" \"SM:$1\" \"PL:$${SEQ_PLATFORM}\" \
                                       --readFilesIn $$(^) \
                                       --readFilesCommand zcat")
                                   
star/$1/$1.Aligned.sortedByCoord.out.bam.bai : star/$1/$1.Aligned.sortedByCoord.out.bam
	$$(call RUN,-n 1 -s 2G -m 4G,"set -o pipefail && \
				      $$(SAMTOOLS) index $$(<)")

bam/$1.bam : star/$1/$1.Aligned.sortedByCoord.out.bam
	$$(call RUN,-n 1 -s 2G -m 4G,"set -o pipefail && \
				      cp $$(<) $$(@)")
                                   
bam/$1.bam.bai : star/$1/$1.Aligned.sortedByCoord.out.bam.bai
	$$(call RUN,-n 1 -s 2G -m 4G,"set -o pipefail && \
                                      cp $$(<) $$(@)")
                                   
endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call align-fastq,$(sample))))
	
	
    
..DUMMY := $(shell mkdir -p version; \
             echo "STAR" > version/star_align.txt; \
	     STAR --version >> version/star_align.txt; \
             $(SAMTOOLS) --version >> version/star_align.txt)
.SECONDARY: 
.DELETE_ON_ERROR:
.PHONY: star
