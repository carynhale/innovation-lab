include innovation-lab/Makefile.inc
include innovation-lab/config/align.inc

LOGDIR = log/star_align.$(NOW)

SEQ_PLATFROM = illumina
STAR_THREADS = 16

STAR_OPTS = --genomeDir $(STAR_REF) \
	    --runThreadN $(STAR_THREADS) \
	    --sjdbGTFfile $(HOME)/share/lib/resource_files/Ensembl/GRCh37/Annotation/Genes/mirna.gtf \
	    --alignEndsType EndToEnd \
	    --outFilterMismatchNmax 1 \
	    --outFilterMultimapScoreRange 0 \
	    --quantMode TranscriptomeSAM GeneCounts \
	    --outReadsUnmapped Fastx \
	    --outSAMtype BAM SortedByCoordinate \
	    --outFilterMultimapNmax 10 \
	    --outSAMunmapped Within \
	    --outFilterScoreMinOverLread 0 \
	    --outFilterMatchNminOverLread 0 \
	    --outFilterMatchNmin 16 \
	    --alignSJDBoverhangMin 1000 \
	    --alignIntronMax 1 \
	    --outWigType wiggle \
	    --outWigStrand Stranded \
	    --outWigNorm RPM \

star : $(foreach sample,$(SAMPLES),star/$(sample).Aligned.sortedByCoord.out.bam \
                                   star/$(sample).Aligned.sortedByCoord.out.bam.bai \
                                   bam/$(sample).bam \
                                   bam/$(sample).bam.bai \
                                   bam/$(sample).bai)

define align-split-fastq
star/$1.Aligned.sortedByCoord.out.bam : $3
	$$(call RUN,-n $(STAR_THREADS) -s 2G -m 4G,"set -o pipefail && \
						    STAR $$(STAR_OPTS) \
						    --outFileNamePrefix star/$1. \
						    --runThreadN $$(STAR_THREADS) \
						    --outSAMattrRGline \"ID:$1\" \"LB:$1\" \"SM:$1\" \"PL:$${SEQ_PLATFORM}\" \
						    --readFilesIn $$^ \
						    --readFilesCommand zcat")
                                   
star/$1.Aligned.sortedByCoord.out.bam.bai : star/$1.Aligned.sortedByCoord.out.bam
	$$(call RUN,-n 1 -s 2G -m 4G,"set -o pipefail && \
                                  $$(SAMTOOLS) index $$(<)")

bam/$1.bam : star/$1.Aligned.sortedByCoord.out.bam
	$$(call RUN,-n 1 -s 2G -m 4G,"set -o pipefail && \
                                  cp $$(<) $$(@)")
                                   
bam/$1.bam.bai : star/$1.Aligned.sortedByCoord.out.bam.bai
	$$(call RUN,-n 1 -s 2G -m 4G,"set -o pipefail && \
                                  cp $$(<) $$(@)")
                                   
bam/$1.bai : star/$1.Aligned.sortedByCoord.out.bam.bai
	$$(call RUN,-n 1 -s 2G -m 4G,"set -o pipefail && \
                                  cp $$(<) $$(@)")

endef
$(foreach ss,$(SPLIT_SAMPLES), \
	$(if $(fq.$(ss)), \
	$(eval $(call align-split-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))
    
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
