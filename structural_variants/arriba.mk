include innovation-lab/Makefile.inc
include innovation-lab/config/arriba.inc

LOGDIR ?= log/arriba.$(NOW)

arriba : $(foreach sample,$(SAMPLES),arriba/$(sample)/$(sample).1.fastq.gz) \
		 $(foreach sample,$(SAMPLES),arriba/$(sample)/$(sample).2.fastq.gz) \
		 $(foreach sample,$(SAMPLES),arriba/$(sample)/$(sample).Aligned.out.bam) \
		 $(foreach sample,$(SAMPLES),arriba/$(sample)/fusions.tsv)


define merged-fastq
arriba/$1/$1.1.fastq.gz : $$(foreach split,$2,$$(word 1, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 2G -m 4G,"set -o pipefail && \
									 mkdir -p arriba/$1 && \
									 cp $$(^) $$(@)")
arriba/$1/$1.2.fastq.gz : $$(foreach split,$2,$$(word 2, $$(fq.$$(split))))
	$$(call RUN,-c -n 1 -s 2G -m 4G,"set -o pipefail && \
									 mkdir -p arriba/$1 && \
									 cp $$(^) $$(@)")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call merged-fastq,$(sample),$(split.$(sample)))))

define star-align
arriba/$1/$1.Aligned.out.bam : arriba/$1/$1.1.fastq.gz arriba/$1/$1.2.fastq.gz
	$$(call RUN,-c -n $(STAR_THREADS) -s 1G -m 2G -v $(ARRIBA_ENV),"set -o pipefail && \
																	STAR \
																	--runThreadN $$(STAR_THREADS) \
																	--genomeDir $$(STAR_INDEX_DIR) \
																	--genomeLoad NoSharedMemory \
																	--readFilesIn $$(<) $$(<<) \
																	--readFilesCommand zcat \
																	--outStd BAM_Unsorted \
																	--outSAMtype BAM Unsorted \
																	--outSAMunmapped Within \
																	--outBAMcompression 0 \
																	--outFilterMultimapNmax 50 \
																	--peOverlapNbasesMin 10 \
																	--alignSplicedMateMapLminOverLmate 0.5 \
																	--alignSJstitchMismatchNmax 5 -1 5 5 \
																	--chimSegmentMin 10 \
																	--chimOutType WithinBAM HardClip \
																	--chimJunctionOverhangMin 10 \
																	--chimScoreDropMax 30 \
																	--chimScoreJunctionNonGTAG 0 \
																	--chimScoreSeparation 1 \
																	--chimSegmentReadGapMax 3 \
																	--chimMultimapNmax 50 \
																	--outFileNamePrefix arriba/$1/$1. > $$(*)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call star-align,$(sample))))
		
define run-arriba
arriba/$1/fusions.tsv : arriba/$1/$1.Aligned.out.bam
	$$(call RUN,-c -n 1 -s 24G -m 36G -v $(ARRIBA_ENV),"set -o pipefail && \
														$$(ARRIBA_EXE) -x $$(<) \
														-o arriba/$1/fusions.tsv -O arriba/$1/discarded.tsv \
														-a $$(ASSEMBLY_FA) \
														-g $$(ANNOTATION_GTF) \
														-b $$(BLACKLIST_TSV) \
														-k $$(KNOWN_FUSIONS_TSV) \
														-t $$(KNOWN_FUSIONS_TSV) \
														-p $$(PROTEIN_DOMAINS_GFF3)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call run-arriba,$(sample))))


..DUMMY := $(shell mkdir -p version; \
			 $(ARRIBA_EXE) -h > version/arriba.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: arriba
