include innovation-lab/Makefile.inc
include innovation-lab/config/arriba.inc

LOGDIR ?= log/arriba.$(NOW)

SEQ_PLATFROM = illumina
STAR_OPTS = --genomeDir $(STAR_INDEX_DIR) \
            --outSAMtype BAM SortedByCoordinate \
	    --twopassMode Basic \
            --outReadsUnmapped None \
            --outSAMunmapped Within \
            --chimSegmentMin 12 \
            --chimJunctionOverhangMin 12 \
            --alignSJDBoverhangMin 10 \
            --alignMatesGapMax 200000 \
            --alignIntronMax 200000 \
            --chimSegmentReadGapMax parameter 3 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --chimOutType WithinBAM \
	    --quantMode GeneCounts

arriba : $(foreach sample,$(SAMPLES),arriba/$(sample)/$(sample).1.fastq.gz) \
	 $(foreach sample,$(SAMPLES),arriba/$(sample)/$(sample).2.fastq.gz) \
	 $(foreach sample,$(SAMPLES),arriba/$(sample)/$(sample).Aligned.sortedByCoord.out.bam) \
	 $(foreach sample,$(SAMPLES),arriba/$(sample)/fusions.tsv) \
	 arriba/summary.txt


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
arriba/$1/$1.Aligned.sortedByCoord.out.bam : arriba/$1/$1.1.fastq.gz arriba/$1/$1.2.fastq.gz
	$$(call RUN,-c -n $(STAR_THREADS) -s 1G -m 2G -v $(ARRIBA_ENV),"set -o pipefail && \
									STAR $$(STAR_OPTS) \
									--outFileNamePrefix arriba/$1/$1. \
									--runThreadN $$(STAR_THREADS) \
									--outSAMattrRGline \"ID:$1\" \"LB:$1\" \"SM:$1\" \"PL:$${SEQ_PLATFORM}\" \
									--readFilesIn $$(<) $$(<<) \
									--readFilesCommand zcat")

arriba/$1/$1.Aligned.sortedByCoord.out.bam.bai : arriba/$1/$1.Aligned.sortedByCoord.out.bam
	$$(call RUN,-n 1 -s 2G -m 4G,"set -o pipefail && \
				      $$(SAMTOOLS) index $$(<)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call star-align,$(sample))))
		
define run-arriba
arriba/$1/fusions.tsv : arriba/$1/$1.Aligned.sortedByCoord.out.bam arriba/$1/$1.Aligned.sortedByCoord.out.bam.bai
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
		
arriba/summary.txt : $(foreach sample,$(SAMPLES),arriba/$(sample)/fusions.tsv)
	$(call RUN, -c -n 1 -s 12G -m 16G,"set -o pipefail && \
					   $(RSCRIPT) $(SCRIPTS_DIR)/rna_seq/summarize_arriba.R --sample_names '$(SAMPLES)'")


..DUMMY := $(shell mkdir -p version; \
	     $(ARRIBA_EXE) -h > version/arriba.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: arriba
