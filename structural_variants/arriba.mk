include innovation-lab/Makefile.inc
include innovation-lab/config/arriba.inc

LOGDIR ?= log/arriba.$(NOW)

arriba : $(foreach sample,$(SAMPLES),arriba/$(sample)/fusions.tsv) \
	 arriba/summary.txt

define arriba
arriba/$1/$1.1.fastq.gz : bam/$1.bam
	$$(call RUN,-n 4 -s 4G -m 9G,"set -o pipefail && \
				      $$(SAMTOOLS) sort -T bam/$1 -O bam -n -@ 4 -m 6G $$(<) | \
				      bedtools bamtofastq -i - -fq >(gzip -c > arriba/$1/$1.1.fastq.gz) -fq2 >(gzip -c > arriba/$1/$1.2.fastq.gz)")

arriba/$1/$1.Aligned.out.bam : arriba/$1/$1.1.fastq.gz
	$$(call RUN,-c -n $(THREADS) -s 3G -m 4G -v $(ARRIBA_ENV) -w 72:00:00,"set -o pipefail && \
									       STAR \
									       --runThreadN $$(THREADS) \
									       --genomeDir $$(STAR_INDEX_DIR) \
									       --genomeLoad NoSharedMemory \
									       --readFilesIn arriba/$1/$1.1.fastq.gz arriba/$1/$1.2.fastq.gz \
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
									       --outFileNamePrefix arriba/$1/$1. > arriba/$1/$1.Aligned.out.bam")


arriba/$1/fusions.tsv : arriba/$1/$1.Aligned.out.bam
	$$(call RUN,-c -n 1 -s 24G -m 36G -v $(ARRIBA_ENV),"set -o pipefail && \
							    $$(ARRIBA_EXE) -x arriba/$1/$1.Aligned.out.bam \
							    -o arriba/$1/fusions.tsv \
							    -O arriba/$1/discarded.tsv \
							    -a $$(ASSEMBLY_FA) \
							    -g $$(ANNOTATION_GTF) \
							    -b $$(BLACKLIST_TSV) \
							    -k $$(KNOWN_FUSIONS_TSV) \
							    -t $$(KNOWN_FUSIONS_TSV) \
							    -p $$(PROTEIN_DOMAINS_GFF3)")

endef
$(foreach sample,$(SAMPLES),\
                $(eval $(call arriba,$(sample))))
				
arriba/summary.txt : $(foreach sample,$(SAMPLES),arriba/$(sample)/fusions.tsv)
	$(call RUN, -c -n 1 -s 16G -m 24G,"set -o pipefail && \
					   $(RSCRIPT) $(SCRIPTS_DIR)/rna_seq/summarize_arriba.R --sample_names '$(SAMPLES)'")
									   
..DUMMY := $(shell mkdir -p version; \
	     $(ARRIBA_EXE) -h > version/arriba.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: arriba
